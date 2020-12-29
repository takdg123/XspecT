import os
import numpy as np
from scipy.interpolate import interp1d
from .lookup.lookup import LookupFile
from gbm.data.phaii import Cspec, Ctime, TTE
from gbm.data.drm import RSP
from gbm.data.scat import DetectorData
from gbm.data.containers import RateHisto, EventList
from .selection import Selection
from .background.background import Background
from .plugin_manager import PluginManager
from . import file
from . import time as gbmtime
from .write_xspec import writeXSPECFiles
import copy

class GspecManager:
    """Class to manager multiple data files for GSPEC.  The data model for GSPEC is a
    centralized manager to coordinate all selection and operations for all
    detectors/data files loaded into a GSPEC for a single observation.  Inherits from 
    GspecLookup which provides an observation-based lookup recording system for all 
    user selections and actions

    Attributes:
    -----------
    
    Public Methods:
    ---------------
    add_data:
        Adds a new data file
    remove_data:
        Removes a data file
    add_response:
        Adds a new response file associated with loaded data
    time_binning:
        Set the time binning
    energy_binning:
        Set the energy binning
    add_background:
        Add a background selection and fitting algorithm
    remove_background:
        Remove a background selection and fit
    source_selection:
        Set the source selection
    remove_source_selections:
        Remove all source selections
    source_by_snr:
        Select signal above a S/N ratio
    energy_selection:
        Set the energy selection
    remove_energy_selections:
        Remove all energy selections
    import_time_binning:
        Import time binning to other data files from a given data file
    import_energy_binning:
        Import energy binning to other data files from a given data file
    import_background:
        Import background to other data files from a given data file
    import_source_selection:
        Import source selection to other data files from a given data file
    import_time_display:
        Import time display to other data files from a given data file
    import_energy_display:
        Import energy display to other data files from a given data file
    clear_time_binning:
        Reset the time binning for a detector to the native resolution
    joint_snr_binning:
        Apply a time binning to all detectors based on the joint SNR for a subset 
        of detectors
    refresh:
        Force a recalculation of the user selections and background
    export_to_xspec:
        Export the source selection to XSPEC-readable files
    display_gui:
        Display the root GSPEC GUI
    display_detector:
        Display the data for detector that has already been loaded
    
    
    Static Methods:
    ---------------
    """

    def __init__(self, *args):
        """Initializes data manager and centralized lookup

        Parameters:
        -----------
        *args:
            Arguments to pass to GspecLookup
        """
        super(GspecManager, self).__init__(*args)
        self.lookup = LookupFile()
        self.data = {}
        self._plugins = PluginManager()
        
        # default energy ranges - cut out under and overflow
        self.ctime_nai_default = [20.0, 900.0]
        self.ctime_bgo_default = [500.0, 38000.0]
        self.cspec_nai_default = [9.0, 900.0]
        self.cspec_bgo_default = [250.0, 38000.0]
        self.cspec_lle_default = [31000.0, 100000.0]
        
        self._openWindows = {}
        self._export_files = None
        self.export_dirty = True
        self.lookup_dirty = False

    def add_data(self, filename, lookup=None, response=None, no_tte_binning=False):
        """Add a data file to the data manager.  This will apply a given lookup file to 
        the data and create all of the relevant selections, background fit, and binning.
        Can also load an associated response file.
        
        Parameters:
        -----------
        filename: str
            The filename of the file to load
        lookup: RmfitLookup or GspecLookup, optional
            The lookup file
        response: str, optional
            The filename of the response file associated with the data
        
        Returns:
        -----------
        dataname: str
            The name of the file added, to be used as a key for all other operations
        """
        # load data
        f = file.GbmFile.from_path(filename)
        dataname = f.basename()
        detector = f.detector
        datatype = f.data_type
        if datatype == 'cspec':
            data = Cspec(filename)
        elif datatype == 'ctime':
            data = Ctime(filename)
        elif datatype == 'tte':
            data = TTE(filename)
        else:
            raise ValueError('{0} is not a valid data type'.format(dataname))

        # intialize the data dictionary
        self.data[dataname] = {'data': data, 'fullview': None, 'timeview': None,
                                'energyview': None, 'bkgdview': None,
                                'bkgdmodel': None, 'sourceview': None,
                                'response': None}
        self.lookup.add_data_file(filename)
        self.lookup_dirty = True


        # set the full data view
        if not(detector.is_bgo()) and not(detector.is_nai()):
            self._full_selection_lle(dataname)

        else:
            self._full_selection(dataname)

        # set default energy range
        if detector.is_nai():
            if datatype == 'ctime':
                self.energy_selection(dataname, self.ctime_nai_default)
            else:
                self.energy_selection(dataname, self.cspec_nai_default)
        elif detector.is_bgo():
            if datatype == 'ctime':
                self.energy_selection(dataname, self.ctime_bgo_default)
            else:
                self.energy_selection(dataname, self.cspec_bgo_default)
        else:
            if datatype == 'ctime':
                self.energy_selection(dataname, self.cspec_lle_default)

        # default binning for TTE data
        if (data.type == 'unbinned') and (self.lookup[dataname].binnings.time is None):
            if not no_tte_binning:
                self.time_binning(dataname, 'Temporal Resolution', 1.024, 
                                  datatype='unbinned', time_ref=0.0, nosave=True)

        # read in lookup if provided
        if lookup is not None:
            self.load_gspec_lookup(dataname, lookup)

        # if lookup is not None:
        #     if isinstance(lookup, RmfitLookup):
        #         self.load_rmfit_lookup(dataname, lookup)
        #     elif isinstance(lookup, GspecLookup):
        #         self.load_gspec_lookup(dataname, lookup)
        #     self.lookup_dirty = False

        # read in response if provided
        if response is not None:
            self.add_response(dataname, response)
        
        self.export_dirty = True
        
        return dataname

    def remove_data(self, dataname):
        """Remove a data file from data manager.
        
        Parameters:
        -----------
        dataname: str
            The data filename to remove
        """
        self.lookup.assert_has_datafile(dataname)
        del self.lookup[dataname]
        del self.data[dataname]
        self.export_dirty = True
        self.lookup_dirty = True

    def add_response(self, dataname, rsp_filename):
        """Add a response file and associate it with loaded data
        
        Parameters:
        -----------
        dataname: str
            The data filename that the response is to be associated with
        rsp_filename: str
            The filename of the response file associated with the data
        """
        self.lookup.assert_has_datafile(dataname)
        # ensure that response matches the detector and data type
        rsp = RSP(rsp_filename)
        det = self.data[dataname]['data'].detector
        if rsp.detector != det:
            print('Warning: Requested response ({0}) is not valid for this '
                  'detector ({1})'.format(rsp.detector, det))
        datatype = self.data[dataname]['data'].header['PRIMARY']['DATATYPE']
        if (rsp.data_type == 'CSPEC') & ((datatype != 'CSPEC') & (datatype != 'TTE')):
            raise ValueError('Requested response type (CSPEC) is not valid for '
                             'this datatype ({0})'.format(datatype))
        if (rsp.data_type == 'CTIME') & (datatype != 'CTIME'):
            raise ValueError('Requested response type (CTIME) is not valid for '
                             'this datatype ({0})'.format(datatype))

        self.data[dataname]['response'] = rsp
        self.lookup[dataname].add_response(rsp_filename)
        self.export_dirty = True
        self.lookup_dirty = True

    def time_binning(self, dataname, binning_name, *args, datatype='binned', 
                     start=None, stop=None, nosave=False, **kwargs):
        """Apply a time binning function to the data file
        
        Parameters:
        -----------
        dataname: str
            The data filename to which the binning is applied
        binning_name: str
            The name of the binning function
        *args:
            Arguments to be passed to the binning function
        datatype: str, optional
            Either 'binned' (default) or 'unbinned'
        start: float, optional
            The start of the data range to be rebinned. The default is to start at the 
            beginning of the data.
        stop: float, optional
            The end of the data range to be rebinned. The default is to stop at
            the end of the data.
        nosave: bool, optional
            If True, do not save to lookup.  Default is False.
        **kwargs:
            Keyword arguments to be passed to the binning function
        """
        # make sure function exists and import
        binning_function = self._plugins.get_binning_method(binning_name,
                                                            datatype=datatype)
        # update lookup
        if not nosave:
            self.lookup[dataname].add_time_binning(binning_name, datatype, *args, start=start, stop=stop, **kwargs)

        # pass-through for single detector SNR binning
        if binning_name == 'Signal-to-Noise':
            if start is None:
                start = self.data[dataname]['data'].time_range[0]
            if stop is None:
                stop = self.data[dataname]['data'].time_range[1]
            self.joint_snr_binning([dataname], [start, stop], *args, datatype=datatype)
            return
        
        # FIXME: This is a temporary fix to get this binning function to work
        # on TTE data.  This needs to be handled better.
        if binning_name=='Combine By Factor' and datatype=='unbinned':
            edges = self.data[dataname]['fullview'].lightcurve[0].edges
            args = list(args)
            args.insert(0, edges)
            args = tuple(args)

        self.data[dataname]['fullview'].set_time_binning(binning_function, *args,
                                                     start=start, stop=stop,
                                                     **kwargs)
        timeview = self.data[dataname]['fullview'].slice(energy_selections=self.lookup[dataname].selections.energy)
        self.data[dataname]['timeview'] = timeview
        energyview = self.data[dataname]['fullview'].slice(time_selections=self.lookup[dataname].selections.source)
        self.data[dataname]['energyview'] = energyview
        
        self.export_dirty = True
        self.lookup_dirty = True

    def energy_binning(self, dataname, binning_name, *args, start=None,
                       stop=None, nosave=False, **kwargs):
        """Apply an energy binning function to the data file
        
        Parameters:
        -----------
        dataname: str
            The data filename to which the binning is applied
        binning_name: str
            The name of the binning function
        *args:
            Arguments to be passed to the binning function
        start: float, optional
            The start of the data range to be rebinned. The default is to start at the 
            beginning of the data.
        stop: float, optional
            The end of the data range to be rebinned. The default is to stop at
            the end of the data.
        nosave: bool, optional
            If True, do not save to lookup.  Default is False.
        **kwargs:
            Keyword arguments to be passed to the binning function
        """
        # make sure the function exists and import
        binning_function = self._plugins.get_binning_method(binning_name)

        # update lookup
        if not nosave:
            self.lookup[dataname].add_energy_binning(binning_name, *args, start=start, stop=stop, **kwargs)

        self.data[dataname]['fullview'].set_energy_binning(binning_function, 
                                                            *args, start=start, 
                                                            stop=stop, **kwargs)
        timeview = self.data[dataname]['fullview'].slice(energy_selections=self.lookup[dataname].selections.energy)
        self.data[dataname]['timeview'] = timeview
        energyview = self.data[dataname]['fullview'].slice(time_selections=self.lookup[dataname].selections.source)
        self.data[dataname]['energyview'] = energyview
        
        self.export_dirty = True
        self.lookup_dirty = True

    def add_background(self, dataname, background_intervals, background_name,
                       *args, datatype='binned', **kwargs):
        """Apply a background selection and fit to the data
        
        Parameters:
        -----------
        dataname: str
            The data filename to which the background is fit/applied
        background_intervals: list
            The background selection intervals
        background_name: str
            The name of background fitting/estimation method
        datatype: str, optional
            Either 'binned' (default) or 'unbinned'
        *args:
            Arguments to be passed to the background fitting/estimation class
        **kwargs:
            Keyword arguments to be passed to the background fitting/estimation class
        """
        # make sure function exists and import
        background_class = self._plugins.get_background_method(background_name,
                                                               datatype=datatype)

        self.lookup[dataname].add_background_selection(background_intervals)
        self.lookup[dataname].add_background_model(background_name, datatype, *args, **kwargs)

        # create the background selection
        bkgd_select = self.data[dataname]['fullview'].slice(time_selections=background_intervals,
                                                            energy_selections=self.lookup[dataname].selections.energy)
        self.data[dataname]['bkgdview'] = bkgd_select

        # perform the background fit/estimate
        bkgd = Background(bkgd_select)
        bkgd.fit(background_class, *args, datatype=datatype)
        self.data[dataname]['bkgdmodel'] = bkgd
        
        self.export_dirty = True
        self.lookup_dirty = True

    def remove_background(self, dataname):
        """Remove the background selections and fit
        
        Parameters:
        -----------
        dataname: str
            The data filename for which the background is to be removed
        """
        self.lookup.assert_has_datafile(dataname)
        self.lookup[dataname].selections.background = None
        self.lookup[dataname].background = None
        self.data[dataname]['bkgdmodel'] = None
        self.data[dataname]['bkgdview'] = None
        self.lookup_dirty = True

    def source_selection(self, dataname, source_intervals, append=False):
        """Apply a source selection to the data file
        
        Parameters:
        -----------
        dataname: str
            The data filename for which the source selection is made
        source_intervals: list
            A list of source interval selections
        append: bool, optional
            If True, then append the source_intervals to the existing intervals.
            Default is False.
        """
        if len(source_intervals) == 0:
            return
        
        if append:
            selections = self.lookup[dataname].selections.source
            if selections is None:
                pass
            elif isinstance(selections[0], (tuple, list)):
                selections.append(source_intervals)
                source_intervals = selections
            else:
                selections = [selections]
                selections.append(source_intervals)
                source_intervals = selections

        self.lookup[dataname].add_source_selection(source_intervals)
        source_select = self.data[dataname]['fullview'].slice(time_selections=source_intervals,
                                                              energy_selections=self.lookup[dataname].selections.energy)
        self.data[dataname]['sourceview'] = source_select
        self.data[dataname]['energyview'] = self.data[dataname]['fullview'].slice(time_selections=source_intervals)
        
        self.export_dirty = True
        self.lookup_dirty = True

    def remove_source_selections(self, dataname):
        """Remove all source selections
        
        Parameters:
        -----------
        dataname: str
            The data filename for which the source selection is made
        """
        self.lookup.assert_has_datafile(dataname)
        self.lookup[dataname].selections.source = None
        self.data[dataname]['sourceview'] = None
        self.lookup_dirty = True

    def source_by_snr(self, dataname, time_intervals, snr):
        """Select signal above a defined signal-to-noise ratio.  This requires
        a background to have been fit.
        
        Parameters:
        -----------
        dataname: str
            The data filename for which the source selection is made
        time_intervals: list
            A list of source interval selections
        snr: float
            The signal-to-noise ratio threshold
        """
        self.lookup.assert_has_datafile(dataname)

        # ensure background exists
        bkgd = self.data[dataname]['bkgdmodel']
        if bkgd is None:
            raise ValueError('Background has not been defined')

        # get the observed rates
        fullview = self.data[dataname]['fullview']
        selection = fullview.slice(time_selections=time_intervals,
                                   energy_selections=self.lookup[dataname].selections.energy)
        rates = RateHisto.merge(selection.lightcurve)

        # background rate
        bkgd_rates, _ = self.data[dataname]['bkgdmodel'].integrate_energy(selection)

        # convert to background counts
        bkgd_counts = bkgd_rates * rates.exposure

        # calculate SNR and mask where it exceeds the SNR threshold
        lc_snr = (rates.values - bkgd_counts) / np.sqrt(bkgd_counts)
        mask = (lc_snr >= snr)

        if np.sum(mask) == 0:
            print('No bins meet a S/N of {0}'.format(snr))
            return

        # reduce the bounds for the source selection
        bounds = (rates.bounds[0][mask], rates.bounds[1][mask])
        reduced_bounds = []
        lo = bounds[0][0]
        for i in range(len(bounds[0]) - 1):
            if bounds[0][i + 1] == bounds[1][i]:
                continue
            else:
                reduced_bounds.append((lo, bounds[1][i]))
                lo = bounds[0][i + 1]
        reduced_bounds.append((lo, bounds[1][-1]))

        # pass to source selection
        self.source_selection(dataname, reduced_bounds)

    def energy_selection(self, dataname, energy_intervals, append=False):
        """Apply an energy selection to the data file
        
        Parameters:
        -----------
        dataname: str
            The data filename for which the energy selection is made
        energy_intervals: list
            A list of energy interval selections
        append: bool, optional
            If True, then append the energy_intervals to the existing intervals.
            Default is False.
        """
        if append:
            selections = self.lookup[dataname].selections.energy
            if selections is None:
                pass
            elif isinstance(selections[0], (tuple, list)):
                selections.append(energy_intervals)
                energy_intervals = selections
            else:
                selections = [selections]
                selections.append(energy_intervals)
                energy_intervals = selections

        self.lookup[dataname].selections.energy = energy_intervals
        timeview = self.data[dataname]['fullview'].slice(energy_selections=energy_intervals)
        self.data[dataname]['timeview'] = timeview
        
        self.export_dirty = True
        self.lookup_dirty = True

    def remove_energy_selections(self, dataname):
        """Remove all energy selections
        
        Parameters:
        -----------
        dataname: str
            The data filename for which the source selection is made
        """
        self.lookup.assert_has_datafile(dataname)
        self.lookup[dataname].selections.energy = None
        timeview = self.data[dataname]['fullview'].slice()
        self.data[dataname]['timeview'] = timeview
        
        self.export_dirty = True
        self.lookup_dirty = True

    # TODO: Rename parameters
    def import_time_binning(self, export_dataname, import_datanames):
        """Import time binning to other data files from a given data file
        
        Parameters:
        -----------
        export_dataname: str
            The data file from which to export the time binning
        import_detanames: list
            The data filename(s) to which the time binning is being imported
        """
        self.lookup.assert_has_datafile(export_dataname)
        binning = self.lookup[export_dataname].binnings.time
        if binning is None:
            return
        for dataname in import_datanames:
            self.lookup.assert_has_datafile(dataname)
            for abinning in binning:
                self.time_binning(dataname, abinning.method, *abinning.args, datatype=abinning.datatype,
                                  **abinning.kwargs)

    def import_energy_binning(self, export_dataname, import_datanames):
        """Import energy binning to other data files from a given data file
        
        Parameters:
        -----------
        export_dataname: str
            The data filename from which to export the energy binning
        import_datanames: list
            The data filename(s) to which the energy binning is being imported
        """
        self.lookup.assert_has_datafile(export_dataname)
        binning = self.lookup[export_dataname].binnings.energy
        if binning is None:
            return
        for dataname in import_datanames:
            self.lookup.assert_has_datafile(dataname)
            for abinning in binning:
                self.energy_binning(dataname, abinning['method'], *abinning['args'],
                                    **abinning['kwargs'])

    def import_background(self, export_dataname, import_datanames):
        """Import time background selection and fit to other data files from a given 
        data file
        
        Parameters:
        -----------
        export_dataname: str
            The data filename from which to export the background
        import_datanames: list
            The data filename(s) to which the background is being imported
        """
        self.lookup.assert_has_datafile(export_dataname)
        bkgd = self.data[export_dataname]['bkgdview']
        if bkgd is None:
            return
        bkgd_intervals = self.lookup[export_dataname].selections.background
        bkgd = self.lookup[export_dataname].background
        if bkgd is not None:
            for dataname in import_datanames:
                self.lookup.assert_has_datafile(dataname)

                self.add_background(dataname, bkgd_intervals, bkgd.method, *bkgd.args, datatype=bkgd.datatype,
                                    **bkgd.kwargs)
        else:
            print("Warning: no available backround to push.")

    def import_source_selection(self, export_dataname, import_datanames):
        """Import source selection to other data files from a given data file
        
        Parameters:
        -----------
        export_dataname: str
            The data filename from which to export the source selection
        import_datanames: list
            The data filename(s) to which the source selection is being imported
        """
        self.lookup.assert_has_datafile(export_dataname)
        source = self.data[export_dataname]['sourceview']
        if source is None:
            return
        source_intervals = self.lookup[export_dataname].selections.source
        for dataname in import_datanames:
            self.lookup.assert_has_datafile(dataname)
            self.source_selection(dataname, source_intervals)

    def import_time_display_view(self, export_dataname, import_datanames):
        """Import time display to other data files from a given data file
        
        Parameters:
        -----------
        export_dataname: str
            The data filename from which to export the time display
        import_datanames: list
            The data filename(s) to which the time display is being imported
        """
        self.lookup.assert_has_datafile(export_dataname)
        time_view = self.lookup[export_dataname].views.time
        for dataname in import_datanames:
            self.lookup.assert_has_datafile(dataname)
            self.lookup[dataname].add_time_display_view(time_view.to_list())
        self.lookup_dirty = True

    def import_energy_display_view(self, export_dataname, import_datanames):
        """Import energy display to other data files from a given data file
        
        Parameters:
        -----------
        export_dataname: str
            The data filename from which to export the energy display
        import_datanames: list
            The data filename(s) to which the energy display is being imported
        """
        self.lookup.assert_has_datafile(export_dataname)
        energy_view = self.lookup[export_dataname].views.energy
        for dataname in import_datanames:
            self.lookup.assert_has_datafile(dataname)
            self.lookup[dataname].add_energy_display_view(energy_view)
        self.lookup_dirty = True

    def clear_time_binning(self, dataname):
        """Clear the time binning for a data file. 
        Results in the native binning resolution
        
        Parameters:
        -----------
        dataname: str
            The data filename
        """
        self.lookup.assert_has_datafile(dataname)
        self.lookup[dataname].binnings.time = None
        self.data[dataname]['fullview'].set_time_binning(None)

        # if TTE, default to the 64 ms time binning
        if self.data[dataname]['data'].type == 'unbinned':
            self.time_binning(dataname, 'Temporal Resolution', 0.064, 
                              datatype='unbinned', time_ref=0.0, nosave=True)
        else:
            # otherwise, show the raw binning
            timeview = self.data[dataname]['fullview'].slice(energy_selections=self.lookup[dataname].selections.energy)
            self.data[dataname]['timeview'] = timeview
            energyview = self.data[dataname]['fullview'].slice(time_selections=self.lookup[dataname].selections.source)
            self.data[dataname]['energyview'] = energyview
        
        self.export_dirty = True
        self.lookup_dirty = True

    def clear_energy_binning(self, dataname):
        """Clear the energy binning for a data file. 
        Results in the native binning resolution
        
        Parameters:
        -----------
        dataname: str
            The data filename
        """
        self.lookup.assert_has_datafile(dataname)
        self.lookup[dataname].binnings.energy = None
        self.data[dataname]['fullview'].set_energy_binning(None)

        timeview = self.data[dataname]['fullview'].slice(energy_selections=self.lookup[dataname].selections.energy)
        self.data[dataname]['timeview'] = timeview
        energyview = self.data[dataname]['fullview'].slice(time_selections=self.lookup[dataname].selections.source)
        self.data[dataname]['energyview'] = energyview
        
        self.export_dirty = True
        self.lookup_dirty = True

    def joint_snr_binning(self, datanames, time_intervals, snr, datatype='binned'):
        """Bin all data files based on the joint SNR of the selected detectors.  A subset 
        of all detectors may be provided for the SNR binning, but the binning will be 
        applied to all detectors that are currently loaded.
        
        Parameters:
        -----------
        datanames: list
            The data files to be used to estimate the joint SNR
        time_intervals: list
            The time intervals over which to bin
        snr: float
            The SNR threshold
        """
        if datatype == 'binned':
            rates = []
            bkgd_rates = []
            for dataname in datanames:
                self.lookup.assert_has_datafile(dataname)

                # ensure background exists
                bkgd = self.data[dataname]['bkgdmodel']
                if bkgd is None:
                    raise ValueError('Background has not been defined')

                # get the observed rates
                selection = Selection(self.data[dataname]['data'], 
                                      time_selections=time_intervals,
                                      energy_selections=self.lookup[dataname].selections.energy)
                rates.append(RateHisto.merge(selection.lightcurve))

                # get the background rates
                bkgd_rate, _ = self.data[dataname]['bkgdmodel'].integrate_energy(selection)
                bkgd_rates.append(bkgd_rate)

            # convert to background counts and sum over all detectors
            bkgd_counts = np.zeros(bkgd_rates[0].shape)
            for i in range(len(rates)):
                bkgd_counts += bkgd_rates[i] * rates[i].exposure
            rates = RateHisto.sum_time(rates)

            # the snr binning
            snr_binning = self._plugins.get_binning_method('Signal-to-Noise',
                                                           datatype=datatype)
            binned_rates = rates.rebin(snr_binning, bkgd_counts, snr)

            # set the binning for each detector
            all_datanames = self.data.keys()
            for dataname in all_datanames:
                self.time_binning(dataname, 'By Time Edge', binned_rates.edges,
                                  datatype=datatype, start=time_intervals[0],
                                  stop=time_intervals[-1])

        elif datatype == 'unbinned':
            events = []
            bkgd_rates = []
            bkgd_times = []
            for dataname in datanames:
                self.lookup.assert_has_datafile(dataname)
                # get the events
                selection = Selection(self.data[dataname]['data'], 
                                      time_selections=time_intervals,
                                      energy_selections=self.lookup[dataname].selections.energy)
                events.append(EventList.merge(selection.lightcurve))

                # get the background rates
                fullview = self.data[dataname]['fullview']
                selection = fullview.slice(time_selections=time_intervals,
                                           energy_selections=self.lookup[dataname].selections.energy)
                bkgd_rate, _ = self.data[dataname]['bkgdmodel'].integrate_energy(selection)
                bkgd_rates.append(bkgd_rate)
                bkgd_times.append(RateHisto.merge(selection.lightcurve).bin_centers())

            # merge and time sort all of the events
            events = EventList.merge(events, sort_axis='TIME')

            # get the background rates summed over all detectors
            bkgd_rate = np.zeros(events.size)
            for idata in range(len(datanames)):
                rate_interpolator = interp1d(bkgd_times[idata], bkgd_rates[idata],
                                             fill_value='extrapolate')
                bkgd_rate += rate_interpolator(events._events['TIME'])

            # the snr binning
            snr_binning = self._plugins.get_binning_method('Signal-to-Noise',
                                                           datatype=datatype)
            binned_rates = events.bin(snr_binning, bkgd_rate, snr)

            # set the binning for each detector
            all_datanames = self.data.keys()
            for dataname in all_datanames:
                self.time_binning(dataname, 'By Time Edge', binned_rates.edges,
                                  datatype=datatype, start=time_intervals[0],
                                  stop=time_intervals[-1])

    def refresh(self, dataname, background=True):
        """Force a recalculation of the user selections and background
        
        Parameters:
        -----------
        dataname: str
            The data filename to refresh
        background: bool, optional
            If True (default) then recalculate background
        """
        self.lookup.assert_has_datafile(dataname)

        # if we have a background, refresh the model fit
        if (self.data[dataname]['bkgdview'] is not None) & background:
            self.add_background(dataname, self.lookup[dataname].selections.background,
                                self.lookup[dataname].background.method,
                                *self.lookup[dataname].background.args,
                                datatype=self.lookup[dataname].background.datatype,
                                **self.lookup[dataname].background.kwargs)

        # if we have a source selection, refresh the source view
        if self.data[dataname]['sourceview'] is not None:
            self.source_selection(dataname, self.lookup[dataname].selections.source)

    def export_to_xspec(self, dataname, directory=None, MULTI=False, verbose=False, **kwargs):
        """Write out (PHAI or PHAII), BAK, and (RSP or RSP2) files for XSPEC
        
        Parameters:
        -----------
        dataname: str
            The data filename to export
        directory: str, optional
            The directory to write the files to.  If not set, then will write
            to the directory where the original files exist
        
        Returns:
        -----------
        phafilename: str
            The output PHAI filename
        """
        self.lookup.assert_has_datafile(dataname)

        # make sure we have all of the pieces we need
        sourceview = self.data[dataname]['sourceview']
        if sourceview is None:
            raise ValueError('Source selection has not been defined')
        rsp = self.data[dataname]['response']
        if rsp is None:
            raise ValueError('Response file has not been defined')
        bkgd = self.data[dataname]['bkgdmodel']
        if bkgd is None:
            raise ValueError('Background has not been defined')

        # create the filenames
        datafile = file.GbmFile.from_path(self.data[dataname]['data'].filename)
        if directory is None:
            subdir = gbmtime.Met.now().datetime.strftime("%Y-%m-%dT%H-%M-%S")
            directory = os.path.join(datafile.directory, subdir)

            # Check to see if the directory already exists
            if os.path.isdir(directory) == False:
                os.mkdir(directory)


        bakfile = file.GbmFile.create(uid=datafile.uid,
                                      data_type=datafile.data_type,
                                      detector=datafile.detector,
                                      trigger=datafile.trigger,
                                      version=datafile.version,
                                      extension='bak', meta='_xspec',
                                      directory=directory)
        bakfilename = bakfile.path()
        phafile = bakfile
        phafile.extension = 'pha'
        phafilename = phafile.path()
        rspfile = phafile
        rspfile.extension = 'rsp'
        rspfilename = rspfile.path()

        # and here we need to hook to writeXSPEC
        writeXSPECFiles(sourceview, bkgd, rsp, phafilename, bakfilename, rspfilename, MULTI=MULTI, **kwargs)

        if verbose == True:
            print(phafilename)
            print(rspfilename)        
            print(bakfilename)

        return phafilename

    def save_lookup(self, filename=None):
        """Write out a GSPEC lookup file
        
        Parameters:
        -----------
        filename: str, optional
            The filename of the lookup file.  If not set, then a standard
            name will be produced based on the original input filenames.
        
        Returns:
        -----------
        filename: str
            The full path + filename of the written lookup file        
        """
        if filename is None:
            datanames = list(self.data.keys())
            datafile = file.GbmFile.from_path(self.data[datanames[0]]['data'].filename)
            directory = datafile.directory
            uid = datafile.uid
            data_type = datafile.data_type
            trigger = datafile.trigger
            meta = '_lookup'
            extension = 'json'
            lufile = file.GbmFile.create(directory=directory, uid=uid,
                                         data_type=data_type, trigger=trigger,
                                         meta=meta, extension=extension)
            filename = lufile.path()

        self.lookup.write_to(filename)
        self.lookup_dirty = False
        return filename

    def load_rmfit_lookup(self, dataname, lookup):
        """Apply an RMFit lookup file to the data file
        
        Parameters:
        -----------
        dataname: str
            The data filename to which the lookup is applied
        lookup: RmfitLookup
            The RMFit lookup object
        """
        data_type = type(self.data[dataname]['data']).__name__

        # apply energy selection
        if lookup.has_energy_selection:
            energy_ranges = lookup.energy_selection()
            energy_ranges = [e.tolist() for e in energy_ranges]
            self.energy_selection(dataname, energy_ranges)

        # apply time binning
        if lookup.is_time_rebinned:
            if (data_type == 'Ctime') | (data_type == 'Cspec'):
                bin_method = 'By Edge Index'
                arg = lookup.rebinned_time_edges()
                datatype = 'binned'
            elif data_type == 'TTE':
                bin_method = 'By Time Edge'
                arg = lookup.tte_edges()
                datatype = 'unbinned'
            self.time_binning(dataname, bin_method, arg, datatype=datatype)

        # apply energy binning
        if lookup.is_energy_rebinned:
            bin_method = 'By Edge Index'
            self.energy_binning(dataname, bin_method, lookup.rebinned_energy_edges())

        # apply background fit
        if lookup.has_background_selection:
            background_ranges = lookup.background_selection()
            background_ranges = [t.tolist() for t in background_ranges]
            self.add_background(dataname, background_ranges, 'Polynomial',
                                lookup.background_poly_order)

        # apply source selection
        if lookup.has_time_selection:
            time_ranges = lookup.time_selection()
            time_ranges = [t.tolist() for t in time_ranges]
            self.source_selection(dataname, time_ranges)

        # apply window views
        self.lookup[dataname].add_time_display_view(lookup.time_display_range)
        self.lookup[dataname].add_energy_display_view(lookup.energy_display_range)

    def load_gspec_lookup(self, dataname, lookup):
        """Apply an Gspec lookup file to the data file
        
        Parameters:
        -----------
        dataname: str
            The data filename to which the lookup is applied
        lookup: GspecLookup
            The Gspec lookup object
        """
        # the gspec lookup may contain more than one detector, so only load that
        # detector's lookup
        self.lookup.merge_lookup(lookup.split_off_dataname(dataname), overwrite=True)

        # add and load response
        responsefile = self.lookup[dataname].response
        if responsefile is not None:
            dirname = os.path.dirname(self.data[dataname]['data'].filename)
            self.add_response(dataname, os.path.join(dirname, responsefile))
        
        # apply energy selection
        energy_ranges = self.lookup[dataname].selections.energy
        if energy_ranges is not None:
            self.energy_selection(dataname, energy_ranges)        
        
        # apply the background model
        bkgd = self.lookup[dataname].background
        bkgd_select = self.lookup[dataname].selections.background
        if bkgd is not None:
            self.add_background(dataname, bkgd_select, bkgd.method, *bkgd.args, datatype=bkgd.datatype, **bkgd.kwargs)

        # apply any special time binning
        tb = copy.deepcopy(self.lookup[dataname].binnings.time)
        if tb is not None:
            for onetb in tb:
                self.time_binning(dataname, onetb.method, *onetb.args, datatype=onetb.datatype, start=onetb.start,
                                  stop=onetb.stop, **onetb.kwargs, nosave=True)

        # apply any special energy binning
        eb = copy.deepcopy(self.lookup[dataname].binnings.energy)
        if eb is not None:
            for oneeb in eb:
                self.energy_binning(dataname, oneeb.method, *oneeb.args, start=oneeb.start, stop=oneeb.stop,
                                    **oneeb.kwargs)

        # apply source selection
        source = self.lookup[dataname].selections.source
        if source is not None:
            self.source_selection(dataname, source)

    def get_detector_data(self, dataname):
        """Get the detector info for the given data file

        Parameters:
        -----------
        dataname: str
            The data filename to which the lookup is applied

        Returns:
        -----------
        detdata: DetectorData
            The detector info
        """
        instrument = self.data[dataname]['data'].header['PRIMARY']['INSTRUME']
        detector = self.data[dataname]['data'].detector
        datatype = self.data[dataname]['data'].data_type
        filename = os.path.basename(self.data[dataname]['data'].filename)
        response = self.lookup[dataname].response
        time_range = self.data[dataname]['sourceview'].time_bound
        energy_range = self.data[dataname]['sourceview'].energy_bound
        ebounds = self.data[dataname]['data'].ebounds
        mask = (ebounds['E_MIN'] <= energy_range[0]) & (ebounds['E_MAX'] > energy_range[0])
        startchan = ebounds['CHANNEL'][mask][0]
        mask = (ebounds['E_MIN'] <= energy_range[1]) & (ebounds['E_MAX'] > energy_range[1])
        endchan = ebounds['CHANNEL'][mask][0]
        channel_range = (startchan, endchan)
        energy_edges = np.append(ebounds['E_MIN'], ebounds['E_MAX'][-1])
        detdata = DetectorData.create(instrument=instrument, detector=detector,
                                      datatype=datatype, filename=filename,
                                      response=response, time_range=time_range,
                                      energy_range=energy_range,
                                      channel_range=channel_range,
                                      energy_edges=energy_edges)
        return detdata

    def _full_selection(self, dataname):
        """The full view of the data file
        
        Parameters:
        -----------
        dataname: str
            The data filename from which the view is made
        """
        fullview = Selection(self.data[dataname]['data'])
        timeview = fullview.slice(energy_selections=self.lookup[dataname].selections.energy)
        energyview = fullview.slice(time_selections=self.lookup[dataname].selections.source)
        self.data[dataname]['fullview'] = fullview
        self.data[dataname]['timeview'] = timeview
        self.data[dataname]['energyview'] = energyview
   
    def _full_selection_lle(self, dataname):
        """The full view of the data file
        
        Parameters:
        -----------
        dataname: str
            The data filename from which the view is made
        """
        fullview = Selection(self.data[dataname]['data'], energy_selections=self.cspec_lle_default)
        timeview = fullview.slice(energy_selections=self.lookup[dataname].selections.energy)
        energyview = fullview.slice(time_selections=self.lookup[dataname].selections.source)
        self.data[dataname]['fullview'] = fullview
        self.data[dataname]['timeview'] = timeview
        self.data[dataname]['energyview'] = energyview
