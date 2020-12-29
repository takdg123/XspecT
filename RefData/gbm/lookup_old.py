import numpy as np
import os
import warnings
import json
from collections import OrderedDict
from .file import GbmFile

class GspecLookup(object):
    """Class for an Gspec lookup file
    
    The lookup file contains the time selection, energy selection, background model and 
    parameters, time and energy binning function and parameters, and window display ranges

    Attributes:
    -----------
    
    Public Methods:
    -----------
    add_data_file:
        Add a data file to the lookup
    add_response_filename:
        Add a response filename to the lookup
    add_background_model:
        Add a background model
    add_time_binning:
        Add a time binning function
    add_energy_binning:
        Add an energy binning function
    add_source_selection:
        Add source selection(s)
    add_energy_selection:
        Add energy selection(s)
    add_display_view:
        Add the range of the default display
    read_lookup:
        Read and load an existing lookup
    display_lookup:
        Pretty print a lookup for display (in json format)
    write_lookup:
        Write the current lookup to disk
    split_off_detector:
        Create a new lookup for a single detector within the current lookup
    merge_lookup:
        Merge a new lookup into the existing lookup
    """ 
    def __init__(self, filename=None):
        """
        Parameters:
        --------------
        filename: str, optional
            If provided, the filename to load
        """
        self._lookup_dict = OrderedDict()
        self._export_options = ['npy', 'json']
        
        # If given an existing filename, read and load the lookup
        if filename is not None:
            filetype = filename.split('.')[-1]
            self.read_lookup(filename, filetype=filetype)
    
    def add_data_file(self, fulldataname):
        """Add a detector to the lookup
        
        Parameters:
        --------------
        fulldataname: str
            The full data filename
        """
        file = GbmFile.from_path(fulldataname)
        one_lookup = [('detector', file.detector.short_name),
                      ('response', None),
                      ('background',  None), 
                      ('time_binning', None), 
                      ('energy_binning', None), 
                      ('source_selection', None), 
                      ('energy_selection', None), 
                      ('background_selection', None),
                      ('time_display_view', None), 
                      ('energy_display_view', None)] 
        self._lookup_dict[file.basename()] = OrderedDict(one_lookup)
    
    def add_response_filename(self, dataname, rsp_filename):
        """Add a response file for the data
        
        Parameters:
        --------------
        dataname: str
            The the data filename
        rsp_filename: str
            The filename of the response file
        """
        self._check_dataname(dataname)
        self._lookup_dict[dataname]['response'] = os.path.basename(rsp_filename)
    
    def add_background_model(self, dataname, background_name, datatype, *args, **kwargs):
        """Add a new background model for the data file
        
        Parameters:
        --------------
        dataname: str
            The data filename
        background_class: str
            The background fitting/estimation name
        datatype: str
            The datatype the background is applied to. Either 'binned' or 'unbinned'
        *args: 
            Additional arguments used by the background class
        **kwargs:
            Additional keywords used by the background class
        """
        self._check_dataname(dataname)
        self._check_datatype(datatype)
        bkgd_dict = [('method', background_name),
                     ('datatype', datatype),
                     ('args', args),
                     ('kwargs', kwargs)]
        self._lookup_dict[dataname]['background'] = OrderedDict(bkgd_dict)
        
    def add_time_binning(self, dataname, binning_name, datatype, *args, start=None, 
                         stop=None, **kwargs):
        """Add a new time binning function for the data file
        
        Parameters:
        --------------
        dataname: str
            The data filename
        binning_function: str
            The binning function name
        datatype: str
            The datatype the binning is applied to. Either 'binned' or 'unbinned'
        *args: 
            Additional arguments used by the binning function
        start: float, optional
            The start of the data range to be rebinned. The default is to start at the 
            beginning of the RateHisto segment.
        stop: float, optional
            The end of the data range to be rebinned. The default is to stop at
            the end of the RateHisto segment.
        **kwargs:
            Additional keywords used by the binning function
        """
        self._check_dataname(dataname)
        self._check_datatype(datatype)
        bin_dict = [('method', binning_name),
                    ('datatype', datatype),
                    ('args', args), 
                    ('start', start), 
                    ('stop', stop),
                    ('kwargs', kwargs)]
        if self._lookup_dict[dataname]['time_binning'] is None:
            self._lookup_dict[dataname]['time_binning'] = [OrderedDict(bin_dict)]
        else:
            self._lookup_dict[dataname]['time_binning'].append(OrderedDict(bin_dict))
    
    def add_energy_binning(self, dataname, binning_function, *args, start=None, 
                           stop=None, **kwargs):
        """Add a new energy binning function for the data file
        
        Parameters:
        --------------
        dataname: str
            The data filename
        binning_function: function
            The binning function
        *args: 
            Additional arguments used by the binning function
        start: float, optional
            The start of the data range to be rebinned. The default is to start at the 
            beginning of the RateHisto segment.
        stop: float, optional
            The end of the data range to be rebinned. The default is to stop at
            the end of the RateHisto segment.
        **kwargs:
            Additional keywords used by the binning function
        """
        self._check_dataname(dataname)
        bin_dict = [('method', binning_function),
                    ('args', args), 
                    ('start', start), 
                    ('stop', stop),
                    ('kwargs', kwargs)]
        if self._lookup_dict[dataname]['energy_binning'] is None:
            self._lookup_dict[dataname]['energy_binning'] = [OrderedDict(bin_dict)]
        else:
            self._lookup_dict[dataname]['energy_binning'].append(OrderedDict(bin_dict))
            
    def add_source_selection(self, dataname, source_intervals):
        """Add source selection(s) for the data file
        
        Parameters:
        --------------
        dataname: str
            The data filename
        source_intervals: list
            A list of source selection intervals, each item of the list being a tuple 
            of the format (low, high)
        """
        self._check_dataname(dataname)
        source_intervals = self._assert_selections(source_intervals)
        self._lookup_dict[dataname]['source_selection'] = source_intervals
    
    def add_energy_selection(self, dataname, energy_intervals):
        """Add energy selection(s) for the data file
        
        Parameters:
        --------------
        dataname: str
            The data filename
        energy_intervals: list
            A list of energy selection intervals, each item of the list being a tuple 
            of the format (low, high)
        """
        self._check_dataname(dataname)
        energy_intervals = self._assert_selections(energy_intervals)
        self._lookup_dict[dataname]['energy_selection'] = energy_intervals
    
    def add_background_selection(self, dataname, background_intervals):
        """Add background selection(s) for the data file
        
        Parameters:
        --------------
        dataname: str
            The data filename
        background_intervals: list
            A list of background selection intervals, each item of the list being a tuple 
            of the format (low, high)
        """
        self._check_dataname(dataname)
        self._lookup_dict[dataname]['background_selection'] = background_intervals
    
    def add_time_display_view(self, dataname, display_range):
        """Add the display range of the lightcurve for the data file
        
        Parameters:
        --------------
        dataname: str
            The data filename
        display_range: list
            The values of the lightcurve display window in the format
            [xmin, xmax, ymin, ymax]
        """
        self._check_dataname(dataname)
        if isinstance(display_range, np.ndarray):
            display_range = display_range.tolist()
        self._lookup_dict[dataname]['time_display_view'] = display_range

    def add_energy_display_view(self, dataname, display_range):
        """Add the display range of the count spectrum for the data file
        
        Parameters:
        --------------
        dataname: str
            The data filename
        display_range: list
            The values of the count spectrum display window in the format
            [xmin, xmax, ymin, ymax]
        """
        self._check_dataname(dataname)
        if isinstance(display_range, np.ndarray):
            display_range = display_range.tolist()
        self._lookup_dict[dataname]['energy_display_view'] = display_range
    
    def read_lookup(self, filename, filetype='json'):
        """Read a lookup file
        
        Parameters:
        --------------
        filename: str
            The file name
        filetype: str, optional
            The filetype, default is json
        """
        ext = filename.split('.')[-1]
        if filetype not in self._export_options:
            raise ValueError("File format '{0}' is not available for lookup filetype".format(filetype)) 
        
        if not os.path.isfile(filename):
            raise IOError("Gspec lookup file does not exist")
            
        if filetype == 'npy':
            self._lookup_dict = np.load(filename)[0]
        elif filetype == 'json':
            with open(filename, 'r') as f:
                json_str = f.read()
            thedict = json.loads(json_str, cls=NumpyArrayDecoder)
            self._lookup_dict = self._deserialize_lookup(thedict)
        else:
            pass        
    
    def write_lookup(self, filename, export='json'):
        """Write a lookup file
        
        Parameters:
        --------------
        filename: str
            The filename to write out 
        export: str, optional
            The export filetype, default is json
        """
        if export not in self._export_options:
            raise ValueError("File format '{0}' is not available for lookup export".format(export)) 
        
        if export == 'npy':
            np.save(filename, np.array([self._lookup_dict]))
        elif export == 'json':
            with open(filename, 'w') as f:
                json.dump(self._lookup_dict, f, indent=4, separators=(',', ': '), 
                          cls=NumpyArrayEncoder)
        else:
            pass
    
    def _check_dataname(self, dataname):
        """Check to see if the data file has been added to the lookup 
        
        Parameters:
        --------------
        dataname: str
            The data file name
        """
        if dataname not in self._lookup_dict.keys():
            raise KeyError('File {0} not currently tracked. Add this file to '
                           'the lookup and try again.'.format(dataname))
    
    def split_off_dataname(self, dataname):
        """Return a new lookup object containing only the requested data file
        
        Parameters:
        --------------
        dataname: str
            The requested data filename
        
        Returns:
        -----------
        new_lookup: GspecLookup
            The new lookup object
        """
        self._check_dataname(dataname)
        new_lookup = GspecLookup()
        new_lookup._lookup_dict[dataname] = self._lookup_dict[dataname]
        return new_lookup
        
    def _assert_selections(self, selections):
        """Check to ensure the selections are of the correct form.

        Parameters:
        -----------
        selections: tuple or list of tuples
            The selection(s) to check
        
        Returns:
        --------
        selections: list
        """
        if (all(isinstance(selection, list) for selection in selections)) | \
           (all(isinstance(selection, tuple) for selection in selections)):
            if any(len(selection) != 2 for selection in selections):
                raise ValueError('Each range in selections must be of the '
                                 'form (lo, hi)')
            else:
                return selections
        else:
            if len(selections) != 2:
                raise ValueError('Selections must either be a range of '
                                   'the form (lo, hi) or a list of ranges')
            else:
                return [selections]
    
    def merge_lookup(self, lookup, overwrite=False):
        """Merge an existing lookup into this lookup
        
        Parameters:
        --------------
        lookup: GspecLookup
            The lookup object to be merged into this lookup
        overwrite: bool, optional
            If set to True, then any datanames in the current lookup will be overwritten
            if those same datanames are in the input lookup.  Default is False
        """
        # get datanames of the input lookup
        datanames = lookup._lookup_dict.keys()
        for dataname in datanames:
            # if dataname is already in this lookup and we don't want to overwrite
            if (dataname in self._lookup_dict) & (not overwrite):
                continue
            self._lookup_dict[dataname] = lookup._lookup_dict[dataname]
    
    def display_lookup(self):
        """Pretty print a lookup for display (in json format)
        """
        lu = json.dumps(self._lookup_dict, indent=4, separators=(',', ': '), 
                        cls=NumpyArrayEncoder)
        return lu
            
    def _deserialize_lookup(self, thedict):
        """Deserializes the lookup dictionary so that it can be used.
        Specifically, this converts the background and binning classes/functions
        described by the tuple of their module and name strings to actual
        compiled classes/functions.
        """
        for d in thedict.keys():
            if thedict[d]['source_selection'] is not None:
                thedict[d]['source_selection'] = thedict[d]['source_selection'].tolist()
            if thedict[d]['energy_selection'] is not None:
                thedict[d]['energy_selection'] = thedict[d]['energy_selection'].tolist()
            if thedict[d]['background_selection'] is not None:
                thedict[d]['background_selection'] = thedict[d]['background_selection'].tolist()
            if thedict[d]['time_binning'] is not None:
                thedict[d]['time_binning'] = thedict[d]['time_binning'].tolist()
            if thedict[d]['energy_binning'] is not None:
                thedict[d]['energy_binning'] = thedict[d]['energy_binning'].tolist()
            
        return thedict
    
    def _check_datatype(self, datatype):
        """Check validity of datatype string
        """
        if (datatype != 'binned') & (datatype != 'unbinned'):
            raise ValueError("datatype must be either 'binned' or 'unbinned'")

        
class RmfitLookup(object):
    """Container class for an RMfit lookup file
    
    The lookup file contains the time selection, energy selection, background selection, 
    rebinned time edges, rebinned energy edges, window display ranges, and polynomial
    background fit order.

    Attributes:
    -----------
    is_energy_rebinned: bool
        True if the lookup contains rebinned energy edges
    is_time_rebinned: bool
        True if the lookup contains rebinned time edges
    has_energy_selection: bool
        True if the lookup contains energy selection(s)
    has_time_selection: bool
        True if the lookup contains time selection(s)
    has_background_selection: bool
        True if the lookup contains background selection(s)
    background_poly_order: int
        The background polynomial fit order
    time_display_range: np.array
        The time window display range in the order of (xmin, xmax, ymin, ymax)
    energy_display_range: np.array
        The energy window display range in the order of (xmin, xmax, ymin, ymax)
    
    Public Methods:
    -----------
    rebinned_energy_edges:
        Return the rebinned energy edges, if available
    energy_selection:
        Return the energy selection(s), if available
    rebinned_time_edges:
        Return the rebinned time edges, if available
    time_selection:
        Return the time selection(s), if available
    background_selection:
        Return the background selection(s), if available
    """ 
    def __init__(self, filename, ti_file=None):
        """
        Parameters:
        --------------
        filename: str
            The filename of the lookup file to load
        ti_file: str, optional
            The time index file for TTE files
        """
        # set default attributes
        self.is_energy_rebinned = False
        self.is_time_rebinned = False
        self.has_energy_selection = False
        self.has_time_selection = False
        self.has_background_selection = False
        self.background_poly_order = None
        self.time_display_range = None
        self.energy_display_range = None
        
        self._rebinned_energy_edges = None
        self._energy_selection = None
        self._rebinned_time_edges = None
        self._time_selection = None
        self._background_selection = None
        self._tte_edges = None
        
        # read the lookup file
        self._read_file(filename)
        
        # read the ti file (if given)
        if ti_file is not None:
            self._read_ti(ti_file)
    
    def _read_file(self, filename):
        """
        Read and parse the lookup file
        
        Parameters:
        --------------
        filename: str
            The filename of the lookup file to load
        """
        if not os.path.isfile(filename):
            raise IOError("RMFit lookup file does not exist")
        
        try:
            with open(filename, 'r') as f:
                txt = list(f)
        except:
            raise IOError('Incorrect format for an RMFit lookup file')
        txt = [t.strip() for t in txt]
        txt = ' '.join(txt).split()
        try:
            self._parse_file(txt)
        except:
            raise IOError('Incorrect format for an RMFit lookup file')
    
    def _parse_file(self, txt):
        """
        Parse the lookup file array
        
        Parameters:
        --------------
        txt: list
            The contents of the lookup file, converted to a list
        """
        
        # first number is number of rebinned energy edges.  Even if it's zero, we still 
        # have to read the next number (which should also be zero).
        num_energy_rebinned = int(txt[0])
        if num_energy_rebinned == 0:
            num_energy_rebinned = 1
        else:
            self.is_energy_rebinned = True
        idx1 = num_energy_rebinned+1
        self._rebinned_energy_edges = np.array(txt[1:idx1], dtype=int)
        
        # energy selections
        num_energy_selections = int(txt[idx1])
        if num_energy_selections == 0:
            num_energy_selections = 1
        else:
            self.has_energy_selection = True
        idx1 += 1
        idx2 = idx1+num_energy_selections
        self._energy_selection = np.array(txt[idx1:idx2], dtype=float)
        
        # rebinned time edges
        idx1 = idx2
        num_time_rebinned = int(txt[idx1])
        if num_time_rebinned == 0:
            num_time_rebinned = 1
        else:
            self.is_time_rebinned = True
        idx1 += 1
        idx2 = idx1+num_time_rebinned
        self._rebinned_time_edges = np.array(txt[idx1:idx2], dtype=int)
        
        # time selections
        idx1 = idx2
        num_time_selections = int(txt[idx1])
        if num_time_selections == 0:
            num_time_selections = 1
        else:
            self.has_time_selection = True
        idx1 += 1
        idx2 = idx1+num_time_selections
        self._time_selection = np.array(txt[idx1:idx2], dtype=float)
        
        # background selections
        idx1 = idx2
        num_bkgd_selections = int(txt[idx1])
        if num_bkgd_selections == 0:
            num_bkgd_selections = 1
        else:
            self.has_background_selection = True
        idx1 += 1
        idx2 = idx1+num_bkgd_selections
        self._background_selection = np.array(txt[idx1:idx2], dtype=float)
        
        # time and energy window ranges: (xmin, xmax, ymin, ymax)
        idx1 = idx2+4
        idx2 = idx1+4
        self.time_display_range = np.array(txt[idx1:idx2], dtype=float) 
        idx1 = idx2
        idx2 += 4
        self.energy_display_range = np.array(txt[idx1:idx2], dtype=float)
        # polynomial background order
        self.background_poly_order = int(txt[idx2])
    
    def _read_ti(self, ti_file):
        """
        Read the time index file for TTE binning
        
        Parameters:
        --------------
        ti_file: str
            The filename of the ti file to load
        """
        if not os.path.isfile(ti_file):
            raise IOError("Rmfit ti file does not exist")
        with open(ti_file, 'r') as f:
            txt = list(f)
        txt = txt[1:]
        self._tte_edges = np.array([t.strip() for t in txt], dtype=float)
        self.is_time_rebinned = True
    
    def rebinned_energy_edges(self):
        """Return the rebinned energy edges, if available
        
        Returns:
        -----------
        np.array
            The indices of the rebinned energy edges
        """
        if not self.is_energy_rebinned:
            warnings.warn("No rebinned energy edges in file")
            return None
        
        return self._rebinned_energy_edges

    def rebinned_time_edges(self):
        """Return the rebinned time edges, if available
        
        Returns:
        -----------
        np.array
            The indices of the rebinned time edges
        """
        if not self.is_time_rebinned:
            warnings.warn("No rebinned time edges in file")
            return None
        
        return self._rebinned_time_edges
    
    def energy_selection(self):
        """Return the energy selections, if available
        
        Returns:
        -----------
        selections: tuple
            Tuple of arrays that conain the range of each selection
        """
        if not self.has_energy_selection:
            warnings.warn("No energy selection in file")
            return None
        
        arr = self._energy_selection.reshape(2, -1).T
        # we don't care about the convex hull of the selections
        arr = arr[1:]
        selections = tuple([selection for selection in arr])
        return selections

    def time_selection(self):
        """Return the time selections, if available
        
        Returns:
        -----------
        selections: tuple
            Tuple of arrays that conain the range of each selection
        """
        if not self.has_time_selection:
            warnings.warn("No time selection in file")
            return None
        
        arr = self._time_selection.reshape(2, -1).T
        # we don't care about the convex hull of the selections
        arr = arr[1:]
        selections = tuple([selection for selection in arr])
        return selections
        
    def background_selection(self):
        """Return the background selections, if available
        
        Returns:
        -----------
        selections: tuple
            Tuple of arrays that conain the range of each selection
        """
        if not self.has_background_selection:
            warnings.warn("No background selection in file")
            return None
        
        arr = self._background_selection.reshape(2, -1).T
        selections = tuple([selection for selection in arr])
        return selections
    
    def tte_edges(self):
        """Return the TTE time edges, if available
        
        Returns:
        -----------
        tte_edges: np.array
            The time bin edges for a TTE file
        """
        if self._tte_edges is None:
            warnings.warn("No TTE edges found.  Need '.ti' file")
            return None
        
        return self._tte_edges


class NumpyArrayEncoder(json.JSONEncoder):
    """Custom JSON encoder for numpy arrays. Converts them to a list.
    """
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

class NumpyArrayDecoder(json.JSONDecoder):
    """Custom JSON decoder to turn JSON lists into numpy arrays
    """
    def __init__(self, *args, **kwargs):
        json.JSONDecoder.__init__(self, object_hook=self.object_hook, 
                                  *args, **kwargs)
    def object_hook(self, obj):
        # if object is a dictionary
        if type(obj) == dict:
            for key in obj.keys():
                # and if the value is a list, change to numpy array
                if type(obj[key]) == list:
                    obj[key] = np.array(obj[key], dtype=type(obj[key]))
            
        return obj
