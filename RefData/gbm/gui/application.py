import os
import shutil
from tkinter import N, E, W, S, BROWSE, SUNKEN, END, ACTIVE
from tkinter import Tk, Frame, Listbox, Button, filedialog, messagebox

import astropy.io.fits as fits
from numpy import ceil, sqrt, vstack

import gbm
from gbm.data.scat import SCat, ModelFit
from gbm.gspec import GspecManager, LookupFile
from gbm.xspec.model_manager import XspecModelManager
from gbm.xspec.xspec import XSPEC
from .dataviewer import DataViewer
from .dialogs import Splash, FitOptionsDialog
from .fitplotter import FitPlotter
from .globals import BGCOLOR, FRAMERELIEF, GUIFONTCOLOR
from .logger import Logger
from .menubar import MenuBar
from .tooltip import ListboxToolTip
from .window import Window


class Application(Window):
    """Primary GSPEC application class"""

    def __init__(self):

        super(Application, self).__init__()
        self.name = 'GSpec v{:}'.format(gbm.__version__)
        self.gspec = GspecManager()
        self._open_windows = {}
        self.filenames = []

        self._export_files = None
        self._last_fit_batch = None
        self._current_fit_batch = None
        self.xspec = None
        self.logger = None
        self.scat = None
        self.models = XspecModelManager()
        self.xspec_version = None
        self.temp_dirs = []

        # valid files for GSpec
        self.datafiletypes = (("PHA Files", "*.pha"), ("FITS Files", "*.fit"))
        self.lookupfiletypes = (("GSpec Lookup", "*.json"),)
        self.combinedfiletypes = (("PHA Files", "*.pha"), ("FITS Files", "*.fit"), ("GSpec Lookup", "*.json"))

        # Create the tkinter interface
        self.build_gui()

    def build_gui(self):

        # Create the root window
        root = Tk()
        root.title(self.name)
        self.width = self._geometry.width_fraction_to_npix(0.2)
        self.height = self._geometry.height_fraction_to_npix(0.25)
        self.aspectratio = self.width / self.height

        # Set a soft maximum window size
        if self.width > 500:
            self.width = 500
        elif self.width < 450:
            self.width = 450
        if self.height > 236:
            self.height = 236

        root.geometry("{0}x{1}+0+0".format(self.width, self.height))
        root.rowconfigure(0, weight=1)
        root.columnconfigure(0, weight=1)

        # splash = Splash(root)

        # Set a minimum window size
        root.minsize(400, 200)

        # Make a reference to the root window so that we can kill it later
        self.root = root

        # Override the frame inilization so as to pass the root window
        Frame.__init__(self, root)
        self.config(bg=BGCOLOR)

        # Tell the frame to expand and fill the root window
        self.grid(sticky=N + W + E + S)

        # self.about = AboutGspec(parent=self.root)

        # Create the menu bar
        self.menubar = MenuBar(root, self.name)
        items = ["About...", "Quit"]
        commands = [lambda: Splash(root), self.on_window_close]
        self.gspecmenu = self.menubar.register_menu('GSpec', items, commands)

        # file menu
        items = ["Load Data...", "Load From Lookup..."]
        commands = [self.load_file, self.load_lookup]
        self.filemenu = self.menubar.register_menu('File', items, commands)

        # task menu
        items = ["Display", "Hide", "Delete"]
        commands = [self.display, self.hide, self.delete_file]
        self.optionmenu = self.menubar.register_menu('Options', items, commands)

        # export menu
        items = ["Selected Data to XSPEC PHAI Format...", "Selected Data to XSPEC PHAII Format...", "Fit Results to FITS File...", "Fit Log to Text File..."]
        commands = [self.prepare_export_to_xspec_pha1, self.prepare_export_to_xspec_pha2, self.export_fit_results, self.export_fit_log]
        self.exportmenu = self.menubar.register_menu('Export', items, commands)
        # self.exportmenu.entryconfig(0, state='disable')
        # self.exportmenu.entryconfig(1, state='disable')
        # self.exportmenu.entryconfig(2, state='disable')

        # windows menu
        items = ["Cascade Windows", "Tile Windows"]
        commands = [self.cascade_windows, self.tile_windows]
        self.windowmenu = self.menubar.register_menu('Windows', items, commands,
                                                     postcommand=self.refresh_window_list)
        # self.windowmenu.add_separator()
        self.menubar.add_separator('Windows')

        # Add the menu bar to the root window
        root.config(menu=self.menubar.top)

        # Create the listbox
        self.listbox = Listbox(self, selectmode=BROWSE,
                               borderwidth=2,
                               selectbackground=BGCOLOR,
                               selectborderwidth=0,
                               highlightthickness=0,
                               fg=GUIFONTCOLOR,
                               relief=SUNKEN)

        # Setup the listbox
        self.tooltip = ListboxToolTip(self.listbox, self.filenames)
        self.listbox.bind("<Leave>", self.reset_tooltip)

        # Create the buttons
        frame_center = Frame(self, bg=BGCOLOR, relief=FRAMERELIEF)
        loadbutton = Button(frame_center, text=u"Load", relief=FRAMERELIEF, fg=GUIFONTCOLOR,
                            highlightbackground=BGCOLOR, command=self.load_file, width=7)
        displaybutton = Button(frame_center, text=u"Display", relief=FRAMERELIEF, fg=GUIFONTCOLOR,
                               highlightbackground=BGCOLOR, command=self.display, width=7)
        hidebutton = Button(frame_center, text=u"Hide", relief=FRAMERELIEF, fg=GUIFONTCOLOR,
                            highlightbackground=BGCOLOR, command=self.hide, width=7)
        deletebutton = Button(frame_center, text=u"Delete", relief=FRAMERELIEF, fg=GUIFONTCOLOR,
                              highlightbackground=BGCOLOR, command=self.delete_file, width=7)

        # Place up the list box
        self.listbox.grid(column=0, row=0, padx=10, pady=(10, 0), sticky=N + E + W + S)

        # Pack all the buttons together
        loadbutton.grid(row=0, column=0, sticky=E + W, pady=10)
        displaybutton.grid(row=0, column=1)
        hidebutton.grid(row=0, column=2)
        deletebutton.grid(row=0, column=3)

        # Place the frame within the root window
        frame_center.grid(column=0, row=1, sticky=N + S)

        # Configure the weights of the columns
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=0)
        self.columnconfigure(0, weight=1)

        # Bind the onWindowClose function to the window close action
        root.protocol("WM_DELETE_WINDOW", self.on_window_close)

        def scrollwheel(event):
            print('hello')

        root.bind('<MouseWheel>', scrollwheel)

        # root.protocol("WM_DELETE_WINDOW", self.update_window_status)

        self.register_window(self.name, self)

        # Warn the user if xspec could not be found in their path
        # or it is older than v12.10
        root.update()
        self.test_xspec()
        root.update()

        self.window_is_open = True

        root.mainloop()

        # Resorting to using update loop rather than mainloop to avoid a known issue with tk
        # https://bugs.python.org/issue10731
        # while self.window_is_open == True:
            # Tk.update(self)

        return
        

    def test_xspec(self):
        # load XSPEC and get the version string
        try:
            test_xspec = XSPEC()
            v = test_xspec.do_command('version')[-1]
            v = v.split(':')[-1].strip()
            major, minor, build = v.split('.')
            self.xspec_version = v
            if (int(major) < 12) or (int(major) == 12 and int(minor) < 10):
                msg = 'You are using XSPEC version {0}.\n'.format(v)
                msg += 'Some functionality will not be available for versions < 12.10. Please consult the GSPEC documentation for more information.'
                # messagebox.showerror('Old XSPEC Version', msg, parent=self.root)
                messagebox.showerror('Warning', msg)
                # remove models that depend on v12.10
                self.models.remove_model('Band (Epeak)')
                self.models.remove_model('Comptonized')

        except:
            messagebox.showerror('XSPEC not found', 'GSpec cannot find XSPEC.\n'
                                                    'Make sure it is in your path.', parent=self.root)

    def reset_tooltip(self, event=None):
        self.tooltip.index = -1
        self.tooltip.deactivate()

    def register_window(self, name, window):

        # # Enable fit log exporting
        # if 'Fit Log' in name:
        #     print('enabling fit log exporting')
        # self.exportmenu.entryconfig(2, state='normal')

        # Enable fit result exporting
        # if 'Spectral Fit Display' in name:
        #     print('enabling fit result exporting')
        # self.exportmenu.entryconfig(1, state='normal')

        self._open_windows[name] = window

    def unregister_window(self, name):

        # Disable fit log exporting
        # if 'Fit Log' in name:
        #     print('disabling fit log exporting')
        # self.exportmenu.entryconfig(2, state='disable')

        # Disable fit result exporting
        # if 'Spectral Fit Display' in name:
        #     print('disabling fit result exporting')
        # self.exportmenu.entryconfig(1, state='disable')

        if name in self._open_windows.keys():
            del self._open_windows[name]

    def load_file(self):
        selected_filename = filedialog.askopenfilenames(filetypes=self.datafiletypes,
                                                        parent=self.root)
        if selected_filename == '':
            return

        num_current = len(self.filenames)
        selected_filename = self.root.tk.splitlist(selected_filename)
        for filename in selected_filename:
            # if the file is already open
            if filename in self.filenames:
                index = self.filenames.index(filename)
                self.listbox.activate(index)
                self.display()
                continue

            self.filenames.append(filename)
            self.listbox.insert(END, os.path.basename(filename))
            self.gspec.add_data(filename)

        self.display(index_range=(num_current, END))

    def load_lookup(self):
        # only allow gspec lookups right now
        selected_filename = filedialog.askopenfilename(filetypes=self.lookupfiletypes,
                                                       parent=self.root)
        if selected_filename == '':
            return

        # load lookup, register with Gspec, and launch window for each detector
        # in lookup
        num_current = len(self.filenames)
        dir = os.path.dirname(selected_filename)
        lu = LookupFile.read_from(selected_filename)
        for dataname in lu.files():
            datafilename = os.path.join(dir, dataname)
            try:
                # self.gspec.add_data(datafilename, lookup=lu)
                self.gspec.add_data(datafilename)
                self.gspec.load_gspec_lookup(dataname, lu)

            except:
                pass
            self.filenames.append(datafilename)
            self.listbox.insert(END, dataname)
        self.display(index_range=(num_current, END))

    def display(self, index_range=None):

        # Get the selected filename
        if index_range is None:
            selected_filename = self.listbox.get(ACTIVE)
        else:
            selected_filename = self.listbox.get(index_range[0], index_range[1])
        if selected_filename == '':
            return

        # Remove any whitespaces that the listbox may have added to the string
        selected_filename = self.root.tk.splitlist(selected_filename)
        for filename in selected_filename:
            dataname = os.path.basename(filename)
            if dataname in self._open_windows.keys():
                print("Window already open")
                self._open_windows[dataname].window_focus()
                continue

            # Display the data
            viewer = DataViewer(filename, gspec_root=self)

            # register with the window menu
            self.register_window(dataname, viewer)

            # Enable data exporting to xspec
            # print('enabling data exporting')
            # self.exportmenu.entryconfig(0, state='normal')

    def hide(self):
        # Get the selected filename
        selected_filename = self.listbox.get(ACTIVE)
        if selected_filename == '':
            return

        # Remove any whitespaces that the listBox may have added to the string
        selected_filename = "".join(selected_filename.split())
        selected_filename = os.path.basename(selected_filename)

        # Hide the data
        try:
            window = self._open_windows[selected_filename]
            self.unregister_window(window.name)
            window.on_window_close()
        except:
            print('Window is already hidden')

    def delete_file(self):
        selected_filename = self.listbox.get(ACTIVE)
        if selected_filename == '':
            return
        selected_filename = [f for f in self.filenames if selected_filename in f][0]
        self.hide()
        self.filenames.remove(selected_filename)
        self.listbox.delete(ACTIVE)
        self.gspec.remove_data(os.path.basename(selected_filename))

    def cascade_windows(self):
        # calculate number of columns and rows
        num_windows = len(self._open_windows)
        if num_windows == 0:
            return

        iwindow = 0
        for _, window in self._open_windows.items():
            xpos = int(iwindow * 45.0)
            ypos = int(iwindow * 45.0)
            window.root.geometry("{0}x{1}+{2}+{3}".format(window.width,
                                                          window.height,
                                                          xpos, ypos))
            window.root.lift()
            iwindow += 1

    def tile_windows(self):
        # calculate number of columns and rows
        num_windows = len(self._open_windows) - 1
        if num_windows <= 0:
            return

        # window grid and sizing
        num_cols = ceil(sqrt(num_windows))
        num_rows = ceil(num_windows / num_cols)
        window_width = self._geometry.ncolumns_to_npix(num_cols)
        window_height = self._geometry.nrows_to_npix(num_rows) - 45

        # go through each open window, move and resize
        icolumn = 0
        irow = 0
        for window_name, window in self._open_windows.items():
            # ignore root window
            if window_name == self.name:
                continue
            if icolumn > num_cols - 1:
                icolumn = 0
                irow += 1
            xpos = icolumn * window_width
            ypos = irow * (window_height + 45)
            window.root.geometry("{0}x{1}+{2}+{3}".format(window_width,
                                                          window_height,
                                                          xpos, ypos))
            icolumn += 1

    def on_window_close(self):
        if self.gspec.lookup_dirty:
            title = 'Save Lookup?'
            msg = 'Settings have changed for this data set.\n'
            msg += 'Would you like to save a lookup file?'
            a = messagebox.askquestion(title, msg, parent=self.root)
            if a == 'yes':
                lufile = self.write_lookup(self.root)

                # if save is canceled, also cancel the exit
                if lufile == '':
                    return

        # if messagebox.askokcancel("Quit", "Are you sure you want to quit?"):

        # pop the application window because we are going to directly
        # kill it after all other windows are killed
        self._open_windows.pop(self.name)

        openwindows = list(self._open_windows.values())
        num_windows = len(openwindows)
        for i in range(num_windows):
            # try:
            openwindows[i].on_window_close()
            # except:
            #    pass

        # clean up temporary directories
        if len(self.temp_dirs) > 0:
            for dir in self.temp_dirs:
                shutil.rmtree(dir)

        self.window_is_open = False

        # bye bye!
        self.root.destroy()

    def write_lookup(self, parent):
        datanames = list(self.gspec.data.keys())
        datafile = gbm.file.GbmFile.from_path(self.gspec.data[datanames[0]]['data'].filename)
        directory = datafile.directory
        uid = datafile.uid
        data_type = datafile.data_type
        trigger = datafile.trigger
        meta = '_lookup'
        extension = 'json'
        lufile = gbm.file.GbmFile.create(uid=uid, data_type=data_type,
                                         trigger=trigger, meta=meta,
                                         extension=extension)
        initialfile = lufile.basename()
        filename = filedialog.asksaveasfilename(initialdir=directory,
                                                initialfile=initialfile,
                                                title="Save Lookup File",
                                                filetypes=(('JSON file', '*.json'),),
                                                defaultextension='.json',
                                                parent=parent)
        if filename != '':
            filename = self.gspec.save_lookup(filename=filename)
        return filename

    def refresh_window_list(self):

        numwin = len(self._open_windows)
        names = list(self._open_windows.keys())
        windows = list(self._open_windows.values())
        for i in range(numwin):
            windows[i].menubar.remove_menu_item_range("Windows", 2, numwin + 3)
            for j in range(numwin):
                windows[i].menubar.add_menu_item('Windows', names[j], windows[j].window_focus)

    def prepare_spectral_fit(self, batch=False):

        # Specify if the current fit is a batch fit
        self._current_fit_batch = batch

        # Check to make sure all of the needed selections have been made
        datanames = self.gspec.data.keys()
        directory = None
        for dataname in datanames:
            try:
                parent = self._open_windows[dataname]
            except:
                parent = self

            # Create an error message if no background fit exists
            if self.gspec.data[dataname]['bkgdmodel'] is None:
                message = 'No background model exists for the file:\n\n ' \
                          '%s\n\nA background fit must exist to perform' \
                          ' a spectral fit.' % os.path.basename(self.gspec.data[dataname]['data'].filename)
                response = messagebox.showerror('Error', message, parent=parent.root)
                return

            # Create an error message if no response file exists
            if self.gspec.data[dataname]['response'] is None:

                title = 'Warning'
                message = 'No response file has been loaded for the file:\n\n  ' \
                          '%s\n\nWould you like to select a response file now?' \
                          % os.path.basename(self.gspec.data[dataname]['data'].filename)
                response = messagebox.askyesno(title, message, parent=parent.root)

                if response == True:
                    datafilename = self.gspec.data[dataname]['data'].filename.split('/')[-1]
                    selected_filename = parent.load_response()
                    if selected_filename != '':
                        self.gspec.lookup[dataname].add_response(selected_filename)
                    else:
                        return
                else:
                    return

        # Check to see if the user changed from single to batch fitting (or vice versa)        
        if self._last_fit_batch is not None:
            if self._last_fit_batch != batch:
                self.gspec.export_dirty = True

        # Something has changed (selection, fitting type, etc), so write out new data
        if self.gspec.export_dirty:
            directory = None
            phafilenames = []
            for dataname in datanames:

                # Create the pha and bak files required by xspec
                phafilename = self.gspec.export_to_xspec(dataname,
                                                         directory=directory,
                                                         MULTI=batch)

                # write multiple detector files to one directory
                if directory is None:
                    directory = os.path.dirname(phafilename)
                phafilenames.append(phafilename)

            # Double check that the number of specta match in each file before proceeding
            print("\nChecking for matching number of spectra...")

            # Get the number of spectra in each file
            num_specs = []
            for phafilename in phafilenames:
                with fits.open(phafilename) as hdulist:
                    num_spec = hdulist['SPECTRUM'].header['NAXIS2']
                    num_specs.append(num_spec)

            # Check to see if all the files have the same number of spectra
            if all(x == num_specs[0] for x in num_specs) == False:
                print("*** Error ***\nNon-matching number of spectra. Fit canceled.")

                title = 'Batch Fit Error'
                message = 'Non-matching number of spectra found.\n\n' \
                          'Time resolved fits require the same number of time bin selections for each file.\n\n' \
                          'Please adjust your temporal selections and try again.'
                response = messagebox.showerror(title, message)

                return


            else:
                print("Done.\n")

            # print('Working files written to {0}'.format(directory))
            self.temp_dirs.append(directory)
            self.export_files = phafilenames
            self.gspec.export_dirty = False
            self._last_fit_batch = batch

        # # Now do the model selection, etc.
        # models = XspecModelManager()
        # dialog = FitOptionsDialog(self.root, models, self.perform_spectral_fit)

        # now do the model selection, etc.
        dialog = FitOptionsDialog(self.root, self.models, self.perform_spectral_fit)

    def perform_spectral_fit(self, model, statistic, weight, undetermined, batch=False):

        # print("Current batch fit: %s" % self._current_fit_batch)
        batch = self._current_fit_batch

        # print(model)
        # print(statistic)
        # print(weight)
        # print(undetermined)
        # print(batch)

        # print("\nFitting data...")

        # create the SCAT object and store the detector data
        with fits.open(self.export_files[0]) as hdulist:
            prihdr = hdulist['PRIMARY'].header
        self.scat = SCat(tstart=prihdr['TSTART'], tstop=prihdr['TSTOP'],
                         trigtime=prihdr['TSTOP'])

        for filename in self.filenames:
            dataname = os.path.basename(filename)
            self.scat.add_detector_data(self.gspec.get_detector_data(dataname))

        if 'Fit Log' in self._open_windows.keys():
            self.logger = self._open_windows['Fit Log']
        else:
            self.xspec = XSPEC(models=self.models)
            self.logger = Logger(self.xspec, gspec_root=self)

        # Time integrated fit
        if not batch:

            print("\nPerforming time integrated fit...")

            # # Get the number of spectra in the file
            # with fits.open(self.export_files[0]) as hdulist:
            #     spectrum = hdulist['SPECTRUM'].data
            #     counts = spectrum['COUNTS']
            #     time_start = spectrum['TIME']
            #     time_end = spectrum['ENDTIME']
            #     exposure = spectrum['EXPOSURE']

            for i in range(len(self.export_files)):
                self.xspec.load_pha(self.export_files[i], slot=1, extension=i + 1)

            self.xspec.set_energy_range(10.0, 1000.0)
            self.xspec.set_fit_statistic(statistic)
            self.xspec.set_fit_weight(weight)
            self.xspec.set_model(model)
            self.xspec.do_fit()

            num_groups = len(self.export_files)
            num_params = self.xspec.num_params
            data = self.xspec.extract_scat_data(num_groups, sum(num_params),
                                                self.scat.detectors[0].time_range)

            # store fit results
            model_results = ModelFit.from_xspec_data_dict(data)
            model_results.name = model
            self.scat.add_model_fit(model_results)

            # update the detector info with the counts model
            for i in range(num_groups):
                self.scat.detectors[i].photon_counts = data['counts'][i]
                self.scat.detectors[i].photon_errors = data['counts_error'][i]
                self.scat.detectors[i].photon_model = data['counts_model'][i]

            print("Done.")

        # Time resolved fits
        else:

            print("\nPerforming time resolved fit...\n")

            # Get the number of spectra in the file
            with fits.open(self.export_files[0]) as hdulist:
                spectrum = hdulist['SPECTRUM'].data
                counts = spectrum['COUNTS']
                time_start = spectrum['TIME']
                time_end = spectrum['ENDTIME']
                exposure = spectrum['EXPOSURE']

            num_extensions = len(counts)

            data = {}
            data['parameter_values'] = []
            data['parameter_sigmas'] = []

            # progressBar = ProgressBar(num_extensions, 0)

            # Load the pha files. Testing - Adam
            # self.xspec.load_phaii_files(self.export_files)
            # self.xspec.set_fit_statistic(statistic)
            # self.xspec.set_fit_weight(weight)
            # self.xspec.set_model(model, numgroups=num_extensions)
            # num_params = self.xspec.num_params
            # self.xspec.do_fit()

            # Select the fitting statistic, weight, and model
            self.xspec.set_fit_statistic(statistic)
            self.xspec.set_fit_weight(weight)
            self.xspec.set_model(model)
            num_params = self.xspec.num_params

            for extension in range(num_extensions):

                # The first extension should always be one, not zero
                extension = extension + 1

                # Load the pha file.
                for i in range(len(self.export_files)):
                    self.xspec.load_pha(self.export_files[i], extension=extension, slot=i + 1)

                print("Fitting spectrum %s" % extension)

                # Fit the data
                self.xspec.do_fit()

                # # Update the progress bar
                # progressBar.updateValue(extension)
                # progressBar.lift()

                # Extract the fit results for this extension
                num_groups = len(self.export_files)
                time_range = (time_start[extension - 1], time_end[extension - 1])
                data_single_fit = self.xspec.extract_scat_data(num_groups,
                                                               sum(num_params),
                                                               time_range)

                # store fit results
                model_results = ModelFit.from_xspec_data_dict(data_single_fit)
                model_results.name = model
                self.scat.add_model_fit(model_results)

                # update the detector info with the counts model
                # for i in range(num_groups):
                #    self.scat.detectors[i].photon_counts = data_single_fit['counts'][i]
                #    self.scat.detectors[i].photon_errors = data_single_fit['counts_error'][i]
                #    self.scat.detectors[i].photon_model = data_single_fit['counts_model'][i]

                if extension == 1:
                    data['time_start'] = time_start
                    data['time_end'] = time_end
                    data['exposure'] = exposure
                    data['parameter_names'] = data_single_fit['parameter_names']
                    data['parameter_values'] = data_single_fit['parameter_values']
                    data['parameter_sigmas'] = data_single_fit['parameter_sigmas']

                else:
                    data['parameter_values'] = vstack((data['parameter_values'], data_single_fit['parameter_values']))
                    data['parameter_sigmas'] = vstack((data['parameter_sigmas'], data_single_fit['parameter_sigmas']))

            # # Destroy the progress bar
            # progressBar.destroy()

            print("\nDone.")

        # Specify the model associated with this data
        data['model'] = model

        # Specify if the data is associated with a batch fit or single fit
        data['batch'] = batch

        if 'Spectral Fit Display' not in self._open_windows.keys():

            # fit_plotter = FitPlotter(self, data, gspec_root=self.gspec)
            fit_plotter = FitPlotter(data, gspec_root=self)

            # # Send the window the new fit results
            # fit_plotter.data = data

            # # Bring the window to the foreground
            # fit_plotter.root.lift()

            # # Replot the default spectrum
            # if fit_plotter.data['batch'] == False:
            #     fit_plotter.xlog = True
            #     fit_plotter.ylog = True
            #     fit_plotter.plotSpectrum()

        else:

            # Get the reference to existing fit plotter window
            fit_plotter = self._open_windows['Spectral Fit Display']

            # Send the window the new fit results
            fit_plotter.data = data

            # Bring the window to the foreground
            fit_plotter.root.lift()

            # Replot the default spectrum
            if fit_plotter.data['batch'] == False:
                fit_plotter.xlog = True
                fit_plotter.ylog = True
                fit_plotter.plotSpectrum()

            else:
                fit_plotter.plotBatchFitResults()


    def prepare_export_to_xspec_pha1(self, clobber=False):
        self.prepare_export_to_xspec(MULTI=False)

    def prepare_export_to_xspec_pha2(self, clobber=False):
        self.prepare_export_to_xspec(MULTI=True)

    def prepare_export_to_xspec(self, clobber=False, MULTI=False):

        # Check to make sure all of the needed selections have been made
        datanames = self.gspec.data.keys()

        if len(datanames) == 0:
            response = messagebox.showerror('Export Error', 'No data selections available for export')
            return

        for dataname in datanames:

            try:
                parent = self._open_windows[dataname]
            except:
                parent = self

            # Create an error message if no background fit exists
            if self.gspec.data[dataname]['bkgdmodel'] is None:
                message = 'No background model exists for the file:\n\n ' \
                          '%s\n\nA background fit must exist to perform' \
                          ' a spectral fit.' % os.path.basename(self.gspec.data[dataname]['data'].filename)
                response = messagebox.showerror('Error', message, parent=parent.root)
                return

            # Create an error message if no response file exists
            if self.gspec.data[dataname]['response'] is None:

                title = 'Warning'
                message = 'No response file has been loaded for the file:\n\n  ' \
                          '%s\n\nWould you like to select a response file now?' \
                          % os.path.basename(self.gspec.data[dataname]['data'].filename)
                response = messagebox.askyesno(title, message, parent=parent.root)

                if response == True:
                    datafilename = self.gspec.data[dataname]['data'].filename.split('/')[-1]
                    selected_filename = parent.load_response()
                    if selected_filename != '':
                        self.gspec.add_response(dataname, selected_filename)
                    else:
                        return
                else:
                    return

        # Ask the user to select a directory
        directory = filedialog.askdirectory()

        if directory == '':
            print('\nExport canceled.')
            return

        # Write out new data
        # directory=None
        # directory = "/Users/Kocevski/Desktop/GBM/"
        phafilenames = []

        print("\nExporting selected data to XSpec format.")

        for dataname in datanames:

            # Check to see if an existing file with the final xspec suffix already exists
            dataname_xspec = dataname.replace('_v', '_xspec_v')
            does_file_exist = os.path.isfile(directory + '/' + dataname_xspec)

            # Ask the user if they want to replace existing files.  This is asked only once and applied to all existing files
            if does_file_exist == True and clobber == False:

                title = "Warning"
                message = '“%s” already exists. Do you want to replace all existing files at this location?' % dataname_xspec
                response = messagebox.askyesno(title, message)

                if response == False:
                    print("\nExport canceled.")
                    return
                else:
                    clobber = True

            # Remove any existing files
            if does_file_exist == True and clobber == True:
                # print("\nRemoving %s" % (directory + '/' + dataname_xspec))
                os.remove(directory + '/' + dataname_xspec)
                # print("\nRemoving %s" % (directory + '/' + dataname_xspec.replace('.pha','.bak')))
                os.remove(directory + '/' + dataname_xspec.replace('.pha', '.bak'))
                # print("\nRemoving %s" % (directory + '/' + dataname_xspec.replace('.pha','.rsp')))
                os.remove(directory + '/' + dataname_xspec.replace('.pha', '.rsp'))

            # Create the pha and bak files required by xspec
            phafilename = self.gspec.export_to_xspec(dataname,
                                                     directory=directory,
                                                     MULTI=MULTI,
                                                     verbose=True)

            # write multiple detector files to one directory
            if directory is None:
                directory = os.path.dirname(phafilename)

            phafilenames.append(phafilename)

        print("Done.")

        # print('Working files written to {0}'.format(directory))
        # self.temp_dirs.append(directory)
        # self.export_files = phafilenames

    def export_fit_results(self):

        if 'Spectral Fit Display' in self._open_windows.keys():
            fit_plotter = self._open_windows['Spectral Fit Display']
            fit_plotter.commandSelection('Write Results to File')
        else:
            response = messagebox.showerror('Export Error', 'No current fit available for export')

        return

    def export_fit_log(self):

        if 'Fit Log' in self._open_windows.keys():
            logger = self._open_windows['Fit Log']
            logger.save_text()
        else:
            response = messagebox.showerror('Export Error', 'No current fit log available for export')

        return

    def command(self):
        pass
