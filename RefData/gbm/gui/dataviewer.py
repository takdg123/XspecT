import tkinter
from tkinter import Tk, Frame, Button, filedialog, StringVar, BooleanVar, OptionMenu, Checkbutton, Label
from tkinter import N, E, W, S, LEFT, RIGHT, TOP, BOTTOM, X, Y, BOTH, CENTER, ACTIVE
from tkinter import messagebox
import os
import copy

import matplotlib
from numpy import where

matplotlib.rcParams.update({'font.size': 10})
matplotlib.use("TkAgg")
from matplotlib import rcParams
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.patches import Rectangle
import warnings
from matplotlib.cbook import MatplotlibDeprecationWarning
warnings.filterwarnings("ignore", category=MatplotlibDeprecationWarning)

from .globals import BGCOLOR, FRAMERELIEF, OM_JUSTIFY, GUIFONTCOLOR, PLOTFONTSIZE
from .window import Window
from .dialogs import TextDisplayDialog, OptionDialog, PlotDialog, TextOptionDialog, ManualInputDialog
from gbm.plot.lightcurve import Lightcurve
from gbm.plot.spectrum import Spectrum
from gbm.plot.chisq import BackgroundChisq
from gbm import file

# Use the classic matplotlib style sheet
# if float(matplotlib.__version__[0]) >= 2.0:
#     import matplotlib.style
#     matplotlib.style.use('classic')

from gbm.lookup.lookup import LookupFile, View

class DataViewer(Window):
    """Class for displaying pha files"""

    def __init__(self, filename, gspec_root=None):
        
        super(DataViewer, self).__init__()
        self.filename = filename
        self.dataname = os.path.basename(filename)
        self.name = self.dataname
        
        # get reference to the root application window, register this window
        # with the application, and get the reference to the GspecManager
        if gspec_root is not None:
            self.gspec_root = gspec_root
            self.gspec_root.register_window(self.dataname, self)
            self.gspec = self.gspec_root.gspec

        # Extract the detector name for the data file
        self.detector = self.gspec.lookup[self.dataname].detector
        
        # Set the current binning
        if type(self.gspec.data[self.dataname]['data']).__name__ == 'Cspec':
            self.current_binning = 1.024
        else:
            self.current_binning = 0.064
        self.current_binning_function = None
        self.currentSNR = None

        self.lightcurve_display = True
        # Setting default state parameters
        self.source_selection_active = False
        self.snr_selection_active = False
        self.background_selection_active = False
        self.binning_selection_active = False
        #self.xZoomSelectionActive = False
        self.ylog_lightcurve = False
        self.xlog_spectrum = True
        self.ylog_spectrum = True

        self.lookupfiletypes = (("GSPEC Lookups", "*.json"), 
                                ("RMfit Lookups", "*.lu"))
        self.responsefiletypes = (("RSP", "*.rsp"), ("RSP2", "*.rsp2"))

        # Create the tkinter interface
        self.createTKinterInterface()

    ############################################

    def createTKinterInterface(self):

        # Create the root window
        # root = Tk()
        # root.title(self.name)
        # self.width = self._geometry.width_fraction_to_npix(0.67)
        # self.height = self._geometry.height_fraction_to_npix(0.50)
        # self.aspectratio = self.width/self.height
        # root.geometry("{0}x{1}+0+0".format(self.width, self.height))

        # Create the root window
        root = Tk()
        root.title(self.name)
        self.width = self._geometry.width_fraction_to_npix(0.67)
        self.height = self._geometry.height_fraction_to_npix(0.50)

        # Set a soft maximum window size of 1024x550
        if self.width > 1024: self.width = 1024
        if self.height > 550: self.height = 550

        self.aspectratio = self.width/self.height
        root.geometry("{0}x{1}+0+0".format(self.width, self.height))

        # Set a minimum window size
        root.minsize(600, 400)

        # Make a reference to the root window so that we can kill it later
        self.root = root

        # Initilize a frame instance with the root as its parent
        Frame.__init__(self, root)

        # Configure the frame
        self.config(bg=BGCOLOR)
        self.grid(sticky=N+W+E+S)

        # reference the menubar
        self.menubar = copy.copy(self.gspec_root.menubar)
        self.menubar.update_root(self.root, self.name)
        
        # update the list of available headers
        items = ['Header...', 'Lookup...', 'Response...', 'Screenshot', 'Dismiss']
        commands = [self.command, self.command, self.command, self.command, self.on_window_close]
        self.filemenu = self.menubar.replace_menu('File', items, commands)
        
        # add a header submenu
        headers = self.gspec.data[self.dataname]['data'].header
        items = list(headers.keys())
        commands = [lambda ext=ext: self.display_header(ext) for ext in headers.keys()]
        self.menubar.register_submenu('File', 'Header...', items, commands)
        
        # add a lookup submenu
        items = ["Save Lookup", "Read Lookup", "File Content", "Erase Current"]
        commands = [self.write_lookup, self.read_lookup, self.display_lookup, self.erase_lookup]
        self.menubar.register_submenu('File', 'Lookup...', items, commands)
        
        # add a response submenu
        items = ["Load Response", "Remove Response", "Display Response"]
        commands = [self.load_response, self.remove_response, self.command]
        respmenu = self.menubar.register_submenu('File', 'Response...', items, commands)
        respmenu.entryconfig(2, state='disable')
        
        items = ["Show Current Selection", "Colors", "Plot Configuration"]
        commands = [self.command, self.command, self.command]
        self.options_menu = self.menubar.replace_menu('Options', items, commands)
        self.options_menu.entryconfig(0, state='disable')
        self.options_menu.entryconfig(1, state='disable')
        self.options_menu.entryconfig(2, state='disable')
        
        #items = ["Refresh Plot"]
        #commands = [self.redraw]
        #self.window_menu = self.gspec_root.register_menu('Window', items, commands)

        root.config(menu=self.menubar.top)
        
        ########### Buttons ###########
        self.frame_controls = Frame(self, bg=BGCOLOR, relief=FRAMERELIEF, pady=10, padx=10, height=100)

        self.toggleButton = Button(self.frame_controls, text=u"Show Spectrum", highlightbackground=BGCOLOR,
                                   command=self.toggle, width=14, fg=GUIFONTCOLOR)
        self.fit_background_button = Button(self.frame_controls, text=u"Fit Background",
                                          highlightbackground=BGCOLOR, width=14, command=self.select_background, fg=GUIFONTCOLOR)

        # Rebin menu
        self.rebin_selection = StringVar(self.frame_controls)
        self.rebin_selection.set('Rebin Data')
        choices = ['Full Resolution', 'Temporal Resolution...', 'Signal to Noise...', 'Combine Source Intervals',
                   'Combine Into Single Bin', 'Combine by Factor']
        self.rebin_menu = OptionMenu(self.frame_controls, self.rebin_selection, *choices, command=self.rebin)
        self.rebin_menu.config(bg=BGCOLOR, width=14, justify=OM_JUSTIFY, fg=GUIFONTCOLOR, highlightcolor=GUIFONTCOLOR)

        # Source selection menu
        self.source_selection = StringVar(self.frame_controls)
        self.source_selection.set('Source Selection')
        choices = ['Source Interactive...', 'Source By S/N...', 'Clear Selections']
        self.source_selection_menu = OptionMenu(self.frame_controls, self.source_selection,
                                                *choices, command=self.select_source)
        self.source_selection_menu.config(bg=BGCOLOR, width=14, justify=OM_JUSTIFY, 
                                          fg=GUIFONTCOLOR)

        # Adjust source selection menu
        self.adjustSourceSelection = StringVar(self.frame_controls)
        self.adjustSourceSelection.set('Adjust Source')
        choices = ['< Shift Selection', '> Shift Selection', '< Left Selection', '> Left Selection',
                   '< Right Selection', '> Right Selection']
        self.adjustSourceSelectionMenu = OptionMenu(self.frame_controls, self.adjustSourceSelection, *choices,
                                                    command=self.adjust_source)
        self.adjustSourceSelectionMenu.config(bg=BGCOLOR, width=14, justify=OM_JUSTIFY, fg=GUIFONTCOLOR)

        # Export Settings menu
        self.exportSettingsSelection = StringVar(self.frame_controls)
        self.exportSettingsSelection.set('Push Selections')
        choices = ['Temporal Binning', 'Background', 'Source Selection', 'View Range', 'All of the Above']
        self.exportSettingsMenu = OptionMenu(self.frame_controls, self.exportSettingsSelection, *choices,
                                             command=self.export_settings)
        self.exportSettingsMenu.config(bg=BGCOLOR, width=14, justify=OM_JUSTIFY, fg=GUIFONTCOLOR)

        # Spectral fitting menu
        self.spectralFittingSelection = StringVar(self.frame_controls)
        self.spectralFittingSelection.set('Spectral Fitting')
        choices = ['Fit Selection...', 'Batch Fit Selections...']
        self.spectralFittingSelectionMenu = OptionMenu(self.frame_controls, self.spectralFittingSelection, *choices,
                                                       command=self.spectral_fit_selection)
        self.spectralFittingSelectionMenu.config(bg=BGCOLOR, width=14, justify=OM_JUSTIFY, fg=GUIFONTCOLOR)
        if self.gspec_root.xspec_version == None:
            self.spectralFittingSelectionMenu['state'] = 'disabled'

        # Pack up the buttons
        self.toggleButton.pack(padx=0, side=TOP, fill=X)
        self.rebin_menu.pack(padx=0, side=TOP, fill=X)
        self.fit_background_button.pack(padx=0, side=TOP, fill=X)
        self.source_selection_menu.pack(padx=0, side=TOP, fill=X)
        self.adjustSourceSelectionMenu.pack(padx=0, side=TOP, fill=X)
        self.exportSettingsMenu.pack(padx=0, side=TOP, fill=X)
        self.spectralFittingSelectionMenu.pack(padx=0, side=TOP, fill=X)

        ########### Check boxes ###########

        self.xlogTK = BooleanVar(self.frame_controls)
        self.ylogTK = BooleanVar(self.frame_controls)

        self.checkbutton_xscale = Checkbutton(self.frame_controls, bg=BGCOLOR, text='X Log', variable=self.xlogTK,
                                             command=lambda: self.on_select_logscale('x'), fg=GUIFONTCOLOR)
        self.checkbutton_xscale['state'] = 'disabled'
        self.checkbutton_xscale.pack(padx=(10, 0), pady=10, side=LEFT)
        self.checkbutton_yscale = Checkbutton(self.frame_controls, bg=BGCOLOR, text='Y Log', variable=self.ylogTK,
                                             command=lambda: self.on_select_logscale('y'), fg=GUIFONTCOLOR)
        self.checkbutton_yscale.pack(padx=(5, 0), pady=10, side=LEFT)

        # Pack up the frame control
        self.frame_controls.grid(column=0, row=0, sticky=W+N+E)

        ########### Labels ###########

        # The label frame
        header = {key: None for key in ['OBJECT', 'TELESCOP', 'INSTRUME', 'DATATYPE']}
        header = self.gspec.data[self.dataname]['data'].header['PRIMARY']
        labels = Frame(self, bg=BGCOLOR, relief=FRAMERELIEF, pady=0, padx=10)

        self.eventLabel = Label(labels, bg=BGCOLOR, text="Event: %s" % header['OBJECT'], fg=GUIFONTCOLOR)
        self.telescopeLabel = Label(labels, bg=BGCOLOR, text="Telescope: %s" % header['TELESCOP'], fg=GUIFONTCOLOR)
        self.instrumentLabel = Label(labels, bg=BGCOLOR, text="Instrument: %s" % header['INSTRUME'], fg=GUIFONTCOLOR)
        self.detectorLabel = Label(labels, bg=BGCOLOR, text="Detector: %s" % self.detector, fg=GUIFONTCOLOR)
        self.dataTypes = Label(labels, bg=BGCOLOR, text="Data types: %s" % header['DATATYPE'], fg=GUIFONTCOLOR)
        self.blankSpace = Label(labels, bg=BGCOLOR, text="").grid(row=6, column=0, sticky=N+W)
        self.selectionLabel = Label(labels, bg=BGCOLOR, text="Selections:", fg=GUIFONTCOLOR)

        self.time_selection_text = StringVar(labels)
        # self.set_time_selection_label(self.gspec.data[self.dataname]['energyview'].time_bound)
        self.time_selection_text.set(" -- : -- s")
        self.time_selection = Label(labels, bg=BGCOLOR, 
                                    textvariable=self.time_selection_text, 
                                    fg=GUIFONTCOLOR)

        self.energy_selection_text = StringVar(labels)
        self.set_energy_selection_label(self.gspec.data[self.dataname]['timeview'].energy_bound)
        #self.energy_selection_text.set(" -- - -- keV")
        self.energy_selection = Label(labels, bg=BGCOLOR, 
                                      textvariable=self.energy_selection_text, 
                                      fg=GUIFONTCOLOR)


        # Pack up the labels
        self.eventLabel.grid(row=1, column=0, sticky=N+W)
        self.telescopeLabel.grid(row=2, column=0, sticky=N+W)
        self.instrumentLabel.grid(row=3, column=0, sticky=N+W)
        self.detectorLabel.grid(row=4, column=0, sticky=N+W)
        self.dataTypes.grid(row=5, column=0, sticky=N+W)
        self.selectionLabel.grid(row=7, column=0, sticky=N+W)
        self.time_selection.grid(row=8, column=0, sticky=N + W)
        self.energy_selection.grid(row=9, column=0, sticky=N + W)


        # Add the labels to the grid
        labels.grid(column=0, row=1, sticky=N+W)


        ########### Keyboard bindings ###########

        # Bind the return key to the rebin function
        # RDP 6/27/2018: Added some adjust_source key bindings and removed 
        # (redundant) x and y bindings.
        root.bind('<x>', self.keyboard_event)
        root.bind('<y>', self.keyboard_event)
        root.bind('<z>', self.keyboard_event)
        root.bind('<s>', self.keyboard_event)
        root.bind('<b>', self.keyboard_event)
        root.bind('<t>', self.keyboard_event)
        root.bind('<R>', self.keyboard_event)
        root.bind('<Key-4>', self.keyboard_event)
        root.bind('<Key-5>', self.keyboard_event)
        root.bind('<6>', self.keyboard_event)
        root.bind('<7>', self.keyboard_event)
        root.bind('<8>', self.keyboard_event)
        root.bind('<9>', self.keyboard_event)
        root.bind('<Escape>', self.keyboard_event)

        ########### Plot Window ###########

        # Create a tk frame that will contain the plotting window
        # and a child that contains the matplotlib toolbar
        self.plot_frame = Frame(root, bg=BGCOLOR)
        self.lightcurve_toolbar_frame = Frame(self.plot_frame, bg=BGCOLOR)
        self.spectrum_toolbar_frame = Frame(self.plot_frame, bg=BGCOLOR)
        
        # the lightcurve canvas/toolbar
        self.lightcurve, self.lightcurve_toolbar = self._construct_plotter(plottype='lightcurve')
        self.current_plotter = self.lightcurve
        self.current_toolbar = self.lightcurve_toolbar
        self.lightcurve.plot()
        
        # the spectrum canvas/toolbar
        self.spectrum, self.spectrum_toolbar = self._construct_plotter(plottype='spectrum')
        self.spectrum.plot()

        # display the lightcurve on initialization
        self.lightcurve.canvas.get_tk_widget().grid(column=0, row=0, padx=0, 
                                                    pady=0, sticky=N+W+S+E)
        self.lightcurve_toolbar_frame.grid(column=0, row=1, sticky=N+W+E)
        self.current_toolbar_frame = self.lightcurve_toolbar_frame

        self.filemenu.entryconfig('Screenshot', command=self.lightcurve_toolbar.save_figure)



        # self.lightcurve_toolbar.rowconfigure(0, weight=1)
        # self.lightcurve_toolbar.columnconfigure(1, weight=1)

        # self.selection_mode_buttons.rowconfigure(0, weight=1)
        # self.selection_mode_buttons.columnconfigure(0, weight=1)

        # self.selection_mode_buttons(1, weight=1)

        # self.lightcurve_toolbar.update()

        # self.selection_mode_buttons.pack(expand=1, fill=X)
        # self.selection_mode_buttons.rowconfigure(0, weight=1)
        # self.selection_mode_buttons.columnconfigure(0, weight=1)
        # self.selection_mode_buttons.columnconfigure(1, weight=1)
        # self.selection_mode_buttons.columnconfigure(2, weight=1)

        # self.lightcurve_toolbar_frame.columnconfigure(3, weight=1)


        # self.button1.grid(column=0, row=0, sticky=E)
        # self.button2.grid(column=1, row=0, sticky=E)
        # self.button3.grid(column=2, row=0, sticky=E)

        # self.buttonFrame.grid(column=0, row=1, sticky=E, padx=25)

        # self.lightcurve_toolbar_frame.grid(column=0, row=1, sticky=N+W)

        ########### Row/Column Configurations ###########

        # Configuration needed with matplotlib toolbar
        # self.rowconfigure(0, weight=1)
        # self.rowconfigure(1, weight=1)
        # self.columnconfigure(0, weight=0)
        # self.columnconfigure(1, weight=1)

        root.rowconfigure(0, weight=1)
        root.columnconfigure(0, weight=0)
        root.columnconfigure(1, weight=2)

        # clicking the close window icon
        root.protocol("WM_DELETE_WINDOW", self.on_window_close)

        return
    
    def on_window_close(self):
        if self.gspec_root is not None:
            self.gspec_root.unregister_window(self.name)
        self.root.destroy()            
    
    def _construct_plotter(self, plottype=None):
        # construct the plotter and toolbar for the lightcurve/spectrum
        if plottype == 'lightcurve':
            plotter = Lightcurve(self.gspec, self.dataname,
                                    canvas=(FigureCanvasTkAgg, self.plot_frame))
            toolbar = NavigationToolbar2TkAgg(plotter.canvas, self.lightcurve_toolbar_frame)
            # toolbar = CustomToolbar(plotter.canvas, self.lightcurve_toolbar_frame)


        elif plottype == 'spectrum':
            plotter = Spectrum(self.gspec, self.dataname,
                                     canvas=(FigureCanvasTkAgg, self.plot_frame))
            toolbar = NavigationToolbar2TkAgg(plotter.canvas, self.spectrum_toolbar_frame)
        else:
            return

        plotter.canvas.show()
        # Register callback function for button press and motion
        plotter.canvas.callbacks.connect('button_press_event', self.on_mouse_click)
        plotter.canvas.callbacks.connect('motion_notify_event', self.on_mouse_move)
        # Register callback functions to catch when the user zooms or pans
        plotter.ax.callbacks.connect('xlim_changed', self.on_view_change)
        plotter.ax.callbacks.connect('ylim_changed', self.on_view_change)

        toolbar.config(bg=BGCOLOR)
        toolbar._message_label.config(background=BGCOLOR)
        toolbar.update()
        self.plot_frame.rowconfigure(0, weight=1)
        self.plot_frame.rowconfigure(1, weight=0)
        self.plot_frame.columnconfigure(0, weight=1)
        self.plot_frame.columnconfigure(1, weight=0)
        self.plot_frame.grid(column=1, row=0, padx=0, pady=0, sticky=N+W+S+E, rowspan=2)

        return plotter, toolbar

    def disable_controls(self):
        for child in self.frame_controls.winfo_children():

            child.configure(state='disable')

            if type(child) == OptionMenu:
                menu = child["menu"]
                last = menu.index("end")
                for i in range(last+1):
                    menu.entryconfigure(i, state="disabled")

    def enable_controls(self):
        for child in self.frame_controls.winfo_children():
            child.configure(state='normal')
            if type(child) == OptionMenu:
                menu = child["menu"]
                last = menu.index("end")
                for i in range(last+1):
                    menu.entryconfigure(i, state="normal")

        if self.lightcurve_display:
            self.checkbutton_xscale['state'] = 'disabled'
        else:
            self.fit_background_button['state'] = 'disabled'
        if self.gspec_root.xspec_version == None:
            self.spectralFittingSelectionMenu['state'] = 'disabled'

    def disable_toolbar_actions(self):
        if self.current_toolbar.mode == 'zoom rect':
            self.current_toolbar.zoom()
        if self.current_toolbar.mode == 'pan/zoom':
            self.current_toolbar.pan()
    
    def disable_toolbar_buttons(self):
        for child in self.current_toolbar.winfo_children():
            try:
                child.configure(state='disable')
            except:
                pass

    def disable_interactve_mode_buttons(self):
        for child in self.selection_mode_buttons.winfo_children():
            try:
                child.configure(state='disable')
            except:
                pass

    def enable_interactve_mode_buttons(self):
        for child in self.selection_mode_buttons.winfo_children():
            try:
                child.configure(state='normal')
            except:
                pass


    def enable_toolbar_buttons(self):
        for child in self.current_toolbar.winfo_children():
            try:
                child.configure(state='normal')
            except:
                pass

    def set_time_selection_label(self, bounds):
        text = " %0.2f : %0.2f s" % bounds
        self.time_selection_text.set(text)

    def set_energy_selection_label(self, bounds):
        text = " %0.2f : %0.2f keV" % bounds
        self.energy_selection_text.set(text)
    
    def toggle(self):

        # Disable any ongoing toolbar actions
        self.disable_toolbar_actions()

        if self.source_selection_active:
            self.accept_selections()

        # toggle between the lightcurve and spectrum display
        # the canvases and toolbars are already constructed and applied to the
        # frame, but we are toggling of the display of the widget frame itself
        if self.lightcurve_display:
            self.lightcurve.canvas.get_tk_widget().grid_forget()
            self.lightcurve_toolbar_frame.grid_forget()
            self.spectrum.canvas.get_tk_widget().grid(column=0, row=0, padx=0, 
                                                      pady=0, sticky=N+W+S+E)
            self.spectrum_toolbar_frame.grid(column=0, row=1, sticky=N+W+E)
            self.current_plotter = self.spectrum
            self.current_toolbar = self.spectrum_toolbar
            self.current_toolbar_frame = self.spectrum_toolbar_frame
            self.filemenu.entryconfig('Screenshot', command=self.spectrum_toolbar.save_figure)
            self.checkbutton_xscale['state'] = 'normal'
            self.xlogTK.set(self.xlog_spectrum)
            self.ylogTK.set(self.ylog_spectrum)
            self.fit_background_button['state'] = 'disabled'
            # update the binning menu
            # TODO: Fix Signal to Noise and add back in
            binchoices = ['Full Resolution',
                          'Combine Into Single Bin', 'Combine by Factor']
            self.rebin_menu['menu'].delete(0, 'end')
            for item in binchoices:
                self.rebin_menu['menu'].add_command(label=item, 
                                     command=lambda item=item: self.rebin(item))
            # update the selection menu
            selectchoices = ['Source Interactive', 'Clear Selections']
            self.source_selection_menu['menu'].delete(0, 'end')
            for item in selectchoices:
                self.source_selection_menu['menu'].add_command(label=item, 
                              command=lambda item=item: self.select_source(item))
            self.toggleButton.config(text=u"Show Lightcurve")
        else:
        
            self.spectrum.canvas.get_tk_widget().grid_forget()
            self.spectrum_toolbar_frame.grid_forget()
            self.lightcurve.canvas.get_tk_widget().grid(column=0, row=0, padx=0, 
                                                        pady=0, sticky=N+W+S+E)
            self.lightcurve_toolbar_frame.grid(column=0, row=1, sticky=N+W+E)
            self.current_plotter = self.lightcurve
            self.current_toolbar = self.lightcurve_toolbar
            self.current_toolbar_frame = self.lightcurve_toolbar_frame
            self.filemenu.entryconfig('Screenshot', command=self.lightcurve_toolbar.save_figure)
            self.checkbutton_xscale['state'] = 'disabled'
            self.xlogTK.set(False)
            self.ylogTK.set(self.ylog_lightcurve)
            self.fit_background_button['state'] = 'normal'
            # update the binning menu
            binchoices = ['Full Resolution', 'Temporal Resolution', 
                          'Signal to Noise', 'Combine Source Intervals', 
                          'Combine Into Single Bin', 'Combine by Factor']
            self.rebin_menu['menu'].delete(0, 'end')
            for item in binchoices:
                self.rebin_menu['menu'].add_command(label=item, 
                                     command=lambda item=item: self.rebin(item))
            # update the selection menu
            selectchoices = ['Source Interactive', 'Source By S/N', 
                             'Clear Selections']
            self.source_selection_menu['menu'].delete(0, 'end')
            for item in selectchoices:
                self.source_selection_menu['menu'].add_command(label=item, 
                              command=lambda item=item: self.select_source(item))
            self.toggleButton.config(text=u"Show Spectrum")

        self.lightcurve_display = not self.lightcurve_display
        
    def on_view_change(self, ax):
        # update the view stored in the lookup
        xlim = self.current_plotter.ax.get_xlim()
        ylim = self.current_plotter.ax.get_ylim()
        view = View(xlim[0], xlim[1], ylim[0], ylim[1])
        if self.lightcurve_display:
            self.gspec.lookup[self.dataname].views.time = view
        else:
            self.gspec.lookup[self.dataname].views.energy = view
            
    def display_clickzones(self):

        figbox = self.current_plotter.ax.figbox.bounds
        bounds = (figbox[0], figbox[0]+figbox[2], figbox[1], figbox[1]+figbox[3])
        # self.accept_zone = Rectangle((0,0), bounds[0], 1.0, fill=True, color='#afd7b4', alpha=0.1,
        #                  transform=self.current_plotter.figure.transFigure,
        #                  figure=self.current_plotter.figure)
        # self.current_plotter.figure.patches.extend([self.accept_zone])
        # self.accept_text = self.current_plotter.ax.annotate('Accept', xy=(0.02, 0.95), zorder=1000,
                                         # color='#afd7b4', xycoords='figure fraction', alpha=0)

        self.clear_zone = Rectangle((bounds[1],0), 1.0, 1.0, fill=True, color='#e39f9f', alpha=0.1,
                         transform=self.current_plotter.figure.transFigure,
                         figure=self.current_plotter.figure)
        self.current_plotter.figure.patches.extend([self.clear_zone])
        self.clear_text = self.current_plotter.ax.annotate('Clear', xy=(0.94, 0.95), zorder=1000,
                                         color='#e39f9f', xycoords='figure fraction', alpha=0, size=PLOTFONTSIZE)

        self.input_zone = Rectangle((bounds[0],0), bounds[1]-bounds[0], bounds[2], fill=True, color='#e4e56d', alpha=0.1,
                         transform=self.current_plotter.figure.transFigure,
                         figure=self.current_plotter.figure)
        self.current_plotter.figure.patches.extend([self.input_zone])
        self.input_text = self.current_plotter.ax.annotate('Manual Input', xy=(0.75, 0.02), zorder=1000,
                                         color='#e4e56d', xycoords='figure fraction', alpha=0, size=PLOTFONTSIZE)

    def hide_clickzones(self):
        self.accept_text.remove()
        self.clear_text.remove()
        self.input_text.remove()
        self.current_plotter.figure.patches = []
        self.accept_text = None
        self.clear_text = None
        self.input_text = None
        self.accept_zone = None
        self.clear_zone = None
        self.input_zone = None

    def on_mouse_click(self, event):

        if not self.source_selection_active and \
           not self.background_selection_active and \
           not self.binning_selection_active and \
           not self.snr_selection_active:
            return
        
        # within plot axes
        if event.inaxes is not None:
            if self.source_selection_active:
                self.current_plotter.plot_source_selection_line(event.xdata)

            if self.binning_selection_active or self.snr_selection_active:
                self.current_plotter.plot_binning_selection_line(event.xdata)
            
            elif self.background_selection_active:
                self.current_plotter.plot_background_selection_line(event.xdata)

            self.accept_cancel_text = 'Accept Selections'
            self.exit_selection_button["text"] = u"Accept Selections"

            self.current_plotter.canvas.draw()
        
        # outside plot axes
        elif event.inaxes is None:
            fig = self.current_plotter.figure
            bbox = self.current_plotter.ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            
            # exit and register selections and hide any tooltips
            if event.x < (bbox.x0 * fig.dpi):
                self.accept_selections()

            # clear selection
            elif event.x > ((bbox.x0 + bbox.width) * fig.dpi):        
                self.clear_selections()
                
            # manual input
            elif event.y < (bbox.y0 * fig.dpi):
                self.manual_entry()

    def on_mouse_move(self, event):

        if (not self.source_selection_active) and \
           (not self.background_selection_active) and \
           (not self.binning_selection_active) and \
           (not self.snr_selection_active):

            return
        
        # Within plot axes
        if event.inaxes is not None: 

            text = '(x: %.2f, y: %.2f)' % (event.xdata, event.ydata)
        
            self.current_plotter.plot_text(text, (0.02, 0.02),
                                           xycoords='figure fraction', size=PLOTFONTSIZE)

            self.current_plotter.canvas.draw()

        # Outside plot axes
        elif event.inaxes is None: 
            fig = self.current_plotter.figure
            bbox = self.current_plotter.ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            if event.x < (bbox.x0 * fig.dpi):

                # Update the tooltip text
                text = self.accept_cancel_text

            elif event.x > ((bbox.x0 + bbox.width) * fig.dpi):
                text = 'Clear Selections'
            elif event.y < (bbox.y0 * fig.dpi):
                text = 'Manual Entry'
            else:
                text = ''
                self.current_plotter.figure.patches = []

            self.current_plotter.plot_text(text, (0.02, 0.02),
                                           xycoords='figure fraction', size=PLOTFONTSIZE)

            self.current_plotter.canvas.draw()



    def _get_selection_intervals(self):
        selection_intervals = []
        temp_list = []
        artists = self.current_plotter.annotations['selections']['artists']
        # Make sure there are an even number of selections
        if len(artists) % 2  == 1:
            artists = artists[:-1]
        
        # Get the x coordinates of the user selections
        xlimit_selections = [artist[0].get_xdata()[0] for artist in artists]
        
        while len(xlimit_selections) != 0:
            # Get the upper and lower selections
            xmax = xlimit_selections.pop()
            xmin = xlimit_selections.pop()

            # Swap the xmin and xmax values if the selections were made in 
            # reverse order
            if (xmax > xmin) == False:
                xmin, xmax = xmax, xmin
            
            # Add the selections to a temporary list
            temp_list.append([xmin, xmax])
        
        # now sort all of the selection and merge any overlapped selections
        temp_list.sort()
        num_selects = len(temp_list)
        if num_selects == 0:
            return selection_intervals
        
        selection_intervals.append(temp_list[0])
        for i in range(1,num_selects):
            if temp_list[i][0] <= temp_list[i-1][1]:
                if temp_list[i][1] > temp_list[i-1][1]:
                    selection_intervals[-1][1] = temp_list[i][1]
            else:
                selection_intervals.append(temp_list[i])
            
        return selection_intervals
    
    def manual_selection(self, x0, x1, y0, y1):

        if x0 is None and x1 is None and y0 is None and y1 is None:
            self.enable_interactve_mode_buttons()
            return

        if y0 is None and y1 is None:
            if self.background_selection_active:
                self.current_plotter.plot_background_selection_line(x0)
                self.current_plotter.plot_background_selection_line(x1)
            elif self.source_selection_active or self.binning_selection_active:
                self.current_plotter.plot_source_selection_line(x0)
                self.current_plotter.plot_source_selection_line(x1)
            else:
                pass
        else:
            pass # zoom

        self.current_plotter.canvas.draw()
        self.accept_cancel_text = 'Accept Selections'
        self.exit_selection_button["text"] = u"Accept Selections"

        self.enable_interactve_mode_buttons()
    
    def write_lookup(self):
        datanames = list(self.gspec.data.keys())
        datafile = file.GbmFile.from_path(self.gspec.data[datanames[0]]['data'].filename)
        directory = datafile.directory
        uid = datafile.uid
        data_type = datafile.data_type
        trigger = datafile.trigger
        meta = '_lookup'
        extension = 'json'
        lufile = file.GbmFile.create(uid=uid, data_type=data_type,
                                     trigger=trigger, meta=meta,
                                     extension=extension)
        initialfile = lufile.basename()
        filename = filedialog.asksaveasfilename(initialdir=directory,
                                                initialfile=initialfile,
                                                title="Save Lookup File",
                                                filetypes=(('JSON file', '*.json'),),
                                                defaultextension='.json',
                                                parent=self.root)

        filename = self.gspec.save_lookup(filename=filename)
        print('Lookup file saved to {0}'.format(filename))

    def read_lookup(self):
        filename = filedialog.askopenfilename(filetypes=self.lookupfiletypes, 
                                              parent=self.root)
        if filename == '':
            return

        # load in the appropriate lookups
        if filename.split('.')[-1] == 'lu':
            tifile = None
            if 'tte' in filename:
                test_tifile = filename.split('.')[0] + '.ti'
                if os.path.isfile(test_tifile):
                    tifile = test_tifile
            lu = LookupFile.read_from_rmfit(filename, ti_file=tifile, dataname=self.dataname)
        else:
            lu = LookupFile.read_from(filename)
        self.gspec.load_gspec_lookup(self.dataname, lu)

        # RDP: 6/27/2018 Selection Labels weren't updating. Fixed.
        self.set_time_selection_label(self.gspec.data[self.dataname]['energyview'].time_bound)
        self.set_energy_selection_label(self.gspec.data[self.dataname]['timeview'].energy_bound)

        self.lightcurve.clear()
        self.lightcurve.display_lookup_view()
        self.lightcurve.plot()
        self.lightcurve.canvas.draw()
        self.spectrum.clear()
        self.spectrum.display_lookup_view()
        self.spectrum.plot()
        self.spectrum.canvas.draw()
    
    def display_header(self, ext):
        headers = self.gspec.data[self.dataname]['data'].header
        header = headers[ext].tostring(sep='\n', padding=False)
        title = '{0} header for {1}'.format(ext, self.dataname)
        font = ("Andale Mono", 12)
        height = len(headers[ext])+1
        if (height*16.0) > 608:
            height=38
        TextDisplayDialog(title, header, font=font, height=height)

    def display_lookup(self):
        lu = self.gspec.lookup.split_off_dataname(self.dataname).display_lookup()
        title = 'Lookup for {0}'.format(self.dataname)
        font = ("Roboto Mono", 12)
        height = lu.count('\n')+1
        if (height*16.0) > 608:
            height=38
        TextDisplayDialog(title, lu, font=font, width=40, height=height)
    
    def erase_lookup(self):
        self.gspec.remove_data(self.dataname)
        self.gspec.add_data(self.dataname)
        self.lighcurve.clear()
        self.lightcurve.plot()
        self.lightcurve.canvas.draw()
        self.spectrum.clear()
        self.spectrum.plot()
        self.spectrum.canvas.draw()
    
    def load_response(self):
        # allowed filetypes
        filename = filedialog.askopenfilename(filetypes=self.responsefiletypes, 
                                              parent=self.root)
        if filename == '':
            return

        self.gspec.add_response(self.dataname, filename)
        return filename
    
    def remove_response(self):
        self.gspec.data[self.dataname]['response'] = None
        self.gspec.lookup[self.dataname].add_response(None)
            
    def rebin(self, menu_selection):

        # Reset the menu
        self.rebin_selection.set('Rebin Data')

        # Disable any ongoing toolbar actions (zoom/pan)
        self.disable_toolbar_actions()
        
        title = 'Rebin'
        command = lambda value, option: self.rebin_command(value, option)

        if 'Full Resolution' in menu_selection:
            if self.lightcurve_display:
                self.gspec.clear_time_binning(self.dataname)
            else:
                self.gspec.clear_energy_binning(self.dataname)
            self.gspec.refresh(self.dataname)
            # remove old lightcurve and selection artists
            self.current_plotter.remove_data()
            self.current_plotter.remove_source()
            # plot the new lightcurve and selection artists
            self.current_plotter.plot_data()
            self.current_plotter.plot_source()
            self.current_plotter.canvas.draw()
            
        elif 'Temporal Resolution' in menu_selection:

            # Disable the dataviewer controls while the dialog widget is open
            self.disable_controls()
            self.disable_toolbar_buttons()

            self.current_binning_function = 'Temporal Resolution'
            message = ['Enter new data resolution in seconds',
                       'NOTE: Resolution must be an integer\nmultiple of the independent axis data']
            TextOptionDialog(title, master=self.root, width=275, height=195, 
                             message=message, default_value=self.current_binning,
                             button_label='Rebin', button_command=command)

        elif 'Signal to Noise' in menu_selection:

            # Disable the dataviewer controls while the dialog widget is open
            self.disable_controls()
            self.disable_toolbar_buttons()

            # check for required background
            if self.gspec.data[self.dataname]['bkgdmodel'] is None:
                message = 'No background model exists.  Background fit must ' \
                          'exist to perform S/N binning.'
                response = messagebox.showerror('Error', message, parent=self.root)
                return

            self.current_binning_function = 'Signal-to-Noise'
            message = ['Enter the desired SNR per bin']
            TextOptionDialog(title, master=self.root, width=275, height=145, 
                             message=message, default_value=5.0,
                             button_label='Rebin',
                             button_options=['Current Selection','Custom Range'],
                             button_command=command)
        
        elif 'Combine Source Intervals' in menu_selection:
            if self.gspec.data[self.dataname]['sourceview'] is None:
                raise ValueError('No Source Intervals Selected')

            times = self.gspec.data[self.dataname]['sourceview']._time_selections
            for selection in times:
                self.current_binning_function = 'Combine Into Single Bin'
                self.current_binning = None              
                self.apply_binning(start=selection[0], stop=selection[1])

        elif 'Combine Into Single Bin' in menu_selection:
            self.current_binning_function = 'Combine Into Single Bin'
            self.current_binning = None
            self.binning_selection_active = True
            # self.display_clickzones()

        
        elif 'Combine by Factor' in menu_selection:

            # Disable the dataviewer controls while the dialog widget is open
            self.disable_controls()
            self.disable_toolbar_buttons()

            self.current_binning_function = 'Combine By Factor'
            message = ['Enter the desired rebinning factor']
            TextOptionDialog(title, master=self.root, width=275, height=150, 
                             message=message, default_value=2,
                             button_label='Rebin', button_command=command)
       
    def rebin_command(self, value, fullrange):

        # Re-enable controls and exit if the user closed the window without specifying anything
        if value is None and fullrange is None:
            self.enable_controls()
            self.enable_toolbar_buttons()            
            return

        try:
            value_float = float(value.get())
        except:
            value_float = value

        if value_float is not None:
            self.current_binning = value_float
        
        fullrange = fullrange.get()
        if fullrange == 'Full':

            # Apply the binning over the full range
            self.apply_binning()
            
            # Re-enable controls
            self.enable_controls()
            self.enable_toolbar_buttons()

        else:
            # Activeate binning selection. This state change will be picked up by on_mouse_move()
            self.binning_selection_active = True

            # Activate interactive mode
            self.activate_interactive_mode()



    def apply_binning(self, start=None, stop=None):

        if self.lightcurve_display:
            datatype = self.gspec.data[self.dataname]['data'].type
            if self.current_binning is None:
                self.gspec.time_binning(self.dataname, 
                                        self.current_binning_function,
                                        datatype=datatype, start=start, stop=stop)
            else:
                self.gspec.time_binning(self.dataname, 
                                        self.current_binning_function,
                                        self.current_binning, datatype=datatype,
                                        start=start, stop=stop)
            self.gspec.refresh(self.dataname, background=False)
        # spectral binning
        else:
            if self.current_binning is None:
                self.gspec.energy_binning(self.dataname, 
                                          self.current_binning_function, 
                                          start=start, stop=stop)
            else:
                self.gspec.energy_binning(self.dataname, 
                                          self.current_binning_function, 
                                          self.current_binning, start=start, 
                                          stop=stop)
            self.gspec.refresh(self.dataname)

        self.current_plotter.remove_data()
        self.current_plotter.remove_source()
        self.current_plotter.plot_data()
        self.current_plotter.plot_source()
        self.current_plotter.canvas.draw()

    def select_background(self):

        # # Disable the control buttons
        # self.disable_controls() 

        # # Disable the toolbar actions and gray out the buttons
        # self.disable_toolbar_actions()
        # self.disable_toolbar_buttons()        
        
        # # Add the guidance tool tip
        # self.tooltip = self.current_plotter.ax.annotate('Click on the plot to make selections. Cancel or accept selections to exit selection mode.', (0.16, 0.94), xycoords='figure fraction', size=PLOTFONTSIZE)

        # # Create a blue border to indicate that the user is in selection mode
        # self.current_plotter.ax.spines['bottom'].set_color('blue')
        # self.current_plotter.ax.spines['top'].set_color('blue')
        # self.current_plotter.ax.spines['left'].set_color('blue')
        # self.current_plotter.ax.spines['right'].set_color('blue')

        # # Color the ticks as well as the border
        # for tick in self.current_plotter.ax.xaxis.get_majorticklines():
        #     tick.set_color('blue')
        # for tick in self.current_plotter.ax.yaxis.get_majorticklines():
        #     tick.set_color('blue')
        # for tick in self.current_plotter.ax.xaxis.get_minorticklines():
        #     tick.set_color('blue')
        # for tick in self.current_plotter.ax.yaxis.get_minorticklines():
        #     tick.set_color('blue')


        # if self.gspec.data[self.dataname]['bkgdview'] is not None:
        #     for selection in self.gspec.data[self.dataname]['bkgdview']._time_selections:
        #         self.lightcurve.plot_background_selection_line(selection[0])
        #         self.lightcurve.plot_background_selection_line(selection[1])
                
        # # Update the plot
        # self.current_plotter.canvas.draw() 

        # # Add the selection mode buttons
        # self.selection_mode_buttons_super.pack()

        # Activate interactive selection mode
        self.activate_interactive_mode()

        # Activate selection mode       
        self.background_selection_active = True



    def background_fit_options(self):
        
        # Setup the message
        message = ['Select Order of', 'Background Polynomial']
        options = ['Cancel', '0', '1', '2', '3', '4']
        commands = [self.cancel_background_fit, 
                    lambda: self.fit_background_selections(0),
                    lambda: self.fit_background_selections(1),
                    lambda: self.fit_background_selections(2),
                    lambda: self.fit_background_selections(3),
                    lambda: self.fit_background_selections(4)]
        
        # Display the dialog
        self.background_options_dialog = OptionDialog('Background', options, 
                                                      commands, master=self.root,
                                                      message=message)
         
    def fit_background_selections(self, order):

        # Close the dialog window
        self.background_options_dialog.root.destroy()
        self.background_options_dialog = None
        selections = self.background_selections
        
        # Add a background model to the detector
        self.gspec.add_background(self.dataname, selections, 'Polynomial', order)
        self.gspec.refresh(self.dataname)
        # Redraw the plotting window
        self.lightcurve.remove_background()
        self.lightcurve.plot_background()
        self.lightcurve.canvas.show()
        self.spectrum.remove_background()
        self.spectrum.plot_background()
        self.spectrum.canvas.show()
        
        # create the reduced chisq plot dialog
        height = self.height
        width = int(round(1.33*height))
        plotdialog = PlotDialog('Background Reduced Chi-Square', master=self.root, 
                                width=width, height=height)
        plotter = BackgroundChisq(self.gspec, self.dataname, 
                                  canvas=(FigureCanvasTkAgg, plotdialog.root))
        plotter.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        plotter.plot_reduced_chisq()
        plotter.canvas.show()

        # Re-enable the toolbar
        self.enable_toolbar_buttons()        


    def cancel_background_fit(self):
        self.background_options_dialog.root.destroy()
        self.background_options_dialog = None

    def zoom_blocker(self, *args, **kwargs):
        print('zoom pressed')
        return "break"

    def activate_interactive_mode(self):
    
        # Disable controls
        self.disable_controls()

        # Disable the toolbar actions and gray out the buttons
        self.disable_toolbar_actions()
        self.disable_toolbar_buttons()

        # Create the frames that will contain the interactive mode buttons
        self.selection_mode_buttons_super = Frame(self.current_toolbar, bg=BGCOLOR)
        self.selection_mode_buttons = Frame(self.selection_mode_buttons_super, bg=BGCOLOR)
        self.selection_mode_buttons.pack(expand=1, fill=X)

        # Create the interactive mode buttons
        self.manual_entry_button = Button(self.selection_mode_buttons, text=u"Manual Entry", highlightbackground=BGCOLOR,
                                   command=self.manual_entry, fg=GUIFONTCOLOR, width=14)
        self.clear_selection_button = Button(self.selection_mode_buttons, text=u"Clear Selections", highlightbackground=BGCOLOR,
                                   command=self.clear_selections, fg=GUIFONTCOLOR, width=14)

        self.exit_selection_button = Button(self.selection_mode_buttons, text=u"Cancel", highlightbackground=BGCOLOR,
                                   command=self.accept_selections, fg=GUIFONTCOLOR, width=14, state=ACTIVE)

        # Set the default tooltip text 
        self.accept_cancel_text = 'Cancel'        

        # Packup the buttons
        self.manual_entry_button.grid(row=0, column=0, sticky=N+E+W+S, pady=5)
        self.clear_selection_button.grid(row=0, column=1, sticky=N+E+W+S, pady=5)
        self.exit_selection_button.grid(row=0, column=2, sticky=N+E+W+S, pady=5, padx=20)

        # Create a blue border to indicate that the user is in selection mode
        self.current_plotter.ax.spines['bottom'].set_color('blue')
        self.current_plotter.ax.spines['top'].set_color('blue')
        self.current_plotter.ax.spines['left'].set_color('blue')
        self.current_plotter.ax.spines['right'].set_color('blue')

        # Color the ticks as well as the border
        for tick in self.current_plotter.ax.xaxis.get_majorticklines():
            tick.set_color('blue')
        for tick in self.current_plotter.ax.yaxis.get_majorticklines():
            tick.set_color('blue')
        for tick in self.current_plotter.ax.xaxis.get_minorticklines():
            tick.set_color('blue')
        for tick in self.current_plotter.ax.yaxis.get_minorticklines():
            tick.set_color('blue')

        # Add the guidance tool tip
        self.tooltip = self.current_plotter.ax.annotate('Click on the plot to make selections. Cancel or accept selections to exit selection mode.', (0.16, 0.94), xycoords='figure fraction', size=PLOTFONTSIZE)

        # Make the source selection buttons visible
        self.selection_mode_buttons_super.pack()
        self.current_toolbar_frame.update()

        # Update the plot
        self.current_plotter.canvas.draw()   

    def deactivate_interactive_mode(self):

        # Re-enable controls
        self.enable_controls()

        # Re-enable the toolbar
        self.enable_toolbar_buttons()

        # Create a blue border to indicate that the user is in selection mode
        self.current_plotter.ax.spines['bottom'].set_color('black')
        self.current_plotter.ax.spines['top'].set_color('black')
        self.current_plotter.ax.spines['left'].set_color('black')
        self.current_plotter.ax.spines['right'].set_color('black')
        self.current_plotter.ax.tick_params(axis='x', colors='black')

        # Color the ticks as well as the border
        for tick in self.current_plotter.ax.xaxis.get_majorticklines():
            tick.set_color('black')
        for tick in self.current_plotter.ax.yaxis.get_majorticklines():
            tick.set_color('black')
        for tick in self.current_plotter.ax.xaxis.get_minorticklines():
            tick.set_color('black')
        for tick in self.current_plotter.ax.yaxis.get_minorticklines():
            tick.set_color('black')

        # Remove any selection lines
        self.current_plotter.remove_selection_lines()

        # Remove the tooltip text
        self.current_plotter.remove_text()

        # Destroy all of the frame's children
        for child in self.selection_mode_buttons_super.winfo_children():
            child.destroy()

        # Hide the selection mode buttons
        self.selection_mode_buttons_super.destroy()

        # Redraw the plot
        self.lightcurve.canvas.draw()
        self.spectrum.canvas.draw()

    def select_source(self, menu_selection):

        # Reset the menu selection
        self.source_selection.set('Source Selection')

        # Do nothing if the user is already in the source selection mode
        if self.source_selection_active:
            return

        # interactive source selection
        if 'Source Interactive' in menu_selection:

            # Active the interactive mode
            self.activate_interactive_mode()

            # self.display_clickzones()
            self.source_selection_active = True


        # select bins above a S/N within the interactive source selection
        elif 'Source By S/N' in menu_selection:
    
            # check for required background
            if self.gspec.data[self.dataname]['bkgdmodel'] is None:
                message = 'No background model exists.  Background fit must ' \
                          'exist to perform S/N selection.'
                response = messagebox.showerror('Error', message, parent=self.root)
                return
            
            # Disable the dataviewer controls
            self.disable_controls()
            self.disable_toolbar_buttons()

            title = 'Source by S/N'
            message = ['Enter Signal-to-Noise Ratio']
            command = lambda value, option: self.snr_select(value, option)            
            TextOptionDialog(title, master=self.root, width=275, height=150, 
                             message=message, default_value=5.0,
                             button_label='Apply', button_command=command)
        
        # clear the time selections:
        elif 'Clear Selections' in menu_selection:
            self.clear_selections()
            return

        Tk.update_idletasks(self)

        # Disable controls
        self.disable_controls()

    def manual_entry(self):

        self.disable_interactve_mode_buttons()

        title = 'Manual Selection'
        message = 'Manual Selection Input'
        command = lambda x0, x1, y0, y1: self.manual_selection(x0,x1,y0,y1)
        ManualInputDialog(title, command, master=self.root, 
                          message=message, xinput=True, yinput=False)

    def clear_selections(self):

        # source selection
        if self.lightcurve_display:
            self.gspec.remove_source_selections(self.dataname)
            self.set_time_selection_label(self.gspec.data[self.dataname]['energyview'].time_bound)
            self.time_selection_text.set(" -- : -- s")
            self.gspec.refresh(self.dataname)
            self.current_plotter.remove_source()   
        else:
            self.gspec.remove_energy_selections(self.dataname)
            self.set_energy_selection_label(self.gspec.data[self.dataname]['timeview'].energy_bound)
            self.energy_selection_text.set(" -- : -- keV")
            self.gspec.refresh(self.dataname)
            self.current_plotter.remove_source()   
        
        # background selection
        if self.background_selection_active:
            self.gspec.remove_background(self.dataname)
            self.gspec.refresh(self.dataname)
            self.current_plotter.remove_background()
        
        self.current_plotter.remove_selection_lines()
        self.current_plotter.canvas.draw()

        self.accept_cancel_text = 'Cancel'
        self.exit_selection_button["text"] = u"Cancel"

    def accept_selections(self):

        try:
            # Remove any tooltips
            self.tooltip.remove()
        except:
            pass

        # self.hide_clickzones()
        selections = self._get_selection_intervals()

        # source selection
        if self.source_selection_active:
            self.source_selection_active = False
            if len(selections) > 0:
                self.spectrum.clear()
                self.lightcurve.clear()
                if self.lightcurve_display:
                    self.gspec.source_selection(self.dataname, selections)
                    self.gspec.refresh(self.dataname)
                    self.set_time_selection_label((selections[0][0], 
                                                   selections[-1][1]))
                    self.lightcurve.plot(yscaling=False)
                    self.spectrum.plot(yscaling=True)
                else:
                    self.gspec.energy_selection(self.dataname, selections)
                    self.gspec.refresh(self.dataname)
                    self.set_energy_selection_label((selections[0][0], 
                                                     selections[-1][1]))
                    self.lightcurve.plot(yscaling=True)
                    self.spectrum.plot(yscaling=False)

        # background selection
        elif self.background_selection_active:
            self.background_selection_active = False
            self.background_selections = selections
            if len(self.background_selections) > 0:
                self.background_fit_options()
            
        # binning selection
        elif self.binning_selection_active:
            self.binning_selection_active = False
            if len(selections) > 0:
                for selection in selections:
                    self.apply_binning(start=selection[0], stop=selection[1])

        # snr selection
        elif self.snr_selection_active:
            self.snr_selection_active = False
            self.gspec.source_by_snr(self.dataname, selections, 
                                     self.current_snr)
            self.gspec.refresh(self.dataname, background=False)
            self.lightcurve.remove_source()
            self.lightcurve.plot_source()
            self.spectrum.clear()
            self.spectrum.plot()
                    
        
        self.deactivate_interactive_mode()

    def snr_select(self, value, fullrange):

        # Re-enable controls and exit if the user closed the window without specifying anything
        if value is None and fullrange is None:
            self.enable_controls()
            self.enable_toolbar_buttons()            
            return

        try:
            snr = float(value.get())
        except:
            snr = value

        if fullrange.get() == 'Full':
            time_interval = self.gspec.data[self.dataname]['data'].time_range
            self.gspec.source_by_snr(self.dataname, time_interval, snr)
            self.gspec.refresh(self.dataname, background=False)
            self.current_plotter.remove_source()
            self.current_plotter.plot_source()
        else:
            self.snr_selection_active = True
            self.current_snr = snr
            # self.display_clickzones()

    def adjust_source(self, menu_selection):
        self.disable_toolbar_actions()

        if self.lightcurve_display:
            if self.gspec.data[self.dataname]['sourceview'] is None:
                self.adjustSourceSelection.set('Adjust Source')
                return
                
            # Get the currently displayed light curve
            lc = self.gspec.data[self.dataname]['timeview'].lightcurve[0].bin_centers()

            # Get the lower and upper bounds of the current source selection
            selection = self.gspec.data[self.dataname]['sourceview'].lightcurve[0].bin_centers()
            selection_lower = selection[0]
            selection_upper = selection[-1]

            # Find the index of the lower and upper bounds of the current source selection
            index_lower = where(lc == selection_lower)[0]
            index_upper = where(lc == selection_upper)[0]

            # Get the maximum allowable index
            index_max = len(lc) - 1

        else:

            # Get the currently displayed spectrum
            spec = self.gspec.data[self.dataname]['energyview'].spectrum[0].bin_centers()
            
            # Get the lower and upper bounds of the current source selection
            selection = self.gspec.data[self.dataname]['timeview'].spectrum[0].bin_centers()
            selection_lower = selection[0]
            selection_upper = selection[-1]

            # Find the index of the lower and upper bounds of the current source selection
            index_lower = where(spec == selection_lower)[0]
            index_upper = where(spec == selection_upper)[0]

            # Get the maximum allowable index
            index_max = len(spec) - 1

        # Shift the index of the lower and upper bounds of the current source 
        # selection accordingly
        if '< Shift Selection' in menu_selection:
            if index_lower != 0:
                index_lower = index_lower - 1
                index_upper = index_upper - 1

        elif '> Shift Selection' in menu_selection:
            if index_upper != index_max:
                index_lower = index_lower + 1
                index_upper = index_upper + 1

        elif '< Left Selection' in menu_selection:
            if index_lower != 0:
                index_lower = index_lower - 1

        elif '> Left Selection' in menu_selection:
            if index_lower != index_max - 1:
                index_lower = index_lower + 1

        elif '< Right Selection' in menu_selection:
            if index_upper != 1:
                index_upper = index_upper - 1

        elif '> Right Selection' in menu_selection:
            if index_upper != index_max - 1:
                index_upper = index_upper + 1
        else:
            pass
        
        # RDP 6/27/2018: Added refresh and plot sequences to each branch to address the 
        # bug that prevented the selections from being updated.
        if self.lightcurve_display:
            selection_new = (lc[index_lower][0], lc[index_upper][0])
            self.gspec.source_selection(self.dataname, selection_new)
            self.set_time_selection_label(self.gspec.data[self.dataname]['energyview'].time_bound)            
        else:
            selection_new = (spec[index_lower][0], spec[index_upper][0])
            self.gspec.energy_selection(self.dataname, selection_new)
            self.set_energy_selection_label(self.gspec.data[self.dataname]['timeview'].energy_bound)
            
        # Refresh the detector and redraw the plot
        # RDP 6/27/2018: Removed refresh and plotting commands to fix selection update issue
        self.gspec.refresh(self.dataname)
        self.lightcurve.clear()
        self.lightcurve.plot()
        self.lightcurve.canvas.draw()
        self.spectrum.clear()
        self.spectrum.plot()
        self.spectrum.canvas.draw()

        # Reset the menu button
        self.adjustSourceSelection.set('Adjust Source')
    
    def export_settings(self, menu_selection):
        self.disable_toolbar_actions()

        # get other detectors
        datafiles = []
        windows = []
        for key, value in self.gspec_root._open_windows.items():
            if (self.dataname not in key) and isinstance(value, DataViewer):
                datafiles.append(key)
                windows.append(value)

        # the settings to export
        if 'Temporal Binning' in menu_selection:
            self.gspec.import_time_binning(self.dataname, datafiles)
            self.gspec.import_time_display_view(self.dataname, datafiles)
        elif 'Background' in menu_selection:
            self.gspec.import_background(self.dataname, datafiles)
            # self.gspec.import_time_display_view(self.dataname, datafiles)
        elif 'Source Selection' in menu_selection:
            self.gspec.import_source_selection(self.dataname, datafiles)
            self.gspec.import_time_display_view(self.dataname, datafiles)
        elif 'View Range' in menu_selection:
            self.gspec.import_time_display_view(self.dataname, datafiles)
        elif 'All of the Above' in menu_selection:
            self.gspec.import_time_binning(self.dataname, datafiles)
            self.gspec.import_background(self.dataname, datafiles)
            self.gspec.import_source_selection(self.dataname, datafiles)
            self.gspec.import_time_display_view(self.dataname, datafiles)

        # redraw the windows
        for window in windows:
            window.lightcurve.clear()
            window.lightcurve.plot()
            window.lightcurve.canvas.draw()
            window.spectrum.clear()
            window.spectrum.plot()
            window.spectrum.canvas.draw()

        self.exportSettingsSelection.set('Push Selections')
       
    def on_select_logscale(self, axis):
        if axis == 'x':
            self.current_plotter.toggle_xscale()
        elif axis == 'y':
            self.current_plotter.toggle_yscale()
        self.current_plotter.canvas.draw()
        
        if self.lightcurve_display:
            self.ylog_lightcurve = not self.ylog_lightcurve
        else:
            if axis == 'x':
                self.xlog_spectrum = not self.xlog_spectrum
            else:
                self.ylog_spectrum = not self.ylog_spectrum                
    
    def keyboard_event(self, event):
        # RDP 6/27/2018: Added some adjust_source key bindings and removed 
        # (redundant) x and y bindings.
        """
        if 'x' in event.char:
            self.current_toolbar.zoom()

        if 'y' in event.char:
            self.current_toolbar.zoom()
        """
        
        if 'z' in event.char:
            self.current_toolbar.zoom()

        if 's' in event.char:
            self.select_source(['Source Interactive'])

        if 't' in event.char:
            self.toggle()
        
        if 'R' in event.char:
            self.read_lookup()

        if '6' in event.char:
            self.adjust_source(['< Shift Selection'])

        if '7' in event.char:
            self.adjust_source(['> Shift Selection'])

        if '4' in event.char:
            self.adjust_source(['< Left Selection'])

        if '5' in event.char:
            self.adjust_source(['> Left Selection'])

        if '8' in event.char:
            self.adjust_source(['< Right Selection'])

        if '9' in event.char:
            self.adjust_source(['> Right Selection'])

        # immediately cancel out of all interactive selections, do not pass go
        if event.char == '\x1b':
            self.source_selection_active = False
            self.background_selection_active = False
            self.binning_selection_active = False
            self.snr_selection_active = False

            try:
                # self.hide_clickzones()
                pass
            except:
                pass
            self.current_plotter.remove_text()
            self.current_plotter.remove_selection_lines()
            self.current_plotter.canvas.draw()
            self.enable_controls()

    def spectral_fit_selection(self, menu_selection):

        self.disable_toolbar_actions()
        
        if 'Fit Plotter' in menu_selection: 
            # Launch the fit plotter and pass the current xspec session
            #if './Spectral Fit Display' not in self.gspec_root._openWindows.keys():
            #    FitPlotter(self.gspec, None)
            #else:
            #    fitPlotter = self.gspec_root._openWindows['./Spectral Fit Display']
            #    fitPlotter.root.lift()
            return
        elif 'Batch Fit Selections' in menu_selection: 
            batch = True            
        elif 'Fit Selection' in menu_selection:
            batch = False
        
        self.gspec_root.prepare_spectral_fit(batch=batch)
        
    def command(self):
        pass


# # Custom toolbar
# class CustomToolbar(NavigationToolbar2TkAgg):
#     def plot_axes(self):
#         # This function currently makes it so that the 'original view' is lost
#         # TODO Fix the above bug
#         self.canvas.figure.axes[0].set_xlim([10,60])
#         self.canvas.draw()

#     def __init__(self,canvas_,parent_):

#         NavigationToolbar2TkAgg.__init__(self,canvas_,parent_)

#     def _Button(self, text, file, command, extension='.png'):
#         img_file = os.path.join(
#             rcParams['datapath'], 'images', file + extension)
#         im = tkinter.PhotoImage(master=self, file=img_file)
#         b = tkinter.Button(
#             master=self, text=text, padx=2, pady=2, image=im, command=command)
#         b._ntimage = im
#         b.pack(side=Tk.LEFT)
#         return b
            
