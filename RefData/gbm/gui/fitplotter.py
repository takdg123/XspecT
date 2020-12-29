#! /usr/bin/env python3

from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
from tkinter import ttk

import matplotlib
matplotlib.rcParams.update({'font.size': 10})
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure#
from matplotlib.ticker import AutoMinorLocator#
import matplotlib.gridspec as gridspec#
import copy
import os

import numpy
from .globals import BGCOLOR, FRAMERELIEF, OM_JUSTIFY, GUIFONTCOLOR, PLOTFONTSIZE
from .window import Window
from gbm import data

from .dialogs import PlotModelParametersDialog

import random

# Use the classic matplotlib style sheet
if float(matplotlib.__version__[0]) >= 2.0:
    import matplotlib.style
    matplotlib.style.use('classic')

##########################################################################################

class FitPlotter(Window):
    """Class for displaying pha files"""
    
    name = 'Spectral Fit Display'
    
    # def __init__(self, root, filename):
    def __init__(self, data, gspec_root=None):
        super(FitPlotter, self).__init__()
        # Create the root window
        root = Tk()
        root.title(self.name)
        self.root = root

        # register this window with the root application window
        if gspec_root is not None:
            self.gspec_root = gspec_root
            self.gspec_root.register_window(self.name, self)
                
        # Create an instance variable to store the reference to the gspecManager 
        self.gspecManager = self.gspec_root.gspec

        self.data = data

        # Create the tkinter interface
        self.createTKinterInterface()
        

    ############################################

    def createTKinterInterface(self):

        # root = Tk()
        # root.title(self.name)
        self.width = self._geometry.width_fraction_to_npix(0.55)
        self.height = self._geometry.height_fraction_to_npix(0.75)
        #self.width  = 1024
        #self.height = 768

        # Set a soft maximum window size of 1024x550
        if self.width > 1024: self.width = 1024
        if self.height > 768: self.height = 768
        
        x_pos = self._geometry.width_fraction_to_npix(1.0)-self.width
        self.aspectratio = self.width/self.height
        self.root.geometry("{0}x{1}+{2}+0".format(self.width, self.height, x_pos))

        self.root.rowconfigure(0, weight=1)

        # Initialize a frame instance with the root as its parent
        Frame.__init__(self, self.root)

        self.config(bg=BGCOLOR)
        self.grid(sticky=N + W + E + S)

        # Set a minimum window size
        self.root.minsize(600, 500)

        # # Create the menu bar
        # menubar = Menu(root)

        # fileMenu = Menu(menubar)
        # fileMenu.add_cascade(label="Print Fit Info on Plot", command=self.command, state=DISABLED)
        # fileMenu.add_command(label="Screenshot", command=self.command)
        # fileMenu.add_command(label="Dismiss", command=self.on_window_close)

        # datasetMenu = Menu(menubar)
        # datasetMenu.add_command(label="All", command=self.command, state=DISABLED)
        # datasetMenu.add_command(label="<filename>", command=self.command, state=DISABLED)

        # optionsMenu = Menu(menubar)
        # optionsMenu.add_command(label="Colors", command=self.command, state=DISABLED)
        # optionsMenu.add_command(label="Plot Configuration", command=self.command, state=DISABLED)
        # optionsMenu.add_command(label="Set Fit Repeat Count", command=self.command, state=DISABLED)
        # optionsMenu.add_command(label="Set Fluence Energies", command=self.command, state=DISABLED)

        # menubar.add_cascade(label="File", menu=fileMenu)
        # menubar.add_cascade(label="Dataset", menu=datasetMenu)
        # menubar.add_cascade(label="Options", menu=optionsMenu)

        # root.config(menu=menubar)

        # reference the menubar
        self.menubar = copy.copy(self.gspec_root.menubar)
        self.menubar.update_root(self.root, self.name)
        
        # update the file menu
        items = ['Dismiss']
        commands = [self.on_window_close]
        self.filemenu = self.menubar.replace_menu('File', items, commands)
        
        # update the options menu
        items = ['Colors', 'Plot Configuration', 'Set Fit Repeat Count', 
                 'Set Fluence Energies']
        commands = [self.command, self.command, self.command, self.command] 
        self.options_menu = self.menubar.replace_menu('Options', items, commands)
        self.options_menu.entryconfig(0, state='disable')
        self.options_menu.entryconfig(1, state='disable')
        self.options_menu.entryconfig(2, state='disable')
        self.options_menu.entryconfig(3, state='disable')
        
        self.root.config(menu=self.menubar.top)


        ########### Buttons ###########

        frame_controls = Frame(self, bg=BGCOLOR, relief=FRAMERELIEF, pady=10, padx=0, height=100)

        # Adjust source selection menu
        self.spectrumSelectionTK = StringVar(frame_controls)
        self.spectrumSelectionTK.set('Fit Display Options')
        choices = ['Cumulative', 'Raw Counts', 'Counts Spectrum', 'Photon Spectrum', 'Energy Spectrum',
                   'Nu Fnu Spectrum']
        self.fitDisplayOptionseMenu = OptionMenu(frame_controls, self.spectrumSelectionTK, *choices,
                                                 command=self.onFitDisplaySelection)
        self.fitDisplayOptionseMenu.config(bg=BGCOLOR, width=16, justify=OM_JUSTIFY)

        # Fit results menu
        self.fitResultsSelection = StringVar(frame_controls)
        self.fitResultsSelection.set('Export Fit Results')
        choices = ['Write Results to File']#, 'Write Fit Params Only', 'Read Fit Results File']
        self.fitResultsMenu = OptionMenu(frame_controls, self.fitResultsSelection, *choices,
                                         command=self.commandSelection)
        self.fitResultsMenu.config(bg=BGCOLOR, width=16, justify=OM_JUSTIFY)
        #self.fitResultsMenu['state'] = 'disabled'

        # Redo fit button
        self.redoLastFitButton = Button(frame_controls, text=u"Redo Last Fit",
                                        highlightbackground=BGCOLOR, 
                                        command=self.gspec_root.prepare_spectral_fit, width=14)

        # Spectral fitting menu
        self.spectralFittingSelectionTK = StringVar(frame_controls)
        self.spectralFittingSelectionTK.set('Spectral Fitting')
        choices = ['Redo Last Fit', 'Synthesize Burst', 'Error Interval', 'ChiSqr 1D Plot', 'ChiSqr 2D Plot', 'Batch Fit Selections']
        self.spectralFittingMenu = OptionMenu(frame_controls, self.spectralFittingSelectionTK, *choices,
                                              command=self.spectralFittingSelection)
        self.spectralFittingMenu.config(bg=BGCOLOR, width=16, justify=OM_JUSTIFY)
        # self.spectralFittingMenu['state'] = 'disabled'

        self.spectralFittingMenu['menu'].entryconfigure(choices[1], state = "disabled")
        self.spectralFittingMenu['menu'].entryconfigure(choices[2], state = "disabled")
        self.spectralFittingMenu['menu'].entryconfigure(choices[3], state = "disabled")
        self.spectralFittingMenu['menu'].entryconfigure(choices[4], state = "disabled")
        # self.spectralFittingMenu['menu'].entryconfigure(choices[5], state = "disabled")

        # Batch fit options menu
        self.batchFitOptionsSelectionTK = StringVar(frame_controls)
        self.batchFitOptionsSelectionTK.set('Batch Fit Options')
        choices = ['Plot Model Parameters...', 'Residual Contours', 'NuFnu Contours', 'Stack Spectra']
        self.batchFitOptionsMenu = OptionMenu(frame_controls, self.batchFitOptionsSelectionTK, *choices,
                                              command=self.onBatchFitOptionsSelection)
        self.batchFitOptionsMenu.config(bg="#e1e1e1", width=16, justify=CENTER)

        self.batchFitOptionsMenu['menu'].entryconfigure(choices[1], state = "disabled")
        self.batchFitOptionsMenu['menu'].entryconfigure(choices[2], state = "disabled")
        self.batchFitOptionsMenu['menu'].entryconfigure(choices[3], state = "disabled")

        if self.data['batch'] == False:
            self.batchFitOptionsMenu['state'] = 'disabled'


        # Pack up the buttons
        self.fitDisplayOptionseMenu.pack(padx=0, side=TOP)        
        # self.redoLastFitButton.pack(padx=0, side=TOP)
        self.spectralFittingMenu.pack(padx=0, side=TOP)
        self.batchFitOptionsMenu.pack(padx=0, side=TOP)
        self.fitResultsMenu.pack(padx=0, side=TOP)

        ########### Check boxes ###########

        # Create the frame that will contain the toggle buttons and associated labels
        self.frame_toggleOptions = Frame(self, bg=BGCOLOR, relief=FRAMERELIEF, pady=0, padx=0)

        # # Create the count display options text label
        label_countDisplayOptions = Label(self.frame_toggleOptions, text="Spectral Data Options:", bg=BGCOLOR)
        label_countDisplayOptions.config(justify=OM_JUSTIFY)
        label_countDisplayOptions.grid(column=0, row=0, columnspan=2, sticky=N + W, pady=(8, 0), padx=5)

        # Create the count display options checkboxes
        self.showCountDisplayTK = BooleanVar(self.frame_toggleOptions)
        self.checkbutton_showCountDisplay = Radiobutton(self.frame_toggleOptions, bg=BGCOLOR, text='Show', value=True,
                                                        variable=self.showCountDisplayTK,
                                                        command=self.onCountDisplaySelection)
        self.checkbutton_hideCountDisplay = Radiobutton(self.frame_toggleOptions, bg=BGCOLOR, text='Hide', value=False,
                                                        variable=self.showCountDisplayTK,
                                                        command=self.onCountDisplaySelection)
        self.checkbutton_showCountDisplay.grid(column=0, row=1, sticky=W + N, padx=(15, 0))
        self.checkbutton_hideCountDisplay.grid(column=1, row=1, sticky=W + N)

        # # Create the model display options text label
        label_modelDisplayOptions = Label(self.frame_toggleOptions, text="Spectral Model Options:", bg=BGCOLOR)
        label_modelDisplayOptions.config(justify=OM_JUSTIFY)
        label_modelDisplayOptions.grid(column=0, row=2, columnspan=2, sticky=N + W, pady=(8, 0), padx=5)

        # Create the model display options checkboxes
        self.showModelDisplayTK = BooleanVar(self.frame_toggleOptions)
        self.checkbutton_showModelDisplay = Radiobutton(self.frame_toggleOptions, bg=BGCOLOR, text='Show', value=True,
                                                        variable=self.showModelDisplayTK,
                                                        command=self.onModelDisplaySelection)
        self.checkbutton_hideModelDisplay = Radiobutton(self.frame_toggleOptions, bg=BGCOLOR, text='Hide', value=False,
                                                        variable=self.showModelDisplayTK,
                                                        command=self.onModelDisplaySelection)

        self.checkbutton_showModelDisplay['state'] = 'disabled'
        self.checkbutton_hideModelDisplay['state'] = 'disabled'

        self.checkbutton_showModelDisplay.grid(column=0, row=3, sticky=W + N, padx=(15, 0))
        self.checkbutton_hideModelDisplay.grid(column=1, row=3, sticky=W + N)


        # # Create the spectral model components text label
        # label_spectralModelComponents = Label(self.frame_toggleOptions, text='Spectral Model Components:', bg=BGCOLOR)
        # label_spectralModelComponents.config(justify=OM_JUSTIFY)
        # label_spectralModelComponents.grid(column=0, row=4, columnspan=2, sticky=N + W, pady=(8, 0), padx=5)

        # # Create the spectral model components checkboxes
        # self.showSpectralModelComponentsTK = BooleanVar(self.frame_toggleOptions)
        # self.checkbutton_showSpectralModelComponents = Radiobutton(self.frame_toggleOptions, bg=BGCOLOR, text='Show',
        #                                                            value=True,
        #                                                            variable=self.showSpectralModelComponentsTK,
        #                                                            command=self.onModelComponentsDisplaySelection)
        # self.checkbutton_hideSpectralModelComponents = Radiobutton(self.frame_toggleOptions, bg=BGCOLOR, text='Hide',
        #                                                            value=False,
        #                                                            variable=self.showSpectralModelComponentsTK,
        #                                                            command=self.onModelComponentsDisplaySelection)

        # self.checkbutton_showSpectralModelComponents.grid(column=0, row=5, sticky=W + N, padx=(15, 0))
        # self.checkbutton_hideSpectralModelComponents.grid(column=1, row=5, sticky=W + N)


        # Residual Display Options  
        self.label_residualDisplayOptions = Label(self.frame_toggleOptions, text='Residual Display Options:', bg=BGCOLOR)
        self.label_residualDisplayOptions.config(justify=OM_JUSTIFY)
        self.label_residualDisplayOptions.grid(column=0, row=6, columnspan=2, sticky=N + W, pady=(8, 0), padx=5)

        self.residualSelectionTK = StringVar(self.frame_toggleOptions)
        self.checkbutton_SigmaResiduals = Radiobutton(self.frame_toggleOptions, bg=BGCOLOR, text='Sigma Residuals',
                                                      value='Sigma Residuals', variable=self.residualSelectionTK,
                                                      command=self.onResidualSelection, state = DISABLED)
        self.checkbutton_countResiduals = Radiobutton(self.frame_toggleOptions, bg=BGCOLOR, text='Count Residuals',
                                                      value='Count Residuals', variable=self.residualSelectionTK,
                                                      command=self.onResidualSelection)
        self.checkbutton_noResiduals = Radiobutton(self.frame_toggleOptions, bg=BGCOLOR, text='No Residuals',
                                                   value='No Residuals', variable=self.residualSelectionTK,
                                                   command=self.onResidualSelection)

        # self.checkbutton_SigmaResiduals['state'] = 'disabled'
        # self.checkbutton_SigmaResiduals.configure(state = DISABLED)

        self.checkbutton_SigmaResiduals.grid(column=0, row=7, columnspan=2, sticky=W + N, padx=(15, 0))
        self.checkbutton_countResiduals.grid(column=0, row=8, columnspan=2, sticky=W + N, padx=(15, 0))
        self.checkbutton_noResiduals.grid(column=0, row=9, columnspan=2, sticky=W + N, padx=(15, 0))


        # Legend Options
        label_legendDisplayOptions = Label(self.frame_toggleOptions, text="Legend Options:", bg=BGCOLOR)
        label_legendDisplayOptions.config(justify=OM_JUSTIFY)
        label_legendDisplayOptions.grid(column=0, row=10, columnspan=2, sticky=N + W, pady=(8, 0), padx=5)

        # Create the model display options checkboxes
        self.showLegendDisplayTK = BooleanVar(self.frame_toggleOptions)
        self.checkbutton_showLegendDisplay = Radiobutton(self.frame_toggleOptions, bg=BGCOLOR, text='Show', value=True,
                                                        variable=self.showLegendDisplayTK,
                                                        command=self.onLegendDisplaySelection)
        self.checkbutton_hideLegendDisplay = Radiobutton(self.frame_toggleOptions, bg=BGCOLOR, text='Hide', value=False,
                                                        variable=self.showLegendDisplayTK,
                                                        command=self.onLegendDisplaySelection)
        self.checkbutton_showLegendDisplay.grid(column=0, row=11, sticky=W + N, padx=(15, 0))
        self.checkbutton_hideLegendDisplay.grid(column=1, row=11, sticky=W + N)



       # Create the frame that will contain the toggle buttons and associated labels
        # self.frame_toggleOption_lower = Frame(self, bg=BGCOLOR, relief=FRAMERELIEF, pady=0, padx=0)

        label_scaleDisplayOptions = Label(self.frame_toggleOptions, text="Plot Scale Options:", bg=BGCOLOR)
        label_scaleDisplayOptions.config(justify=OM_JUSTIFY)
        label_scaleDisplayOptions.grid(column=0, row=12, columnspan=2, sticky=N + W, pady=(8, 0), padx=5)

        # self.xlog = BooleanVar()     
        self.checkbutton_xlog = Checkbutton(self.frame_toggleOptions, bg=BGCOLOR, text='X Log', command=self.setXLogScale)
        self.checkbutton_xlog.grid(column=0, row=13, sticky=W + N, pady=(8, 0), padx=(15, 0))

        # self.ylog = BooleanVar()     
        self.checkbutton_ylog = Checkbutton(self.frame_toggleOptions, bg=BGCOLOR, text='Y Log', command=self.setYLogScale)
        self.checkbutton_ylog.grid(column=1, row=13, sticky=W + N, pady=(8, 0))

        # Set the default values
        self.showCountDisplay = True
        self.countDataVisible = True
        self.showModelDisplay = True
        self.showSpectralModelComponents = True
        self.modelDisplayVisible = True
        self.showLegendDisplay = True
        self.residualSelection = 'No Residuals'
        self.spectrumSelectionTK.set('Counts Spectrum')
        self.spectrumSelection = 'Counts Spectrum'
        self.checkbutton_showCountDisplay.select()
        self.checkbutton_showModelDisplay.select()
        # self.checkbutton_hideSpectralModelComponents.select()
        self.checkbutton_showLegendDisplay.select()
        self.checkbutton_noResiduals.select()


        if self.data is not None and self.data['batch'] == False:
            self.checkbutton_xlog.select()
            self.checkbutton_ylog.select()
            self.xlog = True
            self.ylog = True
        else:
            self.xlog = False
            self.ylog = False

            # else:
            # self.checkbutton_xlog.deselect()
            # self.checkbutton_ylog.deselect()


        ########### Filename Readout ###########

        # # Create the frame that will contain the toggle buttons and associated labels
        # frame_filaname  = Frame(self, bg=BGCOLOR, relief=FRAMERELIEF, pady=0, padx=0)
        # label_filename = Label(frame_filaname, text='File:', bg=BGCOLOR)
        # label_filename.config(justify=LEFT)
        # label_filename.grid(column=0, row=1, columnspan=3, sticky=W, pady=(8,0), padx=5) 

        ########### Plot Window ###########

        # Create a tk frame that will contain the plotting window
        self.plot_frame = Frame(self.root, bg=BGCOLOR)

        # Create a figure instance and attach it to self for future reference
        self.figure = Figure(figsize=(5, 4), dpi=100)

        # # Add the subplot axis to self for future reference
        self.ax = self.figure.add_subplot(111)
        self.ax2 = None

        # Set the format of the coordinate readout
        self.ax.format_coord = lambda x, y: ""

        # Set the x minor tick frequency
        minorLocator = AutoMinorLocator()
        self.ax.xaxis.set_minor_locator(minorLocator)

        # Set the y minor tick frequency
        minorLocator = AutoMinorLocator()
        self.ax.yaxis.set_minor_locator(minorLocator)

        # With matplotlib toolbar
        self.canvas = FigureCanvasTkAgg(self.figure, self.plot_frame)
        self.canvas.show()
        self.canvas.callbacks.connect('button_press_event', self.onMouseClick)
        self.canvas.callbacks.connect('motion_notify_event', self.onMouseMove)

        # Add the matplotlib toolbar
        toolbar_frame = Frame(self.plot_frame, bg=BGCOLOR)
        toolbar_frame.grid(column=0, row=1, sticky=N + W)
        toolbar_frame.config(bg=BGCOLOR)
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, toolbar_frame)
        self.toolbar.config(bg=BGCOLOR)

        # Change the message label and button background colors
        self.toolbar._message_label.config(background=BGCOLOR)
        for button in self.toolbar.winfo_children():
            button.config(background=BGCOLOR)

        self.toolbar.update()
        self.canvas.get_tk_widget().grid(column=0, row=0, padx=0, pady=0, sticky=N + W + S + E)

        self.plot_frame.rowconfigure(0, weight=1)
        # self.plot_frame.rowconfigure(1, weight=0)
        self.plot_frame.columnconfigure(0, weight=1)
        # self.plot_frame.columnconfigure(1, weight=0)

        
        ########### Pack it all up ###########

        frame_controls.grid(column=0, row=0, sticky=N)
        self.frame_toggleOptions.grid(column=0, row=1, sticky=W + N)
        # self.frame_toggleOption_lower.grid(column=0, row=2, sticky=W + N)
        self.plot_frame.grid(column=1, row=0, padx=0, pady=0, sticky=N + W + S + E, rowspan=2)

        # frame_countDisplayOptions.grid(column=0, row=1, sticky=N)
        # frame_spectralModelComponents.grid(column=0, row=2, sticky=N)
        # frame_modelColorOptions.grid(column=0, row=3, sticky=N)

        # Configuration needed with matplotlib toolbar
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=1)
        self.rowconfigure(2, weight=1)
        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)

        self.root.rowconfigure(0, weight=1)
        self.root.columnconfigure(0, weight=0)
        self.root.columnconfigure(1, weight=2)
        
        # Close the window
        self.root.protocol("WM_DELETE_WINDOW", self.on_window_close)

        if self.data is not None and self.data['batch'] == False:
            self.plotSpectrum()
        else:
            self.plotBatchFitResults()

    ############################################

    def onFitDisplaySelection(self, menuSelection):

        # Convert the tkinter variable to a proper string
        self.spectrumSelection = menuSelection

        # Redraw the current spectrum
        self.plotSpectrum()

    ############################################

    def onResidualSelection(self):

        # Convert the tkinter variable to a proper string
        self.residualSelection = self.residualSelectionTK.get()

        # Redraw the current spectrum
        self.plotSpectrum()

    ############################################

    def onCountDisplaySelection(self):

        self.showCountDisplay = self.showCountDisplayTK.get()

        data_successfully_removed = False
        data_successfully_added = False

        if self.showCountDisplay == False:
            for data_container in self.data_containers:

                try:

                    # Hide the line collections only if they are already visible
                    if self.countDataVisible == True:

                        # Hide the upper limit cap
                        data_container[2][0].set_visible(False)

                        # Hide the upper limit bar
                        data_container[2][1].set_visible(False)
 
                        # Hide the upper limit arrow
                        if len(data_container[1]) != 0:
                            data_container[1][0].set_visible(False)
                        
                        data_successfully_removed = True
                except:
                    pass


            # Toggle the variable indicating whether count data is currently displayed
            if data_successfully_removed == True and self.countDataVisible == True:
                self.countDataVisible = False

        else:

            for data_container in self.data_containers:
            
                try:
                    # Add the line collections only if they are already hidden
                    if self.countDataVisible == False:

                        # Show the upper limit cap
                        data_container[2][0].set_visible(True)

                        # Show the upper limit bar
                        data_container[2][1].set_visible(True)

                        # Show the upper limit arrow
                        if len(data_container[1]) != 0:
                            data_container[1][0].set_visible(True)

                        data_successfully_added = True
                except:
                    pass

            # Toggle the variable indicating whether count data is currently displayed
            if data_successfully_added == True and self.countDataVisible == False:
                self.countDataVisible = True

        # Redraw the canvas
        self.canvas.draw()
            

    ############################################

    def onModelDisplaySelection(self):

        self.showModelDisplay = self.showModelDisplayTK.get()

        if self.showModelDisplay == False:
            for model_artist in self.model_artists:

                # Only hide it if it's already visible
                if model_artist.get_visible() == True:
                    model_artist.set_visible(False)
        else:
            for model_artist in self.model_artists:

                # Only show it if it's already hidden
                if model_artist.get_visible() == False:
                    model_artist.set_visible(True)

        self.canvas.draw()

        # Redraw the current spectrum
        # self.plotSpectrum()

    ############################################

    def onLegendDisplaySelection(self):

        self.showLegendDisplay = self.showLegendDisplayTK.get()

        if self.showLegendDisplay == False:

            # Only hide it if it's already visible
            if self.legend.get_visible() == True:
                self.legend.set_visible(False)

        else:

            # Only show it if it's already hidden
            if self.legend.get_visible() == False:
                self.legend.set_visible(True)


        self.canvas.draw()


        # Redraw the current spectrum
        # self.plotSpectrum()


    ############################################

    def onModelComponentsDisplaySelection(self):

        self.showSpectralModelComponents = self.showSpectralModelComponentsTK.get()

        print("Function not yet implemented")
        # Redraw the current spectrum
        # self.plotSpectrum(self.spectrumSelection)

    ############################################

    def onMouseClick(self, event):
        # print("Mouse button clicked")
        pass

    ############################################

    def onMouseMove(self, event):
        pass

    ############################################

    def command(self):
        print('function not yet implmented')

    ############################################

    def commandSelection(self, selection):
        if 'Write Results to File' in selection:
            initialfile = self.gspec_root.scat.filename_obj.basename()
            filename = filedialog.asksaveasfilename(initialfile=initialfile,  
                                                title="Save Fit Results",
                                                filetypes=(('SCAT File', '*.fit'),),
                                                defaultextension='.fit')
            if len(filename) == 0:
                return

            else:
                self.gspec_root.scat.write(os.path.dirname(filename), 
                                           filename=os.path.basename(filename))
        else:
            print('function not yet implmented')


    ############################################

    def spectralFittingSelection(self, selection):
        self.spectralFittingSelectionTK.set('Spectral Fitting')
        if 'Redo Last Fit' in selection:

            # Re-enable the fit display options menu
            self.spectrumSelectionTK.set('Cumulative')
            self.fitDisplayOptionseMenu['state'] = 'normal'

            # self.fitDisplayOptionseMenu.configure(state='normal')
            # menu = self.fitDisplayOptionseMenu["menu"]
            # last = menu.index("end")
            # for i in range(last+1):
            #     menu.entryconfigure(i, state="normal")

            # Disable the batch fit menu
            self.batchFitOptionsMenu['state'] = 'disabled'

            self.gspec_root.prepare_spectral_fit()

        if 'Batch Fit Selections' in selection:

            # Re-enable the fit display options menu
            self.spectrumSelectionTK.set('Fit Display Options')
            self.fitDisplayOptionseMenu['state'] = 'disabled'

            self.gspec_root.prepare_spectral_fit(batch=True)

    ############################################

    def setXLogScale(self):

        if self.xlog == True:
            self.xlog = False
        else:
            self.xlog = True

        # Redraw the current spectrum
        # self.plotSpectrum()

        if data is not None and self.data['batch'] == False:
            self.plotSpectrum()
        else:
            self.plotBatchFitResults()

    ############################################

    def setYLogScale(self):

        if self.ylog == True:
            self.ylog = False
        else:
            self.ylog = True

        # Redraw the current spectrum
        # self.plotSpectrum()

        if data is not None and self.data['batch'] == False:
            self.plotSpectrum()
        else:
            self.plotBatchFitResults()

    ############################################

    def plotSpectrum(self):

        # Show the widgets that may have been hidden during batch fitting
        self.fitDisplayOptionseMenu.configure(state='normal')
        
        # Re-enable all buttons that may have been disabled during batch fitting
        for child in self.frame_toggleOptions.winfo_children():

            # Make sure the sigma residuals are disabled if the statistic is anything other than chisqr
            if 'Sigma Residuals' in child['text']:
                if 'Chi-Squared' in self.data['stat_name']:
                    child.configure(state='normal')
                else:
                    # Deselect the sigma residuals if they were selected during a previous fit
                    if 'normal' in child['state']:
                        self.residualSelection = 'No Residuals'
                        self.residualSelectionTK.set('No Residuals')
        
                    # Disable the sigma residual option 
                    child.configure(state='disabled')
            else:
                child.configure(state='normal')


        # Get the number of datagroups
        numberOfDataGroups = len(self.data['x_ufspec'])

        # Extract the data names to be used as labels
        data_labels = list(self.gspecManager.data.keys())

        # Define a list of plotting colors to cycle throug
        colors = ['#394264', '#6E8846', '#7B3F5B', '#aa611d', '#463965', '#93894B', '#396F45', '#93524B', '#553462', '#8F924A', '#2F545B', '#946F4B']
        # colors_model = ['#8b0000', '#720000', '#580000', '#3f0000', '#a50000', '#be0000', '#d80000', '#9f0000', '#c60000', '#d90000']
        colors_model = ['#8b0000', '#8b0000', '#8b0000', '#8b0000', '#8b0000', '#8b0000', '#8b0000', '#8b0000', '#8b0000', '#8b0000']

        # Shuffle the colors
        # random.shuffle(colors)

        # Clear any existing plotting data
        self.ax.clear()
        if self.ax2 is not None:
            self.ax2.clear()

        # Make sure autoscaling is turned on
        self.ax.set_autoscale_on(True)

        # Do not display any residuals
        if 'No Residuals' in self.residualSelection:
            if self.ax2 is not None:
                self.figure.delaxes(self.ax)
                self.figure.delaxes(self.ax2)
                self.ax = self.figure.add_subplot(111)

                self.ax2 = None

        # Create a second axes to display fit residuals
        if 'No Residuals' not in self.residualSelection:

            if self.ax2 is None:
                gs = gridspec.GridSpec(4, 1)
                self.ax.set_position(gs[0:3].get_position(self.figure))
                self.ax.set_subplotspec(gs[0:3])
                self.ax.xaxis.set_ticklabels([])

                self.ax2 = self.figure.add_subplot(gs[3])

                # self.ax2.set_ylabel('Count Rate')
                self.ax2.set_xlabel('Time (s)', fontsize=PLOTFONTSIZE)

                # Set the x minor tick frequency
                minorLocator = AutoMinorLocator()
                self.ax2.xaxis.set_minor_locator(minorLocator)

                # Set the y minor tick frequency
                minorLocator = AutoMinorLocator()
                self.ax2.yaxis.set_minor_locator(minorLocator)

                # Link the axes
                self.ax.get_shared_x_axes().join(self.ax, self.ax2)

                # Create zero space between the two plots
                self.figure.subplots_adjust(hspace=0)

                # Hide the axis labels
                self.ax.tick_params(labelbottom='off')

        
        # Get the energy bins
        energy_data = self.data['x_ufspec']
        energy_error_data = self.data['xerr_ufspec']
        # energy_summed_model = self.data['energies'][0:-1]  # this bin selection needs to be re-examined
        energy_summed_model = self.data['energies_model']

        # Get the counts
        counts = self.data['counts']
        counts_error = self.data['counts_error']
        counts_model = self.data['counts_model']

        # Get the integrated counts
        counts_integrated = self.data['counts_integrated']
        counts_integrated_model = self.data['counts_integrated_model']

        # Get the photon flux values
        flux_data = self.data['y_ufspec']
        flux_error_data = self.data['yerr_ufspec']

        # Get the model photon flux
        flux_model = self.data['y_ufspec_model']

        # Get the summer model photon flux
        flux_summed_model = self.data['ufspec_model']

        # Get the sigma residuals
        residuals_sigma = self.data['residuals_sigma']

        # Set the x-axis range to match the data across all data groups
        energy_data_flattened = numpy.hstack(energy_data)
        xmin = numpy.floor(numpy.min(energy_data_flattened))
        xmax = numpy.ceil(numpy.max(energy_data_flattened))


        # Set the data alpha
        if self.showCountDisplay == True:
            alpha_data = 1
        else:
            alpha_data = 0

        # Set the model alpha
        if self.showModelDisplay == True:
            alpha_model = 1
        else:
            alpha_model = 0

        # Set the model alpha
        if self.showLegendDisplay == True:
            alpha_legend = 1
        else:
            alpha_legend = 0


        # Determine what kind of data to display
        if self.spectrumSelection == 'Cumulative':
            flux_data = counts_integrated
            flux_error_data = []
            for counts_integrated_group in counts_integrated:
                flux_error_data.append(numpy.zeros(len(counts_integrated_group)))
            flux_error_data = numpy.array(flux_error_data)
            flux_model = counts_integrated_model
            residuals_data = flux_data - flux_model

            self.ax.set_ylabel(r'Normalized integrated counts', fontsize=PLOTFONTSIZE)

        if self.spectrumSelection == 'Raw Counts':
            flux_data = counts
            flux_error_data = counts_error
            flux_model = counts_model
            residuals_data = flux_data - flux_model

            self.ax.set_ylabel(r'Counts bin$^{-1}$', fontsize=PLOTFONTSIZE) 

        # Prepare the count spectrum
        if self.spectrumSelection == 'Counts Spectrum':
            flux_data = counts/energy_error_data
            flux_error_data = (counts_error / counts) * flux_data
            flux_model = counts_model/energy_error_data
            residuals_data = flux_data - flux_model

            self.ax.set_ylabel(r'Counts s$^{-1}$ keV$^{-1}$', fontsize=PLOTFONTSIZE)

        # Prepare the flux density spectrum (Photons cm-2 s-1 keV-1)
        if self.spectrumSelection == 'Photon Spectrum':
            residuals_data = flux_data - flux_model

            self.ax.set_ylabel(r'Flux (Photons cm$^{-2}$ s$^{-1}$ keV$^{-1}$)', fontsize=PLOTFONTSIZE)

        # Prepare the Fv spectrum (Photons cm-2 s-1)
        if self.spectrumSelection == 'Energy Spectrum':
            flux_error_data = (flux_error_data / flux_data) * (flux_data * energy_data)  # photons cm-2
            flux_data = flux_data * energy_data
            flux_model = flux_model * energy_data
            flux_summed_model = flux_summed_model * energy_summed_model
            residuals_data = flux_data - flux_model

            self.ax.set_ylabel(r'Energy (Photons cm$^{-2}$ s$^{-1}$)', fontsize=PLOTFONTSIZE)

        # Prepare the vFv spectrum (Photons keV cm-2 s-1)
        if self.spectrumSelection == 'Nu Fnu Spectrum':
            flux_error_data = (flux_error_data / flux_data) * (flux_data * (energy_data ** 2))  # photons keV
            flux_data = flux_data * (energy_data ** 2)
            flux_model = flux_model * (energy_data ** 2)
            flux_summed_model = flux_summed_model * (energy_summed_model ** 2)
            residuals_data = flux_data - flux_model

            self.ax.set_ylabel(r'$\nu$ F$_{\nu}$ (Photons keV cm$^{-2}$ s$^{-1}$)', fontsize=PLOTFONTSIZE)



        # Determine which data points have error bars that go negative
        index_detections = []
        index_upperLimits = []
        for dataGroup in range(numberOfDataGroups):
            index_detections.append(numpy.where((flux_data[dataGroup]-flux_error_data[dataGroup]) > 0))
            index_upperLimits.append(numpy.where((flux_data[dataGroup]-flux_error_data[dataGroup]) <= 0))

        # Create a list to store data artists
        self.data_containers = []

        # Plot the data points
        for dataGroup, detections, upperLimits in zip(range(numberOfDataGroups), index_detections, index_upperLimits):

            index_random = random.random()
            self.data_detections_container = self.ax.errorbar(energy_data[dataGroup][detections], flux_data[dataGroup][detections], xerr=energy_error_data[dataGroup][detections], capsize=0,
                             fmt='none', zorder=1, ecolor=colors[dataGroup % len(colors)], alpha=alpha_data, yerr=flux_error_data[dataGroup][detections], label=data_labels[dataGroup])

            self.data_limits_container = self.ax.errorbar(energy_data[dataGroup][upperLimits], flux_data[dataGroup][upperLimits], xerr=energy_error_data[dataGroup][upperLimits], capsize=0,
                             fmt='none', zorder=1, ecolor=colors[dataGroup % len(colors)], alpha=alpha_data, yerr=flux_data[dataGroup][upperLimits]/2.0, uplims=True)

            # Add the artists to the list of data artists
            self.data_containers.append(self.data_detections_container)
            self.data_containers.append(self.data_limits_container)

        # Plot the residuals
        for dataGroup in range(numberOfDataGroups):

            if 'Count Residuals' in self.residualSelection and self.ax2 is not None:
                self.ax2.plot([xmin, xmax], [0,0], linestyle='-', linewidth=0.8, c=colors_model[dataGroup % len(colors)], zorder=2, alpha=alpha_model)
                self.ax2.errorbar(energy_data[dataGroup], residuals_data[dataGroup], xerr=energy_error_data[dataGroup],
                                  capsize=0, fmt='none', ecolor=colors[dataGroup % len(colors)], zorder=1, alpha=alpha_data, yerr=flux_error_data[dataGroup])

            if 'Sigma Residuals' in self.residualSelection and self.ax2 is not None:
                self.ax2.plot([xmin, xmax], [0,0], linestyle='-', linewidth=0.8, c=colors_model[dataGroup % len(colors)], zorder=2, alpha=alpha_model)
                self.ax2.errorbar(energy_data[dataGroup], residuals_sigma[dataGroup], xerr=energy_error_data[dataGroup],
                                  capsize=0, fmt='none', ecolor=colors[dataGroup % len(colors)], zorder=1, alpha=alpha_data, yerr=numpy.ones(len(residuals_sigma[dataGroup])))

                self.ax2.set_ylabel('Sigma', fontsize=PLOTFONTSIZE)
                
        # Set the x-axis labels
        self.ax.set_xlabel('Energy (keV)', fontsize=PLOTFONTSIZE) 

        # Set the tick label size
        self.ax.xaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        self.ax.yaxis.set_tick_params(labelsize=PLOTFONTSIZE)

        # Set the x-axis minor tick frequency
        minorLocator = AutoMinorLocator()
        self.ax.xaxis.set_minor_locator(minorLocator)

        # Set the y-axis minor tick frequency
        minorLocator = AutoMinorLocator()
        self.ax.yaxis.set_minor_locator(minorLocator)

        # Add a bit of buffer between the data max and min and the plot edges
        self.ax.set_xlim(xmin/2, xmax*2)

        # Set the plot scale
        if self.xlog == True:
            self.ax.set_xscale('log')

            if self.ax2 is not None:
                self.ax2.set_xscale('log')

        if self.ylog == True:
            self.ax.set_yscale('log')

        # Set the y-axis range to match the data across all data groups
        flux_data_model_flattened = numpy.hstack([numpy.hstack(flux_data), numpy.hstack(flux_model)])
        ymin = numpy.floor(numpy.min(flux_data_model_flattened))
        ymax = numpy.ceil(numpy.max(flux_data_model_flattened))

        # Add a bit of buffer between the data/model max and min and the plot edges
        if ymin == 0:
            ymin = numpy.min(flux_data_model_flattened)
            print(ymin)
        try:
            self.ax.set_ylim(ymin/2, ymax*2)
        except:
            pass

        # # Set a minimum y value
        # dflux = flux_data-flux_error_data
        # for dataGroup in range(numberOfDataGroups):
        #     if ((len(numpy.where(dflux[dataGroup] <= 0)[0]) > 0) == True) and self.spectrumSelection == 'Counts Spectrum':
        #         self.ax.set_ylim(bottom=0.1)

        # Don't let the model affect the plot scale
        self.ax.set_autoscale_on(False)

        # Create a list to store model artists
        self.model_artists = []

        # Plot the individual detector model
        if self.spectrumSelection == 'Raw Counts' or self.spectrumSelection == 'Counts Spectrum' or self.spectrumSelection == 'Cumulative':

            for dataGroup in range(numberOfDataGroups):
                for i in range(len(energy_data[dataGroup])):

                    # Plot the individual model values for the count and cumulutive spectrum
                    energy_low = energy_data[dataGroup][i] - energy_error_data[dataGroup][i]
                    energy_high = energy_data[dataGroup][i] + energy_error_data[dataGroup][i]

                    self.model_artist, =  self.ax.plot([energy_low, energy_high], [flux_model[dataGroup][i], flux_model[dataGroup][i]],
                                 c=colors_model[dataGroup % len(colors_model)], zorder=2, linewidth=0.8, alpha=alpha_model)

                    self.model_artists.append(self.model_artist)
                    try:
                        self.model_artist, =  self.ax.plot([energy_high, energy_high],
                                     [flux_model[dataGroup][i], flux_model[dataGroup][i + 1]],
                                     c=colors_model[dataGroup % len(colors_model)], zorder=2, linewidth=0.8, alpha=alpha_model)

                        self.model_artists.append(self.model_artist)
                    except:
                        pass

            # if 'Count Residuals' in self.residualSelection and self.ax2 is not None: 
            #     self.ax2.plot(energy_data[dataGroup], numpy.zeros(len(energy_data[dataGroup])), linestyle='--', c=colors_model[dataGroup], zorder=2)
            #     self.ax2.errorbar(energy_data[dataGroup], flux_data[dataGroup]-flux_model[dataGroup], xerr=energy_error_data[dataGroup], capsize=0, fmt='none', ecolor=colors[dataGroup], zorder=1)
            #     self.ax2.errorbar(energy_data[dataGroup], flux_data[dataGroup]-flux_model[dataGroup], yerr=flux_error_data[dataGroup], capsize=0, fmt='none', ecolor=colors[dataGroup], zorder=1)

            # if 'Sigma Residuals' in self.residualSelection and self.ax2 is not None: 
            #     self.ax2.plot(energy_data[dataGroup], numpy.zeros(len(energy_data)), linestyle='--', c=colors_model[dataGroup], zorder=2)
            #     self.ax2.errorbar(energy_data[dataGroup], flux_data[dataGroup]-flux_model[dataGroup], xerr=energy_error_data[dataGroup], capsize=0, fmt='none', ecolor=colors[dataGroup], zorder=1)
            #     self.ax2.errorbar(energy_data[dataGroup], flux_data[dataGroup]-flux_model[dataGroup], yerr=flux_error_data[dataGroup], capsize=0, fmt='none', ecolor=colors[dataGroup], zorder=1)
            #     self.ax2.set_ylabel('Sigma')


        else:

            # Plot the summed model
            self.model_artist, = self.ax.plot(energy_summed_model, flux_summed_model, c=colors_model[dataGroup % len(colors_model)], zorder=2, alpha=alpha_model)
            # self.ax.plot(energy_summed_model, flux_summed_model, c='darkred', zorder=2, alpha=1)

            self.model_artists.append(self.model_artist)


        # Hide the last tick label on the second axes, if it exists
        if self.ax2 is not None:
            yticks = self.ax2.get_yticks()
            yticks = yticks[1:-1]
            self.ax2.set_yticks(yticks)

            # Set the tick label size
            self.ax2.xaxis.set_tick_params(labelsize=PLOTFONTSIZE)
            self.ax2.yaxis.set_tick_params(labelsize=PLOTFONTSIZE)

        # Add a legend
        if self.showLegendDisplay == True:
            if self.spectrumSelection == 'Cumulative':
                legend_location = 'lower right'
            if self.spectrumSelection == 'Raw Counts' or self.spectrumSelection == 'Counts Spectrum' or self.spectrumSelection == 'Energy Spectrum':
                legend_location = 'lower left'
            if self.spectrumSelection == 'Photon Spectrum':
                legend_location = 'upper right'
            if self.spectrumSelection == 'Nu Fnu Spectrum':
                legend_location = 'upper left'

            self.legend = self.ax.legend(numpoints=1, scatterpoints=1, fontsize='x-small', frameon=False, loc=legend_location)
            # legend.get_frame().set_linewidth(0.5)


        self.canvas.draw()


    ############################################

    def plotBatchFitResults(self, parameter_index=(0,)):

        # Enable the batch fit menu
        self.batchFitOptionsMenu['state'] = 'normal'

        # Turn off log scaling by default
        if self.xlog == False:
            self.checkbutton_xlog.deselect()

        if self.ylog == False:
            self.checkbutton_ylog.deselect()

        # Clear any existing plotting data
        self.ax.clear()

        # Remove the residual plot if it exists
        if self.ax2 is not None:
            self.figure.delaxes(self.ax)
            self.figure.delaxes(self.ax2)
            self.ax = self.figure.add_subplot(111)
            self.ax2 = None

        
        # Hide the widgets that we don't need for batch fitting
        # self.fitDisplayOptionseMenu.configure(state='disable')

        for child in self.frame_toggleOptions.winfo_children():
            text = child.cget("text")
            if 'X Log' not in text and 'Y Log' not in text:
                child.configure(state='disable')


        # Extract the fit results
        model = self.data['model'][0]
        parameter_names = self.data['parameter_names']
        parameter_values = self.data['parameter_values']
        parameter_sigmas = self.data['parameter_sigmas']

        if len(parameter_index) == 1:
            
            # Extract the extension information
            time_start = self.data['time_start']
            time_end = self.data['time_end']
            time_exposure = self.data['exposure']

            time = time_start + (time_end-time_start)/2.0
            dtime = (time_end-time_start)/2.0

            x_name = 'Time (sec)'
            x_value = time
            x_sigma = dtime

            y_name = parameter_names[parameter_index[0]]
            y_value = parameter_values[:,parameter_index[0]]
            y_sigma = parameter_sigmas[:,parameter_index[0]]

        else:

            x_name = "%s - %s" % (model, parameter_names[parameter_index[0]])
            x_value = parameter_values[:,parameter_index[0]]
            x_sigma = parameter_sigmas[:,parameter_index[0]]


            y_name = parameter_names[parameter_index[1]]
            y_value = parameter_values[:,parameter_index[1]]
            y_sigma = parameter_sigmas[:,parameter_index[1]]


        # Plot the parameter history
        self.ax.errorbar(x_value, y_value, xerr=x_sigma, yerr=y_sigma, capsize=0, fmt='none', zorder=1, ecolor='#394264', alpha=1)

        # Set the x-axis labels
        self.ax.set_xlabel(x_name, fontsize=PLOTFONTSIZE)

        # Set the x-axis labels
        self.ax.set_ylabel("%s - %s" % (model, y_name), fontsize=PLOTFONTSIZE)


        # Set the x-axis minor tick frequency
        minorLocator = AutoMinorLocator()
        self.ax.xaxis.set_minor_locator(minorLocator)

        # Set the y-axis minor tick frequency
        minorLocator = AutoMinorLocator()
        self.ax.yaxis.set_minor_locator(minorLocator)

        # Set the plot scale
        if self.xlog == True:
            self.ax.set_xscale('log')

            if self.ax2 is not None:
                self.ax2.set_xscale('log')

        if self.ylog == True:
            self.ax.set_yscale('log')

        # Set the tick label size
        self.ax.xaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        self.ax.yaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        
        self.canvas.draw()



    ############################################

    def onBatchFitOptionsSelection(self, menuSelection):

        if 'Plot Model Parameters...' in menuSelection:

            # Get the model associated with the current data
            selected_model = self.data['model']

            # Reset the menu selection
            self.batchFitOptionsSelectionTK.set('Batch Fit Options') 

            # Display the model parameter selection dialog
            dialog = PlotModelParametersDialog(self.root, selected_model, self.plotBatchFitResults)

        else:

            print("function not yet implemented")

            # Reset the menu selection
            self.batchFitOptionsSelectionTK.set('Batch Fit Options') 


    ############################################

    def quitConfirmation(self):
        exit()

    ############################################

    # def on_window_close(self):
    #     try:
    #         del self.gspecManager._open_windows['Spectral Fit Display']
    #     except:
    #         pass
    #     self.root.destroy()

    def on_window_close(self):
        if self.gspec_root is not None:
            self.gspec_root.unregister_window(self.name)
        self.root.destroy()   


