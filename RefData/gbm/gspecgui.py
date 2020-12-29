#! /usr/bin/env python3

from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
from tkinter import ttk

import matplotlib

from .xspec import xspec

matplotlib.rcParams.update({'font.size': 10})
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import numpy
from gbm import data
import os
from matplotlib.ticker import AutoMinorLocator
from .lookup_old import RmfitLookup, GspecLookup
import matplotlib.gridspec as gridspec

import resource
import time

import astropy.io.fits as fits


##########################################################################################

class Splash(Toplevel):
    def __init__(self, parent):

        Toplevel.__init__(self, parent)

        xsize = 300
        ysize = 435

        # Center the root box
        width = self.winfo_screenwidth()
        height = self.winfo_screenheight()
        size = [250, 400]
        x = width / 2 - size[0] / 2
        y = height / 2 - size[1] / 2

        # Set the window size
        self.geometry("%sx%s+%s+%s" % (xsize, ysize, int(x), int(y)))
        self.config(bg="#e1e1e1")

        # Remove the window border
        self.overrideredirect(True)
        self.wm_attributes("-topmost", True)
        
        canvas = Canvas(self, height=90, highlightthickness=0, bg="#e1e1e1")   
        canvas.pack(pady=(20,0))                 
        self.image = PhotoImage(file="spaceship_75x75.gif")      
        canvas.create_image(xsize/2 - 9, 50, anchor=CENTER, image=self.image)  

        # Setup the text
        lines = ['\nGSpec',
                         'GBM Spectral Analysis Package', 
                         'Version 0.4\n',
                         'Adam Goldstein',
                         'Daniel Kocevski',
                         'Rob Preece',
                         'William Cleveland\n',
                         # 'Courtesy of:',
                         'Unversities Space Research Association', 
                         'NASA Marshall Space Flight Center', 
                         'University of Alabama in Huntsville', '']

        # Loop through each line and add it to the widget
        for line in lines:
            label = Label(self, text=line, bg="#e1e1e1", pady=2)
            label.pack() 

        # Bind the any key or mouse click to the destroy function
        self.bind('<Key>', self.dismiss)
        self.bind('<Button-1>', self.dismiss)

        ## required to make window show before the program gets to the mainloop
        self.update()


    def dismiss(self, event):

        # Close the window
        self.destroy()


##########################################################################################

class Application(Frame):
    """Primary GSPEC application class"""

    def __init__(self, gspecManager):

        # Create an instance variable to store the reference to the gspecManager 
        self.gspecManager = gspecManager

        # Create the tkinter interface
        self.createTKinterInterface()

    ############################################

    def createTKinterInterface(self):

        # Create the root window
        root = Tk()
        root.title('GSPEC v0.4')
        root.geometry("650x250+0+0")  # +300+300")
        root.minsize(650, 250)
        root.rowconfigure(0, weight=1)
        root.columnconfigure(0, weight=1)
        # root.config(bg="#e1e1e1")

        # Hide the window while the splash screen is up
        # root.withdraw()

        # # Indicate that the splash screen is currently closed
        # self.splashScreenOpen = False 

        # # Create a splash screen 
        # self.splash = Splash(root)

        # # Close the splash screen after 5 seconds 
        # self.splash.after(5000, lambda: self.splash.destroy())

        # Show the window
        # root.deiconify()

        # Make a reference to the root window so that we can kill it later
        self.root = root

        # Override the frame inilization so as to pass the root window
        Frame.__init__(self, root)

        # Configure the background
        self.config(bg="#e1e1e1")

        # Tell the frame to expand and fill the root window
        self.grid(sticky=N + W + E + S)
        
        self.about = AboutGspec(parent=self.root)
        
        # Create the menu bar
        menubar = Menu(root)

        filemenu = Menu(menubar, tearoff=0)
        filemenu.add_command(label="Load Data...", command=self.loadFile)
        filemenu.add_command(label="Load From Lookup...", command=self.loadLookup) 
        filemenu.add_command(label="Download Data...", command=self.command)
        filemenu.add_command(label="About...", command=self.about.show)
        filemenu.add_command(label="Quit", command=self.quitConfirmation)

        self.tasksmenu = Menu(root, tearoff=0)
        self.tasksmenu.add_command(label="Display", command=self.display)
        self.tasksmenu.add_command(label="Hide", command=self.hide)
        self.tasksmenu.add_command(label="Delete", command=self.deleteFile)

        # self.exportMenu = Menu(root, tearoff=0)
        # self.exportMenu.add_command(label="Selected Data to XSpec Format...", command=self.command, state=DISABLED)
        # self.exportMenu.add_command(label="Fit Results to File...", command=self.command, state=DISABLED)

        self.windowmenu = Menu(root, tearoff=0, postcommand=self.refresh_window_list)
        self.windowmenu.add_command(label='Cascade Windows', command=self.cascade_windows, state=DISABLED)
        self.windowmenu.add_command(label='Tile Windows', command=self.tile_windows, state=DISABLED)
        self.windowmenu.add_separator()

        menubar.add_cascade(label="Gspec", menu=filemenu)
        menubar.add_cascade(label="Options", menu=self.tasksmenu)
        # menubar.add_cascade(label="Export", menu=self.exportMenu)
        menubar.add_cascade(label="Windows", menu=self.windowmenu)

        # Add the menu bar to the root window
        root.config(menu=menubar)


        # Create the listbox
        self.filenames = []
        self.listBox = Listbox(self, selectmode=SINGLE,
                               borderwidth=2,
                               selectbackground='#116cd6',
                               selectborderwidth=0,
                               highlightthickness=0,
                               relief=SUNKEN)

        for i in self.filenames:
            self.listBox.insert(END, i)

        self.listBox.bind("<<ListboxSelect>>", self.onSelect)

        self.listBox.grid(column=0, row=0, padx=10, pady=10, sticky=N + E + W + S)
        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)

        # Create the buttons
        frame_left = Frame(self, bg="#e1e1e1", relief=FLAT)
        frame_center = Frame(self, bg="#e1e1e1", relief=FLAT)
        frame_right = Frame(self, bg="#e1e1e1", relief=FLAT)

        loadButton = Button(frame_center, text=u"Load", relief=FLAT, highlightbackground="#e1e1e1",
                            command=self.loadFile, width=5)
        loadButton.pack(padx=0, side=LEFT)
        # loadButton.grid(row=3, column=0, pady=0)

        displayButton = Button(frame_center, text=u"Display", relief=FLAT, highlightbackground="#e1e1e1",
                               command=self.display, width=5)
        displayButton.pack(padx=0, side=LEFT)
        # displayButton.grid(row=3, column=1, pady=0)

        hideButton = Button(frame_center, text=u"Hide", relief=FLAT, highlightbackground="#e1e1e1", command=self.hide, width=5)
        hideButton.pack(padx=0, side=LEFT)
        # hideButton.grid(row=3, column=2, pady=0)

        button = Button(frame_center, text=u"Delete", relief=FLAT, highlightbackground="#e1e1e1",
                        command=self.deleteFile, width=5)
        button.pack(padx=0, side=LEFT)
        # button.grid(row=3, column=3, pady=0)

        # Place the frame within the root window
        frame_center.grid(column=0, row=1, sticky=N + S)

        # Configure the weights of the columns
        self.rowconfigure(1, weight=1)
        self.columnconfigure(0, weight=1)

        # Fill the list box with any open data files
        self.populateListBox()

        # Bind the onWindowClose function to the window close action
        root.protocol("WM_DELETE_WINDOW", self.onWindowClose)

        # Create a binding to perform an action when the window comes into focus
        # root.bind("<FocusIn>", self.handle_focus)

        root.mainloop()


    ##########################################

    def handle_focus(self, event):

        # Close the splash window
        self.splash.destroy()        


    ##########################################

    def download(self):
        pass
        #self.download_root = Tk()
        #gui = TriggerFinder(self.download_root)
        #self.download_root.mainloop()

    
    ##########################################

    def onSelect(self, selected):

        # sender = selected.widget
        # idx = sender.curselection()
        # value = sender.get(idx)   

        # self.var.set(value)
        pass

    ##########################################

    def quitConfirmation(self):
        exit()

    ##########################################

    def loadFile(self):
        filetypes = (("PHA Files", "*.pha"), ("FITS Files", "*.fit"))
        selectedFilename = filedialog.askopenfilenames(filetypes=filetypes, parent=self.root)
        if selectedFilename == '':
            return

        num_current = len(self.filenames)
        selectedFilename = self.root.tk.splitlist(selectedFilename)
        for filename in selectedFilename:
            self.filenames.append(filename)
            self.listBox.insert(END, filename)
            self.gspecManager.add_data(filename)
        self.display(index_range=(num_current, END))

        # DataViewer(selectedFilename)

    ##########################################

    def loadLookup(self):
        # only allow gspec lookups right now
        filetypes = (("GSpec Lookup", "*.json"),)
        selectedFilename = filedialog.askopenfilename(filetypes=filetypes, parent=self.root)
        if selectedFilename == '':
            return
        
        # load lookup, register with Gspec, and launch window for each detector
        # in lookup
        num_current = len(self.filenames)
        dir = os.path.dirname(selectedFilename)
        lu = GspecLookup(filename=selectedFilename)
        for detector, info in lu._lookup_dict.items():
            datafile = os.path.join(dir, info['filename'])
            try:
                self.gspecManager.add_data(datafile, lookup=lu)
            except:
                pass
            self.filenames.append(datafile)
            self.listBox.insert(END, datafile)
        self.display(index_range=(num_current, END))

    ##########################################

    def deleteFile(self):

        selectedFilename = self.listBox.get(ACTIVE)
        if selectedFilename == '':
            return
        
        self.hide()
        self.filenames.remove(selectedFilename)
        self.listBox.delete(ACTIVE)        

        self.gspecManager.lookup.remove_data(selectedFilename)

    ##########################################

    def display(self, index_range=None):

        # Get the selected filename
        if index_range is None:
            selectedFilename = self.listBox.get(ACTIVE)
        else:
            selectedFilename = self.listBox.get(index_range[0], index_range[1])
        if selectedFilename == '':
            return

        # Remove any whitespaces that the listBox may have added to the string
        selectedFilename = self.root.tk.splitlist(selectedFilename)
        for filename in selectedFilename:
            # register with the window menu
            self.windowmenu.add_command(label=os.path.basename(filename),
                                        command=lambda fname=filename: self.window_focus(fname))

            # Display the data
            self.gspecManager.display_file(filename)

    ##########################################

    def hide(self):

        # Get the selected filename
        selectedFilename = self.listBox.get(ACTIVE)
        if selectedFilename == '':
            return

        # Remove any whitespaces that the listBox may have added to the string
        selectedFilename = "".join(selectedFilename.split())
        selectedFilename = os.path.basename(selectedFilename)

        # Hide the data
        try:
            dataviewer = self.gspecManager._openWindows[selectedFilename]
            dataviewer.onWindowClose()
        except:
            print('Window is already hidden')

    ##########################################

    def populateListBox(self):

        loadedFiles = self.gspecManager._data

        for detector in loadedFiles.keys():
            filename = loadedFiles[detector]['data'].filename
            self.filenames.append(filename)
            self.listBox.insert(END, filename)

    ##########################################

    def command(self):
        print('function not yet implmented')

    ##########################################

    def onWindowClose(self):
        if messagebox.askokcancel("Quit", "Are you sure you want to quit?"):
            openwindows = self.gspecManager._openWindows
            openwindows = [openwindows[key] for key in openwindows.keys()]
            num_windows = len(openwindows)
            for i in range(num_windows):
                try:
                    openwindows[i].onWindowClose()
                except:
                    pass
            self.root.destroy()

    ##########################################

    def tile_windows(self):
        # calculate number of columns and rows
        num_windows = len(self.gspecManager._openWindows)
        if num_windows == 0:
            return
        num_cols = numpy.ceil(numpy.sqrt(num_windows))
        num_rows = numpy.ceil(num_windows / num_cols)

        # screen dimensions       
        width = self.root.winfo_screenwidth()
        height = self.root.winfo_screenheight()

        # window tile dimensions keeping the same aspect ratio
        window_width = int(numpy.floor(width / num_cols))
        aspect_ratio = 0.5
        window_height = int(numpy.floor(aspect_ratio * float(window_width)))
        if (float(window_height + 45.0) * num_rows) > height:
            window_height = int(numpy.floor((height - 45.0 * num_rows) / num_rows))
            window_width = int(numpy.floor(float(window_height) / aspect_ratio))

        # go through each open window, move and resize
        icolumn = 0
        irow = 0
        for _, window in self.gspecManager._openWindows.items():
            if icolumn > num_cols - 1:
                icolumn = 0
                irow += 1
            xpos = icolumn * window_width
            ypos = irow * (window_height + 45)
            window.root.geometry("{0}x{1}+{2}+{3}".format(window_width,
                                                          window_height,
                                                          xpos, ypos))
            icolumn += 1

    ##########################################

    def cascade_windows(self):
        # calculate number of columns and rows
        num_windows = len(self.gspecManager._openWindows)
        if num_windows == 0:
            return

        # set up window dimensions
        width = self.root.winfo_screenwidth()
        height = self.root.winfo_screenheight()
        window_width = int(numpy.floor(width * 0.67))
        window_height = int(numpy.floor(window_width * 0.5))

        iwindow = 0
        for _, window in self.gspecManager._openWindows.items():
            xpos = int(iwindow * 45.0)
            ypos = int(iwindow * 45.0)
            window.root.geometry("{0}x{1}+{2}+{3}".format(window_width,
                                                          window_height,
                                                          xpos, ypos))
            window.root.lift()
            iwindow += 1

    ##########################################

    def window_focus(self, key):

        # set up window dimensions and position    
        width = self.root.winfo_screenwidth()
        height = self.root.winfo_screenheight()
        window_width = int(numpy.floor(width * 0.67))
        window_height = int(numpy.floor(window_width * 0.5))
        xpos = int(numpy.floor((width - window_width) / 2.0))
        ypos = int(numpy.floor((height - window_height) / 2.0))

        window = self.gspecManager._openWindows[key]
        window.root.geometry("{0}x{1}+{2}+{3}".format(window_width,
                                                      window_height, xpos, ypos))
       
        window.root.lift()

    ##########################################

    def refresh_window_list(self):
        
        numwin = len(self.gspecManager._openWindows)
        if numwin > 0:
            self.windowmenu.entryconfig(0, state=NORMAL)
            self.windowmenu.entryconfig(1, state=NORMAL)
        else:
            self.windowmenu.entryconfig(0, state=DISABLED)
            self.windowmenu.entryconfig(1, state=DISABLED)            
        
        # clear the list of data viewers and repopulate
        self.windowmenu.delete(3, index2=self.windowmenu.index("end"))
        for file in self.gspecManager._openWindows.keys():
            try:
                self.gspecManager._openWindows[file].root.winfo_screenwidth()
                self.windowmenu.add_command(label=os.path.basename(file),
                                            command=lambda fname=file: self.window_focus(fname))
            except:
                pass


##########################################################################################

class DataViewer(Frame):
    """Class for displaying pha files"""

    # def __init__(self, root, filename):
    # def __init__(self, filename, gspecManager, root):
    def __init__(self, filename, gspecManager):

        # Create a instance variable to reference the filename to be displayed
        self.filename = filename
        self.basename = os.path.basename(filename)

        # Create an instance variable to store the reference to the gspecManager 
        self.gspecManager = gspecManager
        self.detector = None

        # Extract the detector name for the supplied filename
        for detector in gspecManager._data.keys():
            if gspecManager._data[detector]['data'].filename == filename:
                self.detector = detector
                break

        # Add a pointer to the window to the list of open windows
        self.gspecManager._openWindows[self.basename] = self

        # Warn the user if the filename no longer exists
        if self.detector == None:
            print("\nError: data file not found:\n %s" % filename)

        # Set the current binning
        if type(self.gspecManager._data[self.detector]['data']).__name__ == 'Cspec':
            self.currentBinning = 1.024
        else:
            self.currentBinning = 0.064
        self.currentBinningFunction = None
        self.currentSNR = None

        # Setting default state parameters
        self.userSelectionArtists = []
        self.userSelectionPolyArtists = []
        self.coordinateAnnotationArtist = None
        self.sourceSelectionActive = False
        self.snrSelectionActive = False
        self.backgroundSelectionActive = False
        self.binningSelectionActive = False
        self.xZoomSelectionActive = False
        self.xlog_lightCurve = False
        self.ylog_lightCurve = False
        self.xlog_spectrum = True
        self.ylog_spectrum = True

        # Create the tkinter interface
        self.createTKinterInterface()

    ############################################

    def createTKinterInterface(self):

        # Create the root window
        # root = Toplevel()
        root = Tk()
        root.title(self.basename)
        root.geometry("1050x525+300+0")  # +850+300")
        root.minsize(500, 250)

        # Make a reference to the root window so that we can kill it later
        self.root = root

        # Initilize a frame instance with the root as its parent
        Frame.__init__(self, root)

        # Configure the frame
        self.config(bg="#e1e1e1")
        self.grid(sticky=N + W + E + S)

        # Create the menu bar
        menubar = Menu(root)

        # update the list of available headers
        headers = self.gspecManager._data[self.detector]['data'].header
        headerMenu = Menu(root)
        for ext in headers.keys():
            headerMenu.add_command(label=ext, command=lambda ext=ext: self.displayHeader(ext))

        lookupMenu = Menu(root)
        lookupMenu.add_command(label="Save Lookup", command=self.writeLookup)
        lookupMenu.add_command(label="Read Lookup", command=self.readLookup)
        lookupMenu.add_command(label="File Content", command=self.displayLookup)
        lookupMenu.add_command(label="Erase Current", command=self.eraseLookup)

        responseMenu = Menu(root)
        responseMenu.add_command(label="Load Response", command=self.loadResponse)
        responseMenu.add_command(label="Remove Response", command=self.removeResponse)
        responseMenu.add_command(label="Display Response", command=self.command, state=DISABLED)

        filemenu = Menu(menubar)
        filemenu.add_cascade(label="Header...", menu=headerMenu)
        filemenu.add_cascade(label="Lookup...", menu=lookupMenu)
        filemenu.add_cascade(label="Response...", menu=responseMenu)
        filemenu.add_command(label="Screenshot", command=self.command)
        filemenu.add_command(label="Dismiss", command=self.onWindowClose)

        optionsMenu = Menu(menubar)
        optionsMenu.add_command(label="Show Current Selection", command=self.command, state=DISABLED)
        optionsMenu.add_command(label="Colors", command=self.command, state=DISABLED)
        optionsMenu.add_command(label="Plot Configuration", command=self.command, state=DISABLED)

        windowMenu = Menu(menubar)
        windowMenu.add_command(label="Show GSPEC Window", command=self.command, state=DISABLED)
        windowMenu.add_command(label="Refresh Plot", command=self.redraw)

        menubar.add_cascade(label="File", menu=filemenu)
        menubar.add_cascade(label="Options", menu=optionsMenu)
        menubar.add_cascade(label="Window", menu=windowMenu)

        # Add the menu bar to the root window
        root.config(menu=menubar)

        ########### Buttons ###########

        frame_controls = Frame(self, bg="#e1e1e1", relief=FLAT, pady=10, padx=10, height=100)

        self.toggleButton = Button(frame_controls, text=u"Toggle", relief=FLAT, highlightbackground="#e1e1e1",
                                   command=self.toggle, width=14)
        self.rebinButton = Button(frame_controls, text=u"Rebin:", relief=FLAT, highlightbackground="#e1e1e1", width=14,
                                  command=self.rebin)
        self.fitBackgroundButton = Button(frame_controls, text=u"Fit Background", relief=FLAT,
                                          highlightbackground="#e1e1e1", width=14, command=self.selectBackground)
        self.sourceSelectButton = Button(frame_controls, text=u"Select Source", relief=FLAT,
                                         highlightbackground="#e1e1e1", command=self.selectSource, width=14)
        self.adjustSourceButton = Button(frame_controls, text=u"Adjust Source:", relief=FLAT,
                                         highlightbackground="#e1e1e1", command=self.adjustSource, width=14)
        # self.spectralFittingButton = Button(frame_controls, text=u"Spectral Fitting:", relief=FLAT, highlightbackground="#e1e1e1", width=12, command=self.fitData) 

        # # Zoom menu
        # self.zoomSelection = StringVar(frame_controls)
        # self.zoomSelection.set('Zoom')
        # choices = ['X Zoom','Y Zoom','Zoom Out: Full Range']
        # self.zoomMenu = OptionMenu(frame_controls, self.zoomSelection, *choices, command=self.zoom)
        # self.zoomMenu.config(bg="#e1e1e1", width=14, justify=CENTER)

        # Rebin menu
        self.rebinSelection = StringVar(frame_controls)
        self.rebinSelection.set('Rebin')
        choices = ['Full Resolution', 'Temporal Resolution', 'Signal to Noise', 'Combine Source Intervals',
                   'Combine Into Single Bin', 'Combine by Factor']
        self.rebinMenu = OptionMenu(frame_controls, self.rebinSelection, *choices, command=self.rebin)
        self.rebinMenu.config(bg="#e1e1e1", width=14, justify=CENTER)

        # Source selection menu
        self.sourceSelection = StringVar(frame_controls)
        self.sourceSelection.set('Source Selection')
        choices = ['Source Interactive', 'Source By S/N', 'Clear Selections']
        self.sourceSelectionMenu = OptionMenu(frame_controls, self.sourceSelection, *choices, command=self.selectSource)
        self.sourceSelectionMenu.config(bg="#e1e1e1", width=14, justify=CENTER)

        # Adjust source selection menu
        self.adjustSourceSelection = StringVar(frame_controls)
        self.adjustSourceSelection.set('Adjust Source')
        choices = ['< Shift Selection', '> Shift Selection', '< Left Selection', '> Left Selection',
                   '< Right Selection', '> Right Selection']
        self.adjustSourceSelectionMenu = OptionMenu(frame_controls, self.adjustSourceSelection, *choices,
                                                    command=self.adjustSource)
        self.adjustSourceSelectionMenu.config(bg="#e1e1e1", width=14, justify=CENTER)

        # Export Settings menu
        self.exportSettingsSelection = StringVar(frame_controls)
        self.exportSettingsSelection.set('Push Selections...')
        choices = ['Temporal Binning', 'Background', 'Source Selection', 'View Range', 'All of the Above']
        self.exportSettingsMenu = OptionMenu(frame_controls, self.exportSettingsSelection, *choices,
                                             command=self.exportSettings)
        self.exportSettingsMenu.config(bg="#e1e1e1", width=14, justify=CENTER)

        # Spectral fitting menu
        self.spectralFittingSelection = StringVar(frame_controls)
        self.spectralFittingSelection.set('Spectral Fitting')
        choices = ['Fit Selection', 'Batch Fit Selections', 'Fit Plotter']
        self.spectralFittingSelectionMenu = OptionMenu(frame_controls, self.spectralFittingSelection, *choices,
                                                       command=self.prepareSpectralFit)
        self.spectralFittingSelectionMenu.config(bg="#e1e1e1", width=14, justify=CENTER)

        if self.gspecManager.xspecExists == False:
            self.spectralFittingSelectionMenu['state'] = 'disabled'

        # Pack up the buttons
        self.toggleButton.pack(padx=0, side=TOP, fill=X)
        # self.zoomMenu.pack(padx=0, side=TOP)
        self.rebinMenu.pack(padx=0, side=TOP, fill=X)
        self.fitBackgroundButton.pack(padx=0, side=TOP, fill=X)
        self.sourceSelectionMenu.pack(padx=0, side=TOP, fill=X)
        self.adjustSourceSelectionMenu.pack(padx=0, side=TOP, fill=X)
        self.exportSettingsMenu.pack(padx=0, side=TOP, fill=X)
        self.spectralFittingSelectionMenu.pack(padx=0, side=TOP, fill=X)

        ########### Check boxes ###########

        self.xlogTK = BooleanVar(frame_controls)
        self.ylogTK = BooleanVar(frame_controls)

        self.checkbuttonXScale = Checkbutton(frame_controls, bg="#e1e1e1", text='X Log', variable=self.xlogTK,
                                             command=self.onSelectLogScale)
        self.checkbuttonXScale['state'] = 'disabled'
        self.checkbuttonXScale.pack(padx=(10, 0), pady=10, side=LEFT)
        self.checkbuttonYScale = Checkbutton(frame_controls, bg="#e1e1e1", text='Y Log', variable=self.ylogTK,
                                             command=self.onSelectLogScale)
        self.checkbuttonYScale.pack(padx=(5, 0), pady=10, side=LEFT)

        frame_controls.grid(column=0, row=0, sticky=W + N + E)

        ########### Labels ###########

        # The label frame
        header = self.gspecManager._data[self.detector]['data'].header['PRIMARY']
        # labels = Frame(self, bg="#e1e1e1", relief=FLAT, pady=0, padx=10, height=100)
        labels = Frame(self, bg="#e1e1e1", relief=FLAT, pady=0, padx=10)

        self.eventLabel = Label(labels, bg="#e1e1e1", text="Event: %s" % header['OBJECT']).grid(row=1, column=0,
                                                                                                sticky=N + W)
        self.telescopeLabel = Label(labels, bg="#e1e1e1", text="Telescope: %s" % header['TELESCOP']).grid(row=2,
                                                                                                          column=0,
                                                                                                          sticky=N + W)
        self.instrumentLabel = Label(labels, bg="#e1e1e1", text="Instrument: %s" % header['INSTRUME']).grid(row=3,
                                                                                                            column=0,
                                                                                                            sticky=N + W)
        self.detectorLabel = Label(labels, bg="#e1e1e1", text="Detector: %s" % self.detector).grid(row=4, column=0,
                                                                                                   sticky=N + W)
        self.dataTypes = Label(labels, bg="#e1e1e1", text="Data types: %s" % header['DATATYPE']).grid(row=5, column=0,
                                                                                                      sticky=N + W)
        self.blankSpace = Label(labels, bg="#e1e1e1", text="").grid(row=6, column=0, sticky=N + W)
        self.selectionLabel = Label(labels, bg="#e1e1e1", text="Selections:").grid(row=7, column=0, sticky=N + W)

        self.timeSelectionText = StringVar(labels)
        self.timeSelectionText.set(" -- - -- s")
        self.timeSelection = Label(labels, bg="#e1e1e1", textvariable=self.timeSelectionText).grid(row=8, column=0,
                                                                                                   sticky=N + W)

        self.energySelectionText = StringVar(labels)
        self.energySelectionText.set(" -- - -- keV")
        self.energySelection = Label(labels, bg="#e1e1e1", textvariable=self.energySelectionText).grid(row=9, column=0,
                                                                                                       sticky=N + W)

        # Add the labels to the grid
        labels.grid(column=0, row=1, sticky=N + W)

        ########### Keyboard bindings ###########

        # Bind the return key to the rebin function
        root.bind('<x>', self.keyboardEvent)
        root.bind('<y>', self.keyboardEvent)
        root.bind('<z>', self.keyboardEvent)
        root.bind('<s>', self.keyboardEvent)
        root.bind('<b>', self.keyboardEvent)
        root.bind('<f>', self.keyboardEvent)

        ########### Plot Window ###########

        # matplotlib.projections.register_projection(CustomAxes)

        # Create a tk frame that will contain the plotting window
        self.plot_frame = Frame(root, bg="#e1e1e1")

        # Create a figure instance and attach it to self for future reference
        self.figure = Figure(figsize=(5, 4), dpi=100)
        # self.ax = self.figure.add_subplot(111, projection="CustomAxes")

        # # Add the subplot axis to self for future reference
        self.ax = self.figure.add_subplot(111)
        # self.figure.tight_layout()

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

        # Plot the light curve by default
        self.photonLightCurveDisplayed = True
        self.plotLightcurve()
        # self.plotEnergySpectrum()

        # Add the matplotlib toolbar
        toolbar_frame = Frame(self.plot_frame, bg="#e1e1e1")
        toolbar_frame.grid(column=0, row=1, sticky=N + W)
        toolbar_frame.config(bg="#e1e1e1")
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, toolbar_frame)
        self.toolbar.config(bg="#e1e1e1")

        # Change the message label and button background colors
        self.toolbar._message_label.config(background="#e1e1e1")
        for button in self.toolbar.winfo_children():
            button.config(background="#e1e1e1")

        self.toolbar.update()

        self.canvas.get_tk_widget().grid(column=0, row=0, padx=0, pady=0, sticky=N + W + S + E)

        self.plot_frame.rowconfigure(0, weight=1)
        self.plot_frame.rowconfigure(1, weight=0)
        self.plot_frame.columnconfigure(0, weight=1)
        self.plot_frame.columnconfigure(1, weight=0)
        self.plot_frame.grid(column=1, row=0, padx=0, pady=0, sticky=N + W + S + E, rowspan=2)

        # have to update the the command *after* the toolbar has been instantiated
        filemenu.entryconfig(0, command=self.toolbar.save_figure)

        ########### Status Bar ###########

        # statusBar = Frame(self, bg="#e1e1e1", relief=FLAT, pady=0, padx=10, height=100)
        # self.statusBar_TimeSelection = " --: -- s,"        
        # self.statusBar_EnergySelection = " --: -- keV"
        # self.statusBarText=StringVar()  
        # self.statusBarText.set(self.statusBar_TimeSelection + self.statusBar_EnergySelection) 
        # self.statusBarLabel = Label(statusBar, bg="#e1e1e1", textvariable=self.statusBarText).grid(row=1, column=0, sticky=W)

        # statusBar.grid(column=0, row=3, sticky=N+W, columnspan=2)

        ########### Row and Column Configurations ###########

        # Configuration needed without matplotlib toolbar
        # self.rowconfigure(0, weight=0)
        # self.rowconfigure(1, weight=1)
        # self.columnconfigure(0, weight=0)        
        # self.columnconfigure(1, weight=1)

        # root.rowconfigure(0, weight=1)
        # root.columnconfigure(0, weight=1)

        # Configuration needed with matplotlib toolbar
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)

        root.rowconfigure(0, weight=1)
        root.columnconfigure(0, weight=0)
        root.columnconfigure(1, weight=2)

        root.protocol("WM_DELETE_WINDOW", self.onWindowClose)

        #root.mainloop()

        return

    ############################################

    def displayHeader(self, ext):
        headers = self.gspecManager._data[self.detector]['data'].header
        header = headers[ext].tostring(sep='\n', padding=False)
        title = '{0} header for {1}'.format(ext, self.detector)
        font = ("Andale Mono", 12)
        height = len(headers[ext])
        if (height*16.0) > 608:
            height=38
        self.textDisplay(header, title=title, state=DISABLED, font=font, 
                         height=height, bg="#e1e1e1")

    ############################################

    def displayLookup(self):
        lu = self.gspecManager.split_off_detector(self.detector).display_lookup()
        title = 'Lookup file for {0}'.format(self.detector)
        font = ("Roboto Mono", 12)
        height = lu.count('\n')
        if (height*16.0) > 608:
            height=38
        self.textDisplay(lu, width=40, title=title, state=DISABLED, font=font, 
                         height=height, bg="#e1e1e1")

    ############################################

    def textDisplay(self, thetext, title=None, **kwargs):
        top = Toplevel(self.root)
        top.config(bg="#e1e1e1")
        if title is not None:
            top.title(title)
        text = Text(top, height=20, width=80)
        text.insert(INSERT, thetext)
        text.config(**kwargs)
        text.grid()

        yscroll = Scrollbar(top, command=text.yview, bg="#e1e1e1")
        text['yscrollcommand'] = yscroll.set
        yscroll.grid(row=0, column=1, sticky='NSEW')


        #text.mainloop()

    ############################################

    def plotLightcurve(self):

        self.ax.clear()

        # Extract the rate history from the file
        lightcurves = self.gspecManager._data[self.detector]['timeview'].lightcurve

        for lightcurve in lightcurves:
            # Plot the data
            self.ax.step(lightcurve.bounds[0], lightcurve.rate, where='post', c='#394264', zorder=2)
            self.ax.errorbar(lightcurve.bin_centers(), lightcurve.rate, lightcurve.uncertainty, capsize=0, fmt='none',
                             ecolor='dimgrey', alpha=0.5, zorder=1)

            # Set the axis labels
            self.ax.set_ylabel('Count Rate (counts/s)')
            self.ax.set_xlabel('Time (s)')

        # Set the x-axis minor tick frequency
        minorLocator = AutoMinorLocator()
        self.ax.xaxis.set_minor_locator(minorLocator)

        # Set the y-axis minor tick frequency
        minorLocator = AutoMinorLocator()
        self.ax.yaxis.set_minor_locator(minorLocator)

        # Register that the light curve view is currenty displayed
        self.photonLightCurveDisplayed = True

        # Display any existing user selections
        if self.gspecManager._data[self.detector]['sourceview'] is not None:

            lightcurves = self.gspecManager._data[self.detector]['sourceview'].lightcurve
            spectra = self.gspecManager._data[self.detector]['sourceview'].spectrum

            for lightcurve in lightcurves:
                x_selected = lightcurve.bounds[0]
                x_selected = numpy.append(x_selected, lightcurve.bounds[1][-1])

                rate_selected = lightcurve.rate
                rate_selected = numpy.append(rate_selected, rate_selected[-1])

                rate_uncertainty_selected = lightcurve.uncertainty
                # rate_uncertainty_selected = numpy.append(rate_uncertainty_selected, rate_uncertainty_selected[-1])

                # Change the color of the selected data
                # self.ax.step(x_selected, rate_selected, where='post', color='#006400')                  # Green
                self.ax.step(x_selected, rate_selected, where='post', color='#9a4e0e', zorder=3)  # Orange
                # self.ax.step(x_selected, rate_selected, where='post', color='darkblue', zorder=3)
                # self.ax.step(x_selected, rate_selected, where='post', color='#3e4d8b', zorder=3)
                # self.ax.step(x_selected, rate_selected, where='post', color='darkgreen', zorder=3)
                # self.ax.step(x_selected, rate_selected, where='post', color='blueviolet')

                # self.ax.errorbar(lightcurve.bin_centers(), lightcurve.rate, lightcurve.uncertainty, capsize=0, fmt='none', ecolor='#006400', alpha=0.5, zorder=1)           # Green
                self.ax.errorbar(lightcurve.bin_centers(), lightcurve.rate, lightcurve.uncertainty, capsize=0,
                                 fmt='none', ecolor='#9a4e0e', alpha=0.5, zorder=1)  # Orange
                # self.ax.errorbar(lightcurve.bin_centers(), lightcurve.rate, lightcurve.uncertainty, capsize=0, fmt='none', ecolor='darkblue', alpha=0.5, zorder=1)
                # self.ax.errorbar(lightcurve.bin_centers(), lightcurve.rate, lightcurve.uncertainty, capsize=0, fmt='none', ecolor='#3e4d8b', alpha=0.5, zorder=1)
                # self.ax.errorbar(lightcurve.bin_centers(), lightcurve.rate, lightcurve.uncertainty, capsize=0, fmt='none', ecolor='blueviolet', alpha=0.5, zorder=1)

                # Get the current axis limits
                xmin = self.ax.get_xlim()[0]
                xmax = self.ax.get_xlim()[1]
                ymin = self.ax.get_ylim()[0]
                ymax = self.ax.get_ylim()[1]

                # Fill in the space below the data
                zeros = numpy.zeros(len(x_selected))
                zeros[:] = -1e10

                # Change the color of the boundry of the selected data 
                self.ax.plot([x_selected[0], x_selected[0]], [zeros[0], rate_selected[0]], color='#9a4e0e',
                             zorder=3)  # Orange
                self.ax.plot([x_selected[-1], x_selected[-1]], [zeros[-1], rate_selected[-1]], color='#9a4e0e',
                             zorder=3)  # Orange
                # self.ax.plot([x_selected[0], x_selected[0]], [zeros[0], rate_selected[0]], color='#3e4d8b', zorder=3)
                # self.ax.plot([x_selected[-1], x_selected[-1]], [zeros[-1], rate_selected[-1]], color='#3e4d8b', zorder=3)

                # artist = self.ax.fill_between(x_selected, rate_selected, zeros, color='#006400', alpha=0.2, step='post')
                artist = self.ax.fill_between(x_selected, rate_selected, zeros, color='#9a4e0e', alpha=0.2, step='post',
                                              zorder=4)
                # artist = self.ax.fill_between(x_selected, rate_selected, zeros, color='darkgreen', alpha=0.2, step='post', zorder=4)
                # artist = self.ax.fill_between(x_selected, rate_selected, zeros, color='#3e4d8b', alpha=0.2, step='post', zorder=4)          # Blue

                # Set the axis limits to what they were before plotting the selection region
                self.ax.set_xlim([xmin, xmax])
                self.ax.set_ylim([ymin, ymax])

            # Swap the first and last time selection if first > last
            firstTimeSelection = lightcurves[0].bounds[0][0]
            lastTimeSelection = lightcurves[-1].bounds[1][-1]
            if firstTimeSelection > lastTimeSelection: firstTimeSelection, lastTimeSelection = lastTimeSelection, firstTimeSelection

            # Swap the first and last energy selection if first > last
            firstEnergySelecton = spectra[0].bounds[0][0]
            lastEnergySelection = spectra[-1].bounds[1][-1]
            if firstEnergySelecton > lastEnergySelection: firstEnergySelecton, lastEnergySelection = lastEnergySelection, firstEnergySelecton

            self.timeSelectionText.set(
                " %0.2f - %0.2f s" % (lightcurves[0].bounds[0][0], lightcurves[-1].bounds[1][-1]))
            self.energySelectionText.set(" %0.2f - %0.2f keV" % (spectra[0].bounds[0][0], spectra[-1].bounds[1][-1]))

        # Display the background model, if defined
        background_model = self.gspecManager._data[self.detector]['bkgdmodel']
        if background_model is not None:

            b, u = background_model.integrate_energy(self.gspecManager._data[self.detector]['timeview'])
            # if self.gspecManager._data[self.detector]['sourceview'] is not None:
            #     b2, u2 = background_model.integrate_energy(self.gspecManager._data[self.detector]['sourceview'])
            #     background_spectrum = background_model.source_spectrum(self.gspecManager._data[self.detector]['sourceview'])

            times = []
            for lightcurves in self.gspecManager._data[self.detector]['timeview'].lightcurve:
                times.extend(lightcurves.bin_centers())

            # print("Plotting background model")

            self.ax.plot(times, b, linestyle='--', color='firebrick', alpha=0.85, zorder=11, linewidth=0.75)
            self.ax.fill_between(times, b + u, b - u, color='firebrick', alpha=0.5, zorder=10)

            # Plot the fit residuals
            # rate = self.gspecManager._data[self.detector]['sourceview'].lightcurve[0].rate
            # self.ax2.step(times, rate-b)

        # if background_model is not None:
        #     b, u = background_model.integrate_energy(g._data['n9']['fullview'])
        #         b2, u2 = background_model.integrate_energy(g._data['n9']['sourceview'])
        #         background_spectrum = background_model.source_spectrum(g._data['n9']['sourceview'])

        # Set any existing axes range limits
        if self.gspecManager.lookup._lookup_dict[self.detector]['time_display_view'] is not None:
            # Extract the user's last recorded axis range
            xmin = self.gspecManager.lookup._lookup_dict[self.detector]['time_display_view'][0]
            xmax = self.gspecManager.lookup._lookup_dict[self.detector]['time_display_view'][1]
            ymin = self.gspecManager.lookup._lookup_dict[self.detector]['time_display_view'][2]
            ymax = self.gspecManager.lookup._lookup_dict[self.detector]['time_display_view'][3]

            # Set the axis range
            self.ax.set_xlim([xmin, xmax])
            self.ax.set_ylim([ymin, ymax])

        # Register callback functions to catch when the user zooms or pans
        self.ax.callbacks.connect('xlim_changed', self.on_xlims_change)
        self.ax.callbacks.connect('ylim_changed', self.on_ylims_change)

        # Set the x and y scales
        if self.xlog_lightCurve == True:
            self.ax.set_xscale('log')
        if self.ylog_lightCurve == True:
            self.ax.set_yscale('log')

        self.photonLightCurveDisplayed = True
        self.canvas.draw()

    ############################################

    def plotEnergySpectrum(self):

        self.ax.clear()

        # Extract the energy spectrum from the file
        spectra = self.gspecManager._data[self.detector]['energyview'].spectrum

        for spectrum in spectra:
            # Setup the energy
            x = spectrum.bounds[0]
            x = numpy.append(x, spectrum.bounds[1][-1])

            # Set up the count rate
            rate = spectrum.rate
            rate = numpy.append(rate, spectrum.rate[-1])

            # Setup the count rate uncertainty
            uncertainty = spectrum.uncertainty
            uncertainty = numpy.append(uncertainty, spectrum.uncertainty[-1])

            # self.ax.step(x, rate, where='post', c='darkblue')
            self.ax.step(x, rate, where='post', c='#394264', zorder=2)
            self.ax.errorbar((spectrum.bounds[1] * spectrum.bounds[0]) ** 0.5, spectrum.rate, spectrum.uncertainty,
                             capsize=0, fmt='none', ecolor='#394264', alpha=0.5, zorder=1)

        # Set the x and y labels
        self.ax.set_xlabel('Energy (keV)')
        self.ax.set_ylabel('Rate (counts/s-keV)')

        # Set the x and y scales
        if self.xlog_spectrum == True:
            self.ax.set_xscale('log')
        if self.ylog_spectrum == True:
            self.ax.set_yscale('log')

        # Display any existing user selections
        if self.gspecManager._data[self.detector]['timeview'] is not None:

            # Extract the data
            lightcurves = self.gspecManager._data[self.detector]['energyview'].lightcurve
            spectra = self.gspecManager._data[self.detector]['timeview'].spectrum

            # Iterate through each of the spectra
            for spectrum in spectra:
                # Get the upper and lower yaxis limits of the current plot
                ymin, ymax = self.ax.get_ylim()

                # Create the artist to contain the user selection
                self.ax.plot([spectrum.bounds[0][0], spectrum.bounds[0][0]], [1e-99, 1e99], '--', c='#9a4e0e',
                             linewidth=0.75)
                self.ax.plot([spectrum.bounds[1][-1], spectrum.bounds[1][-1]], [1e-99, 1e99], '--', c='#9a4e0e',
                             linewidth=0.75)

                # Set the limits to the way they were before plotting the selection boundries
                self.ax.set_ylim([ymin, ymax])

            # Update the selection labels
            self.timeSelectionText.set(
                " %0.2f - %0.2f s" % (lightcurves[0].bounds[0][0], lightcurves[-1].bounds[1][-1]))
            self.energySelectionText.set(" %0.2f - %0.2f keV" % (spectra[0].bounds[0][0], spectra[-1].bounds[1][-1]))

        # Display the background model, if defined
        background_model = self.gspecManager._data[self.detector]['bkgdmodel']
        if background_model is not None:
            if self.gspecManager._data[self.detector]['sourceview'] is not None:
                bkgd_spectrum = background_model.source_spectrum_integrated(
                    self.gspecManager._data[self.detector]['sourceview'])
            else:
                bkgd_spectrum = background_model.source_spectrum_integrated(
                    self.gspecManager._data[self.detector]['timeview'])

                # have to multiply by binwidth this is because Rob demanded that he be given counts/sec for xspec
            bkgd_spectrum.exposure = bkgd_spectrum.exposure * bkgd_spectrum.bin_widths()
            bkgd_spectrum.rate = data.RateHisto.calc_rate(bkgd_spectrum.values, bkgd_spectrum.exposure)

            self.ax.step(bkgd_spectrum.bounds[0], bkgd_spectrum.rate, where='post', c='firebrick', alpha=0.85,
                         zorder=11, linewidth=0.75)

            t = numpy.array(bkgd_spectrum.bounds).T.flatten().tolist()
            r1 = [bkgd_spectrum.rate - bkgd_spectrum.uncertainty,
                  bkgd_spectrum.rate - bkgd_spectrum.uncertainty]
            r1 = numpy.array(r1).T.flatten().tolist()
            r2 = [bkgd_spectrum.rate + bkgd_spectrum.uncertainty,
                  bkgd_spectrum.rate + bkgd_spectrum.uncertainty]
            r2 = numpy.array(r2).T.flatten().tolist()

            self.ax.fill_between(t, r1, r2, color='firebrick', alpha=0.5, zorder=10)

        # Set any existing axes range limits
        if self.gspecManager.lookup._lookup_dict[self.detector]['energy_display_view'] is not None:
            # Extract the user's last recorded axis range
            xmin = self.gspecManager.lookup._lookup_dict[self.detector]['energy_display_view'][0]
            xmax = self.gspecManager.lookup._lookup_dict[self.detector]['energy_display_view'][1]
            ymin = self.gspecManager.lookup._lookup_dict[self.detector]['energy_display_view'][2]
            ymax = self.gspecManager.lookup._lookup_dict[self.detector]['energy_display_view'][3]

            # Set the axis range
            self.ax.set_xlim([xmin, xmax])
            self.ax.set_ylim([ymin, ymax])

        # Register callback functions to catch when the user zooms or pans
        self.ax.callbacks.connect('xlim_changed', self.on_xlims_change)
        self.ax.callbacks.connect('ylim_changed', self.on_ylims_change)

        self.photonLightCurveDisplayed = False
        self.canvas.draw()

    ############################################

    # Declare a callback function for xlimit changes
    def on_xlims_change(self, axes):

        # Update the lookup table's view display, depending on what plot is currently displayed
        if self.photonLightCurveDisplayed == True:
            self.gspecManager.lookup.add_time_display_view(self.detector, [self.ax.get_xlim()[0], self.ax.get_xlim()[1],
                                                                    self.ax.get_ylim()[0], self.ax.get_ylim()[1]])

        else:
            self.gspecManager.lookup.add_energy_display_view(self.detector, [self.ax.get_xlim()[0], self.ax.get_xlim()[1],
                                                                      self.ax.get_ylim()[0], self.ax.get_ylim()[1]])

        # print("updated xlims: ", self.ax.get_xlim())

    ############################################

    # Declare a callback function for ylimit changes
    def on_ylims_change(self, axes):

        # Update the lookup table's view display, depending on what plot is currently displayed
        if self.photonLightCurveDisplayed == True:
            self.gspecManager.lookup.add_time_display_view(self.detector, [self.ax.get_xlim()[0], self.ax.get_xlim()[1],
                                                                    self.ax.get_ylim()[0], self.ax.get_ylim()[1]])

        else:
            self.gspecManager.lookup.add_energy_display_view(self.detector, [self.ax.get_xlim()[0], self.ax.get_xlim()[1],
                                                                      self.ax.get_ylim()[0], self.ax.get_ylim()[1]])

        # print("updated ylims: ", self.ax.get_ylim())

    ############################################

    def keyboardEvent(self, event):

        # print(event.char)

        if 'x' in event.char:
            self.toolbar.zoom()

        if 'y' in event.char:
            self.toolbar.zoom()

        if 'z' in event.char:
            self.toolbar.zoom()

        if 's' in event.char:
            self.selectSource(['Source Interactive'])

        # if 'b' in event.char:
        #     self.selectBackground()

        # if 'f' in event.char:
        #     self.prepareSpectralFit(['Fit Selection'])


    ############################################

    def selectSource(self, menuSelection):

        # Interactively select the time selection
        if 'Source Interactive' in menuSelection:
            self.sourceSelectionActive = True

            # Disable all of the buttons while selecting the source interval
            self.toggleButton['state'] = 'disabled'
            self.rebinMenu['state'] = 'disabled'
            self.fitBackgroundButton['state'] = 'disabled'
            self.sourceSelectionMenu['state'] = 'disabled'
            self.adjustSourceSelectionMenu['state'] = 'disabled'
            self.exportSettingsMenu['state'] = 'disabled'
            self.spectralFittingSelectionMenu['state'] = 'disabled'

        if 'Source By S/N' in menuSelection:
            txt = ['Enter Signal-to-Noise Ratio']
            self.createTextInputDialog(275, 150, 'Source by S/N', txt, 5.0,
                                       'Source Select', self.snr_select, '')
        # Clear all time selections
        if 'Clear Selections' in menuSelection:

            if self.photonLightCurveDisplayed == True:

                # Clear the temporal selections
                self.gspecManager.remove_source_selections(self.detector)
                self.gspecManager.refresh(self.detector)
                self.redraw()

                self.timeSelectionText.set(" -- - -- s")

                # self.plotLightcurve()

            else:

                # Clear the energy selections
                self.gspecManager.remove_energy_selections(self.detector)
                self.gspecManager.refresh(self.detector)
                self.redraw()

                self.energySelectionText.set(" -- - -- s")

                # self.plotEnergySpectrum()

        # Reset the menu selection
        self.sourceSelection.set('Source Selection')

    ############################################

    def snr_select(self, value, fullrange, dummy):

        try:
            snr = float(value.get())
        except:
            snr = value

        if fullrange.get() == 'Full':
            time_interval = self.gspecManager._data[self.detector]['data'].time_range
            self.gspecManager.source_by_snr(self.detector, time_interval, snr)
            self.gspecManager.refresh(self.detector, background=False)
            self.redraw()
        else:
            self.snrSelectionActive = True
            self.currentSNR = snr

        # Close the dialog window
        self.root_dialog.destroy()

        return

    ############################################

    def adjustSource(self, menuSelection):

        # print("Adjust Source Selected")

        if self.gspecManager._data[self.detector]['sourceview'] is not None:

            if self.photonLightCurveDisplayed == True:

                # Get the currently displayed light curve
                # time_displayed = self.gspecManager._data[self.detector]['timeview'].lightcurve[0].bounds[0]
                time_displayed = self.gspecManager._data[self.detector]['timeview'].lightcurve[0].bin_centers()

                # Get the lower and upper bounds of the current source selection
                # selection_lower = self.gspecManager._data[self.detector]['sourceview'].lightcurve[0].bounds[0][0]
                # selection_upper = self.gspecManager._data[self.detector]['sourceview'].lightcurve[0].bounds[1][-1]
                selection_lower = self.gspecManager._data[self.detector]['sourceview'].lightcurve[0].bin_centers()[0]
                selection_upper = self.gspecManager._data[self.detector]['sourceview'].lightcurve[0].bin_centers()[-1]

                # Find the index of the lower and upper bounds of the current source selection
                index_lower = numpy.where(time_displayed == selection_lower)[0]
                index_upper = numpy.where(time_displayed == selection_upper)[0]

                # Get the maximum allowable index
                index_max = len(time_displayed) - 1

                print("\nMaximum allowable time/index: %s, %s" % (time_displayed[index_max], index_max))


            else:

                # Get the currently displayed spectrum
                # energy_displayed = self.gspecManager._data[self.detector]['energyview'].spectrum[0].bounds[0]
                energy_displayed = self.gspecManager._data[self.detector]['energyview'].spectrum[0].bin_centers()

                index_max = len(energy_displayed)

                # Get the lower and upper bounds of the current source selection
                # selection_lower = self.gspecManager._data[self.detector]['sourceview'].spectrum[0].bounds[0][0]
                # selection_upper = self.gspecManager._data[self.detector]['sourceview'].spectrum[0].bounds[1][-1]
                selection_lower = self.gspecManager._data[self.detector]['sourceview'].spectrum[0].bin_centers()[0]
                selection_upper = self.gspecManager._data[self.detector]['sourceview'].spectrum[0].bin_centers()[-1]

                # Find the index of the lower and upper bounds of the current source selection
                index_lower = numpy.where(energy_displayed == selection_lower)[0]
                index_upper = numpy.where(energy_displayed == selection_upper)[0]

                # Get the maximum allowable index
                index_max = len(energy_displayed) - 1

                print("\nMaximum allowable energy/index: %s, %s" % (energy_displayed[index_max], index_max))

            # Shift the index of the lower and upper bounds of the current source selection accordingly
            if '< Shift Selection' in menuSelection:

                if index_lower != 0:
                    index_lower = index_lower - 1
                    index_upper = index_upper - 1

            if '> Shift Selection' in menuSelection:

                if index_upper != index_max:
                    index_lower = index_lower + 1
                    index_upper = index_upper + 1

            if '< Left Selection' in menuSelection:

                if index_lower != 0:
                    index_lower = index_lower - 1

            if '> Left Selection' in menuSelection:

                if index_lower != index_max - 1:
                    index_lower = index_lower + 1

            if '< Right Selection' in menuSelection:

                if index_upper != 1:
                    index_upper = index_upper - 1

            if '> Right Selection' in menuSelection:

                if index_upper != index_max - 1:
                    index_upper = index_upper + 1

            if self.photonLightCurveDisplayed == True:

                # Create a new source selection
                selection_newLower = time_displayed[index_lower]
                selection_newUpper = time_displayed[index_upper]
                selection_new = (selection_newLower, selection_newUpper)

            else:

                # Create a new source selection
                selection_newLower = energy_displayed[index_lower]
                selection_newUpper = energy_displayed[index_upper]
                selection_new = (selection_newLower, selection_newUpper)

        if self.photonLightCurveDisplayed == True:

            print("New selection")
            print(index_lower, time_displayed[index_lower])
            print(index_upper, time_displayed[index_upper])

            # Apply the new source selection
            self.gspecManager.source_selection(self.detector, [selection_new])

        else:

            print("New selection")
            print(index_lower, energy_displayed[index_lower])
            print(index_upper, energy_displayed[index_upper])

            # Apply the new source selection
            self.gspecManager.energy_selection(self.detector, [selection_new])

        # Refresh the detector and redraw the plot
        self.gspecManager.refresh(self.detector)
        self.redraw()

        # Reset the menu button
        self.adjustSourceSelection.set('Adjust Source')

    ############################################

    def onSelectLogScale(self):

        if self.photonLightCurveDisplayed == True:
            self.xlog_lightCurve = self.xlogTK.get()
            self.ylog_lightCurve = self.ylogTK.get()

            if self.ylog_lightCurve == True:
                self.ax.set_ylim(bottom=1)
        else:
            self.xlog_spectrum = self.xlogTK.get()
            self.ylog_spectrum = self.ylogTK.get()

        self.redraw()

    ############################################

    def onMouseClick(self, event):

        # Plot the user zoom selection
        if event.inaxes is not None and self.xZoomSelectionActive == True:

            # Get the upper and lower yaxis limits of the current plot
            ymin, ymax = self.ax.get_ylim()

            # Plot the selection 
            artist, = self.ax.plot([event.xdata, event.xdata], [ymin, ymax], '--', c='dimgray')
            self.ax.set_ylim(ymin, ymax)
            
            # Decide whether to display the next selection or exit out of the selection mode and redraw the plot using the new limits
            if len(self.userSelectionArtists) == 0:

                # Add the artist to the list of user selection artists and redraw
                self.userSelectionArtists.append(artist)
                self.canvas.draw()

            else:

                # Get the selection values
                xmin = self.userSelectionArtists[0].get_xdata()[0]
                xmax = event.xdata

                # Add the xmin value to the lookup
                self.gspecManager.lookup._lookup_dict[self.detector]['time_display_view'][0] = xmin
                self.gspecManager.lookup._lookup_dict[self.detector]['time_display_view'][1] = xmax

                self.ax.set_xlim(xmin, xmax)

                # Reset the list of user selection artists
                self.userSelectionArtists = []

                # Deactivate the z xoom selection mode
                self.xZoomSelectionActive = False

                print("Deactivating the zoom select")

                # Redraw the plot using the new limits
                # self.canvas.draw()
                self.redraw()

                # Plot the user source or background selections
        if event.inaxes is not None and (self.sourceSelectionActive or
                                         self.backgroundSelectionActive or
                                         self.binningSelectionActive or
                                         self.snrSelectionActive):

            # Get the upper and lower yaxis limits of the current plot
            ymin, ymax = self.ax.get_ylim()

            # Create the artist to contain the user selection
            # artist, = self.ax.plot([event.xdata, event.xdata], [ymin, ymax], '--', c='#de7370')

            if self.backgroundSelectionActive or self.binningSelectionActive:
                artist, = self.ax.plot([event.xdata, event.xdata], [ymin, ymax], '--', c='firebrick', linewidth=0.75)

            if self.sourceSelectionActive or self.snrSelectionActive:
                artist, = self.ax.plot([event.xdata, event.xdata], [ymin, ymax], '--', c='#9a4e0e', linewidth=0.75)
            
            self.ax.set_ylim(ymin, ymax)

            # Add the artist to the list of user selection artists and redraw
            self.userSelectionArtists.append(artist)
            self.canvas.draw()

        # Exit user selection mode
        if event.inaxes is None and (self.sourceSelectionActive or
                                     self.backgroundSelectionActive or
                                     self.binningSelectionActive or
                                     self.snrSelectionActive):

            # Get the upper and lower yaxis limits of the current plot
            xmin, xmax = self.ax.get_xlim()

            # Get the size of the drawing window
            bbox = self.ax.get_window_extent().transformed(self.figure.dpi_scale_trans.inverted())
            width = bbox.width * self.figure.dpi
            halfwidth = width / 2.0

            # User clicked the left margin to exit the selection mode
            if event.x < (bbox.x0 * self.figure.dpi):

                # print("Exit Button Clicked")

                # Hide the buttons while the user is making the selection
                self.toggleButton['state'] = 'normal'
                # self.zoomMenu['state'] = 'normal'
                self.rebinMenu['state'] = 'normal'
                if self.photonLightCurveDisplayed:
                    self.fitBackgroundButton['state'] = 'normal'
                self.sourceSelectionMenu['state'] = 'normal'
                self.adjustSourceSelectionMenu['state'] = 'normal'
                self.exportSettingsMenu['state'] = 'normal'
                self.spectralFittingSelectionMenu['state'] = 'normal'

                # Reset the coordinate annotation artist
                if self.coordinateAnnotationArtist is not None:
                    self.coordinateAnnotationArtist.remove()
                    self.coordinateAnnotationArtist = None

                # Draw the selection area before exiting
                if len(self.userSelectionArtists) != 0:

                    # Create an empty list to contain the source intervals
                    selection_intervals = []

                    # Create a list to contain the xlimits of the user selection
                    xlimit_selections = []

                    # Make sure there are an even number of selections
                    if len(self.userSelectionArtists) % 2 == 1:

                        # Remove the last selection
                        del self.userSelectionArtists[-1]

                    # Get the x coordinates of the user selections
                    for artist in self.userSelectionArtists:
                        xlimit_selections.append(artist.get_xdata()[0])

                    # Populate the list of source intervals with tuples of the format (low, high)
                    while len(xlimit_selections) != 0:

                        # Get the upper and lower selections
                        xmax = xlimit_selections.pop()
                        xmin = xlimit_selections.pop()

                        # Swap the xmin and xmax values if the selections were made in reverse order
                        if (xmax > xmin) == False:
                            xmin, xmax = xmax, xmin

                        # Add the selections to the list
                        selection_intervals.append((xmin, xmax))


                    # Add the selection_intervals
                    if self.sourceSelectionActive and len(selection_intervals) != 0:

                        if self.photonLightCurveDisplayed:
                            self.gspecManager.source_selection(self.detector, selection_intervals)
                            self.gspecManager.refresh(self.detector)

                        if not self.photonLightCurveDisplayed:
                            self.gspecManager.energy_selection(self.detector, selection_intervals)
                            self.gspecManager.refresh(self.detector)

                    # Fit the background
                    if self.backgroundSelectionActive:
                        self.backgroundSelectionIntervals = selection_intervals
                        
                        # Create the polynomial selection window
                        if len(self.userSelectionArtists) > 0:
                            self.createPolyomialSelectionDialog()

                    # binning selection
                    if self.binningSelectionActive:
                        start = selection_intervals[0][0]
                        stop = selection_intervals[-1][-1]
                        # temporal binning
                        if self.photonLightCurveDisplayed:
                            datatype = self.gspecManager._data[self.detector]['data'].type
                            if self.currentBinning is None:
                                self.gspecManager.time_binning(self.detector, self.currentBinningFunction,
                                                               datatype=datatype,
                                                               start=start, stop=stop)
                            else:
                                self.gspecManager.time_binning(self.detector, self.currentBinningFunction,
                                                               self.currentBinning, datatype=datatype,
                                                               start=start, stop=stop)
                            self.gspecManager.refresh(self.detector, background=False)
                        # spectral binning
                        else:
                            if self.currentBinning is None:
                                self.gspecManager.energy_binning(self.detector, self.currentBinningFunction, 
                                                                 start=start, stop=stop)
                            else:
                                self.gspecManager.energy_binning(self.detector, self.currentBinningFunction, 
                                                                 self.currentBinning, start=start, stop=stop)
                                
                            self.gspecManager.refresh(self.detector)

                if self.snrSelectionActive:
                    self.gspecManager.source_by_snr(self.detector,
                                                    selection_intervals,
                                                    self.currentSNR)
                    self.gspecManager.refresh(self.detector, background=False)

                # Exit the selection mode
                self.sourceSelectionActive = False
                self.backgroundSelectionActive = False
                self.binningSelectionActive = False
                self.snrSelectionActive = False

                self.redraw()

                # # Redraw the canvas
                # self.canvas.draw()

                # for artist in self.userSelectionArtists:
                #     artist.remove()
                #     self.canvas.draw()

                self.userSelectionArtists = []

                # # Update the lookup dictionary for this detector
                # self.gspecManager

            # User clicked the right margin to clear the selections
            # if event.x > halfwidth:
            if event.x > ((bbox.x0 + bbox.width) * self.figure.dpi):

                print("Clear Selections Button Clicked")

                if self.sourceSelectionActive:

                    if self.photonLightCurveDisplayed:
                        print("Clearing temporal selection")

                        # Set the temporal selection to None
                        self.gspecManager.remove_source_selections(self.detector)
                        self.gspecManager.refresh(self.detector)

                        # Update the selection labels
                        self.timeSelectionText.set(" -- - -- s")

                        # Disable the source selection mode
                        self.sourceSelectionActive = False

                    if not self.photonLightCurveDisplayed:
                        print("Clearing energy selection")

                        # Set the spectral selection to None
                        self.gspecManager.remove_energy_selections(self.detector)
                        self.gspecManager.refresh(self.detector)

                        # Update the selection labels
                        self.energySelectionText.set(" -- - -- keV")
                        # Disable the source selection mode
                        # self.sourceSelectionActive = False

                if self.backgroundSelectionActive:
                    # Set the temporal selection to None
                    self.gspecManager.remove_background(self.detector)
                    self.gspecManager.refresh(self.detector)
                    self.userSelectionArtists = []

                    # Disable the background selection mode
                    #self.backgroundSelectionActive = False

                # for artist in self.userSelectionArtists:
                #     artist.remove()

                # for artist in self.userSelectionPolyArtists:
                #     artist.remove()

                # self.userSelectionArtists = []
                # self.userSelectionPolyArtists = []

                self.coordinateAnnotationArtist = None

                # Redraw the plot
                self.redraw()
                # self.canvas.draw()

            if event.y < (bbox.y0 * self.figure.dpi):

                print("Entering Manual Input")

                if self.sourceSelectionActive or self.backgroundSelectionActive or \
                        self.binningSelectionActive or self.snrSelectionActive:
                    self.manualSelectionDialog(xinput=True, yinput=False)

                # if self.zoomSelectActive:
                #    self.manualSelectionDialog(xinput=True, yinput=True)

    ############################################

    def onMouseMove(self, event):

        if event.inaxes is not None and (self.sourceSelectionActive or
                                         self.backgroundSelectionActive or
                                         self.binningSelectionActive or
                                         self.snrSelectionActive):

            if self.coordinateAnnotationArtist is not None:
                self.coordinateAnnotationArtist.remove()

            self.coordinateAnnotationArtist = self.ax.annotate('(x: %.2f, y:%.2f)' % (event.xdata, event.ydata),
                                                               xy=(0.02, 0.02), xytext=(0.02, 0.02),
                                                               xycoords='figure fraction')
            self.canvas.draw()

            # outside the plotting area
        elif event.inaxes is None and (self.sourceSelectionActive or
                                       self.backgroundSelectionActive or
                                       self.binningSelectionActive or
                                       self.snrSelectionActive):

            xmin, xmax = self.ax.get_xlim()

            bbox = self.ax.get_window_extent().transformed(self.figure.dpi_scale_trans.inverted())
            width = bbox.width * self.figure.dpi
            halfwidth = width / 2.0

            if event.x < (bbox.x0 * self.figure.dpi):
                if self.coordinateAnnotationArtist is not None:
                    self.coordinateAnnotationArtist.remove()

                self.coordinateAnnotationArtist = self.ax.annotate('Exit', xy=(0.02, 0.02), xytext=(0.02, 0.02),
                                                                   xycoords='figure fraction')
                self.canvas.draw()

            elif event.x > ((bbox.x0 + bbox.width) * self.figure.dpi):

                if self.coordinateAnnotationArtist is not None:
                    self.coordinateAnnotationArtist.remove()

                self.coordinateAnnotationArtist = self.ax.annotate('Clear Selections', xy=(0.02, 0.02),
                                                                   xytext=(0.02, 0.02), xycoords='figure fraction')
                self.canvas.draw()

            elif event.y < (bbox.y0 * self.figure.dpi):
                if self.coordinateAnnotationArtist is not None:
                    self.coordinateAnnotationArtist.remove()

                self.coordinateAnnotationArtist = self.ax.annotate('Manual Input', xy=(0.02, 0.02), xytext=(0.02, 0.02),
                                                                   xycoords='figure fraction')
                self.canvas.draw()

    ############################################

    def zoom(self, menuSelection):

        if 'Zoom Out: Full Range' in menuSelection:

            if self.photonLightCurveDisplayed == True:
                self.gspecManager.lookup._lookup_dict[self.detector]['time_display_view'] = []

            else:
                self.gspecManager.lookup._lookup_dict[self.detector]['energy_display_view'] = []

            self.redraw()

        self.zoomSelection.set('Zoom')

    ############################################

    def toggle(self, pushAction=False):

        print('Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


        if pushAction == True:
            for key in self.gspecManager._openWindows.keys():
                if self.basename not in key:
                    self.gspecManager._openWindows[key].toggle(pushAction=False)

        # Display the light curve plot
        if self.photonLightCurveDisplayed == False:

            # print("Displaying photon spectrum")

            self.checkbuttonXScale.deselect()
            self.checkbuttonXScale['state'] = 'disabled'
            self.fitBackgroundButton['state'] = 'normal'
            # update the binning menu
            binchoices = ['Full Resolution', 'Temporal Resolution', 'Signal to Noise', 
                   'Combine Source Intervals', 'Combine Into Single Bin', 
                   'Combine by Factor']
            self.rebinMenu['menu'].delete(0, 'end')
            for item in binchoices:
                self.rebinMenu['menu'].add_command(label=item, command=lambda item=item: self.rebin(item))
            # update the selection menu
            selectchoices = ['Source Interactive', 'Source By S/N', 'Clear Selections']
            self.sourceSelectionMenu['menu'].delete(0, 'end')
            for item in selectchoices:
                self.sourceSelectionMenu['menu'].add_command(label=item, command=lambda item=item: self.selectSource(item))
           
            if self.ylog_lightCurve == True:
                self.checkbuttonYScale.select()
            else:
                self.checkbuttonYScale.deselect()

            self.plotLightcurve()

            self.photonLightCurveDisplayed = True

            # Display the energy spectrum plot
        else:

            # print("Displaying energy spectrum")
            self.checkbuttonXScale['state'] = 'normal'
            self.fitBackgroundButton['state'] = 'disabled'
            # update the binning menu
            binchoices = ['Full Resolution', 'Signal to Noise',
                          'Combine Into Single Bin', 'Combine by Factor']
            self.rebinMenu['menu'].delete(0, 'end')
            for item in binchoices:
                self.rebinMenu['menu'].add_command(label=item, command=lambda item=item: self.rebin(item))
            # update the selection menu
            selectchoices = ['Source Interactive', 'Clear Selections']
            self.sourceSelectionMenu['menu'].delete(0, 'end')
            for item in selectchoices:
                self.sourceSelectionMenu['menu'].add_command(label=item, command=lambda item=item: self.selectSource(item))
                
            
            if self.xlog_spectrum == True:
                self.checkbuttonXScale.select()
            else:
                self.checkbuttonXScale.deselect()

            if self.ylog_spectrum == True:
                self.checkbuttonYScale.select()
            else:
                self.checkbuttonYScale.deselect()

            self.plotEnergySpectrum()

            self.photonLightCurveDisplayed = False

    ############################################

    def format_coord(self, x, y):

        return "%.2f" % x, "%.2f" % y

    ############################################

    def redraw(self):

        if self.photonLightCurveDisplayed == True:

            self.plotLightcurve()

        else:

            self.plotEnergySpectrum()

    ############################################

    def manualSelectionDialog(self, xinput=True, yinput=True, xdefaults=None, ydefaults=None):

        # Create a new root window
        bgcolor = "#e1e1e1"
        self.root_dialog = Toplevel(self.root)
        self.root_dialog.title('Selection')
        self.root_dialog.config(bg=bgcolor)

        # determine size
        xsize = 165
        ysize = 120
        if xinput:
            ysize += 50
        if yinput:
            ysize += 50

        # Center the root box
        width = self.root_dialog.winfo_screenwidth()
        height = self.root_dialog.winfo_screenheight()
        size = [250, 250]
        x = width / 2 - size[0] / 2
        y = height / 2 - size[1] / 2
        self.root_dialog.geometry("%sx%s+%s+%s" % (xsize, ysize, int(x), int(y)))

        # Create a dialog frame within the new root window 
        dialogframe = Frame(self.root_dialog, bg=bgcolor, height=100)
        dialogframe.grid(column=0, row=0)
        label = Label(dialogframe, text='Manual Selection Input', bg=bgcolor)
        label.grid(row=0, columnspan=2, pady=(5, 10))

        # initialize values
        self._x0default = DoubleVar(dialogframe, value=None)
        self._x1default = DoubleVar(dialogframe, value=None)
        self._y0default = DoubleVar(dialogframe, value=None)
        self._y1default = DoubleVar(dialogframe, value=None)

        # x axis input
        irow = 1
        if xinput:
            if xdefaults is None:
                xdefaults = (0.0, 0.0)
            label1 = Label(dialogframe, text='X Lo: ', bg=bgcolor).grid(row=irow, column=0)
            self._x0default.set(xdefaults[0])
            entry1 = ttk.Entry(dialogframe, textvariable=self._x0default, width=7).grid(row=irow, column=1)
            label2 = Label(dialogframe, text='X Hi: ', bg=bgcolor).grid(row=irow + 1, column=0)
            self._x1default.set(xdefaults[1])
            entry2 = ttk.Entry(dialogframe, textvariable=self._x1default, width=7).grid(row=irow + 1, column=1)
            irow += 2

        # y axis inputs
        if yinput:
            if ydefaults is None:
                ydefaults = (0.0, 0.0)
            label1 = Label(dialogframe, text='Y Lo: ', bg=bgcolor).grid(row=irow, column=0)
            self._y0default.set(ydefaults[0])
            entry1 = ttk.Entry(dialogframe, textvariable=self._y0default, width=7).grid(row=irow, column=1)
            label2 = Label(dialogframe, text='Y Hi: ', bg=bgcolor).grid(row=irow + 1, column=0)
            self._y1default.set(ydefaults[1])
            entry2 = ttk.Entry(dialogframe, textvariable=self._y1default, width=7).grid(row=irow + 1, column=1)
            irow += 2

        # accept and cancel buttons
        acceptButton = Button(dialogframe, text=u"Accept", relief=FLAT, highlightbackground=bgcolor,
                              command=lambda: self.manualSelectionAccept(xinput=xinput, yinput=yinput), width=12)
        acceptButton.grid(row=irow, columnspan=2, sticky=N + W + S + E, pady=(10, 0))
        cancelButton = Button(dialogframe, text=u"Cancel", relief=FLAT, highlightbackground=bgcolor,
                              command=self.root_dialog.destroy, width=12)
        cancelButton.grid(row=irow + 1, columnspan=2, sticky=N + W + S + E, pady=0)

        self.root_dialog.grid_columnconfigure(0, weight=1)
        self.root_dialog.grid_rowconfigure(0, weight=1)

        self.root_dialog.lift()

    ############################################

    # apply the inputs
    def manualSelectionAccept(self, xinput=None, yinput=None):
        if xinput:
            try:
                x0 = self._x0default.get()
                x1 = self._x1default.get()
            except:
                print('Entry must be a number')

        if yinput:
            try:
                y0 = self._y0default.get()
                y1 = self._y1default.get()
            except:
                print('Entry must be a number')

        # some sort of selection
        if xinput and not yinput:
            ymin, ymax = self.ax.get_ylim()
            if self.backgroundSelectionActive or self.binningSelectionActive:
                color = 'firebrick'
            elif self.sourceSelectionActive:
                color = '#9a4e0e'

            artist, = self.ax.plot([x0, x0], [ymin, ymax], '--', c=color, linewidth=0.75)
            self.userSelectionArtists.append(artist)
            artist, = self.ax.plot([x1, x1], [ymin, ymax], '--', c=color, linewidth=0.75)
            self.userSelectionArtists.append(artist)
            self.canvas.draw()

        # zoom
        if xinput and yinput:
            pass

        self.root_dialog.destroy()

    ############################################

    def selectBackground(self):

        # Disable all of the buttons while selecting the bcakground to be fit
        self.toggleButton['state'] = 'disabled'
        # self.zoomMenu['state'] = 'disabled'
        self.rebinMenu['state'] = 'disabled'
        self.fitBackgroundButton['state'] = 'disabled'
        self.sourceSelectionMenu['state'] = 'disabled'
        self.adjustSourceSelectionMenu['state'] = 'disabled'
        self.exportSettingsMenu['state'] = 'disabled'
        self.spectralFittingSelectionMenu['state'] = 'disabled'

        if self.gspecManager._data[self.detector]['bkgdview'] is not None:
            ymin, ymax = self.ax.get_ylim()
            for selection in self.gspecManager._data[self.detector]['bkgdview']._time_selections:
                artist, = self.ax.plot([selection[0], selection[0]], [ymin, ymax], '--', c='firebrick', linewidth=0.75)
                self.userSelectionArtists.append(artist)
                artist, = self.ax.plot([selection[1], selection[1]], [ymin, ymax], '--', c='firebrick', linewidth=0.75)
                self.userSelectionArtists.append(artist)
            self.canvas.draw()
                
        
        # Activate the user selection mode
        self.backgroundSelectionActive = True
        
        # The GUI now goes into user selection mode. Upon exiting the selection mode, the background selections will be fit

    ############################################

    def fitBackground(self, order):

        # Add a background model to the detector
        self.gspecManager.lookup.add_background(self.detector, self.backgroundSelectionIntervals, 'Polynomial', order)

        # Refresh the detector
        self.gspecManager.refresh(self.detector)

        # Close the dialog window
        self.root_dialog.destroy()

        # Redraw the plotting window
        self.redraw()

        # Create the reduced chisq diagnostic window
        top = Toplevel(self.root)
        width = self.root.winfo_screenwidth()
        height = self.root.winfo_screenheight()
        xsize, ysize = (640, 480)
        x = width / 2 - xsize / 2
        y = height / 2 - ysize / 2
        top.geometry("%sx%s+%s+%s" % (xsize, ysize, int(x), int(y)))

        # Create a tk frame that will contain the plotting window
        plot_frame = Frame(top, bg="#e1e1e1")

        # Create a figure and axis
        figure = Figure(figsize=(5, 4), dpi=100)
        ax = figure.add_subplot(111)

        # Set the format of the coordinate readout
        ax.format_coord = lambda x, y: ""

        # Set the minor tick frequency
        minorLocator = AutoMinorLocator()
        ax.xaxis.set_minor_locator(minorLocator)
        minorLocator = AutoMinorLocator()
        ax.yaxis.set_minor_locator(minorLocator)

        # Plot the chisq
        chisq = self.gspecManager._data[self.detector]['bkgdmodel'].chisq
        dof = self.gspecManager._data[self.detector]['bkgdmodel'].dof
        bounds = self.gspecManager._data[self.detector]['bkgdview']._get_energy_channel_bounds()
        edges = numpy.array(bounds)[:, 1]
        ax.step(edges, chisq / dof, where='post', color='#394264')
        ax.plot((edges[0], edges[-1]), (1.0, 1.0), '--', c='#8b0000')
        ax.set_ylabel('Reduced Chi-Square')
        ax.set_xlabel('Energy (keV)')
        ax.set_xscale('log')
        ax.set_xlim((edges[0], edges[-1]))
        ax.set_ylim((0.9 * numpy.min(chisq / dof), 1.1 * numpy.max(chisq / dof)))

        # Create a new plotting canvas and attach it to the widget
        canvas = FigureCanvasTkAgg(figure, master=top)
        canvas.show()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        # Bring the window to the foreground
        top.lift()

    ############################################

    def cancelBackgroundFit(self):

        # Close the dialog window
        self.root_dialog.destroy()

        pass

    ############################################

    def createPolyomialSelectionDialog(self):

        # Create a new root window
        self.root_dialog = Toplevel(self.root)
        self.root_dialog.title('Background')
        self.root_dialog.config(bg="#e1e1e1")

        xsize = 160
        ysize = 220

        # Center the root box
        width = self.root_dialog.winfo_screenwidth()
        height = self.root_dialog.winfo_screenheight()
        size = [250, 250]
        x = width / 2 - size[0] / 2
        y = height / 2 - size[1] / 2
        # self.root_dialog.geometry("275x175+%s+%s" % (int(x), int(y)))
        self.root_dialog.geometry("%sx%s+%s+%s" % (xsize, ysize, int(x), int(y)))

        # Create a dialog frame within the new root window 
        dialogframe = Frame(self.root_dialog, bg="#e1e1e1", height=100)
        dialogframe.config(bg="#e1e1e1")
        dialogframe.grid(column=0, row=0, sticky=(N, W, E, S))
        dialogframe.columnconfigure(0, weight=1)
        dialogframe.rowconfigure(0, weight=1)

        # Create the text label 
        label1 = Label(dialogframe, text='Select Order of', bg="#e1e1e1").pack(padx=0, side=TOP)
        label2 = Label(dialogframe, text='Background Polynomial', bg="#e1e1e1").pack(padx=0, side=TOP)

        # Create the buttons
        self.cancelButton = Button(dialogframe, text=u"Cancel", relief=FLAT, highlightbackground="#e1e1e1",
                                   command=self.cancelBackgroundFit, width=12)
        self.zerothOrderButton = Button(dialogframe, text=u"0", relief=FLAT, highlightbackground="#e1e1e1",
                                        command=lambda: self.fitBackground(0), width=12)
        self.firstOrderButton = Button(dialogframe, text=u"1", relief=FLAT, highlightbackground="#e1e1e1",
                                       command=lambda: self.fitBackground(1), width=12)
        self.secondOrderButton = Button(dialogframe, text=u"2", relief=FLAT, highlightbackground="#e1e1e1",
                                        command=lambda: self.fitBackground(2), width=12)
        self.thirdOrderButton = Button(dialogframe, text=u"3", relief=FLAT, highlightbackground="#e1e1e1",
                                       command=lambda: self.fitBackground(3), width=12)
        self.fourthOrderButton = Button(dialogframe, text=u"4", relief=FLAT, highlightbackground="#e1e1e1",
                                        command=lambda: self.fitBackground(4), width=12)

        # Pack up the buttons
        self.cancelButton.pack(padx=0, side=TOP)
        self.zerothOrderButton.pack(padx=0, side=TOP)
        self.firstOrderButton.pack(padx=0, side=TOP)
        self.secondOrderButton.pack(padx=0, side=TOP)
        self.thirdOrderButton.pack(padx=0, side=TOP)
        self.fourthOrderButton.pack(padx=0, side=TOP)

        self.root_dialog.lift()

    ############################################

    def createTextInputDialog(self, xsize, ysize, title, message, defaultValue, buttonLabel, buttonCommand, function,
                              allowFull=True, allowCustom=True):

        # Create a new root window
        self.root_dialog = Toplevel(self.root)
        self.root_dialog.title('Rebin')
        self.root_dialog.config(bg="#e1e1e1")
        self.root_dialog.resizable(False, False)

        # Center the root box
        width = self.root_dialog.winfo_screenwidth()
        height = self.root_dialog.winfo_screenheight()
        size = [250, 250]
        x = width / 2 - size[0] / 2
        y = height / 2 - size[1] / 2
        self.root_dialog.geometry("%sx%s+%s+%s" % (xsize, ysize, int(x), int(y)))

        # Expand the column to match the window size
        self.root_dialog.columnconfigure(0, weight=1)

        # Create the text label 
        row_root = 0
        for line in message:
            label = Label(self.root_dialog, text=line, bg="#e1e1e1")
            label.grid(column=0, row=row_root, columnspan=2, pady=7)
            row_root = row_root + 1

        # Create a new frame within the root window 
        entryframe = Frame(self.root_dialog, bg="#e1e1e1", height=100)
        entryframe.config(bg="#e1e1e1")
        entryframe.grid(column=0, row=row_root, pady=5)
        row_root = row_root + 1

        # Create a new row index for the entry dialog
        row_entry = 0
        
        # Create the entry form with use the specified default value
        defaultValueVar = DoubleVar(entryframe, value=defaultValue)
        entry = ttk.Entry(entryframe, width=7, textvariable=defaultValueVar)
        entry.grid(column=0, row=row_entry, rowspan=2)


        # create the radio button
        fullrange = StringVar(entryframe, value='Full')
        if allowFull:
            fullrange.set('Full')
        else:
            fullrange.set('Custom')

        if allowFull:
            radio1 = Radiobutton(entryframe, bg="#e1e1e1", text='Full Range',
                                 variable=fullrange, value='Full')
            radio1.grid(row=row_entry, column=1, sticky=W, padx=(20,5))
            row_entry = row_entry + 1

        if allowCustom:
            radio2 = Radiobutton(entryframe, bg="#e1e1e1", text='Custom Range',
                                 variable=fullrange, value='Custom')
            radio2.grid(row=row_entry, column=1, sticky=W, padx=(20,5))
            row_entry = row_entry + 1


        # Create the action button and bind it to the specified button command
        button = ttk.Button(self.root_dialog, text=buttonLabel, width=5,
                   command=lambda: buttonCommand(defaultValueVar, fullrange, function))
        button.grid(column=0, row=row_root, columnspan=2, sticky=S, pady=(10,5))

        # Bind the return key to the specified button command 
        self.root_dialog.bind('<Return>',
                              lambda event, value=defaultValueVar: buttonCommand(value, fullrange, function))
        self.root_dialog.lift()

    ##########################################

    def writeLookup(self):
        filename = self.gspecManager.save_lookup()
        print('Lookup file saved to {0}'.format(filename))

    ##########################################

    def readLookup(self):
        # allowed filetypes
        filetypes = (("GSPEC Lookups", "*.json"), ("RMfit Lookups", "*.lu"))
        selectedFilename = filedialog.askopenfilename(filetypes=filetypes, parent=self.root)
        if selectedFilename == '':
            return

        # load in the appropriate lookups
        if selectedFilename.split('.')[-1] == 'lu':
            tifile = None
            if 'tte' in selectedFilename:
                test_tifile = selectedFilename.split('.')[0] + '.ti'
                if os.path.isfile(test_tifile):
                    tifile = test_tifile
                    print('has ti file')
            lu = RmfitLookup(selectedFilename, ti_file=tifile)
            self.gspecManager.load_rmfit_lookup(self.detector, lu)
        else:
            lu = GspecLookup(selectedFilename)
            self.gspecManager.load_gspec_lookup(self.detector, lu)

        self.redraw()

    ##########################################

    def eraseLookup(self):
        self.gspecManager.remove_data(self.filename)
        self.gspecManager.add_data(self.filename)
        self.redraw()

    ##########################################

    def loadResponse(self):
        # allowed filetypes
        filetypes = (("RSP", "*.rsp"), ("RSP2", "*.rsp2"))
        selectedFilename = filedialog.askopenfilename(filetypes=filetypes, parent=self.root)
        if selectedFilename == '':
            return

        self.gspecManager.lookup.add_response(self.detector, selectedFilename)
        print('Loaded Response')

    ##########################################

    def removeResponse(self):
        self.gspecManager._data[self.detector]['response'] = None
        self.gspecManager.lookup._lookup_dict[self.detector]['response'] = None
        print('Removed Response')

    ############################################

    # def createYesNoDialog(xsize, ysize, title, message, defaultValue, buttonLabel, buttonCommand, function):

    #     # Create a new root window
    #     self.root_dialog = Tk()
    #     self.root_dialog.title('Rebin')
    #     self.root_dialog.config(bg="#e1e1e1")

    #     # Center the root box
    #     width = self.root_dialog.winfo_screenwidth()
    #     height = self.root_dialog.winfo_screenheight()
    #     size = [250,250]
    #     x = width/2 - size[0]/2
    #     y = height/2 - size[1]/2
    #     # self.root_dialog.geometry("275x175+%s+%s" % (int(x), int(y)))
    #     self.root_dialog.geometry("%sx%s+%s+%s" % (xsize, ysize, int(x), int(y)))

    #     # Create a dialog frame within the new root window 
    #     dialogframe = Frame(self.root_dialog, bg="#e1e1e1", height=100)
    #     dialogframe.config(bg="#e1e1e1")
    #     dialogframe.grid(column=0, row=0, sticky=(N, W, E, S))
    #     dialogframe.columnconfigure(0, weight=1)
    #     dialogframe.rowconfigure(0, weight=1)

    #     # Create the text label 
    #     row = 1
    #     for line in message:
    #         label = Label(dialogframe, text=line, bg="#e1e1e1").grid(column=1, row=row, sticky=N)
    #         # label.place(x=xsize/2, y=ysize/2, anchor="center")
    #         row = row + 1

    ############################################

    def exportSettings(self, menuSelection):

        # get other detectors
        detectors = []
        for windowKey in self.gspecManager._openWindows.keys():
            if self.basename not in windowKey:
                detectors.append(self.gspecManager._openWindows[windowKey].detector)

        # the settings to export
        if 'Temporal Binning' in menuSelection:
            self.gspecManager.import_time_binning(self.detector, detectors)
        elif 'Background' in menuSelection:
            self.gspecManager.import_background(self.detector, detectors)
        elif 'Source Selection' in menuSelection:
            self.gspecManager.import_source_selection(self.detector, detectors)
        elif 'View Range' in menuSelection:
            self.gspecManager.import_time_display_view(self.detector, detectors)
        elif 'All of the Above' in menuSelection:
            self.gspecManager.import_time_binning(self.detector, detectors)
            self.gspecManager.import_background(self.detector, detectors)
            self.gspecManager.import_source_selection(self.detector, detectors)
            self.gspecManager.import_time_display_view(self.detector, detectors)

        # redraw the windows
        for _, window in self.gspecManager._openWindows.items():
            window.redraw()

        self.exportSettingsSelection.set('Push Selections...')

    ############################################

    def rebinCommand(self, value, fullrange, binning_function):

        try:
            value_float = float(value.get())
        except:
            value_float = value
        if value_float is not None:
            self.currentBinning = value_float

        if fullrange.get() == 'Full':
            if self.photonLightCurveDisplayed:
                datatype = self.gspecManager._data[self.detector]['data'].type
                self.gspecManager.time_binning(self.detector, binning_function,
                                               self.currentBinning, datatype=datatype)
                self.gspecManager.refresh(self.detector, background=False)
            else:
                self.gspecManager.energy_binning(self.detector, binning_function,
                                                 self.currentBinning)
                self.gspecManager.refresh(self.detector)                
            self.redraw()
        else:
            self.binningSelectionActive = True
            self.currentBinningFunction = binning_function

        self.rebinSelection.set('Rebin')

        # Close the dialog window
        self.root_dialog.destroy()

        return

    ############################################

    def rebin(self, menuSelection):

        if 'Full Resolution' in menuSelection:
            if self.photonLightCurveDisplayed:
                self.gspecManager.clear_time_binning(self.detector)
            else:
                self.gspecManager.clear_energy_binning(self.detector)
            self.gspecManager.refresh(self.detector)
            self.redraw()

        if 'Temporal Resolution' in menuSelection:
            self.currentBinningFunction = 'Temporal Resolution'
            self.createTextInputDialog(275, 195, 'Rebin', [
                'Enter new data resolution in sec\n\nNOTE: Resolution must be an integer\nmultiple of the indepdentant axis data.'],
                                       self.currentBinning, 'Rebin', self.rebinCommand, 'Temporal Resolution')

        if 'Signal to Noise' in menuSelection:
            self.currentBinningFunction = 'Signal-to-Noise'
            self.createTextInputDialog(275, 145, 'Rebin', ['Enter the desired SNR per bin'], '', 'Rebin',
                                       self.rebinCommand, 'Signal-to-Noise')

        # if 'Refine by Half':
        #    self.rebinCommand(self.currentBinning/2.0, self.currentBinningFunction)

        if 'Combine Source Intervals' in menuSelection:
            if self.gspecManager._data[self.detector]['sourceview'] is None:
                raise ValueError('No Source Intervals Selected')

            times = self.gspecManager._data[self.detector]['sourceview']._time_selections
            for selection in times:
                datatype = self.gspecManager._data[self.detector]['data'].type
                self.gspecManager.time_binning(self.detector, 'Combine Into Single Bin',
                                               datatype=datatype, start=selection[0],
                                               stop=selection[1])
            self.gspecManager.refresh(self.detector, background=False)
            self.redraw()

        if 'Combine Into Single Bin' in menuSelection:
            self.currentBinningFunction = 'Combine Into Single Bin'
            self.currentBinning = None
            self.binningSelectionActive = True

        if 'Combine by Factor' in menuSelection:
            self.currentBinningFunction = 'Combine By Factor'
            self.createTextInputDialog(275, 150, 'Rebin', ['Enter the desired rebinning factor'], '', 'Rebin',
                                       self.rebinCommand, 'Combine By Factor')

        self.rebinSelection.set('Rebin')

        return

    ############################################

    def customSave(self, *args, **kwargs):
        print('function not yet implmented')
        return

    ############################################

    def custom_save_figure(self):
        print('function not yet implmented')
        return

    ############################################

    def command(self):
        print('function not yet implmented')

    ############################################

    def prepareSpectralFit(self, menuSelection):

        # print("Menu Selection:")
        # print(menuSelection)

        if 'Batch Fit Selections' in menuSelection: 
            
            print("\nBatch fit selected.")
            self.performBatchFit = True
            self.gspecManager.perform_spectral_fit(batch=True)

        elif 'Fit Selection' in menuSelection:

            self.performBatchFit = False
            self.gspecManager.perform_spectral_fit()

        elif 'Fit Plotter' in menuSelection: 

            # Launch the fit plotter and pass the current xspec session
            if './Spectral Fit Display' not in self.gspecManager._openWindows.keys():
                FitPlotter(self.gspecManager, None)
            else:
                # self.window_focus('./Spectral Fit Display')
                fitPlotter = self.gspecManager._openWindows['./Spectral Fit Display']
                fitPlotter.root.lift()
        else:

            print("Function not yet implemented")

        self.spectralFittingSelection.set('Spectral Fitting')

        return

    ############################################

    def onWindowClose(self):
        try:
            del self.gspecManager._openWindows[self.basename]
        except:
            pass
        self.root.destroy()


##########################################################################################


class SpectralFitPreperationDialog(Frame,):
    """Class for displaying the spectral fit preperation window"""

    # def __init__(self, root, filename):
    def __init__(self, gspecManager, dataFiles, batch=False):

        # Create an instance variable to store the reference to the gspecManager 
        self.gspecManager = gspecManager

        # Create a new xspec session instance if one does not already exist
        try:
            self.gspecManager._openWindows['./XSPEC Log'].root.winfo_screenwidth()
            self.xspec = self.gspecManager._openWindows['./XSPEC Log']
        except:
            self.xspec = xspec.XSPEC()
            print('\nCreating XSPEC interface')            

        # Create an instance variable to store the reference to the data list 
        self.dataFiles = dataFiles

        self.performBatchFit = batch

        # Set some defaults
        self.fitStatistic = 'chi'
        self.fitWeight = 'standard'

        # Create the tkinter interface
        self.createTKinterInterface()

    ############################################

    def createTKinterInterface(self):

        # Create a new root window
        self.root_dialog = Tk()
        self.root_dialog.title('Photon Model')
        self.root_dialog.config(bg="#e1e1e1")

        xsize = 375
        ysize = 610

        # Center the root box
        width = self.root_dialog.winfo_screenwidth()
        height = self.root_dialog.winfo_screenheight()
        size = [250, 250]
        x = width / 2 - size[0] / 2
        y = height / 2 - size[1] / 2
        # self.root_dialog.geometry("275x175+%s+%s" % (int(x), int(y)))
        self.root_dialog.geometry("%sx%s+%s+%s" % (xsize, ysize, int(x), int(y)))

        # Create a dialog frame within the new root window 
        dialogframe = Frame(self.root_dialog, bg="#e1e1e1")
        dialogframe.pack(anchor=N, fill=BOTH, expand=True, side=TOP)

        # Create the text label 
        label_title = Label(dialogframe, text='Select one or more photon models', bg="#e1e1e1")
        label_title.config(justify=CENTER)

        # Create the listbox
        self.filenames = ['Power Law', 'Smoothly Broken Power Law']
        self.photonModelListBox = Listbox(dialogframe, selectmode=SINGLE,
                                          borderwidth=2,
                                          selectbackground='#116cd6',
                                          selectborderwidth=0,
                                          highlightthickness=0,
                                          relief=SUNKEN)
        self.photonModelListBox.bind('<<ListboxSelect>>', self.onSelectModel)
        
        # Populate the list box with available models
        #models = self.xspec.paramStr
        models = self.xspec.models.model_names()

        for model in models:
            self.photonModelListBox.insert(END, model)

        # Bind the list box and add it to the grid
        # self.listBox.bind("<<ListboxSelect>>", self.onSelect)    
        # self.listBox.grid(column=0, row=1, padx=10, pady=10, sticky=N+E+W+S)
        # self.rowconfigure(0, weight=1)
        # self.columnconfigure(0, weight=1)

        # Create a dialog frame within the new root window
        frame_modelParameters = Frame(dialogframe, bg="#e1e1e1")

        # Create the text label 
        label_modelParameters = Label(frame_modelParameters, text='Photon Model Parameters:', bg="#e1e1e1")
        label_modelParameters.config(justify=CENTER)
        label_modelParameters.pack(pady=5, side=TOP)

        self.setParameters = Button(frame_modelParameters, text=u"Set Parameters", relief=FLAT, highlightbackground="#e1e1e1",
                                          command=self.onSelectNewParameters, width=5)
        self.setParameters['state'] = 'disabled'
        self.setParameters.pack(fill=X)
        #self.modelParameterStrategyTK = StringVar(frame_modelParameters)
        #self.modelParameterStrategyTK.set('Keep current')
        #self.checkbutton_keepCurrentParameters = Radiobutton(frame_modelParameters, bg="#e1e1e1", text='Keep current',
        #                                                     value='Keep current',
        #                                                     variable=self.modelParameterStrategyTK,
        #                                                     command=self.onSelectKeepParameters)
        #self.checkbutton_setNewParameters = Radiobutton(frame_modelParameters, bg="#e1e1e1", text='Set parameters',
        #                                                value='Set parameters', variable=self.modelParameterStrategyTK,
        #                                                command=self.onSelectNewParameters)
        
        #self.checkbutton_keepCurrentParameters.config(bg="#e1e1e1", width=12, justify=CENTER)
        #self.checkbutton_setNewParameters.config(bg="#e1e1e1", width=12, justify=CENTER)
        #self.checkbutton_keepCurrentParameters.pack(side=LEFT, padx=5)
        #self.checkbutton_setNewParameters.pack(side=LEFT, padx=5)

        # Create a dialog frame within the new root window
        frame_fitStatistics_upper = Frame(dialogframe, bg="#e1e1e1")
        frame_fitStatistics_lower = Frame(dialogframe, bg="#e1e1e1")

        # Create the text label   
        label_fitStatistic = Label(frame_fitStatistics_upper, text='Fitting Statistic:', bg="#e1e1e1")
        label_fitStatistic.config(justify=CENTER)
        label_fitStatistic.pack(pady=5, side=TOP)

        # Create the fit statistic variables
        self.fitStatisticTK = StringVar(frame_fitStatistics_upper)
        self.checkbutton_chiSqr = Radiobutton(frame_fitStatistics_upper, bg="#e1e1e1", text='Chi2', value='chi',
                                              variable=self.fitStatisticTK, command=self.onSelectFitStatistic)
        self.checkbutton_cstat = Radiobutton(frame_fitStatistics_upper, bg="#e1e1e1", text='C-Stat', value='cstat',
                                             variable=self.fitStatisticTK, command=self.onSelectFitStatistic)
        self.checkbutton_lstat = Radiobutton(frame_fitStatistics_upper, bg="#e1e1e1", text='L-Stat', value='lstat',
                                             variable=self.fitStatisticTK, command=self.onSelectFitStatistic)
        self.checkbutton_pgstat = Radiobutton(frame_fitStatistics_lower, bg="#e1e1e1", text='PG-Stat', value='pgstat',
                                              variable=self.fitStatisticTK, command=self.onSelectFitStatistic)
        self.checkbutton_pstat = Radiobutton(frame_fitStatistics_lower, bg="#e1e1e1", text='P-Stat', value='pstat',
                                             variable=self.fitStatisticTK, command=self.onSelectFitStatistic)
        self.checkbutton_whittle = Radiobutton(frame_fitStatistics_lower, bg="#e1e1e1", text='Whittle', value='whittle',
                                               variable=self.fitStatisticTK, command=self.onSelectFitStatistic)

        self.checkbutton_chiSqr.pack(side=LEFT)
        self.checkbutton_cstat.pack(side=LEFT)
        self.checkbutton_lstat.pack(side=LEFT)
        self.checkbutton_pgstat.pack(side=LEFT)
        self.checkbutton_pstat.pack(side=LEFT)
        self.checkbutton_whittle.pack(side=LEFT)

        # Create a dialog frame within the new root window 
        frame_fitWeights = Frame(dialogframe, bg="#e1e1e1")

        # Create the text label   
        label_fitWeights = Label(frame_fitWeights, text='Fit Weighting:', bg="#e1e1e1")
        label_fitWeights.config(justify=CENTER)
        label_fitWeights.pack(pady=5, side=TOP)

        # Create the fit statistic variables
        self.fitWeightingTK = StringVar(frame_fitWeights)
        self.checkbutton_standard = Radiobutton(frame_fitWeights, bg="#e1e1e1", text='Standard', value='standard',
                                                variable=self.fitWeightingTK, command=self.onSelectFitWeight)
        self.checkbutton_model = Radiobutton(frame_fitWeights, bg="#e1e1e1", text='Model', value='model',
                                             variable=self.fitWeightingTK, command=self.onSelectFitWeight)
        self.checkbutton_gehrels = Radiobutton(frame_fitWeights, bg="#e1e1e1", text='Gehrels', value='gehrels',
                                               variable=self.fitWeightingTK, command=self.onSelectFitWeight)
        self.checkbutton_churazov = Radiobutton(frame_fitWeights, bg="#e1e1e1", text='Churazov', value='churazov',
                                                variable=self.fitWeightingTK, command=self.onSelectFitWeight)

        self.checkbutton_standard.pack(side=LEFT)
        self.checkbutton_model.pack(side=LEFT)
        self.checkbutton_gehrels.pack(side=LEFT)
        self.checkbutton_churazov.pack(side=LEFT)

        # Undetermined Values in Batch Fit 
        frame_undeterminedValues = Frame(dialogframe, bg="#e1e1e1")
        label_undeterminedValues = Label(frame_undeterminedValues, text='Undetermined Values in Batch Fit:',
                                         bg="#e1e1e1")
        label_undeterminedValues.config(justify=CENTER)
        label_undeterminedValues.pack(pady=5, side=TOP)

        self.undeterminedValues = StringVar(frame_undeterminedValues)
        self.checkbuttonLeaveFree = Radiobutton(frame_undeterminedValues, bg="#e1e1e1", text='Leave free',
                                                value='Leave free', variable=self.undeterminedValues)
        self.checkbuttonAutoFix = Radiobutton(frame_undeterminedValues, bg="#e1e1e1", text='Automatically fix',
                                              value='Automatically fix', variable=self.undeterminedValues)
        self.checkbuttonLeaveFree.pack(side=LEFT)
        self.checkbuttonAutoFix.pack(side=LEFT)

        frame_buttons = Frame(dialogframe, bg="#e1e1e1")
        self.acceptFitParameters = Button(frame_buttons, text=u"Accept", relief=FLAT, highlightbackground="#e1e1e1",
                                          command=self.performSpectralFit, width=5)
        self.acceptFitParameters['state'] = 'disabled'
        self.omitFitParameters = Button(frame_buttons, text=u"Omit", relief=FLAT, highlightbackground="#e1e1e1",
                                        command=self.omitFit, width=5)
        self.omitFitParameters['state'] = 'disabled'
        self.restoreFitParameters = Button(frame_buttons, text=u"Restore", relief=FLAT, highlightbackground="#e1e1e1",
                                           command=self.restoreFit, width=5)
        self.restoreFitParameters['state'] = 'disabled'
        self.cancelFitParameters = Button(frame_buttons, text=u"Cancel", relief=FLAT, highlightbackground="#e1e1e1",
                                          command=self.cancelFit, width=5)

        self.acceptFitParameters.pack(side=LEFT)
        self.omitFitParameters.pack(side=LEFT)
        self.restoreFitParameters.pack(side=LEFT)
        self.cancelFitParameters.pack(side=LEFT)

        # Set the default values
        #self.checkbutton_keepCurrentParameters.select()
        self.checkbutton_chiSqr.select()
        self.checkbutton_standard.select()
        self.checkbuttonLeaveFree.select()

        # Pack it all together
        label_title.pack(pady=10, side=TOP)
        self.photonModelListBox.pack(padx=10, side=TOP, fill=X)
        frame_modelParameters.pack(pady=10, side=TOP)
        frame_fitStatistics_upper.pack(pady=(10, 0), side=TOP)
        frame_fitStatistics_lower.pack(pady=(0, 10), side=TOP)
        frame_fitWeights.pack(pady=10, side=TOP)
        frame_undeterminedValues.pack(pady=10, side=TOP)
        frame_buttons.pack(padx=4, pady=15, side=TOP)

    ############################################

    def onSelectModel(self, event):
        self.acceptFitParameters['state'] = 'normal'
        self.setParameters['state'] = 'normal'
        
    def onSelectKeepParameters(self):
        pass

    def onSelectNewParameters(self):
        selectedPhotonModel = self.photonModelListBox.get(ACTIVE)
        params = self.xspec.models.parameter_settings(selectedPhotonModel)
        self.parameterSelectionDialog(selectedPhotonModel, params) 

    def onSelectLeaveFree(self):
        self.checkbuttonAutoFix.toggle()

    def onSelectAutoFix(self):
        self.checkbuttonLeaveFree.toggle()

    def onSelectFitStatistic(self):
        self.fitStatistic = self.fitStatisticTK.get()

    def onSelectFitWeight(self):
        self.fitWeight = self.fitWeightingTK.get()

    def onSelect(self):
        # print(self.undeterminedValues.get())
        pass

    def cancelFit(self):
        # Destroy the dialog window
        self.root_dialog.destroy()

    def restoreFit(self):
        print('function not yet implmented')

    def omitFit(self):
        print('function not yet implmented')

    ############################################

    def performSpectralFit(self):

        # Get the selected model name
        selectedPhotonModel = self.photonModelListBox.get(ACTIVE)

        # Destroy the dialog window that brought us here
        self.root_dialog.destroy()

        # Display the xspec session logger
        try:
            self.gspecManager._openWindows['./XSPEC Log'].root.winfo_screenwidth()
        except:
            self.xspec.display_logger()
            self.gspecManager._openWindows['./XSPEC Log'] = self.xspec

        # Display the current fit parameters
        print("\nSelected fit parameters...")
        print("Model: %s" % selectedPhotonModel)
        print("Statistic: %s" % self.fitStatistic)
        print("Weight: %s" % self.fitWeight)
        print("\nData Files:")
        for phafilename in self.dataFiles:
            print(phafilename)

        print("\nFitting data...")


        if self.performBatchFit == False:

            # Load the pha file. 
            slotNumber = 0
            for phafilename in self.dataFiles:
                slotNumber = slotNumber + 1
                self.xspec.data_cmd(fullfilename=phafilename, slotNumber=slotNumber)

            # Select the fitting statistic
            self.xspec.sel_cmd(picked=self.fitStatistic)

            # Select the model to be fit
            self.xspec.weigh_cmd(picked=self.fitWeight)

            # Select the model
            self.xspec.select_model(selectedPhotonModel)

            # Extract the number of parameters associated with the model
            numberOfParameters = self.xspec.selnparam

            # Fit the data    
            self.xspec.fit_cmd()

            numberOfDataGroups=len(self.dataFiles)
            data = self.xspec.extract_plot_data(numberOfDataGroups, numberOfParameters) 

            data['batch'] = False

        else:

            # Get the HDU list
            hdulist = fits.open(self.dataFiles[0])

            # Get the number of spectra in the file
            spectrum = hdulist['SPECTRUM'].data
            counts = spectrum['COUNTS']
            time_start = spectrum['TIME']
            time_end = spectrum['ENDTIME']
            exposure = spectrum['EXPOSURE']

            numberOfExtensions = len(counts)

            data = {}
            data['parameter_values'] = []
            data['parameter_sigmas'] = []

            progressBar = ProgressBar(numberOfExtensions, 0)

            for extensionNumber in range(numberOfExtensions):

                # The first extension should always be one, not zero
                extensionNumber = extensionNumber + 1

                # Load the pha file. 
                slotNumber = 0
                for phafilename in self.dataFiles:
                    slotNumber = slotNumber + 1
                    self.xspec.data_cmd(fullfilename=phafilename, extensionNumber=extensionNumber, slotNumber=slotNumber)

                if extensionNumber == 1:

                    # Select the fitting statistic
                    self.xspec.sel_cmd(picked=self.fitStatistic)

                    # Select the model to be fit
                    self.xspec.weigh_cmd(picked=self.fitWeight)

                    # Select the model
                    self.xspec.select_model(selectedPhotonModel)

                    # Extract the number of parameters associated with the model
                    numberOfParameters = self.xspec.selnparam

                print("Fitting spectrum %s" % extensionNumber)

                # Fit the data    
                self.xspec.fit_cmd()

                # Update the progress bar
                progressBar.updateValue(extensionNumber)
                progressBar.lift()

                numberOfDataGroups=len(self.dataFiles)
                data_single_fit = self.xspec.extract_plot_data(numberOfDataGroups, numberOfParameters) 

                if extensionNumber == 1:
                    data['time_start'] = time_start
                    data['time_end'] = time_end
                    data['exposure'] = exposure
                    data['parameter_names'] = data_single_fit['parameter_names']
                    data['parameter_values'] = data_single_fit['parameter_values']
                    data['parameter_sigmas'] = data_single_fit['parameter_sigmas']

                else:
                    data['parameter_values'] = numpy.vstack((data['parameter_values'], data_single_fit['parameter_values']))
                    data['parameter_sigmas'] = numpy.vstack((data['parameter_sigmas'], data_single_fit['parameter_sigmas']))

                data['batch'] = True

            # Destroy the progress bar
            progressBar.destroy()

        # print("Done.")

        # Launch the fit plotter and pass the current xspec session
        if './Spectral Fit Display' not in self.gspecManager._openWindows.keys():
            FitPlotter(self.gspecManager, data)
        else:

            fitplotter = self.gspecManager._openWindows['./Spectral Fit Display']

            fitplotter.data = data

            if fitplotter.data['batch'] == False:
                fitplotter.xlog = True
                fitplotter.ylog = True
                fitplotter.plotSpectrum()
            else: 
                fitplotter.xlog = False
                fitplotter.ylog = False
                fitplotter.plotBatchFitResults()


    def parameterSelectionDialog(self, modelname, params):
        
        # Create a new root window
        bgcolor = "#e1e1e1"
        self.root = Toplevel(self.root_dialog)
        self.root.title('Model Parameters')
        self.root.config(bg=bgcolor)

        # determine size
        nparams = len(params)
        xsize = 260
        ysize = 110
        ysize += 50*nparams

        # Center the root box
        width = self.root.winfo_screenwidth()
        height = self.root.winfo_screenheight()
        size = [250, 250]
        x = width / 2 - size[0] / 2
        y = height / 2 - size[1] / 2
        self.root.geometry("%sx%s+%s+%s" % (xsize, ysize, int(x), int(y)))

        # Create a dialog frame within the new root window 
        dialogframe = Frame(self.root, bg=bgcolor, height=100)
        label1 = Label(dialogframe, text=modelname, bg=bgcolor)
        label2 = Label(dialogframe, text='Fixed', bg=bgcolor)
        seperator = ttk.Separator(dialogframe)

        dialogframe.grid(column=0, row=0)
        label1.grid(row=0, column=0, columnspan=2, sticky=W)
        label2.grid(row=0, column=2, columnspan=2, sticky=E)
        seperator.grid(row=1, columnspan=4, sticky=E+W, pady=(4,14))

        # initialize values
        self.values = [DoubleVar(dialogframe, value=param['default']) for param in params]
        self.states = [BooleanVar(dialogframe, value=param['freeze']) for param in params]
        
        row = 2
        for i in range(nparams):
            label1 = Label(dialogframe, text=params[i]['name'], bg=bgcolor)
            entry = ttk.Entry(dialogframe, textvariable=self.values[i], width=7)
            label2 = Label(dialogframe, text=params[i]['units'], bg=bgcolor)
            check = Checkbutton(dialogframe, variable=self.states[i], bg=bgcolor)

            label1.grid(row=row+i, column=0, sticky=W, pady=2)
            entry.grid(row=row+i, column=1, pady=2)
            label2.grid(row=row+i, column=2, pady=2)
            check.grid(row=row+i, column=3, sticky=E, padx=(20,0), pady=2)

        row += nparams    
        
        # accept and cancel buttons
        acceptButton = Button(dialogframe, text=u"Accept", relief=FLAT, highlightbackground=bgcolor,
                              command=self.set_params, width=5)
        acceptButton.grid(row=row, columnspan=4, sticky=N + W + S + E, pady=(12, 0))
        cancelButton = Button(dialogframe, text=u"Cancel", relief=FLAT, highlightbackground=bgcolor,
                              command=self.root.destroy, width=5)
        cancelButton.grid(row=row + 1, columnspan=4, sticky=N + W + S + E, pady=0)



        self.root.grid_columnconfigure(0, weight=1)
        self.root.grid_rowconfigure(0, weight=1)

        self.root.lift()
        #self.root.mainloop()
                
    def set_params(self):
        modelname = self.photonModelListBox.get(ACTIVE)
        params = self.xspec.models.parameter_settings(modelname)
        nparams = len(params)
        for i in range(nparams):
            params[i]['default'] = self.values[i].get()
            params[i]['freeze'] = self.states[i].get()
        self.xspec.models.update_parameter_settings(modelname, params)
        self.root.destroy()
        


##########################################################################################

class FitPlotter(Frame):
    """Class for displaying pha files"""

    # def __init__(self, root, filename):
    def __init__(self, gspecManager, data):

        # Create an instance variable to store the reference to the gspecManager 
        self.gspecManager = gspecManager

        self.data = data

        # Create the tkinter interface
        self.createTKinterInterface()
        
        self.gspecManager._openWindows['./Spectral Fit Display'] = self

    ############################################

    def createTKinterInterface(self):

        # Create the root window
        root = Tk()
        root.title('Fit Display')
        #root.geometry("1150x825+850+300") # this geometry is way too big for a single screen laptop
        root.geometry("765x640+850+0")
        root.minsize(800, 600)
        # root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)

        # Make a reference to the root window so that we can kill it later
        self.root = root

        # Initilize a frame instance with the root as its parent
        Frame.__init__(self, root)

        self.config(bg="#e1e1e1")
        self.grid(sticky=N + W + E + S)

        # Create the menu bar
        menubar = Menu(root)

        fileMenu = Menu(menubar)
        fileMenu.add_cascade(label="Print Fit Info on Plot", command=self.command, state=DISABLED)
        fileMenu.add_command(label="Screenshot", command=self.command)
        fileMenu.add_command(label="Dismiss", command=self.onWindowClose)

        datasetMenu = Menu(menubar)
        datasetMenu.add_command(label="All", command=self.command, state=DISABLED)
        datasetMenu.add_command(label="<filename>", command=self.command, state=DISABLED)

        optionsMenu = Menu(menubar)
        optionsMenu.add_command(label="Colors", command=self.command, state=DISABLED)
        optionsMenu.add_command(label="Plot Configuration", command=self.command, state=DISABLED)
        optionsMenu.add_command(label="Set Fit Repeat Count", command=self.command, state=DISABLED)
        optionsMenu.add_command(label="Set Fluence Energies", command=self.command, state=DISABLED)

        # exportMenu = Menu(menubar)
        # exportMenu.add_command(label="Selected Data to XSpec Format...", command=self.command, state=DISABLED)
        # exportMenu.add_command(label="Fit Results to File...", command=self.command, state=DISABLED)

        menubar.add_cascade(label="File", menu=fileMenu)
        menubar.add_cascade(label="Dataset", menu=datasetMenu)
        menubar.add_cascade(label="Options", menu=optionsMenu)
        # menubar.add_cascade(label="Export", menu=exportMenu)

        root.config(menu=menubar)

        ########### Buttons ###########

        frame_controls = Frame(self, bg="#e1e1e1", relief=FLAT, pady=10, padx=0, height=100)

        # Fit results menu
        self.fitResultsSelection = StringVar(frame_controls)
        self.fitResultsSelection.set('Fit Results')
        choices = ['Write Results to File', 'Write Fit Params Only', 'Read Fit Results File']
        self.fitResultsMenu = OptionMenu(frame_controls, self.fitResultsSelection, *choices,
                                         command=self.commandSelection)
        self.fitResultsMenu.config(bg="#e1e1e1", width=16, justify=CENTER)
        self.fitResultsMenu['state'] = 'disabled'

        # Redo fit button
        self.redoLastFitButton = Button(frame_controls, text=u"Redo Last Fit", relief=FLAT,
                                        highlightbackground="#e1e1e1", 
                                        command=self.gspecManager.perform_spectral_fit, width=14)

        # Spectral fitting menu
        self.spectralFittingSelectionTK = StringVar(frame_controls)
        self.spectralFittingSelectionTK.set('Spectral Fitting')
        choices = ['Synthesize Burst', 'Error Interval', 'ChiSqr 1D Plot', 'ChiSqr 2D Plot', 'Batch Fit Selections']
        self.spectralFittingMenu = OptionMenu(frame_controls, self.spectralFittingSelectionTK, *choices,
                                              command=self.commandSelection)
        self.spectralFittingMenu.config(bg="#e1e1e1", width=16, justify=CENTER)
        self.spectralFittingMenu['state'] = 'disabled'

        # Adjust source selection menu
        self.spectrumSelectionTK = StringVar(frame_controls)
        self.spectrumSelectionTK.set('Fit Display Options')
        choices = ['Cumulative', 'Raw Counts', 'Counts Spectrum', 'Photon Spectrum', 'Energy Spectrum',
                   'Nu Fnu Spectrum']
        self.fitDisplayOptionseMenu = OptionMenu(frame_controls, self.spectrumSelectionTK, *choices,
                                                 command=self.onFitDisplaySelection)
        self.fitDisplayOptionseMenu.config(bg="#e1e1e1", width=16, justify=CENTER)

        # Pack up the buttons
        self.fitResultsMenu.pack(padx=0, side=TOP)
        self.redoLastFitButton.pack(padx=0, side=TOP)
        self.spectralFittingMenu.pack(padx=0, side=TOP)
        self.fitDisplayOptionseMenu.pack(padx=0, side=TOP)

        ########### Check boxes ###########

        # Create the frame that will contain the toggle buttons and associated labels
        self.frame_toggleOptions = Frame(self, bg="#e1e1e1", relief=FLAT, pady=0, padx=0)

        # # Create the count display options text label
        label_countDisplayOptions = Label(self.frame_toggleOptions, text="Spectral Data Options:", bg="#e1e1e1")
        label_countDisplayOptions.config(justify=CENTER)
        label_countDisplayOptions.grid(column=0, row=0, columnspan=2, sticky=N + W, pady=(8, 0), padx=5)

        # Create the count display options checkboxes
        self.showCountDisplayTK = BooleanVar(self.frame_toggleOptions)
        self.checkbutton_showCountDisplay = Radiobutton(self.frame_toggleOptions, bg="#e1e1e1", text='Show', value=True,
                                                        variable=self.showCountDisplayTK,
                                                        command=self.onCountDisplaySelection)
        self.checkbutton_hideCountDisplay = Radiobutton(self.frame_toggleOptions, bg="#e1e1e1", text='Hide', value=False,
                                                        variable=self.showCountDisplayTK,
                                                        command=self.onCountDisplaySelection)
        self.checkbutton_showCountDisplay.grid(column=0, row=1, sticky=W + N, padx=(15, 0))
        self.checkbutton_hideCountDisplay.grid(column=1, row=1, sticky=W + N)

        # # Create the count display options text label
        label_modelDisplayOptions = Label(self.frame_toggleOptions, text="Spectral Model Options:", bg="#e1e1e1")
        label_modelDisplayOptions.config(justify=CENTER)
        label_modelDisplayOptions.grid(column=0, row=2, columnspan=2, sticky=N + W, pady=(8, 0), padx=5)

        # Create the count display options checkboxes
        self.showModelDisplayTK = BooleanVar(self.frame_toggleOptions)
        self.checkbutton_showModelDisplay = Radiobutton(self.frame_toggleOptions, bg="#e1e1e1", text='Show', value=True,
                                                        variable=self.showModelDisplayTK,
                                                        command=self.onModelDisplaySelection)
        self.checkbutton_hideModelDisplay = Radiobutton(self.frame_toggleOptions, bg="#e1e1e1", text='Hide', value=False,
                                                        variable=self.showModelDisplayTK,
                                                        command=self.onModelDisplaySelection)
        self.checkbutton_showModelDisplay.grid(column=0, row=3, sticky=W + N, padx=(15, 0))
        self.checkbutton_hideModelDisplay.grid(column=1, row=3, sticky=W + N)

        # # Create the spectral model components text label
        label_spectralModelComponents = Label(self.frame_toggleOptions, text='Spectral Model Components:', bg="#e1e1e1")
        label_spectralModelComponents.config(justify=CENTER)
        label_spectralModelComponents.grid(column=0, row=4, columnspan=2, sticky=N + W, pady=(8, 0), padx=5)

        # Create the spectral model components checkboxes
        self.showSpectralModelComponentsTK = BooleanVar(self.frame_toggleOptions)
        self.checkbutton_showSpectralModelComponents = Radiobutton(self.frame_toggleOptions, bg="#e1e1e1", text='Show',
                                                                   value=True,
                                                                   variable=self.showSpectralModelComponentsTK,
                                                                   command=self.onModelComponentsDisplaySelection)
        self.checkbutton_hideSpectralModelComponents = Radiobutton(self.frame_toggleOptions, bg="#e1e1e1", text='Hide',
                                                                   value=False,
                                                                   variable=self.showSpectralModelComponentsTK,
                                                                   command=self.onModelComponentsDisplaySelection)
        self.checkbutton_showSpectralModelComponents.grid(column=0, row=5, sticky=W + N, padx=(15, 0))
        self.checkbutton_hideSpectralModelComponents.grid(column=1, row=5, sticky=W + N)

        # # Create the text label   
        # label_modelColorOptions = Label(frame_toggleOptions, text='Model Color Options:', bg="#e1e1e1")
        # label_modelColorOptions.config(justify=CENTER)
        # label_modelColorOptions.grid(column=0, row=4, columnspan=2, sticky=N+W, pady=(8,0), padx=5) 

        # self.defaultModelColors = BooleanVar(frame_toggleOptions)
        # self.checkbutton_defaultModelColors = Checkbutton(frame_toggleOptions, bg="#e1e1e1", text='Show', variable=self.defaultModelColors, command=self.command)
        # self.checkbutton_defaultModelColors.grid(column=0, row=5, sticky=W+N, padx=(15,0))

        # self.matchingModelColorOptions = BooleanVar(root)     
        # self.checkbutton_matchingModelColor = Checkbutton(frame_toggleOptions, bg="#e1e1e1", text='Hide', variable=self.matchingModelColorOptions, command=self.command)
        # self.checkbutton_matchingModelColor.grid(column=1, row=5, sticky=W+N)

        # Create the text label   
        self.label_residualDisplayOptions = Label(self.frame_toggleOptions, text='Residual Display Options:', bg="#e1e1e1")
        self.label_residualDisplayOptions.config(justify=CENTER)
        self.label_residualDisplayOptions.grid(column=0, row=6, columnspan=2, sticky=N + W, pady=(8, 0), padx=5)

        self.residualSelectionTK = StringVar(self.frame_toggleOptions)
        self.checkbutton_SigmaResiduals = Radiobutton(self.frame_toggleOptions, bg="#e1e1e1", text='Sigma Residuals',
                                                      value='Sigma Residuals', variable=self.residualSelectionTK,
                                                      command=self.onResidualSelection)
        self.checkbutton_countResiduals = Radiobutton(self.frame_toggleOptions, bg="#e1e1e1", text='Count Residuals',
                                                      value='Count Residuals', variable=self.residualSelectionTK,
                                                      command=self.onResidualSelection)
        self.checkbutton_noResiduals = Radiobutton(self.frame_toggleOptions, bg="#e1e1e1", text='No Residuals',
                                                   value='No Residuals', variable=self.residualSelectionTK,
                                                   command=self.onResidualSelection)

        self.checkbutton_SigmaResiduals.grid(column=0, row=7, columnspan=2, sticky=W + N, padx=(15, 0))
        self.checkbutton_countResiduals.grid(column=0, row=8, columnspan=2, sticky=W + N, padx=(15, 0))
        self.checkbutton_noResiduals.grid(column=0, row=9, columnspan=2, sticky=W + N, padx=(15, 0))

       # Create the frame that will contain the toggle buttons and associated labels
        # self.frame_toggleOption_lower = Frame(self, bg="#e1e1e1", relief=FLAT, pady=0, padx=0)

        # self.xlog = BooleanVar()     
        self.checkbutton_xlog = Checkbutton(self.frame_toggleOptions, bg="#e1e1e1", text='X Log', command=self.setXLogScale)
        self.checkbutton_xlog.grid(column=0, row=10, sticky=W + N, pady=(8, 0), padx=(15, 0))

        # self.ylog = BooleanVar()     
        self.checkbutton_ylog = Checkbutton(self.frame_toggleOptions, bg="#e1e1e1", text='Y Log', command=self.setYLogScale)
        self.checkbutton_ylog.grid(column=1, row=10, sticky=W + N, pady=(8, 0))

        # Set the default values
        self.showCountDisplay = True
        self.showModelDisplay = True
        self.showSpectralModelComponents = True
        self.residualSelection = 'No Residuals'
        self.spectrumSelectionTK.set('Counts Spectrum')
        self.spectrumSelection = 'Counts Spectrum'
        self.checkbutton_showCountDisplay.select()
        self.checkbutton_showModelDisplay.select()
        self.checkbutton_hideSpectralModelComponents.select()
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
        # frame_filaname  = Frame(self, bg="#e1e1e1", relief=FLAT, pady=0, padx=0)
        # label_filename = Label(frame_filaname, text='File:', bg="#e1e1e1")
        # label_filename.config(justify=LEFT)
        # label_filename.grid(column=0, row=1, columnspan=3, sticky=W, pady=(8,0), padx=5) 

        ########### Plot Window ###########

        # Create a tk frame that will contain the plotting window
        self.plot_frame = Frame(root, bg="#e1e1e1")

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
        toolbar_frame = Frame(self.plot_frame, bg="#e1e1e1")
        toolbar_frame.grid(column=0, row=1, sticky=N + W)
        toolbar_frame.config(bg="#e1e1e1")
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, toolbar_frame)
        self.toolbar.config(bg="#e1e1e1")

        # Change the message label and button background colors
        self.toolbar._message_label.config(background="#e1e1e1")
        for button in self.toolbar.winfo_children():
            button.config(background="#e1e1e1")

        self.toolbar.update()
        self.canvas.get_tk_widget().grid(column=0, row=0, padx=0, pady=0, sticky=N + W + S + E)

        self.plot_frame.rowconfigure(0, weight=1)
        # self.plot_frame.rowconfigure(1, weight=0)
        self.plot_frame.columnconfigure(0, weight=1)
        # self.plot_frame.columnconfigure(1, weight=0)

        # have to update the the command *after* the toolbar has been instantiated
        fileMenu.entryconfig(1, command=self.toolbar.save_figure)

        
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

        root.rowconfigure(0, weight=1)
        root.columnconfigure(0, weight=0)
        root.columnconfigure(1, weight=2)
        
        root.protocol("WM_DELETE_WINDOW", self.onWindowClose)

        if data is not None and self.data['batch'] == False:
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

        # Redraw the current spectrum
        self.plotSpectrum()

    ############################################

    def onModelDisplaySelection(self):

        self.showModelDisplay = self.showModelDisplayTK.get()

        # Redraw the current spectrum
        self.plotSpectrum()

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
        print('function not yet implmented')

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
        
        for child in self.frame_toggleOptions.winfo_children():
            child.configure(state='normal')


        # Get the number of datagroups
        numberOfDataGroups = len(self.data['x_ufspec'])

        # Define a list of plotting colors to cycle through
        colors = ['#394264', '#39645b', '#645b39', '#586439', '#644639', '#5b3964', '#39645b', '#644639', '#395764',
                  '#426439']  # https://www.colorhexa.com/394264
        colors_model = ['#8b0000', '#720000', '#580000', '#3f0000', '#a50000', '#be0000', '#d80000', '#9f0000',
                        '#c60000', '#d90000']  # https://www.colorhexa.com/8b0000

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
                self.ax2.set_xlabel('Time (s)')

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
        energy_summed_model = self.data['energies'][0:-1]  # this bin selection needs to be re-examined

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
        flux_model = self.data['model_ufspec']

        # Get the summer model photon flux
        flux_summed_model = self.data['model']

        # Get the sigma residuals
        residuals_sigma = self.data['residuals_sigma']

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

        # Determine what kind of data to display
        if self.spectrumSelection == 'Cumulative':
            flux_data = counts_integrated
            flux_error_data = numpy.zeros(len(counts_integrated))
            flux_model = counts_integrated_model
            residuals_data = flux_data - flux_model

            self.ax.set_ylabel(r'Normalized integrated counts')

        if self.spectrumSelection == 'Raw Counts':
            flux_data = counts
            flux_error_data = counts_error
            flux_model = counts_model
            residuals_data = flux_data - flux_model

            # energy_data = numpy.arange(len(counts))
            # energy_error_data = numpy.ones(len(counts))

            self.ax.set_ylabel(r'Counts bin$^{-1}$')

        # Prepare the count spectrum
        if self.spectrumSelection == 'Counts Spectrum':
            flux_data = counts/energy_error_data
            flux_error_data = (counts_error / counts) * flux_data
            flux_model = counts_model/energy_error_data
            residuals_data = flux_data - flux_model

            self.ax.set_ylabel(r'Counts s$^{-1}$ keV$^{-1}$')

        # Prepare the flux density spectrum (Photons cm-2 s-1 keV-1)
        if self.spectrumSelection == 'Photon Spectrum':
            residuals_data = flux_data - flux_model

            self.ax.set_ylabel(r'Flux (Photons cm$^{-2}$ s$^{-1}$ keV$^{-1}$)')

        # Prepare the Fv spectrum (Photons cm-2 s-1)
        if self.spectrumSelection == 'Energy Spectrum':
            flux_error_data = (flux_error_data / flux_data) * (flux_data * energy_data)  # photons cm-2
            flux_data = flux_data * energy_data
            flux_model = flux_model * energy_data
            flux_summed_model = flux_summed_model * energy_summed_model
            residuals_data = flux_data - flux_model

            self.ax.set_ylabel(r'Energy (Photons cm$^{-2}$ s$^{-1}$)')

        # Prepare the vFv spectrum (Photons keV cm-2 s-1)
        if self.spectrumSelection == 'Nu Fnu Spectrum':
            flux_error_data = (flux_error_data / flux_data) * (flux_data * (energy_data ** 2))  # photons keV
            flux_data = flux_data * (energy_data ** 2)
            flux_model = flux_model * (energy_data ** 2)
            flux_summed_model = flux_summed_model * (energy_summed_model ** 2)
            residuals_data = flux_data - flux_model

            self.ax.set_ylabel(r'$\nu$ F$_{\nu}$ (Photons keV cm$^{-2}$ s$^{-1}$)')


        # Determine which data points have error bars that go negative
        index_detections = []
        index_upperLimits = []
        for dataGroup in range(numberOfDataGroups):
            index_detections.append(numpy.where((flux_data[dataGroup]-flux_error_data[dataGroup]) > 0))
            index_upperLimits.append(numpy.where((flux_data[dataGroup]-flux_error_data[dataGroup]) <= 0))

        # Plot the data points
        for dataGroup, detections, upperLimits in zip(range(numberOfDataGroups), index_detections, index_upperLimits):

            self.ax.errorbar(energy_data[dataGroup][detections], flux_data[dataGroup][detections], xerr=energy_error_data[dataGroup][detections], capsize=0,
                             fmt='none', zorder=1, ecolor=colors[dataGroup], alpha=alpha_data, yerr=flux_error_data[dataGroup][detections])

            self.ax.errorbar(energy_data[dataGroup][upperLimits], flux_data[dataGroup][upperLimits], xerr=energy_error_data[dataGroup][upperLimits], capsize=0,
                             fmt='none', zorder=1, ecolor=colors[dataGroup], alpha=alpha_data, yerr=flux_data[dataGroup][upperLimits]/2.0, uplims=True)



        # Plot the residuals
        for dataGroup in range(numberOfDataGroups):
            if 'Count Residuals' in self.residualSelection and self.ax2 is not None:
                self.ax2.plot(energy_data[dataGroup], numpy.zeros(len(energy_data[dataGroup])), linestyle='--',
                              c=colors_model[dataGroup], zorder=2, alpha=alpha_model)
                self.ax2.errorbar(energy_data[dataGroup], residuals_data[dataGroup], xerr=energy_error_data[dataGroup],
                                  capsize=0, fmt='none', ecolor=colors[dataGroup], zorder=1, alpha=alpha_data, yerr=flux_error_data[dataGroup])

            if 'Sigma Residuals' in self.residualSelection and self.ax2 is not None:
                self.ax2.plot(energy_data[dataGroup], numpy.zeros(len(energy_data[dataGroup])), linestyle='--',
                              c=colors_model[dataGroup], zorder=2, alpha=alpha_model)
                self.ax2.errorbar(energy_data[dataGroup], residuals_sigma[dataGroup], xerr=energy_error_data[dataGroup],
                                  capsize=0, fmt='none', ecolor=colors[dataGroup], zorder=1, alpha=alpha_data, yerr=numpy.ones(len(residuals_sigma[dataGroup])))

                self.ax2.set_ylabel('Sigma')

        # Set the x-axis labels
        self.ax.set_xlabel('Energy (keV)')

        # Set the x-axis minor tick frequency
        minorLocator = AutoMinorLocator()
        self.ax.xaxis.set_minor_locator(minorLocator)

        # Set the y-axis minor tick frequency
        minorLocator = AutoMinorLocator()
        self.ax.yaxis.set_minor_locator(minorLocator)

        # Set the x-axis range to match the data
        xmin = 1e99
        xmax = 0
        for dataGroup in range(numberOfDataGroups):
            xmin_group = numpy.min(energy_data[dataGroup])
            xmax_group = numpy.min(energy_data[dataGroup])
            if xmin_group < xmin:
                xmin = xmin_group
            if xmax_group > xmax:
                xmax = xmax_group


        self.ax.set_xlim(numpy.min(energy_data[0]), numpy.max(energy_data[-1]))

        # Set the plot scale
        if self.xlog == True:
            self.ax.set_xscale('log')

            if self.ax2 is not None:
                self.ax2.set_xscale('log')

        if self.ylog == True:
            self.ax.set_yscale('log')

        # Set a minimum y value
        dflux = flux_data-flux_error_data
        for dataGroup in range(numberOfDataGroups):
            if ((len(numpy.where(dflux[dataGroup] <= 0)[0]) > 0) == True) and self.spectrumSelection == 'Counts Spectrum':
                self.ax.set_ylim(bottom=0.1)

        # Don't let the model affect the plot scale
        self.ax.set_autoscale_on(False)

        # Plot the individual detector model
        if self.spectrumSelection == 'Raw Counts' or self.spectrumSelection == 'Counts Spectrum' or self.spectrumSelection == 'Cumulative':

            for dataGroup in range(numberOfDataGroups):
                for i in range(len(energy_data[dataGroup])):

                    # Plot the individual model values for the count and cumulutive spectrum
                    energy_low = energy_data[dataGroup][i] - energy_error_data[dataGroup][i]
                    energy_high = energy_data[dataGroup][i] + energy_error_data[dataGroup][i]

                    self.ax.plot([energy_low, energy_high], [flux_model[dataGroup][i], flux_model[dataGroup][i]],
                                 c=colors_model[dataGroup], zorder=2, alpha=alpha_model)

                    try:
                        self.ax.plot([energy_high, energy_high],
                                     [flux_model[dataGroup][i], flux_model[dataGroup][i + 1]],
                                     c=colors_model[dataGroup], zorder=2, alpha=alpha_model)
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
            self.ax.plot(energy_summed_model, flux_summed_model, c=colors_model[dataGroup], zorder=2, alpha=alpha_model)

        # Hide the last tick label on the second axes, if it exists
        if self.ax2 is not None:
            yticks = self.ax2.get_yticks()
            yticks = yticks[1:-1]
            self.ax2.set_yticks(yticks)

        self.canvas.draw()


    ############################################

    def plotBatchFitResults(self):

        # Turn off log scaling by default
        if self.xlog == False:
            self.checkbutton_xlog.deselect()

        if self.ylog == False:
            self.checkbutton_ylog.deselect()

        # Clear any existing plotting data
        self.ax.clear()

        # Hide the widgets that we don't need for batch fitting
        self.fitDisplayOptionseMenu.configure(state='disable')

        for child in self.frame_toggleOptions.winfo_children():
            text = child.cget("text")
            if 'X Log' not in text and 'Y Log' not in text:
                child.configure(state='disable')

        # Extract the extension information
        time_start = self.data['time_start']
        time_end = self.data['time_end']
        time_exposure = self.data['exposure']

        time = time_start + (time_end-time_start)/2.0
        dtime = (time_end-time_start)/2.0

        # Extract the fit results
        parameter_names = self.data['parameter_names']
        parameter_values = self.data['parameter_values']
        parameter_sigmas = self.data['parameter_sigmas']

        # print("parameter_values:")
        # print(parameter_values)

        # print("parameter_sigmas:")
        # print(parameter_sigmas)

        parameter_name = parameter_names[0]
        parameter_value = parameter_values[:,0]
        parameter_sigmas = parameter_sigmas[:,0]

        # print("parameter_value:")
        # print(parameter_value)

        # print("parameter_sigma:")
        # print(parameter_sigmas)


        # Plot the parameter history
        self.ax.errorbar(time, parameter_value, xerr=dtime, yerr=parameter_sigmas, capsize=0, fmt='none', zorder=1, ecolor='#394264', alpha=1)

        # Set the x-axis labels
        self.ax.set_xlabel('Time (sec)')

        # Set the x-axis labels
        self.ax.set_ylabel(parameter_name)


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

        self.canvas.draw()

    ############################################

    def quitConfirmation(self):
        exit()

    ############################################

    def onWindowClose(self):
        try:
            del self.gspecManager._openWindows['./Spectral Fit Display']
        except:
            pass
        self.root.destroy()


##########################################################################################

class AboutGspec(Frame):
    
    def __init__(self, parent):
        self.parent = parent
        
    # def show(self):
    #     self.root = Toplevel(self.parent)
    #     bgcolor = "#e1e1e1"
    #     self.root.config(bg=bgcolor)
    #     # remove standard border
    #     self.root.overrideredirect(True)
    #     # put on top
    #     self.root.wm_attributes("-topmost", True)
    #     # Turn off the window shadow
    #     self.root.wm_attributes("-transparent", True)
        
    #     # determine size
    #     xsize = 300
    #     ysize = 200
    #     # Center the root box
    #     width = self.root.winfo_screenwidth()
    #     height = self.root.winfo_screenheight()
    #     size = [250, 250]
    #     x = width / 2 - size[0] / 2
    #     y = height / 2 - size[1] / 2
    #     self.root.geometry("%sx%s+%s+%s" % (xsize, ysize, int(x), int(y)))
    #     txt = '\n'.join(['GSpec - GBM Spectral Analysis Software', 
    #                      'Version 0.3',
    #                      '', '', 
    #                      'Courtesy of:',
    #                      'Unversities Space Research Association', 
    #                      'NASA', 
    #                      'University of Alabama in Huntsville', ''])
    #     label = Label(self.root, text=txt, bg = bgcolor)
    #     #self.root.image = PhotoImage(file="fermi.gif")
    #     #self.root.image = self.root.image.subsample(2, 2)
    #     #button = Button(self.root, image=self.root.image)
    #     #label.config(bg='systemTransparent')
    #     #button.pack()
    #     dismiss = Button(self.root, text=u"Dismiss", relief=FLAT, height=50,
    #                      highlightbackground=bgcolor, command=self.root.destroy)
    #     label.pack()
    #     dismiss.pack(fill=BOTH)        
        
    #     self.root.bind('<Return>', self.root.destroy)
    #     self.root.lift()


    def show(self):

        self.root = Toplevel(self.parent)
        self.root.config(bg="#e1e1e1")

        # remove standard border
        # self.root.overrideredirect(True)

        # put on top
        self.root.wm_attributes("-topmost", True)

        # Turn off the window shadow
        # self.root.wm_attributes("-transparent", True)


        xsize = 300
        ysize = 470

        # Center the root box
        width = self.root.winfo_screenwidth()
        height = self.root.winfo_screenheight()
        size = [250, 400]
        x = width / 2 - size[0] / 2
        y = height / 2 - size[1] / 2
        
        self.root.geometry("%sx%s+%s+%s" % (xsize, ysize, int(x), int(y)))
        self.root.config(bg="#e1e1e1")

        canvas = Canvas(self.root, height=90, highlightthickness=0, bg="#e1e1e1")   
        canvas.pack(pady=(20,0))                 
        self.image = PhotoImage(file="spaceship_75x75.gif")      
        canvas.create_image(xsize/2 - 9, 50, anchor=CENTER, image=self.image)  

        # self.canvas = Canvas(self, width = 300, height = 300, highlightthickness=2, bg="#e1e1e1")  
        # self.canvas.pack(pady=(20,0)) 
        # self.image = ImageTk.PhotoImage(Image.open("spaceship2.png"))  
        # self.canvas.create_image(20, 20, anchor=NW, image=self.image)
        # self.canvas.image = self.image   


        lines = ['\nGSpec',
                         'GBM Spectral Analysis Package', 
                         'version 0.4\n',
                         'Adam Goldstein',
                         'Daniel Kocevski',
                         'Rob Preece',
                         'William Cleveland\n',
                         # 'Courtesy of:',
                         'Unversities Space Research Association', 
                         'NASA Marshall Space Flight Center', 
                         'University of Alabama in Huntsville', '']

        for line in lines:
            label = Label(self.root, text=line, bg="#e1e1e1", pady=2)
            label.pack()

        # Create a dismiss button
        dismiss = Button(self.root, text=u"Dismiss", relief=FLAT, width=10,
                         highlightbackground="#e1e1e1", command=self.root.destroy)
        dismiss.pack()   

        self.root.bind('<Return>', self.root.destroy)
        self.root.lift()

    def toggle(self):
        if self.state == 'fermi':
            self.root.image = PhotoImage(file="bns-small.png")
            self.button.config(image=self.root.image)
            self.state = 'bns'
            return
        elif self.state == 'bns':
            self.root.image = PhotoImage(file="fermi-small.png")
            self.button.config(image=self.root.image)
            self.state = 'fermi'
            return
        
            



##########################################################################################

def create_root_window(gspecManager):
    print(gspecManager)

    # from tkinter import *
    root = Tk()

    root.title('GSPEC v0.4')

    # root.config(bg="#e1e1e1")
    root.geometry("600x250+300+300")
    root.minsize(600, 250)

    # Application(root)
    Application(root, gspecManager)

    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)

    root.mainloop()


##########################################################################################

class CustomAxes(matplotlib.axes.Axes):
    name = "CustomAxes"

    def drag_pan(self, button, key, x, y):
        previous_ylim = self.get_ylim()
        matplotlib.axes.Axes.drag_pan(self, button, key, x, y)  # pretend key=='x'
        current_ylim = self.get_ylim()
        if current_ylim[0] < -500.0:
            self.set_ylim(previous_ylim)
        # print(self)


##########################################################################################

class ProgressBar(Frame):
    """Class for displaying a custom progress bar"""

    def __init__(self, maximum, value):

        # Create a new root window
        bgcolor = "#e1e1e1"
        self.root = Tk()
        self.root.title('Progress')
        self.root.config(bg=bgcolor)

        xsize = 250
        ysize = 80

        # Center the root box
        width = self.root.winfo_screenwidth()
        height = self.root.winfo_screenheight()
        size = [250, 250]
        x = width / 2 - size[0] / 2
        y = height / 2 - size[1] / 2
        # self.root_dialog.geometry("275x175+%s+%s" % (int(x), int(y)))
        self.root.geometry("%sx%s+%s+%s" % (xsize, ysize, int(x), int(y)))

        self.label = Label(self.root, text="Fitting Progress", bg=bgcolor)
        self.progressbar = ttk.Progressbar(self.root, orient="horizontal", length=250, mode ="determinate")
        # self.progressbar.configure(bg=bgcolor)
        self.progressbar['maximum'] = maximum

        self.label.pack(pady=(10,5))
        self.progressbar.pack(padx=10, pady=(5,15))

        self.root.protocol("WM_DELETE_WINDOW", self.onWindowClose)

        self.root.lift()


    def reset(self):    
        self.progressbar['value'] = 0

    def updateValue(self, value):
        self.progressbar['value'] = value

    def setMaximum(self, maximum):
        self.progressbar['maximum'] = maximum

    def lift(self):
        self.root.lift()

    def onWindowClose(self):
        self.root.destroy()

    def destroy(self):
        self.root.destroy()




