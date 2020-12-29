from tkinter import ttk, Tk, Text, Frame, Menu, Listbox, Button, filedialog, messagebox, Canvas, PhotoImage, Checkbutton
from tkinter import Scrollbar, Toplevel, DoubleVar, StringVar, BooleanVar, Label, Radiobutton
from tkinter import N, E, W, S, DISABLED, BROWSE, EXTENDED, SUNKEN, LEFT, END, ACTIVE, NORMAL, INSERT, TOP, BOTH, CENTER, X, MULTIPLE, ANCHOR
import os

from .globals import BGCOLOR, FRAMERELIEF, OM_JUSTIFY, GUIFONTCOLOR
from .geometry import GeometryManager

from gbm.xspec.model_manager import XspecModelManager
import gbm

class Dialog(object):
    width = 250
    height = 250
    def __init__(self, title, master=None, width=None, height=None):
        self.master = None
        if width is not None:
            self.width = width
        if height is not None:
            self.height = height
        self._geometry = GeometryManager()

        if master is None:
            self.root = Tk()
            self.xpos = 0
            self.ypos = 0
        else:
            self.root = Toplevel(master)
            self.master = master
            self.xpos = self._geometry.child_center_x(self.master, self.width)
            self.ypos = self._geometry.child_center_y(self.master, self.height)
        
        #self.root.geometry("{0}x{1}+{2}+{3}".format(self.width, self.height, 
        #                                            self.xpos, self.ypos))
        self.root.config(bg=BGCOLOR)
        self.root.title(title)

class Splash(Dialog):
    dir = os.path.dirname(os.path.realpath(__file__))
    _IMAGE_PATH = os.path.join(dir, "spaceship_75x75.gif")
    
    def __init__(self, master):
        width = 300
        height = 435
        super(Splash, self).__init__('', master=master, width=width, height=height)
        #Toplevel.__init__(self, parent)

        xpos = self._geometry.center_x(self.width)
        ypos = self._geometry.center_y(self.height)
        # Set the window size
        self.root.geometry("{0}x{1}+{2}+{3}".format(self.width, self.height, 
                                                    xpos, ypos))
        self.root.config(bg=BGCOLOR)

        # Remove the window border
        self.root.overrideredirect(True)
        self.root.wm_attributes("-topmost", True)
        
        canvas = Canvas(self.root, height=90, highlightthickness=0, bg=BGCOLOR)   
        canvas.pack(pady=(20,0))                 
        self.image = PhotoImage(file=self._IMAGE_PATH)      
        canvas.create_image(self.width/2 - 9, 50, anchor=CENTER, image=self.image)  

        # Setup the text
        lines = ['\nGSpec',
                 'GBM Spectral Analysis Package', 
                 'Version {:}\n'.format(gbm.__version__),
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
            label = Label(self.root, text=line, bg=BGCOLOR, pady=2)
            label.pack() 

        # Bind the any key or mouse click to the destroy function
        self.root.bind('<Key>', self.dismiss)
        self.root.bind('<Button-1>', self.dismiss)

        ## required to make window show before the program gets to the mainloop
        #self.update()


    def dismiss(self, event):

        # Close the window
        self.root.destroy()

        
class TextDisplayDialog(Dialog):
    def __init__(self, title, thetext, master=None, **kwargs):
        super(TextDisplayDialog, self).__init__(title, master=master)
        
        text = Text(self.root, height=20, width=80, bg=BGCOLOR, fg=GUIFONTCOLOR)
        text.insert(INSERT, thetext)
        text.config(state=DISABLED, **kwargs)
        text.grid()
        
        yscroll = Scrollbar(self.root, command=text.yview, bg=BGCOLOR)
        text['yscrollcommand'] = yscroll.set
        yscroll.grid(row=0, column=1, sticky=N+S+E+W)


class OptionDialog(Dialog):
    """A dialog that presents a list of options on buttons"""
    def __init__(self, title, options, commands, message=None, master=None):
        width = 160
        height = 220
        super(OptionDialog, self).__init__(title, master=master, width=width,
                                           height=height)
        
        self.root.geometry("{0}x{1}+{2}+{3}".format(self.width, self.height, 
                                                    self.xpos, self.ypos))        

        # Create a dialog frame within the new root window 
        dialogframe = Frame(self.root, bg=BGCOLOR, height=100)
        dialogframe.config(bg=BGCOLOR)
        dialogframe.grid(column=0, row=0, sticky=N+W+E+S)
        dialogframe.columnconfigure(0, weight=1)
        dialogframe.rowconfigure(0, weight=1)

        # Create the text message
        if message is not None:
            for line in message:
                Label(dialogframe, text=line, bg=BGCOLOR).pack(padx=0, side=TOP)

        # Create the option buttons
        numbuttons = len(options)
        for i in range(numbuttons):
            button = Button(dialogframe, text=options[i], command=commands[i], 
                            highlightbackground=BGCOLOR, relief=FRAMERELIEF, width=12)
            button.pack(padx=0, side=TOP)
        
        self.root.lift()


class PlotDialog(Dialog):
    def __init__(self, title, master=None, width=None, height=None):
        super(PlotDialog, self).__init__(title, master=master, width=width,
                                         height=height)

        self.root.geometry("{0}x{1}+{2}+{3}".format(self.width, self.height, 
                                                    self.xpos, self.ypos))        
        self.root.lift()


class TextOptionDialog(Dialog):
    def __init__(self, title, master=None, width=None, height=None, 
                 message=None, default_value=None, button_label=None, 
                 button_options=None, button_command=None):
        super(TextOptionDialog, self).__init__(title, master=master, width=width,
                                                height=height)

        #self.root_dialog.resizable(False, False)
        self.button_command = button_command

        self.root.geometry("{0}x{1}+{2}+{3}".format(self.width, self.height, 
                                                    self.xpos, self.ypos))

        # Make the dialog box not resizable
        self.root.resizable(False, False)

        # Expand the column to match the window size
        self.root.columnconfigure(0, weight=1)

        # Create the text label 
        row_root = 0
        for line in message:
            label = Label(self.root, text=line, bg=BGCOLOR)
            label.grid(column=0, row=row_root, columnspan=2, pady=7)
            row_root = row_root + 1

        # Create a new frame within the root window 
        entryframe = Frame(self.root, bg=BGCOLOR, height=100)
        entryframe.config(bg=BGCOLOR)
        entryframe.grid(column=0, row=row_root, pady=5)
        row_root = row_root + 1

        # Create a new row index for the entry dialog
        row_entry = 0
        
        # Create the entry form with use the specified default value
        self.default_var = DoubleVar(entryframe, value=default_value)
        entry = ttk.Entry(entryframe, width=7, textvariable=self.default_var)
        entry.grid(column=0, row=row_entry, rowspan=2)
        

        # create the radio button
        self.fullrange = StringVar(entryframe, value='Full')
        self.fullrange.set('Full')
        
        if button_options is not None:
            radio_text1 = button_options[0]
            radio_text2 = button_options[1]
        else:
            radio_text1 = 'Full Range'
            radio_text2 = 'Custom Range'
        
        radio1 = Radiobutton(entryframe, bg=BGCOLOR, text=radio_text1,
                             variable=self.fullrange, value='Full')
        radio1.grid(row=row_entry, column=1, sticky=W, padx=(20,5))
        row_entry = row_entry + 1

        radio2 = Radiobutton(entryframe, bg=BGCOLOR, text=radio_text2,
                             variable=self.fullrange, value='Custom')
        radio2.grid(row=row_entry, column=1, sticky=W, padx=(20,5))
        row_entry = row_entry + 1

        # Create the action button and bind it to the specified button command
        button = ttk.Button(self.root, text=button_label, width=5,
                   command=self._evaluate_button_command)
        button.grid(column=0, row=row_root, columnspan=2, sticky=S, pady=(10,5))

        # Bind the return key to the specified button command 
        self.root.bind('<Return>', lambda event: self._evaluate_button_command())
        self.root.lift()

        # Bind the window close action
        self.root.protocol("WM_DELETE_WINDOW", self.on_window_close)

        # self.root.grab_set()
        # self.master.wait_window(self.root)

    def _evaluate_button_command(self):
        self.button_command(self.default_var, self.fullrange)
        self.root.destroy()

    def on_window_close(self):
        self.button_command(None, None)
        self.root.destroy()



class ManualInputDialog(Dialog):
    def __init__(self, title, command, master=None, message=None, xinput=True, 
                 yinput=True, xdefaults=None, ydefaults=None):
        # determine size
        width = 165
        height = 120
        if xinput:
            height += 50
        if yinput:
            height += 50
        
        super(ManualInputDialog, self).__init__(title, master=master, width=width,
                                               height=height)
        
        self.root.geometry("{0}x{1}+{2}+{3}".format(self.width, self.height, 
                                                    self.xpos, self.ypos))

        self.button_command = command

        # Create a dialog frame within the new root window 
        dialogframe = Frame(self.root, bg=BGCOLOR, height=100, relief=FRAMERELIEF)
        dialogframe.grid(column=0, row=0)
        label = Label(dialogframe, text=message, bg=BGCOLOR)
        label.grid(row=0, columnspan=2, pady=(5, 10))

        # initialize values
        self._x0default = DoubleVar(dialogframe, value=None)
        self._x1default = DoubleVar(dialogframe, value=None)
        self._y0default = DoubleVar(dialogframe, value=None)
        self._y1default = DoubleVar(dialogframe, value=None)
        self._xinput = xinput
        self._yinput = yinput

        # x axis input
        irow = 1
        if xinput:
            if xdefaults is None:
                xdefaults = (0.0, 0.0)
            label1 = Label(dialogframe, text='X Lo: ', bg=BGCOLOR)
            label1.grid(row=irow, column=0)
            self._x0default.set(xdefaults[0])
            entry1 = ttk.Entry(dialogframe, textvariable=self._x0default, width=7)
            entry1.grid(row=irow, column=1)
            label2 = Label(dialogframe, text='X Hi: ', bg=BGCOLOR)
            label2.grid(row=irow + 1, column=0)
            self._x1default.set(xdefaults[1])
            entry2 = ttk.Entry(dialogframe, textvariable=self._x1default, width=7)
            entry2.grid(row=irow + 1, column=1)
            irow += 2

        # y axis inputs
        if yinput:
            if ydefaults is None:
                ydefaults = (0.0, 0.0)
            label1 = Label(dialogframe, text='Y Lo: ', bg=BGCOLOR)
            label1.grid(row=irow, column=0)
            self._y0default.set(ydefaults[0])
            entry1 = ttk.Entry(dialogframe, textvariable=self._y0default, width=7)
            entry1.grid(row=irow, column=1)
            label2 = Label(dialogframe, text='Y Hi: ', bg=bgcolor)
            label2.grid(row=irow + 1, column=0)
            self._y1default.set(ydefaults[1])
            entry2 = ttk.Entry(dialogframe, textvariable=self._y1default, width=7)
            entry2.grid(row=irow + 1, column=1)
            irow += 2

        # accept and cancel buttons
        acceptbutton = Button(dialogframe, text=u"Accept", 
                              highlightbackground=BGCOLOR, width=12, 
                              command=self._evaluate_button_command)
        acceptbutton.grid(row=irow, columnspan=2, sticky=N+W+S+E, pady=(10, 0))
        cancelbutton = Button(dialogframe, text=u"Cancel", 
                              highlightbackground=BGCOLOR, 
                              command=self.on_window_close, width=12)
        cancelbutton.grid(row=irow + 1, columnspan=2, sticky=N+W+S+E, pady=0)

        self.root.grid_columnconfigure(0, weight=1)
        self.root.grid_rowconfigure(0, weight=1)

        # Bind the window close action
        self.root.protocol("WM_DELETE_WINDOW", self.on_window_close)

        self.root.lift()

    def _evaluate_button_command(self):
        try:
            x0 = self._x0default.get()
            x1 = self._x1default.get()
            y0 = self._y0default.get()
            y1 = self._y1default.get()
        except:
            raise ValueError('Inputs must be numbers')
            
        if not self._yinput:
            y0 = None
            y1 = None
        if not self._xinput:
            x0 = None
            x1 = None
        self.button_command(x0, x1, y0, y1)
        self.root.destroy()

    def on_window_close(self):
        self.button_command(None, None, None, None)
        self.root.destroy()


class FitOptionsDialog(Dialog):
    def __init__(self, master, models, callback, batch=False):
        title = 'Photon Models'
        width = 3
        height = 700
        super(FitOptionsDialog, self).__init__(title, master=master, width=width,
                                           height=height)
        self.root.geometry("{0}x{1}+{2}+{3}".format(self.width, self.height, 
                                                    self.xpos, self.ypos))        
        
        # Create an instance variable to store the reference to the data list 
        self.batch = batch
        self.models = models
        self.callback = callback

        # Create the tkinter interface
        self._construct()

    def _construct(self):

        # Set the root window weightings
        self.root.rowconfigure(0, weight=1)
        self.root.columnconfigure(0, weight=1)

        # Set a minimum window size
        self.root.minsize(365, 550)

        # Create a dialog frame within the new root window 
        dialogframe = Frame(self.root, bg=BGCOLOR)
        # dialogframe.pack(anchor=N, fill=BOTH, expand=True, side=TOP)

        # Tell the frame to expand and fill the root window
        dialogframe.grid(sticky=N+W+E+S)


        # Create the text label 
        label_title = Label(dialogframe, text='Select one or more photon models', 
                            bg=BGCOLOR)
        label_title.config(justify=CENTER)

        # Create the listbox
        self.model_listbox = Listbox(dialogframe, borderwidth=2, selectbackground='#116cd6',
                                     selectborderwidth=0, highlightthickness=0,
                                     relief=SUNKEN, selectmode=EXTENDED)
        self.model_listbox.bind('<<ListboxSelect>>', self.on_select_model)
        
        # Populate the list box with available models
        models = self.models.model_names()
        for model in models:
            self.model_listbox.insert(END, model)

        # Create a dialog frame within the new root window
        frame_model_parameters = Frame(dialogframe, bg=BGCOLOR)

        # Create the text label 
        label_model_parameters = Label(frame_model_parameters, bg=BGCOLOR,
                                       text='Photon Model Parameters:')
        label_model_parameters.config(justify=CENTER)
        label_model_parameters.pack(pady=5, side=TOP)
        
        # Set parameters button
        self.set_parameters = Button(frame_model_parameters, text=u"Set Parameters", 
                                     highlightbackground=BGCOLOR, width=10,
                                     command=self.on_select_new_parameters)
        self.set_parameters['state'] = 'disabled'
        # self.set_parameters.pack(fill=X)
        self.set_parameters.pack(pady=(0,8))

        # Create a dialog frame within the new root window
        frame_fit_statistics_upper = Frame(dialogframe, bg=BGCOLOR)
        frame_fit_statistics_lower = Frame(dialogframe, bg=BGCOLOR)

        # Create the text label   
        label_fit_statistic = Label(frame_fit_statistics_upper, bg=BGCOLOR, 
                                    text='Fitting Statistic:')
        label_fit_statistic.config(justify=CENTER)
        label_fit_statistic.pack(pady=5, side=TOP)

        # Create the fit statistic variables
        self.fit_statisticTK = StringVar(frame_fit_statistics_upper)
        self.statistic_rb = {}
        texts = ['Chi2', 'C-Stat', 'L-Stat']
        vals = ['chi', 'cstat', 'lstat']
        for i in range(3):
            rb = Radiobutton(frame_fit_statistics_upper, text=texts[i], value=vals[i],   
                             bg=BGCOLOR, variable=self.fit_statisticTK)
            rb.pack(side=LEFT)
            self.statistic_rb[texts[i]] = rb
        texts = ['PG-Stat', 'P-Stat', 'Whittle']
        vals = ['pgstat', 'pstat', 'whittle']
        for i in range(3):
            rb = Radiobutton(frame_fit_statistics_lower, text=texts[i], value=vals[i],   
                             bg=BGCOLOR, variable=self.fit_statisticTK)
            rb.pack(side=LEFT)
            self.statistic_rb[texts[i]] = rb

        # Create a dialog frame within the new root window 
        frame_fit_weights = Frame(dialogframe, bg=BGCOLOR)

        # Create the text label   
        label_fit_weights = Label(frame_fit_weights, text='Fit Weighting:', bg=BGCOLOR)
        label_fit_weights.config(justify=CENTER)
        label_fit_weights.pack(pady=5, side=TOP)

        # Create the fit statistic variables
        self.fit_weightsTK = StringVar(frame_fit_weights)
        self.weight_rb = {}
        texts = ['Standard', 'Model', 'Gehrels', 'Churazov']
        vals = [t.lower() for t in texts]
        for i in range(4):
            rb = Radiobutton(frame_fit_weights, text=texts[i], value=vals[i], 
                             bg=BGCOLOR, variable=self.fit_weightsTK)
            rb.pack(side=LEFT)
            self.weight_rb[texts[i]] = rb
        
        # Undetermined Values in Batch Fit 
        frame_undetermined_values = Frame(dialogframe, bg=BGCOLOR)
        label_undetermined_values = Label(frame_undetermined_values, bg=BGCOLOR, 
                                         text='Undetermined Values in Batch Fit:')
        label_undetermined_values.config(justify=CENTER)
        label_undetermined_values.pack(pady=5, side=TOP)

        self.undetermined_values = StringVar(frame_undetermined_values)
        self.batch_rb = {}
        texts = ['Leave free', 'Automatically fix']
        for i in range(2):
            rb = Radiobutton(frame_undetermined_values, text=texts[i], value=texts[i], 
                             bg=BGCOLOR, variable=self.undetermined_values)
            rb.pack(side=LEFT)
            self.batch_rb[texts[i]] = rb
        
        # Create the button frame
        frame_buttons = Frame(dialogframe, bg=BGCOLOR)

        # Create the buttons
        self.accept_fit_parameters = Button(frame_buttons, text=u"Accept", highlightbackground=BGCOLOR,
                                          command=self.perform_spectral_fit, width=10)
        self.accept_fit_parameters['state'] = 'disabled'
        # self.omit_fit_parameters = Button(frame_buttons, text=u"Omit", highlightbackground=BGCOLOR,
                                        # command=self.command, width=5)
        # self.omit_fit_parameters['state'] = 'disabled'
        # self.restore_fit_parameters = Button(frame_buttons, text=u"Restore", highlightbackground=BGCOLOR,
        #                                    command=self.command, width=5)
        # self.restore_fit_parameters['state'] = 'disabled'
        self.cancel_fit_parameters = Button(frame_buttons, text=u"Cancel", highlightbackground=BGCOLOR,
                                          command=self.cancel_fit, width=10)
        
        # Pack up the buttons
        # self.accept_fit_parameters.pack(side=LEFT)
        # self.omit_fit_parameters.pack(side=LEFT)
        # self.restore_fit_parameters.pack(side=LEFT)
        # self.cancel_fit_parameters.pack(side=LEFT)
        self.accept_fit_parameters.grid(row=0, column=0, sticky=E+W, pady=10)
        # self.omit_fit_parameters.grid(row=0, column=1, sticky=E+W, pady=10)
        # self.restore_fit_parameters.grid(row=0, column=2, sticky=E+W, pady=10)
        self.cancel_fit_parameters.grid(row=0, column=3, sticky=E+W, pady=10)

        # Set the default values
        self.statistic_rb['PG-Stat'].select()
        self.weight_rb['Standard'].select()
        self.batch_rb['Leave free'].select()

        # Pack it all together
        # label_title.pack(pady=10, side=TOP)
        # self.model_listbox.pack(padx=10, side=TOP, fill=X)
        # frame_model_parameters.pack(pady=10, side=TOP)
        # frame_fit_statistics_upper.pack(pady=(10, 0), side=TOP)
        # frame_fit_statistics_lower.pack(pady=(0, 10), side=TOP)
        # frame_fit_weights.pack(pady=10, side=TOP)
        # frame_undetermined_values.pack(pady=10, side=TOP)
        # frame_buttons.pack(padx=4, pady=15, side=TOP)

        label_title.grid(row=0, column=0, sticky=E+W, pady=(10,0))
        self.model_listbox.grid(row=1, column=0, sticky=N+E+W+S, padx=10, pady=10)
        frame_model_parameters.grid(row=2, column=0, sticky=E+W, pady=(10,0))
        frame_fit_statistics_upper.grid(row=3, column=0, sticky=N, pady=(8,2))
        frame_fit_statistics_lower.grid(row=4, column=0, sticky=N, pady=(2,8))
        frame_fit_weights.grid(row=5, column=0, sticky=N, pady=8)
        frame_undetermined_values.grid(row=6, column=0, sticky=N, pady=8)
        frame_buttons.grid(row=7, column=0, sticky=N+S, pady=8)


        # Configure the weights of the rows and columns
        dialogframe.rowconfigure(0, weight=0)
        dialogframe.rowconfigure(1, weight=1)
        dialogframe.rowconfigure(2, weight=0)
        dialogframe.rowconfigure(3, weight=0)
        dialogframe.rowconfigure(4, weight=0)
        dialogframe.rowconfigure(5, weight=0)
        dialogframe.rowconfigure(6, weight=0)
        dialogframe.rowconfigure(7, weight=0)
        dialogframe.rowconfigure(8, weight=0)
        dialogframe.columnconfigure(0, weight=1)



    def on_select_model(self, event):
        self.accept_fit_parameters['state'] = 'normal'
        self.set_parameters['state'] = 'normal'

    def on_select_new_parameters(self):
        selected_model = self.model_listbox.get(ACTIVE)
        model = [self.model_listbox.get(index) for index in 
                 self.model_listbox.curselection()]
        for function in model:
            params = self.models.parameter_settings(function)
            ParameterOptionsDialog(self.root, function, params, self.models) 

    def cancel_fit(self):
        self.root.destroy()
    
    def perform_spectral_fit(self):
        statistic = self.fit_statisticTK.get()
        weights = self.fit_weightsTK.get()
        undetermined = self.undetermined_values.get()
        model = [self.model_listbox.get(index) for index in 
                 self.model_listbox.curselection()]
        self.callback(model, statistic, weights, undetermined, batch=self.batch)
        self.root.destroy()
          
    def command(self):
        pass

class ParameterOptionsDialog(Dialog):
    def __init__(self, master, model_name, params, models):
        title = 'Model Parameters'
        self.model_name = model_name
        self.params = params
        self.models = models
        
        # determine size
        nparams = len(params)
        xsize = 260
        ysize = 110
        ysize += 50*nparams
        
        super(ParameterOptionsDialog, self).__init__(title, master=master, 
                                                     width=xsize, height=ysize)
        self.root.geometry("{0}x{1}+{2}+{3}".format(self.width, self.height, 
                                                    self.xpos, self.ypos))        

        # Create a dialog frame within the new root window 
        dialogframe = Frame(self.root, bg=BGCOLOR, height=100)
        dialogframe.grid(column=0, row=0)
        label1 = Label(dialogframe, text=model_name, bg=BGCOLOR)
        label1.grid(row=0, column=0, columnspan=2, sticky=W)
        label2 = Label(dialogframe, text='Fixed', bg=BGCOLOR)
        label2.grid(row=0, column=2, columnspan=2, sticky=E)
        seperator = ttk.Separator(dialogframe)
        seperator.grid(row=1, columnspan=4, sticky=E+W, pady=(4,14))

        # initialize values
        self.values = [DoubleVar(dialogframe, value=param['default']) for param in params]
        self.states = [BooleanVar(dialogframe, value=param['freeze']) for param in params]
        
        row = 2
        for i in range(nparams):
            label1 = Label(dialogframe, text=params[i]['name'], bg=BGCOLOR)
            label1.grid(row=row+i, column=0, sticky=W, pady=2)
            entry = ttk.Entry(dialogframe, textvariable=self.values[i], width=7)
            entry.grid(row=row+i, column=1, pady=2)
            label2 = Label(dialogframe, text=params[i]['units'], bg=BGCOLOR)
            label2.grid(row=row+i, column=2, pady=2)
            check = Checkbutton(dialogframe, variable=self.states[i], bg=BGCOLOR)
            check.grid(row=row+i, column=3, sticky=E, padx=(20,0), pady=2)

        row += nparams    
        
        # accept and cancel buttons
        accept_button = Button(dialogframe, text=u"Accept", width=5, 
                               highlightbackground=BGCOLOR, command=self.set_params)
        accept_button.grid(row=row, columnspan=4, sticky=N+W+S+E, pady=(12, 0))
        cancel_button = Button(dialogframe, text=u"Cancel", width=5, 
                               highlightbackground=BGCOLOR, command=self.root.destroy)
        cancel_button.grid(row=row + 1, columnspan=4, sticky=N+W+S+E, pady=0)
        self.root.grid_columnconfigure(0, weight=1)
        self.root.grid_rowconfigure(0, weight=1)

        self.root.lift()
        
        self.root.grab_set()
        self.master.wait_window(self.root)

    def set_params(self):
        params = self.models.parameter_settings(self.model_name)
        nparams = len(params)
        for i in range(nparams):
            params[i]['default'] = self.values[i].get()
            params[i]['freeze'] = self.states[i].get()
        self.models.update_parameter_settings(self.model_name, params)
        self.root.destroy()



class PlotModelParametersDialog(Dialog):
    def __init__(self, master, selected_model, callback):
        title = 'Batch Fit Plotter'
        width = 325
        height = 300

        super(PlotModelParametersDialog, self).__init__(title, master=master, width=width,
                                           height=height)
        self.root.geometry("{0}x{1}+{2}+{3}".format(self.width, self.height, 
                                                    self.xpos, self.ypos))        
        
        # Set the callback function
        self.callback = callback

        # Create an instance variable to store the selected model and parameters
        self.selection_history = []
        self.selected_model = selected_model

        # Create the tkinter interface
        self._construct()

        self.root.grab_set()
        self.master.wait_window(self.root)

    def _construct(self):

        # Set the root window weightings
        self.root.rowconfigure(0, weight=1)
        self.root.columnconfigure(0, weight=1)

        # Set a minimum window size
        self.root.minsize(300, 250)

        # Create a dialog frame within the new root window 
        dialogframe = Frame(self.root, bg=BGCOLOR)
        # dialogframe.pack(anchor=N, fill=BOTH, expand=True, side=TOP)

        # Tell the frame to expand and fill the root window
        dialogframe.grid(sticky=N+W+E+S)

        # Create the text label 
        label_title = Label(dialogframe, text='Select one parameter to plot against time\nor two parameter to plot against each other',
            # or two model parameters to plot', 
                            bg=BGCOLOR)
        label_title.config(justify=CENTER)

        # Create the listbox
        self.parameter_listbox = Listbox(dialogframe, borderwidth=2, selectbackground='#116cd6',
                                     selectborderwidth=0, highlightthickness=0,
                                     relief=SUNKEN, selectmode=MULTIPLE)
        self.parameter_listbox.bind('<<ListboxSelect>>', self.on_select_parameter)


        # Get the available models
        models = XspecModelManager()

        # Get the list of parameter names associated with this model
        self.parameter_names = models.parameter_names(self.selected_model[0])
        # params = models.parameter_settings(self.selected_model)

        for parameter_name in self.parameter_names:
            self.parameter_listbox.insert(END, parameter_name)


        
        frame_buttons = Frame(dialogframe, bg=BGCOLOR)
        self.plot_fit_parameter = Button(frame_buttons, text=u"Plot", highlightbackground=BGCOLOR,
                                          command=self.plot_parameter, width=10)
        self.plot_fit_parameter['state'] = 'disabled'
        # self.omit_fit_parameters = Button(frame_buttons, text=u"Omit", highlightbackground=BGCOLOR,
        #                                 command=self.command, width=5)
        # self.omit_fit_parameters['state'] = 'disabled'
        # self.restore_fit_parameters = Button(frame_buttons, text=u"Restore", highlightbackground=BGCOLOR,
        #                                    command=self.command, width=5)
        # self.restore_fit_parameters['state'] = 'disabled'
        self.cancel_fit_parameters = Button(frame_buttons, text=u"Cancel", highlightbackground=BGCOLOR,
                                          command=self.cancel, width=10)
        
        self.plot_fit_parameter.pack(side=LEFT)
        # self.omit_fit_parameters.pack(side=LEFT)
        # self.restore_fit_parameters.pack(side=LEFT)
        self.cancel_fit_parameters.pack(side=LEFT)

        # # Set the default values
        # self.statistic_rb['Chi2'].select()
        # self.weight_rb['Standard'].select()
        # self.batch_rb['Leave free'].select()

        # Pack it all together
        # label_title.pack(pady=10, side=TOP)
        # self.parameter_listbox.pack(padx=10, side=TOP, fill=X)
        # frame_buttons.pack(padx=4, pady=15, side=TOP)

        label_title.grid(row=0, column=0, sticky=E+W, pady=(10,0))
        self.parameter_listbox.grid(row=1, column=0, sticky=N+E+W+S, padx=10, pady=10)
        frame_buttons.grid(row=2, column=0, sticky=N+S, padx=4, pady=(10,20))


        # Configure the weights of the rows and columns
        dialogframe.rowconfigure(0, weight=0)
        dialogframe.rowconfigure(1, weight=1)
        dialogframe.columnconfigure(0, weight=1)

        self.selection = self.parameter_listbox.curselection()

    def on_select_parameter(self, event):

        # # Add the previous selection to the list of selected parameters (ordered by selection history)
        # selected_parameter = self.parameter_listbox.get(ACTIVE)
        # if selected_parameter not in self.selection_history:
        #     self.selection_history.append(self.parameter_listbox.get(ACTIVE))

        # # Get the current name and index of the currently selected parameters (order by index)
        # selected_parameter_idx = self.parameter_listbox.curselection()
        # selected_parameter_names = [self.parameter_listbox.get(idx) for idx in self.parameter_listbox.curselection()]

        # if len(selected_parameter_names) == 0:
        #     self.selection_history = []

        # # Check to see if more than 2 selections exist
        # if len(selected_parameter_names) > 2:

        #     # Get the oldest selected parameter
        #     oldest_selected_parameter = self.selection_history[0]

        #     # Find the index of the oldest selected parameter
        #     for name, idx in zip(selected_parameter_names, selected_parameter_idx):
        #         if name in oldest_selected_parameter:
        #             self.parameter_listbox.selection_clear(idx)

        #     self.selection_history.remove(oldest_selected_parameter)

        # self.plot_fit_parameter['state'] = 'normal'

        # Limit the selection number to two
        if len(self.parameter_listbox.curselection()) > 2:
            for index in self.parameter_listbox.curselection():
                if index not in self.selection:
                    self.parameter_listbox.selection_clear(index)

        # Set the current selections
        self.selection = self.parameter_listbox.curselection()

        # Enable the plot button
        self.plot_fit_parameter['state'] = 'normal'


    def cancel(self):
        self.root.destroy() 

    def plot_parameter(self):
        # selected_parameter = self.parameter_listbox.get(ACTIVE)
        selected_parameter_idx = self.parameter_listbox.curselection()

        # parameter_index = 0
        # for parameter_name in self.parameter_names:
        #     if parameter_name in selected_parameter:
        #         break
        #     else:
        #         parameter_index = parameter_index + 1

        # Call fitplotter.plotBatchFitResults(parameter_index)
        # self.callback(parameter_index)
        self.callback(parameter_index=selected_parameter_idx)

        self.root.destroy()
