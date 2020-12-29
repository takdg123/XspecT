from tkinter import *
from tkinter import ttk
from .globals import BGCOLOR, FRAMERELIEF, OM_JUSTIFY, GUIFONTCOLOR
from .window import Window
import copy

#from gbm.xspec.xspec import XSPEC

class Logger(Window):
    name = 'Fit Log'
    filetypes = ((' ', '*.txt'),)
    down_stack = []
    
    def __init__(self, xspec, gspec_root=None):
        super(Logger, self).__init__()
        root = Tk()
        root.title(self.name)
        self.root = root
        self.width = self._geometry.width_fraction_to_npix(0.405)
        self.height = self._geometry.height_fraction_to_npix(0.68)
        self.aspectratio = self.width/self.height

        # Set soft maximum windows size limits
        if self.width > 585:
            self.width = 585
        if self.height > 618:
            self.height = 618

        self.root.geometry("{0}x{1}+0+0".format(self.width, self.height))
        
        # Set a minimum window size
        root.minsize(560, 200)

        # register this window with the root application window
        self.gspec_root = gspec_root
        if gspec_root is not None:
            self.gspec_root.register_window(self.name, self)

            # reference the menubar
            self.menubar = copy.copy(self.gspec_root.menubar)
            self.menubar.update_root(self.root, self.name)
            root.config(menu=self.menubar.top)
        
        # redirect output to this logger
        self.xspec = xspec
        self.xspec.my_session.logger = self.logger
        

        # Initialize a frame instance with the root as its parent
        Frame.__init__(self, self.root)
        self.config(bg=BGCOLOR)       
        self.grid(sticky=N + W + E + S)         
        '''
        root.columnconfigure(0, weight=1)
        root.columnconfigure(1, weight=1)
        root.columnconfigure(2, weight=1)
        root.columnconfigure(3, weight=1)
        root.columnconfigure(4, weight=1)
        root.columnconfigure(5, weight=1)
        root.columnconfigure(6, weight=1)
        root.columnconfigure(7, weight=1)
        root.rowconfigure(0, weight=1)
        '''
        frame = Frame(self, bg=BGCOLOR, relief=FRAMERELIEF, pady=10, padx=10, height=100)

        text = Text(frame, width=80, height=36, bg=BGCOLOR)
        text.grid(column=0, row=0, columnspan=8, sticky=(N, W, E, S))
        text.columnconfigure(0, weight=1)
        scrollbar = ttk.Scrollbar(frame, orient=VERTICAL, command=text.yview)
        scrollbar.grid(column=8, row=0, sticky=(N, S))
        scrollbar.columnconfigure(0, weight=0)
        text['yscrollcommand'] = scrollbar.set
        self.text = text

        frame_buttons = Frame(self, bg=BGCOLOR, relief=FRAMERELIEF, pady=10, padx=10, height=55)

        label = Label(frame_buttons, text='Enter XSPEC Command:', bg=BGCOLOR)
        label.grid(column=0, row=0)
        command = StringVar()
        my_cmd = ttk.Entry(frame_buttons, textvariable=command)
        my_cmd.grid(column=1, row=0, columnspan = 4, sticky=(N,E,S,W))
        my_cmd.bind('<Return>', (lambda event: self.enter_cmd()))
        my_cmd.bind('<Up>', (lambda event: self.pop_cmd()))
        my_cmd.bind('<Down>', (lambda event: self.put_cmd()))
        self.my_cmd = my_cmd

        self.saveScriptButton = ttk.Button(frame_buttons, text="Save Script", command=self.save_script)
        self.saveScriptButton.grid(column=0, row=1, sticky=(N,E,S,W))
        self.runScriptButton = ttk.Button(frame_buttons, text="Run Script", command=self.run_script)
        self.runScriptButton.grid(column=1, row=1, sticky=(N,E,S,W))

        clear = ttk.Button(frame_buttons, text="Clear Text", command=self.clear_cmd)
        clear.grid(column=2, row=1, sticky=(N,E,S,W))
        save_text = ttk.Button(frame_buttons, text="Save Text", command=self.save_text)
        save_text.grid(column=3, row=1, sticky=(N,E,S,W))
        cancel = ttk.Button(frame_buttons, text="Exit", command=self.on_window_close)
        cancel.grid(column=4, row=1, sticky=(N,E,S,W))
              
        frame.grid(column=0, row=0, sticky=(N,E,S,W))  
        frame_buttons.grid(column=0, row=1, sticky=(W,N))  

        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)

        root.grid_columnconfigure(0, weight=1)
        root.grid_rowconfigure(0, weight=1)

        # Configure the weights of the rows and columns
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)

        root.protocol("WM_DELETE_WINDOW", self.on_window_close)

        # root.grid_columnconfigure(0, weight=1)

        # root.mainloop()

    def on_window_close(self):
        if self.gspec_root is not None:
            self.gspec_root.unregister_window(self.name)
        self.root.destroy()            


    def save_text(self):
        contents = self.text.get(1.0, END)
        fullfilename = filedialog.asksaveasfilename(defaultextension='txt', 
                        title = "Select file for saving Fit Log:", 
                        parent=self.root)
        if fullfilename == '':
            return
        save_file = open(fullfilename, 'w')
        save_file.write(contents)
        save_file.close()
    
    def logger(self, text=None, end=None):
        self.text.insert('end', text)
        self.text.see(END)

    def clear_cmd(self):
        self.text.delete(1.0, 'end')

    def save_script(self):
        # Save the XSPEC commands to a script file
        #initialdir = self.directory
        initialfile = 'xspec_commands.txt'
        filename = filedialog.asksaveasfilename(initialfile=initialfile,  
                                                title="Save XSPEC Commands", 
                                                filetypes=self.filetypes)
        if filename == '':
            return
        self.xspec.save_script(filename)
        
    def run_script(self):
        # Run an XSPEC script file
        filetypes = ((' ', '*.txt'),)
        filename = filedialog.askopenfilename(filetypes=filetypes, 
                                              title='Select an XSPEC Command file')
        if filename == '':
            return
        self.xspec.run_script(filename)

    def enter_cmd(self):
        my_command = self.my_cmd.get()
        self.xspec.do_command(my_command)
        #self.xspec.my_session.command('%s\n' % my_command)
        #self.xspec.my_session.response_nowait()
        self.my_cmd.delete(0, 'end')

    def pop_cmd(self):
        if len(self.xspec.cmd_stack) != 0:
            last_cmd = self.xspec.cmd_stack.pop()
            cur_cmd = self.my_cmd.get()
            self.down_stack.append(cur_cmd)
            self.my_cmd.delete(0, 'end')
            self.my_cmd.insert(0, last_cmd)

    def put_cmd(self):
        if len(self.down_stack) != 0:
            last_cmd = self.down_stack.pop()
            cur_cmd = self.my_cmd.get()
            self.xspec.cmd_stack.append(cur_cmd)
            self.my_cmd.delete(0, 'end')
            self.my_cmd.insert(0, last_cmd)

