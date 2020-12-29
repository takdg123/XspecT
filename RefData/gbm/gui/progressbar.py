from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
from tkinter import ttk

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