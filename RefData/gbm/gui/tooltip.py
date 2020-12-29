from tkinter import *

class ListboxToolTip(object):
    """Class for displaying a Tooltip for tkinter Listbox
    """    
    def __init__(self, master, textlist, delay=500):
        self.delay = delay
        self.textlist = textlist
        self.root = None
        self.callback = None
        self.index = -1
            
        self.master = master
        self.master.bind("<Motion>", self.update)
        self.master.bind("<Enter>", self.update)
        self.master.bind("<Leave>", self.deactivate)
    
    def update(self, event=None):
        """Check if mouse is now hovering over a new list entry.
        Update the tooltip if it is.
        """
        index = self.master.nearest(event.y)
        if index != self.index:
            self.index = index
            self.deactivate()
            self.activate()
    
    def activate(self, event=None):
        """Activate the tooltip display with a hover delay
        """
        self.unschedule()
        self.callback = self.master.after(self.delay, self.display)
    
    def deactivate(self, event=None):
        """Deactivate the tooltip
        """
        self.unschedule()
        root = self.root
        self.root = None
        if root:
            root.destroy()
    
    def unschedule(self):
        """Remove the hover delay callback
        """
        callback = self.callback
        self.callback = None
        if callback:
            self.master.after_cancel(callback)
        
    def display(self, event=None):
        """Show the tooltip for the corresponding Listbox element
        """
        if self.index == -1:
            return
        
        self.root = Toplevel(self.master)
        self.root.wm_overrideredirect(True)

        #x, y, cx, cy = self.master.bbox("insert")
        x, y, cx, cy = self.master.bbox(self.index)
        x += self.master.winfo_rootx() + 25
        y += self.master.winfo_rooty() + 20
        self.root.wm_geometry("+%d+%d" % (x, y))
        
        text = self.textlist[self.index]
        label = Label(self.root, text=text, justify='left',
                       background="#ffffff", relief='solid', borderwidth=1)
        label.pack(ipadx=1)
