from tkinter import Frame
from .geometry import GeometryManager

class Window(Frame):
    def __init__(self):
        self._geometry = GeometryManager()
        self.root = None
        self.width = None
        self.height = None
        self.aspectratio = None
            
    def on_window_close(self):
        try:
            self.root.destroy()
        except:
            pass
    
    def window_focus(self):
        # set up window dimensions and position    
        xpos = self._geometry.center_x(self.width)
        ypos = self._geometry.center_y(self.height)
        self.root.geometry("{0}x{1}+{2}+{3}".format(self.width, self.height, 
                                                    xpos, ypos))
        self.root.lift()
