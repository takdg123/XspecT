from tkinter import Tk
from math import floor


class GeometryManager(object):
    """Manage various window size/placement requests based on display resolution

    Attributes:
    -----------
    width:
        The display width in pixels
    height:
        The display height in pixels
    aspectratio:
        The display aspect ratio
    
    Public Methods:
    ---------------
    center_x:
        Offset in pixels to center a window given its width
    center_y:
        Offset in pixels to center a window given its height
    height_fraction_to_npix:
        Given a fractional height of the display, compute the number of pixels
    ncolumns_to_npix:
        Given n evenly-space columns, calculate the width of a window
    nrows_to_npix:
        Given n evenly-space rows, calculate the height of a window
    width_fraction_to_npix:
        Given a fractional width of the display, compute the number of pixels        
    """
    def __init__(self):    
        root = Tk()
        root.withdraw()
        self.width = root.winfo_screenwidth()
        self.height = root.winfo_screenheight()
        self.aspectratio = self.width/self.height
        root.destroy()
    
    def width_fraction_to_npix(self, width_fraction):
        """Convert fraction if display width to pixels

        Parameters:
        -----------
        width_fraction: float
            The fractional width of the display
            
        Returns:
        -----------
        npix: int
            The requested width in pixels
        """
        npix = self.width * width_fraction
        npix = int(round(npix))
        return npix

    def height_fraction_to_npix(self, height_fraction):
        """Convert fraction if display height to pixels

        Parameters:
        -----------
        height_fraction: float
            The fractional height of the display
            
        Returns:
        -----------
        npix: int
            The requested height in pixels
        """
        npix = self.height * height_fraction
        npix = int(round(npix))
        return npix
        
    def center_x(self, width):
        """The horizontal offset required to center a window given its width.
        The offset is measured relative to the left edge of the display

        Parameters:
        -----------
        width: int
            The width of the window in pixels
            
        Returns:
        -----------
        pixel: int
            The horizontal offset in pixels
        """
        halfwidth = float(width)/2.0
        halfscreenwidth = float(self.width)/2.0
        pixel = int(round(halfscreenwidth-halfwidth))
        return pixel
    
    def center_y(self, height):
        """The vertical offset required to center a window given its height.
        The offset is measured relative to the top edge of the display

        Parameters:
        -----------
        height: int
            The height of the window in pixels
            
        Returns:
        -----------
        pixel: int
            The vertical offset in pixels
        """
        halfheight = float(height)/2.0
        halfscreenheight = float(self.height)/2.0
        pixel = int(round(halfscreenheight-halfheight))
        return pixel
    
    def ncolumns_to_npix(self, ncolumns):
        """Calculate the width of a window if the display is divided into 
        n evenly-spaced columns

        Parameters:
        -----------
        ncolumns: int
            The number of columns into which the display is subdivided
            
        Returns:
        -----------
        npix: int
            The width in pixels
        """
        npix = float(self.width)/ncolumns
        return floor(npix)
    
    def nrows_to_npix(self, nrows):
        """Calculate the height of a window if the display is divided into 
        n evenly-spaced rows

        Parameters:
        -----------
        nrows: int
            The number of rows into which the display is subdivided
            
        Returns:
        -----------
        npix: int
            The height in pixels
        """
        npix = float(self.height)/nrows
        return floor(npix)
        
    def child_center_x(self, parent, child_width):
        parent_x = parent.winfo_x()
        parent_halfwidth = float(parent.winfo_width())/2.0
        child_halfwidth = float(child_width)/2.0
        pixel = int(round(parent_halfwidth-child_halfwidth)+parent_x)
        return pixel
    
    def child_center_y(self, parent, child_height):
        parent_y = parent.winfo_y()
        parent_halfheight = float(parent.winfo_height())/2.0
        child_halfheight = float(child_height)/2.0
        pixel = int(round(parent_halfheight-child_halfheight)+parent_y)
        return pixel
        
