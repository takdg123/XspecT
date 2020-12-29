import os
from collections import OrderedDict
from gbm.detectors import Detector
from gbm.file import GbmFile

class Data(object):
    """Base class for interface to data files

    This class should not be used directly, instead use one that inherits from
     it and is specific to your data type e.g. TTE

    Attributes:
    -----------
    data_type: str
        The GBM data type of the file
    detector: str
        The GBM detector the file is associated with
    filename: str
        file name
    id: str
        The file ID
    is_gbm_file: bool
        True if the file is a GBM standard file, False if it is not
    is_trigger: bool
        True if the file is a trigger file, False if it is not
    
    Public Methods:
    ---------------
    create_filename:
        Create a filename and update Data attributes
    set_filename:
        Set the filename
    """
    def __init__(self, filename):
        """
        Parameters:
        --------------
        fname: str
            name of file
        """
        if filename is not None:
            if not os.path.isfile(filename):
                raise IOError("File {0} does not exist".format(filename))

        self.filename = None
        self.filename_obj = GbmFile()
        self.is_gbm_file = False
        self.detector = None
        self.is_trigger = None
        self.data_type = None
        self.id = None
        
        self.set_filename(filename)
        
        if self.is_gbm_file:
            if self.filename_obj.detector == 'all':
                self.detector = 'all'
            else:
                self.detector = self.filename_obj.detector.short_name
            self.is_trigger = True if self.filename_obj.trigger == 'bn' else False
            self.data_type = self.filename_obj.data_type.upper()
            self.id = self.filename_obj.uid
    
    def create_filename(self, extension, meta=None):
        """Create a filename and update Data attributes

        Parameters:
        -----------
        extension: str
            The filename extension
        meta: str, optional 
            A metadata atttribute to be added to the filename
        
        Returns:
        -----------
        filename: str
            The created filename
        """
        trigger = None
        if self.is_trigger:
            trigger = 'bn'
        filename = GbmFile.create(uid=self.id, data_type=self.data_type,
                                  detector=self.detector, trigger=trigger,
                                  extension=extension, meta=meta)
        self.filename = filename.basename()
        self.filename_obj = filename
        return self.filename
    
    def set_filename(self, filename, directory=None):
        """Set the filename

        Parameters:
        -----------
        filename: str
            The filename
        directory: str, optional 
            The directory
        """
        if directory is not None:
            filename = os.path.join(directory, filename)
        self.filename = filename
        
        try:
            self.filename_obj = GbmFile.from_path(self.filename)
            if self.filename_obj is not None:
                self.is_gbm_file = True
        except:
            pass
    