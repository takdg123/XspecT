import numpy as np
import datetime, sys
from os import path, remove, makedirs
from gbm.detectors import Detector

class TriggerDataFinderFTP(object):
    """Class for finding and downloading GBM trigger data from the FSSC FTP server.
    """
    def __init__(self, download_dir=None):
        """
        Parameters:
        -----------
        download_dir: str, optional
            The path to a directory where file transfers will be downloaded
        """
        ftp_host = 'legacy.gsfc.nasa.gov'
        self._ftp = self._init_ftp(ftp_host)
        self._download_dir = download_dir
        
        self._ftp_url = '/fermi/data/gbm/triggers'
        self._dets = (det.short_name for det in Detector)
        
        self.trigger_num = ''
        self._file_list = []
        self._filename = None

    def __del__(self):
        self._ftp.close()
    
    def _init_ftp(self, host):
        """Initialize the FTP client.
        
        Parameters:
        -----------
        host: str
            The FTP server
            
        Returns:
        --------
        ssh: ftplib.FTP
            The FTP client
        """
        from ftplib import FTP
        ftp = FTP(host)
        ftp.login()
        return ftp
    
    def _construct_path(self, str_trigger_num):
        """Constructs the FTP path
        
        Parameters:
        -----------
        str_trigger_num: str
            The trigger number
        """
        year = '20'+str_trigger_num[0:2]
        dir_path = path.join(self._ftp_url, year, 'bn'+str_trigger_num, 'current')
        return dir_path
    
    def _ls_dir(self, trigger_num):
        """Stores the list of files in the directory
        
        Parameters:
        -----------
        trigger_num: str or int
            The trigger number
        """
        str_trigger_num = str(trigger_num)
        if str_trigger_num == self.trigger_num:
            return
        
        self.trigger_num = str_trigger_num
        dir_path = self._construct_path(str_trigger_num)
        try:
            self._file_list = self._ftp.nlst(dir_path)
        except:
            print('No files for trigger number: {0}'.format(trigger_num))
            self._file_list = []
        
    def _file_filter(self, filetype, extension, dets=None):
        """Filters the directory for the requested filetype, extension, and detectors
        
        Parameters:
        -----------
        filetype: str
            The type of file, e.g. 'cspec'
        extension: str
            The file extension, e.g. '.pha'
        dets: list, optional
            The detectors.  If omitted, then files for all detectors are returned

        Returns:
        --------
        files: list
            The filtered file list
        """
        files = [f for f in self._file_list if (filetype in f) & (f.endswith(extension))]
        if dets is not None:
            if type(dets) == str:
                dets = [dets]
            files = [f for f in files if any('_'+det+'_' in f for det in dets)]
        
        return files
    
    def ctime_files(self, trigger_num, dets=None):
        """Returns the full path of the CTIME files for the requested trigger
        
        Parameters:
        -----------
        trigger_num: str or int
            The trigger number
        dets: list, optional
            The detectors.  If omitted, then files for all detectors are returned

        Returns:
        --------
        ctime_files: list
            The available CTIME files
        """
        self._ls_dir(trigger_num)
        ctime_files = self._file_filter('ctime', '.pha', dets=dets)
        return ctime_files
    
    def cspec_files(self, trigger_num, dets=None):
        """Returns the full path of the CSPEC files for the requested trigger
        
        Parameters:
        -----------
        trigger_num: str or int
            The trigger number
        dets: list, optional
            The detectors.  If omitted, then files for all detectors are returned

        Returns:
        --------
        cspec_files: list
            The available CSPEC files
        """
        self._ls_dir(trigger_num)
        cspec_files = self._file_filter('cspec', '.pha', dets=dets)
        return cspec_files
    
    def tte_files(self, trigger_num, dets=None):
        """Returns the full path of the TTE files for the requested trigger
        
        Parameters:
        -----------
        trigger_num: str or int
            The trigger number
        dets: list, optional
            The detectors.  If omitted, then files for all detectors are returned

        Returns:
        --------
        tte_files: list
            The available TTE files
        """
        self._ls_dir(trigger_num)
        tte_files = self._file_filter('tte', '.fit', dets=dets)
        return tte_files

    def rsp_files(self, trigger_num, ctime=False, cspec=False, rsp2=False, dets=None):
        """Returns the full path of the RSP  or RSP2 files for the requested trigger
        
        Parameters:
        -----------
        trigger_num: str or int
            The trigger number
        ctime: bool, optional
            Set to True to return CTIME response files.  Default is False.
        cspec: bool, optional
            Set to True to return CSPEC response files.  Default is False.
        rsp2: bool, optional
            Set to True to return RSP2 files, otherwise return RSP files.  Defaults is False
        dets: list, optional
            The detectors.  If omitted, then files for all detectors are returned

        Returns:
        --------
        rsp_files: list
            The available RSP or RSP2 files
        """
        self._ls_dir(trigger_num)
        if (ctime == True) & (cspec == True):
            filetype = 'glg_'
        elif (ctime == True) & (cspec == False):
            filetype = 'ctime'
        elif (ctime == False) & (cspec == True):
            filetype = 'cspec'
        else:
            return []
        
        if rsp2 == True:
            extension = 'rsp2'
        else:
            extension = 'rsp'
        rsp_files = self._file_filter(filetype, extension, dets=dets)
        return rsp_files
    
    def trigdat_file(self, trigger_num):
        """Returns the full path of the TrigDat file for the requested trigger
        
        Parameters:
        -----------
        trigger_num: str or int
            The trigger number

        Returns:
        --------
        trigdat_file: list
            The available TrigDat file
        """
        self._ls_dir(trigger_num)
        trigdat_file = self._file_filter('trigdat', '.fit')[0]
        return trigdat_file
    
    def _ftp_status(self, chunk):
        """FTP GET callback function that downloads and reports the percent 
        progress of the download.
        
        Parameters:
        -----------
        chunk: str
            The byte data to be written
        """
        # append to file
        with open(path.join(self._download_dir, self._filename), 'ab') as f:
            f.write(chunk)
        
        self._transferred_bytes += len(chunk)
        percent = float(self._transferred_bytes)/float(self._total_bytes)
        
        # download bar
        bar = ('=' * int(percent * 30)).ljust(30)
        # format percent and print along with download bar
        percent = str("{0:.2f}".format(percent*100.0))
        sys.stdout.write("\r%s [%s] %s%%" % (self._filename, bar, percent))
        # file download is finished
        if self._transferred_bytes == self._total_bytes:
            sys.stdout.write('\n')
        sys.stdout.flush()

    def _ftp_silent(self, chunk):
        """FTP GET callback function that silently downloads a file.
        
        Parameters:
        -----------
        chunk: str
            The byte data to be written
        """
        # append to file
        with open(path.join(self._download_dir, self._filename), 'ab') as f:
            f.write(chunk)        
    
    def ftp_get(self, files, download_dir=None, callback=None):
        """Transfer a file from a FTP server to the local machine
        
        Parameters:
        -----------
        files: str or list
            The filepath of the file or files on the remote server
        download_dir: str, optional
            The directory on the local machine to which the file will be 
            transferred.  If not specified, the directory will default to the
            directory specified in temp_path at initialization.  If temp_path
            was also not specified, then an exception is thrown.
        
        Returns:
        --------
        filepaths: list
           The filepaths of the downloaded file on the local machine.
        """
        
        if type(files) == str:
            files = [files]
        
        # set callback
        if callback is None:
            callback = self._ftp_silent
            
        # set download directory
        if download_dir is None:
            if self._download_dir is None:
                raise ValueError("Must designate download directory")
        else:
            self._download_dir = download_dir

        if path.exists(self._download_dir) == False:
            makedirs(self._download_dir)

        # download each file
        filepaths = []
        for file in files:
            # have to save in self because this can't be passed as an argument
            # in the callback
            self._filename = path.basename(file)
        
            # download file
            self._ftp.voidcmd('TYPE I')
            self._total_bytes = self._ftp.size(file)
            self._transferred_bytes = 0
            self._ftp.retrbinary('RETR '+file, callback=callback)        
            
            filepaths.append(path.join(self._download_dir, self._filename))
        return filepaths

