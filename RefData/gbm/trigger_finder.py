import sys
py_version = sys.version_info[0]
if py_version == 2:
    from urllib2 import urlopen
    from Tkinter import *
    from Tkinter import tkFileDialog as filedialog
    import ttk
elif py_version == 3:
    from urllib.request import urlopen
    from tkinter import *
    from tkinter import ttk
    from tkinter import filedialog
import socket
import datetime
import os
import ssl
from math import floor
import numpy as np
from gbm.time import Met
from .GBMDataFinder import TriggerDataFinderFTP
import astropy.io.fits as fits
from gbm.coords import spacecraft_to_radec, haversine
from gbm.detectors import Detector

class TriggerRecord(object):
    """Container class for trigger record

    Attributes:
    -----------
    bnum: str
        The trigger number
    utc: str
        The time of the trigger, in UTC
    met: str 
        The time of the trigger, in Fermi MET
    type: str
        The type of the trigger
    ra: str
        The RA of the trigger, either in decimal or sexagesimal format
    dec: str
        The Dec of the trigger, either in decimal or sexagesimal format
    loc_err: str
        The 1 sigma circular equivalent localization uncertainty
    loc_instr: str
        The localizing instrument
        
    Public Methods:
    ---------------
    use_decimal:
        self.ra and self.dec will be in decimal degrees
    use_sexagesimal:
        self.ra and self.dec will be in sexagesimal format
    print_line:
        prints a formatted line of the record information
    """    
    def __init__(self, bnum, utc=None, type=None, ra=None, dec=None, loc_err=None, 
                 loc_instr=None):
        """
        Parameters:
        -----------
        bnum: str
            The trigger number
        utc: str, optional
            The time of the trigger
        type: str, optional
            The type of the trigger
        ra: str, optional
            The RA of the trigger
        dec: str, optional
            The Dec of the trigger
        loc_err: str, optional
            The 1 sigma circular equivalent localization uncertainty
        loc_instr: str, optional
            The localizing instrument
        """
        self.bnum = bnum.strip()
        self.utc = utc.strip()
        self.met = self._met()
        self.type = type.strip()
        self._ra = ra.strip()
        self._dec = dec.strip()
        self.use_decimal()
        self.loc_err = loc_err.strip()
        self.loc_instr = loc_instr.strip()
    
    def _met(self):
        """Convert the trigger UTC to MET
        
        Returns:
        --------
        met: float
            The trigger MET
        """
        if self.utc is None:
            return None
        whole, frac = self.utc.split('.')
        utc = datetime.datetime.strptime(whole,'%Y-%m-%dT%H:%M:%S')
        utc = utc.replace(microsecond=int(frac)*10**(6-len(frac)))
        met = '{0:13.3f}'.format(Met.from_datetime(utc).met)
        return met
    
    def _ra_decimal(self):
        """Convert the RA from sexagesimal to decimal
        
        Returns:
        --------
        ra_deci: float
            The RA
        """
        if self._ra is None:
            return None
        hour, min, sec = self._ra.split()
        ra_deci = float(hour)*15.0 + (float(min)*15.0)/60.0 + (float(sec)*15.0)/3600.0
        
        return ra_deci
    
    def _dec_decimal(self):
        """Convert the Dec from sexagesimal to decimal
        
        Returns:
        --------
        dec_deci: float
            The Dec
        """
        if self._dec is None:
            return None
        deg, arcmin, arcsec = self._dec.split()
        sign = np.sign(float(deg))
        dec_deci = np.abs(float(deg)) + float(arcmin)/60.0 + float(arcsec)/3600.0
        return sign*dec_deci
    
    def use_decimal(self):
        self.ra = '{0:6.2f}'.format(self._ra_decimal())
        self.dec = '{0:6.2f}'.format(self._dec_decimal())
    
    def use_sexagesimal(self):
        self.ra = self._ra
        self.dec = self._dec
            
    def print_line(self, met=False, decimal_degrees=True):
        """Print a fixed-width formatted line of the trigger information
        
        Parameters:
        -----------
        met: bool, optional
            If True, print the trigger time in MET. Default is False.
        decimal_degrees: bool, optional
            If True, print the RA and Dec in decimal degrees.  Default is True.
        
        Returns:
        --------
        line: str
            Fixed width formatted line
        """
        if met:
            t = '{0:13.3f}'.format(self.MET())#{0:23.3f}
        else:
            t = self.utc
        
        if decimal_degrees:
            ra = '{0:6.2f}'.format(self.ra_decimal())
            dec = '{0:6.2f}'.format(self.dec_decimal())
        else:
            ra = self.ra
            dec = '{0:>9}'.format(self.dec)
        
        ttype = '{0:<7}'.format(self.type)
        loc_instr = '{0:<14}'.format(self.loc_instr)
        loc_err = '{0:4.1f}'.format(float(self.loc_err))
        line = ' '.join((self.bnum, ttype, t, ra, dec, loc_err, loc_instr))
        return line
        

class TriggerList(object):
    """Container class for the trigger list

    Attributes:
    -----------
    triggers: list
        A list of TriggerRecords
        
    Public Methods:
    ---------------
    trigger_types:
        Returns a list of trigger types in the trigger list
    loc_instruments:
        Returns a list of localizing instruments in the trigger list
    sort_list:
        Sort the trigger list based on a given column
    filter_list:
        Filter the trigger list based on matching values
    """    
    def __init__(self):
        print('Querying HEASARC for trigger list')
        host = 'https://heasarc.gsfc.nasa.gov'
        self._is_connected(host)
        url=host+'/cgi-bin/W3Browse/w3query.pl?tablehead=name%3dBATCHRETRIEVALCATALOG_2%2e0+fermigtrig&Action=Query&Coordinates=%27Equatorial%3a+R%2eA%2e+Dec%27&Equinox=2000&Radius=60&NR=&GIFsize=0&Fields=&varon=trigger_name&varon=trigger_time&varon=trigger_type&varon=ra&varon=dec&varon=error_radius&varon=localization_source&sortvar=trigger_name&ResultMax=1000000&displaymode=BatchDisplay'
        self._read_list(url)
        
    def _is_connected(self, host):
        """Test the host for a connection
        Parameters:
        -----------
        host: str
            The host to test
        """
        try:
            # connect to the host -- tells us if the host is actually
            # reachable
            socket.create_connection((host.split('/')[-1], 80))
            return True
        except OSError:
            raise OSError("Either you are not connected to the internet or "
                          "{0} is down.".format(host))
        return False
    
    def _read_list(self, url):
        """Open and read the trigger list output from the host
        
        Parameters:
        -----------
        url: str
            The url of trigger list request
        """
        context = ssl._create_unverified_context()
        page = urlopen(url, context=context)
        txt = page.read()
        # make sure the output is encoded correctly, using python3 it is not
        txt = txt.decode(encoding='utf-8',errors='ignore')
        txt = txt.split('\n')[3:-3]
        self.numTriggers = len(txt)
        self.triggers = []
        for line in txt:
            s = line.split('|')
            record = TriggerRecord(s[1], utc=s[2], type=s[3], ra=s[4], dec=s[5], 
                                   loc_err=s[6], loc_instr=s[7])
            self.triggers.append(record)
    
    def trigger_types(self):
        """List the trigger types
        
        Returns:
        --------
        types: list
            The list of trigger types
        """
        types = list(set([record.type for record in self.triggers]))
        types = sorted(types)
        return types
    
    def loc_instruments(self):
        """List the localizing instruments
        
        Returns:
        --------
        types: list
            The list of localizing instruments
        """
        loc_instr = list(set([record.loc_instr for record in self.triggers]))
        loc_instr = sorted(loc_instr)
        return loc_instr
          
    def sort_list(self, the_list, column):
        """Sort a list by a chosen column
        
        Parameters:
        -----------
        the_list: list
            A list of TriggerRecords
        column: str
            The name of the 'column' which is just the attribute of the TriggerRecords
        
        Returns:
        --------
        new_list: list
            The sorted list of TriggerRecords
        """
        new_list = sorted(the_list, key=lambda record: getattr(record, column))
        return new_list
        
    def filter_list(self, the_list, filters):
        """Filter a list by one or more specification
        
        Parameters:
        -----------
        the_list: list
            A list of TriggerRecords
        filters: dict
            A dictionary where the keys are the columns (TriggerRecord attributes) and
            the values are the filter values
        
        Returns:
        --------
        new_list: list
            The filtered list of TriggerRecords
        """
        new_list = the_list
        for column in filters.keys():
            filter = filters[column]
            new_list = [rec for rec in new_list if getattr(rec, column).startswith(filter)]
        return new_list 
            

class TriggerFinder(object):
    """Class for displaying all of the GBM triggers in a sortable, filterable list.
    Allows the user to select one or more triggers for downloading.
    """    
    def __init__(self, master):
        """Construct the application
        """
        # initialize values
        self._widgets = {}
        self._filters = {f: StringVar() for f in ('bnum', 'type', 'utc', 'met', 'loc_instr')}
        self._get_trigger_list()
        self._num_selected = 0
        self._selected_ids = []
        self.master = master
        master.title("Fermi GBM Trigger List")
        
        # toggle for UTC or MET display
        timeframe = LabelFrame(master, relief=SUNKEN)
        timeframe.grid(row=0, column=0, sticky=E+W)
        label = Label(timeframe, text='Time: ')
        label.grid(row=0, column=0, rowspan=2)
        self._timetype=StringVar()
        self._timetype.set('utc')
        radio1 = Radiobutton(timeframe, text='UTC', variable=self._timetype, 
                             value='utc', command=self._toggle_timetype)
        radio1.grid(row=0, column=1)
        radio2 = Radiobutton(timeframe, text='MET', variable=self._timetype, 
                             value='met', command=self._toggle_timetype)
        radio2.grid(row=1, column=1)
        
        # toggle for sexagesimal or degree location format
        locframe = LabelFrame(master, relief=SUNKEN)
        locframe.grid(row=0, column=1, sticky=E+W)
        label = Label(locframe, text='Coordinate\n Format: ')
        label.grid(row=0, column=0, rowspan=2)
        self._loctype=StringVar()
        self._loctype.set('dec')
        radio3 = Radiobutton(locframe, text='Sexagesimal', variable=self._loctype, 
                             value='sex', command=self._toggle_loctype)
        radio3.grid(row=0, column=1, sticky=W)
        radio4 = Radiobutton(locframe, text='Decimal', variable=self._loctype, 
                             value='dec', command=self._toggle_loctype)
        radio4.grid(row=1, column=1, sticky=W)
        
        # button to accept selection and get data
        acceptframe = LabelFrame(master, relief=FLAT)
        acceptframe.grid(row=0, column=2, sticky=E+W)
        emptylabel = Label(acceptframe, width=41)
        emptylabel.grid(row=0, column=0, sticky=W+E)
        button = Button(acceptframe, text='Get Data', command=self._get_data, 
                        highlightbackground='dark gray')
        button.grid(row=0, column=1, sticky=E+W)
        self._widgets['get_data'] = button
        label = Label(acceptframe, text='0 of {0} triggers selected'.format(self._num_triggers))
        label.grid(row=1, column=1, sticky=E)
        self._widgets['num_triggers'] = label
        
        # filters for bnum
        filterframe = LabelFrame(master, relief=SUNKEN, text='Filters', 
                                 labelanchor=N)
        filterframe.grid(row=2, column=0, columnspan=7, sticky=E+W)
        self._filters['bnum'].set('BNUM filter')
        self._filters['bnum'].trace("w", self._apply_filters)
        bnumfilter = Entry(filterframe, foreground='gray', width=13, 
                           textvariable=self._filters['bnum'])
        bnumfilter.grid(row=0, column=0)
        bnumfilter.bind("<Button-1>", self._entry_click)
        self._widgets['bnumfilter'] = bnumfilter
        
        # filter for time
        self._filters['utc'].set('Time filter')
        self._filters['utc'].trace("w", self._apply_filters)
        self._filters['met'].set('Time filter')
        self._filters['met'].trace("w", self._apply_filters)
        timefilter = Entry(filterframe, foreground='gray', width=21, 
                           textvariable=self._filters['utc'])
        timefilter.grid(row=0, column=2)
        timefilter.bind("<Button-1>", self._entry_click)
        self._widgets['timefilter'] = timefilter
        emptylabel = Label(filterframe, width=33)
        emptylabel.grid(row=0, column=3, columnspan=3)

        # option menu filter for trigger type
        type_list = ['All Types']
        type_list.extend(self._triggerlist.trigger_types())
        self._filters['type'].set(type_list[0])
        self._filters['type'].trace("w", self._apply_filters)
        typemenu = OptionMenu(filterframe, self._filters['type'], *type_list)
        typemenu.config(width=9)
        typemenu.grid(row=0, column=1)
        
        # option menu filter for localizing instrument      
        instr_list = ['All Instruments']
        instr_list.extend(self._triggerlist.loc_instruments())
        self._filters['loc_instr'].set(instr_list[0])
        self._filters['loc_instr'].trace("w", self._apply_filters)
        instmenu = OptionMenu(filterframe, self._filters['loc_instr'], *instr_list)
        instmenu.config(width=13)
        instmenu.grid(row=0, column=6)
        
        # the trigger table - cycle over each column and bind it with the sort function
        # so that clicking on the header will sort/reverse sort the column
        columns = ['BNUM', 'TYPE', 'UTC', 'RA', 'DEC', 'LOC_ERR', 'LOC_INSTR']
        widths = (130, 100, 200, 100, 100, 100, 130)
        table = ttk.Treeview(master, columns=columns[1:], height=20, padding=0)
        for i in range(len(columns)):
            colnum = '#'+str(i)
            colname = columns[i].lower()
            table.column(colnum, width=widths[i])
            table.heading(colnum, text=columns[i], 
                          command=lambda num=colnum, name=colname: self._sort(num, name))
        table.grid(row=3, columnspan=7)
        table.bind("<<TreeviewSelect>>", self._set_num_selected)
        # bind scrollbar
        scrollbar = ttk.Scrollbar(master)
        scrollbar.grid(row=3, column=7, sticky=N+S)
        scrollbar.config(command=table.yview)
        table.config(yscrollcommand=scrollbar.set)
        self._widgets['table']=table
        
        # fill the trigger table
        self._fill_table()
        
        self.master.mainloop()
            
    def _get_trigger_list(self):
        """Get the GBM Trigger List - instantiates a TriggerList class
        """
        #import numpy as np
        #self._triggerlist = np.load('test_list.npy')[0]
        self._triggerlist = TriggerList()
        self._triggers = self._triggerlist.triggers
        self._num_triggers = len(self._triggers)
    
    def _fill_table(self):
        """Fill the trigger table
        """
        table = self._widgets['table']
        # first delete everything in the table
        table.delete(*table.get_children())
        oddrow=True
        for trigger in self._triggers:
            table.insert('', 'end', text=trigger.bnum, 
                         values=(trigger.type, getattr(trigger, self._timetype.get()), 
                                 trigger.ra, trigger.dec, trigger.loc_err, 
                                 trigger.loc_instr), tags=(str(oddrow)))
            oddrow = not oddrow
        # alternating row colors
        table.tag_configure('True', background='antique white')
        table.tag_configure('False', background='white')
                                 
    def _sort(self, colnum, colname, reverse=False):
        """Sort the selected column
        
        Parameters:
        -----------
        colnum: str
            The number of the column in the TreeView object
        colname: str
            The heading of the column
        reverse: bool, optional
            If True, then reverse sort.  Default is False.
        """
        table = self._widgets['table']
        self._triggers = self._triggerlist.sort_list(self._triggers, colname)
        # reverse sort
        if reverse:
            self._triggers = list(reversed(self._triggers))
        # need to set the reverse kw to opposite so the next time we click on the header
        # it reverses the sort
        table.heading(colnum, text=colname.upper(), 
                      command=lambda: self._sort(colnum, colname, reverse=not reverse))    
        self._fill_table()
    
    def _apply_filters(self,*args):
        """Apply the column filters to the table
        
        Parameters:
        -----------
        *args
            Arguments for the callback.  We don't use these.
        """
        
        # get the values and ignore the default values
        filters = {}
        filters['bnum'] = self._filters['bnum'].get()
        if filters['bnum'] == 'BNUM filter':
            filters['bnum'] = ''
        if self._timetype.get() == 'utc':
            filters['utc'] = self._filters['utc'].get()
            if filters['utc'] == 'Time filter':
                filters['utc'] = ''
        elif self._timetype.get() == 'met':
            filters['met'] = self._filters['met'].get()
            if filters['met'] == 'Time filter':
                filters['met'] = ''
        filters['type'] = self._filters['type'].get()
        if filters['type'] == 'All Types':
            filters['type'] = ''
        filters['loc_instr'] = self._filters['loc_instr'].get()
        if filters['loc_instr'] == 'All Instruments':
            filters['loc_instr'] = ''
        
        # send in the filters
        self._triggers = self._triggerlist.filter_list(self._triggerlist.triggers, filters)
        # re-fill the table
        self._fill_table()
        # update number of triggers
        self._set_num_selected({})
        
    def _entry_click(self, event):
        """Clear the context labels from the text entry when clicked
        
        Parameters:
        -----------
        event: dict
            The event dictionary
        """
        if str(self._widgets['bnumfilter']) == event.widget._w:
            entry = self._filters['bnum'].get()
            if entry == 'BNUM filter':
                self._filters['bnum'].set('bn')
                event.widget.config(foreground='black')
        elif str(self._widgets['timefilter']) == event.widget._w:
            entry = self._filters['utc'].get()
            if entry == 'Time filter':
                self._filters['utc'].set('')
                event.widget.config(foreground='black')
            entry = self._filters['met'].get()
            if entry == 'Time filter':
                self._filters['met'].set('')
                event.widget.config(foreground='black')
         
    def _set_num_triggers(self):
        """Update the number of triggers
        """
        self._num_triggers = len(self._triggers)
        num_triggers = '{0:5}'.format(self._num_triggers)
        num_selected = '{0:5}'.format(self._num_selected)
        text = '{0} of {1} triggers selected'.format(num_selected, num_triggers)
        self._widgets['num_triggers'].config(text=text)
    
    def _set_num_selected(self, event):
        """Get the number of selected rows and update
        """
        table = self._widgets['table']
        self._selected_ids = table.selection()
        self._num_selected = len(self._selected_ids)
        self._set_num_triggers()
    
    def _toggle_loctype(self):
        """Toggle the location type from decimal to sexagesimal and vice versa
        """
        table = self._widgets['table']
        rows = table.get_children()
        if self._loctype.get() == 'dec':
            for i in range(self._num_triggers):
                self._triggers[i].use_decimal()
                table.set(rows[i], column='#3', value=self._triggers[i].ra)
                table.set(rows[i], column='#4', value=self._triggers[i].dec)
        elif self._loctype.get() == 'sex':
            for i in range(self._num_triggers):
                self._triggers[i].use_sexagesimal()
                table.set(rows[i], column='#3', value=self._triggers[i].ra)
                table.set(rows[i], column='#4', value=self._triggers[i].dec)
    
    def _toggle_timetype(self):
        """Toggle the time type from UTC to MET and vise versa.  Updates the column
        header and the time filter as well.
        """
        table = self._widgets['table']
        rows = table.get_children()
        if self._timetype.get() == 'met':
            for i in range(self._num_triggers):
                table.set(rows[i], column='#2', value=self._triggers[i].met)
            table.heading('#2', text='MET', command=lambda: self._sort('#2', 'met'))
            self._widgets['timefilter'].config(textvariable=self._filters['met'])
        elif self._timetype.get() == 'utc':
            for i in range(self._num_triggers):
                table.set(rows[i], column='#2', value=self._triggers[i].utc)
            table.heading('#2', text='UTC', command=lambda: self._sort('#2', 'utc'))
            self._widgets['timefilter'].config(textvariable=self._filters['utc'])
    
    def _get_data(self):
        """Get the data for the selected triggers.  This brings up a dialog box for each
        trigger selected showing the available data and the detector angles.  Once the 
        data is selected, runs the download manager
        """
        # disable the button so that errant clicks aren't collected
        self._widgets['get_data'].config(state=DISABLED)
        
        self._downloads = {}
        table = self._widgets['table']
        for id in self._selected_ids:
            index = table.index(id)
            ra = self._triggers[index]._ra_decimal()
            dec = self._triggers[index]._dec_decimal()
            bnum = table.item(id)['text'][2:]
            # the data selector
            downloader = DataSelector(bnum, ra, dec, parent = self.master, 
                                      callback=self._selection_from_popup)
            # collect the selected data and metadata
            try:
                self._downloads[bnum] = self._download
                metadata = {'ra': self._triggers[index]._ra, 
                            'dec': self._triggers[index]._dec, 
                            'utc': self._triggers[index].utc, 
                            'type': self._triggers[index].type, 
                            'loc_err': self._triggers[index].loc_err, 
                            'loc_instr': self._triggers[index].loc_instr}
                self._downloads[bnum]['metadata'] = metadata
            except:
                # if this fails, then the user has canceled the data selection
                pass
        if len(self._downloads) == 0:
            self._widgets['get_data'].config(state=NORMAL)
            return
        # once all of the selections have been made, run the download manager
        dir = filedialog.askdirectory()
        if dir == '':
            self._widgets['get_data'].config(state=NORMAL)
            return
        downloader = DownloadManager(self._downloads, dir=dir, parent=self.master)
        
        # re-enable the button 
        self._widgets['get_data'].config(state=NORMAL)
    
    def _selection_from_popup(self, selection):
        """Callback for DataSelector
        """       
        self._download = selection
    
 
class DataSelector(object):
    """Class for a data selection widget which displays the data available for
    a trigger and allows the user to select which data to download
    """      
    def __init__(self, bnum, ra, dec, parent=None, callback=None):
        """Construct the widget
        
        Parameters:
        -----------
        bnum: str
            The burst number
        ra: str
            The RA of the source
        dec: str
            The Dec of the source
        parent: widgetID, optional
            The parent widget
        callback: function, optional
            A callback function which accepts the data selection from this 
            class once this widget is destroyed.
        """
        self._parent = parent
        self._bnum = bnum
        self._callback = callback
        self._datavars = {}
        self._detvars = {}
        # get detector angles to the source location
        self._datafinder = TriggerDataFinderFTP(download_dir='.')
        self._angles = self._calc_det_angles(ra, dec)
        # number of each datatype available
        num_data = self._num_files()
        
        # toplevel widget - acts as a popup
        self.top = Toplevel(parent)
        self.top.grab_set()
        # center the popup
        self.top.title('Select data for {0}'.format(bnum))
               
        
        # data type selection
        dataframe = LabelFrame(self.top, relief=SUNKEN, text='Data Types')
        dataframe.grid(row=0, column=0, sticky=E+W+N+S)
        self._databuttons = {}
        dtypes = ('CTIME', 'CSPEC', 'TTE', 'RSP', 'RSP2')
        for i in range(len(dtypes)):
            self._data_button(dtypes[i], dataframe, num_data[i], i)
        # all available data types
        alldataframe = LabelFrame(dataframe, relief=RIDGE)
        alldataframe.grid(row=5, sticky=N+S+W+E)
        self._datavars['ALL'] = IntVar(value=0)
        alldatabutton = Checkbutton(alldataframe, text='All Data Types', 
                                    variable=self._datavars['ALL'], 
                                    command=self._all_data)
        alldatabutton.grid(row=0, column=0, sticky=E+W+N+S)
        self._datavars['ALL'].var = alldatabutton
        
        # detector selection
        detframe = LabelFrame(self.top, relief=SUNKEN, text='Detectors')
        detframe.grid(row=0, column=1, sticky=W)
        dets = ('n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'b0')
        for i in range(len(dets)):
            self._det_button(dets[i], detframe, i, 0)
        spacer = Label(detframe, width=3)
        spacer.grid(row=0, column=1, rowspan=7)       
        dets = ('n6', 'n7', 'n8', 'n9', 'na', 'nb', 'b1')
        for i in range(len(dets)):
            self._det_button(dets[i], detframe, i, 2)
        # all detectors
        alldetframe = LabelFrame(detframe, relief=RIDGE)
        alldetframe.grid(row=7, column=0, columnspan=3, sticky=E+W)
        self._detvars['ALL'] = IntVar(value=0)
        alldetbutton = Checkbutton(alldetframe, text='All Detectors', 
                                   variable=self._detvars['ALL'], 
                                   command=self._all_dets)
        alldetbutton.grid(row=0, column=0, columnspan=1, sticky=E+W+N+S)
        self._detvars['ALL'].var = alldetbutton
        
        # get the data
        donebutton = Button(self.top, text='Download', height=2, 
                            highlightbackground='dark gray', 
                            command=self._get_data)
        donebutton.grid(row=1, column=0, columnspan=2, sticky=E+W)
        

        # force parent to wait
        self._parent.wait_window(self.top)
            
    def _data_button(self, datatype, parent, num_data, row):
        """Create and insert a datatype selector button
        
        Parameters:
        -----------
        datatype: str
            The datatype associated with the button
        parent: widgetID
            The ID of the parent widget in which to place the button
        num_data: int
            The number of files available for the datatype
        row: int
            The row number at which to place the button
        """
        text = '{0} ({1} files)'.format(datatype, num_data)
        state = DISABLED if num_data == 0 else NORMAL
        self._datavars[datatype] = IntVar()
        self._datavars[datatype].set(0)       
        b = Checkbutton(parent, text=text, state=state, 
                        variable=self._datavars[datatype])
        b.grid(row=row, sticky=W)
        self._datavars[datatype].var = b
     
    def _det_button(self, det, parent, row, column):
        """Create and insert a detctor selector button
        
        Parameters:
        -----------
        dete: str
            The detector associated with the button
        parent: widgetID
            The ID of the parent widget in which to place the button
        row: int
            The row number at which to place the button
        column: int
            The column number at which to place the button
        """
        text = '{0} ({1} deg)'.format(det.upper(), self._angles[det])
        self._detvars[det] = IntVar()
        # default select button if NaI angle <= 60 or BGO angle <= 90
        if det[0] == 'n':
            self._detvars[det].set(1) if float(self._angles[det]) <= 60.0 \
            else self._detvars[det].set(0)
        else:
            self._detvars[det].set(1) if float(self._angles[det]) <= 90.0 \
            else self._detvars[det].set(0)          
        b = Checkbutton(parent, text=text, variable=self._detvars[det])
        b.grid(row=row, column=column, sticky=W)
        self._detvars[det].var = b
    
    def _all_dets(self):
        """Select/deselect all detectors
        """       
        for det in self._detvars.keys():
            self._detvars[det].set(self._detvars['ALL'].get())
        
    def _all_data(self):
        """Select/deselect all data
        """       
        for datatype in self._datavars.keys():
            state = self._datavars[datatype].var.cget('state')
            if state != 'disabled':
                self._datavars[datatype].set(self._datavars['ALL'].get())
             
    def _get_data(self):
        """Collects the data selections and sends it to the callback.  
        Then destroys this popup.
        """
        # get selected datatypes and detectors  
        datatypes = [key for key in self._datavars.keys() if 
                     self._datavars[key].get() == 1]
        datatypes = [d for d in datatypes if d != 'ALL']
        dets = [key for key in self._detvars.keys() if 
               self._detvars[key].get() == 1]
        dets = [d for d in dets if d != 'ALL']
        
        # call the callback with the selection and destroy
        download = {'datatypes': datatypes, 'detectors': dets}
        if self._callback is not None:
            self._callback(download)
        self.top.destroy()
    
    def _calc_det_angles(self, ra, dec):
        """Calculate the detector angles to the source location
        
        This requires the trigdat file to be downloaded and the trigger
        quaternion is read.
        
        Parameters:
        -----------
        ra: float
            The source RA
        dec: float
            The source Dec
        
        Returns:
        --------
        det_angles: dict
            The string-formatted detector angles
        """
        # retrieve and download trigdat to get quaternion
        trigdat_file = self._datafinder.trigdat_file(self._bnum)
        trigdat_file = self._datafinder.ftp_get(trigdat_file)[0]
        with fits.open(trigdat_file) as hdulist:
            quaternion = hdulist['TRIGRATE'].data['SCATTITD'][0]
        
        # calculate angle between the source location and each detector normal
        det_angles = {}
        for det in Detector:
            det_ra, det_dec = spacecraft_to_radec(det.azimuth, det.zenith, 
                                                  quaternion)
            angle = haversine(ra, dec, det_ra, det_dec)
            det_angles[det.short_name] = '{0:5.1f}'.format(angle)
        
        os.remove(trigdat_file)
        return det_angles
    
    def _num_files(self):
        """Get the number of files available for each datatype.

        Returns:
        --------
        num_ctimes: int
            The number of available CTIME files
        num_cspecs: int
            The number of available CSPEC files
        num_ttes: int
            The number of available TTE files
        num_rsps: int
            The number of available RSP files
        num_rsp2s: int
            The number of available RSP2 files
        """
        num_ctimes = len(self._datafinder.ctime_files(self._bnum))
        num_cspecs = len(self._datafinder.cspec_files(self._bnum))
        num_ttes = len(self._datafinder.tte_files(self._bnum))
        num_rsps = len(self._datafinder.rsp_files(self._bnum, ctime=True, 
                                                  cspec=True))
        num_rsp2s = len(self._datafinder.rsp_files(self._bnum, ctime=True, 
                                                   cspec=True, rsp2=True))
        return (num_ctimes, num_cspecs, num_ttes, num_rsps, num_rsp2s)


class DownloadManager(object):
    """Class for a trigger data download manager.  Manages the data downloads of one or
    many triggers and displays the progress in a widget.
    """          
    def __init__(self, downloads, dir=None, parent=None):
        """Construct the widget
        
        Parameters:
        -----------
        downloads: dict
            A dictionary of data to download
        parent: widgetID, optional
            The parent widget
        """
        #s = ttk.Style()
        #s.theme_use("clam")
        #s.configure("TProgressbar", thickness=50)
        self._parent = parent
        self._datafinder = TriggerDataFinderFTP()
        #self._rootdir = '../test'
        self._rootdir = dir
        self._multitrig = False
        
        # initialize the download mananger
        self.top = Toplevel(self._parent)
        self.top.title('GBM Trigger Download Manager')
        self.top.grab_set()
        
        # if we have more than one trigger that we are downloading, we want to keep 
        # track of the total progress
        irow = 0
        if len(downloads) > 1:
            alllabel = Label(self.top, text='Overall Progress:')
            alllabel.grid(row=0, sticky=E+W)
            self.allprogress = ttk.Progressbar(self.top, orient='horizontal', length=500, 
                                               mode='determinate')
            self.allprogress.grid(row=1, sticky=E+W)
            self._multitrig = True
            self.allprogress['value'] = 0.0
            self.allprogress['maximum'] = float(len(downloads))
            irow = 2
        
        # the progress of each trigger   
        self.triggerlabel = Label(self.top, text='')
        self.triggerlabel.grid(row=irow, sticky=E+W)
        self.trigprogress = ttk.Progressbar(self.top, orient='horizontal', length=500, 
                                            mode='determinate')
        self.trigprogress.grid(row=irow+1, sticky=E+W)
        # the progress of each file
        self.filelabel = Label(self.top, text='')
        self.filelabel.grid(row=irow+2, sticky=W)
        self.fileprogress = ttk.Progressbar(self.top, orient='horizontal', length=500, 
                                            mode='determinate')
        self.fileprogress.grid(row=irow+3, sticky=E+W)
        self.top.update_idletasks()
        
        # cycle through the download list
        i = 0
        for bnum in downloads.keys():
            self._create_metadata_file(bnum, downloads[bnum])
            self._download_set(bnum, downloads[bnum])
            if self._multitrig:
                self.allprogress['value'] = float(i+1.0)
            i += 1
        
        # die mf
        self.top.destroy()
        
    def _construct_path(self, bnum):
        """Constructs the download path. Creates the path if it doesn't exist
        
        Parameters:
        -----------
        bnum: str
            The burst number
        """
        #dir = os.path.join(self._rootdir, 'triggers', bnum, 'data')
        dir = os.path.join(self._rootdir, bnum)
        if not os.path.exists(dir):
            os.makedirs(dir)
        return dir
    
    def _create_metadata_file(self, bnum, download):
        """Create the metadata file for a trigger download.
        
        Parameters:
        -----------
        bnum: str
            The burst number
        download: dict:
            The download dictionary
        """
        dir = self._construct_path(bnum)
        metadata = download['metadata']
        np.save(os.path.join(dir, '.metadata'), np.array([metadata]))
    
    def _download_set(self, bnum, download):
        """Download all the requested data for a single trigger
        
        Parameters:
        -----------
        bnum: str
            The burst number
        download: dict:
            The download dictionary
        """
        # update with the current trigger
        self.triggerlabel.config(text='Progress: {0}'.format(bnum))
        dir = self._construct_path(bnum)
        
        datatypes = download['datatypes']
        detectors = download['detectors']
        
        # collect the file paths of all of the files on the FTP server
        allfiles = []
        for datatype in datatypes:
            if datatype == 'CTIME':
                files = self._datafinder.ctime_files(bnum, dets=detectors)
            elif datatype == 'CSPEC':
                files = self._datafinder.cspec_files(bnum, dets=detectors)
            elif datatype == 'TTE':
                files = self._datafinder.tte_files(bnum, dets=detectors)
            elif (datatype == 'RSP') | (datatype == 'RSP2'):
                ctime = False
                cspec = False
                rsp2 = False
                if any([d == 'CTIME' for d in datatypes]):
                    ctime = True
                if any([(d == 'CSPEC') | (d == 'TTE') for d in datatypes]):
                    cspec = True
                if datatype=='RSP2':
                    rsp2 = True
                files = self._datafinder.rsp_files(bnum, dets=detectors, ctime=ctime, 
                                                   cspec=cspec, rsp2=rsp2)
            allfiles.extend(files)
            
        # download each file
        num_files = len(allfiles)
        self.trigprogress['value'] = 0.0
        self.trigprogress['maximum'] = float(num_files)
        for i in range(num_files):
            self._download_file(dir, allfiles[i])
            self.trigprogress['value'] = float(i+1.0)
            
    def _download_file(self, dir, file):
        """Download a data file
        
        Parameters:
        -----------
        dir: str
            The destination directory
        download: str
            The remote file to be downloaded
        """
        filename = os.path.basename(file)
        self._file = os.path.join(dir, filename)
        if os.path.exists(self._file):
            return
        # update with the current filename
        self.filelabel.config(text='')
        self.filelabel.config(text=filename)
        # initialize the file progress bar            
        self._datafinder._ftp.voidcmd('TYPE I')
        self.fileprogress['value'] = 0
        self.fileprogress['maximum'] = self._datafinder._ftp.size(file)
        # download the file
        self._datafinder._ftp.retrbinary('RETR '+file, self._callback_file)
    
    def _callback_file(self, chunk):
        """Callback function for the FTP download.  Writes the file to disk and updates
        the file progress bar
        
        Parameters:
        -----------
        chunk: str
            The byte data to be written
        """
        # write a chunk of the file
        with open(self._file, 'ab') as f:
            f.write(chunk)
        # update the file progress bar
        self.fileprogress['value'] += len(chunk)
        self.fileprogress.update_idletasks()
        # update the trigger progress
        self._callback_trigger()
    
    def _callback_trigger(self):
        """Updates the trigger progress bar
        """
        frac = float(self.fileprogress['value'])/float(self.fileprogress['maximum'])
        self.trigprogress['value'] = floor(self.trigprogress['value'])+frac
        self.trigprogress.update_idletasks()
        # if we have more than one trigger to download, update the total progress bar
        if self._multitrig:
            self._callback_all()

    def _callback_all(self):
        """Updates the total data download progress bar
        """
        frac = float(self.trigprogress['value'])/float(self.trigprogress['maximum'])
        self.allprogress['value'] = floor(self.allprogress['value'])+frac
        self.allprogress.update_idletasks()

def main():
    root = Tk()
    gui = TriggerFinder(root)
    root.mainloop()

if __name__ == '__main__':
	main()
