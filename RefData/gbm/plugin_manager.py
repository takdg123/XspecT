import os
from collections import OrderedDict

class PluginManager(object):
    """Class to manage plugins for GSpec
    
    Attributes:
    -----------
    
    Public Methods:
    ---------------
    background_plugins:
        List all available background fitting/estimation plugins
    binning_plugins:
        List all available binning plugins
    get_background_method:
        Get the background method for the background plugin
    get_binning_method:
        Get the binning method for the binning plugin
    """    
    def __init__(self):
        dir = os.path.dirname(os.path.realpath(__file__))
        self._BACKGROUND_PATH = os.path.join(dir, 'background')
        self._BINNING_PATH = os.path.join(dir, 'binning')
        
        # read the plugin registries
        path = os.path.join(self._BACKGROUND_PATH, 'plugins_for_binned_data.txt')
        self._binned_bkgd_plugins = self._read_registry(path)
        path = os.path.join(self._BACKGROUND_PATH, 'plugins_for_unbinned_data.txt')
        self._unbinned_bkgd_plugins = self._read_registry(path)
        path = os.path.join(self._BINNING_PATH, 'plugins_for_binned_data.txt')
        self._binned_binning_plugins = self._read_registry(path)
        path = os.path.join(self._BINNING_PATH, 'plugins_for_unbinned_data.txt')
        self._unbinned_binning_plugins = self._read_registry(path)
    
    
    def background_plugins(self, datatype='binned'):
        """List available background fitting/estimation plugins
    
        Parameters:
        -----------
        datatype: str
            Either 'binned' (default) or 'unbinned'
       
        Returns:
        -----------
        list: list
            The registered plugins
        """
        self._check_datatype(datatype)
        if datatype == 'binned':
            return list(self._binned_bkgd_plugins.keys())
        elif datatype == 'unbinned':
            return list(self._unbinned_bkgd_plugins.keys())

    def binning_plugins(self, datatype='binned'):
        """List available binning fitting/estimation plugins
    
        Parameters:
        -----------
        datatype: str
            Either 'binned' (default) or 'unbinned'
       
        Returns:
        -----------
        list: list
            The registered plugins
        """
        self._check_datatype(datatype)
        if datatype == 'binned':
            return list(self._binned_binning_plugins.keys())
        elif datatype == 'unbinned':
            return list(self._unbinned_binning_plugins.keys())
    
    def get_background_method(self, name, datatype='binned'):
        """Import the background method for the given background plugin name
    
        Parameters:
        -----------
        name: str
            The name of the plugin
        datatype: str
            Either 'binned' (default) or 'unbinned'
       
        Returns:
        -----------
        method_ref: class reference
            The reference to the imported background class
        """
        self._check_datatype(datatype)
        if datatype == 'binned':
            d = self._binned_bkgd_plugins
        elif datatype == 'unbinned':
            d = self._unbinned_bkgd_plugins
        
        try:
            method = os.path.join(self._BACKGROUND_PATH, d[name])
        except:
            raise KeyError("'{0}' is not a registered plugin".format(name))
        
        method_ref = self._import_method(method)
        return method_ref

    def get_binning_method(self, name, datatype='binned'):
        """Import the binning method for the given binning plugin name
    
        Parameters:
        -----------
        name: str
            The name of the plugin
        datatype: str
            Either 'binned' (default) or 'unbinned'
       
        Returns:
        -----------
        method_ref: function reference
            The reference to the import binning function
        """
        if datatype == 'binned':
            d = self._binned_binning_plugins
        elif datatype == 'unbinned':
            d = self._unbinned_binning_plugins
        
        try:
            method = os.path.join(self._BINNING_PATH, d[name])
        except:
            raise KeyError("'{0}' is not a registered plugin".format(name))

        method_ref = self._import_method(method)
        return method_ref
        
    def _check_datatype(self, datatype):
        """Check validity of datatype string
        """
        if (datatype != 'binned') & (datatype != 'unbinned'):
            raise ValueError("datatype must be either 'binned' or 'unbinned'")

    def _import_method(self, method):
        """Takes the method path, and imports the method
    
        Parameters:
        -----------
        method: str
            The method to import
       
        Returns:
        -----------
        method_ref:
            The reference to the imported method
        """        
        method = method.replace('/', '.')
        package = '.'.join(method.split('.')[-4:-1])
        function = method.split('.')[-1]
        _temp = __import__(package, globals(), locals(), [function], 0)
        method_ref = getattr(_temp, function)
        return method_ref
        
    def _read_registry(self, file):
        """Read a plugin registry file
    
        Parameters:
        -----------
        file: str
            The file name
       
        Returns:
        -----------
        d: OrderedDict
            A dictionary mapping the registered plugin names to their methods
        """
        with open(file, 'r') as f:
            txt = list(f)
        
        names = []
        methods = []
        for line in txt:
            s = line.split(',')
            if (s[0].strip() == '') | (s[1].strip() == ''):
                continue
            names.append(s[0].strip())
            methods.append(s[1].strip())
        
        d = OrderedDict([(name, method) for (name, method) in zip(*(names, methods))])
        return d
