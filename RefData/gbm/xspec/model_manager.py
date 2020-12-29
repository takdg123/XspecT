import os
from collections import OrderedDict

class XspecModelManager(object):
    """Class to manage the XSPEC models in GSpec
    
    Attributes:
    -----------
    num_models:
        The number of models
    
    Public Methods:
    ---------------
    model_names:
        List all the available models
    parameter_names:
        List the parameters for a given model
    parameter_settings:
        List the parameter settings for a given model
    remove_model:
        Remove a model from the manager
    set_parameter_default:
        Set the default value of a parameter for a given model
    set_parameter_frozen:
        Freeze or thaw a parameter for a given model
    update_parameter_settings:
        Update the parameter settings for a given model
    xspec_name:
        The XSPEC name for the model
    """    
    # know where the models list is located
    _dir = os.path.dirname(os.path.realpath(__file__))
    _modelfile = os.path.join(_dir, 'gspec_models.dat')
    
    # initialize models list and model dictionary 
    models = []
    _param_dict = [('name', None),
                  ('units', None),
                  ('default', None),
                  ('hard min', None),
                  ('soft min', None),
                  ('soft max', None),
                  ('hard max', None),
                  ('fit delta', None),
                  ('freeze', None)]
    _param_dict = OrderedDict(_param_dict)
    
    _model_dict = [('name', None),
                  ('xspec name', None),
                  ('nparams', None), 
                  ('energy lo', None),
                  ('energy hi', None),
                  ('subroutine', None),
                  ('type', None), 
                  ('model variances', None),
                  ('calculate all', None), 
                  ('params', None)]
    _model_dict = OrderedDict(_model_dict)
    
    def __init__(self):
        self._read_gspec_models()
        self.num_models = len(self.models)
            
    def model_names(self):
        """Return a list of all model names
        
        Returns:
        -----------
        models: list
            The model names
        """
        models = [model['name'] for model in self.models]
        return models
    
    def parameter_names(self, model_name):
        """Return a list of all parameter names for a model
        
        Parameters:
        -----------
        model_name: str
            A valid model name

        Returns:
        -----------
        parameters: list
            The parameter names for the model
        """
        index = self._model_index(model_name)
        params = [param['name'] for param in self.models[index]['params']]
        return params
    
    def parameter_settings(self, model_name):
        """Returns the full parameter settings for a model
        
        Parameters:
        -----------
        model_name: str
            A valid model name

        Returns:
        -----------
        settings: list of OrderedDict
            The parameter settings for the model
        """
        index = self._model_index(model_name)
        settings = self.models[index]['params']
        return settings
    
    def update_parameter_settings(self, model_name, settings):
        """Update the parameter settings
        
        Parameters:
        -----------
        model_name: str
            A valid model name
        settings: list of OrderedDict
            The list of parameter settings
        """
        index = self._model_index(model_name)
        self.models[index]['params'] = settings        
    
    def remove_model(self, model_name):
        """Removes a model from the manager
        
        Parameters:
        -----------
        model_name: str
            The model to be removed
        """
        index = self._model_index(model_name)
        self.models.pop(index)
        self.num_models -= 1
    
    def set_parameter_default(self, model_name, parameter, value):
        """Set the default value for a model parameter
        
        Parameters:
        -----------
        model_name: str
            A valid model name
        parameter: str
            The model parameter name
        value: float
            The default value
        """
        index = self._model_index(model_name)
        for param in self.models[index]['params']:
            if param['name'] == parameter:
                param['default'] = value
                break

    def set_parameter_frozen(self, model_name, parameter, frozen):
        """Set the frozen state for the model parameter
        
        Parameters:
        -----------
        model_name: str
            A valid model name
        parameter: str
            The model parameter name
        frozen: bool
            Set to True to freeze the parameter, False to let it be free
        """
        index = self._model_index(model_name)
        for param in self.models[index]['params']:
            if param['name'] == parameter:
                param['freeze'] = int(frozen)
                break
    
    def xspec_name(self, model_name):
        """Return the XSPEC name of the model
        
        Parameters:
        -----------
        model_name: str
            A valid model name

        Returns:
        -----------
        xspec_name: str
            The XSPEC name of the model
        """
        index = self._model_index(model_name)
        xspec_name = self.models[index]['xspec name']
        return xspec_name
    
    def _read_gspec_models(self):
        """Read the model file
        """
        with open(self._modelfile, 'r') as f:
            txt = list(f)
        
        numlines = len(txt)
        for i in range(numlines):
            if txt[i].startswith('#'):
                num = int(txt[i+1].split()[1])
                lines = txt[i:(i+num+2)]
                self.models.append(self._build_model_dict(lines))
    
    def _build_model_dict(self, lines):
        """Construct the model dictionary for a model
        
        Parameters:
        -----------
        lines: list
            Lines of text defining the model
       
        Returns:
        -----------
        new_model: OrderedDict
            The model dictionary
        """
        # fill up the model dictionary
        new_model = self._model_dict.copy()
        new_model['name'] = ' '.join(lines[0].split())[1:].strip()
        model_line = lines[1].split()
        model_line = [m.strip() for m in model_line]
        new_model['xspec name'] = model_line[0]
        new_model['nparams'] = int(model_line[1])
        new_model['energy lo'] = float(model_line[2])
        new_model['energy hi'] = float(model_line[3])
        new_model['subroutine'] = model_line[4]
        new_model['type'] = model_line[5]
        new_model['model variances'] = int(model_line[6])
        try:
            new_model['calculate all'] = int(model_line[7])
        except:
            new_model['calculate all'] = 0
        
        # fill up the parameter dictionary for each parameter
        params = []
        lines = lines[2:]
        for line in lines:
            new_param = self._param_dict.copy()
            param_line = line.split()
            param_line = [p.strip() for p in param_line]
            numfields = len(param_line)
            new_param['name'] = param_line[0]
            if param_line[1] == '""':
                param_line[1] = ''
            new_param['units'] = param_line[1]
            new_param['default'] = float(param_line[2])
            new_param['hard min'] = float(param_line[3])
            new_param['soft min'] = float(param_line[4])
            new_param['soft max'] = float(param_line[5])
            new_param['hard max'] = float(param_line[6])
            new_param['fit delta'] = float(param_line[7])
            if new_param['fit delta'] <= 0.0:
                new_param['freeze'] = 1
            else:
                new_param['freeze'] = 0
            params.append(new_param)
        
        new_model['params'] = params
        return new_model
        
    def _model_index(self, model_name):
        """Return the index of the model in the list
        
        Parameters:
        -----------
        model_name: str
            A valid model name

        Returns:
        -----------
        index: int
            The index
        """
        index = -1
        for i in range(self.num_models):
            if self.models[i]['name'] == model_name:
                index = i
                break
        
        if index == -1:
            raise KeyError("'{0}' not a valid model name".format(model_name))
        
        return index
    