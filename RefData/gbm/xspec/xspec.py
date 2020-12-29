from .model_manager import XspecModelManager
import gbm.xspec.xspec_interface as Xi
import astropy.io.fits as fits
import time
import os
import numpy
import pickle

def parse_fits(FITSfile):
    my_file = fits.open(FITSfile)
    num_spec = 1
    if my_file[0].header['FILETYPE'] == 'PHAII':
        num_spec = my_file['SPECTRUM'].header['NAXIS2']
    return num_spec


class XSPEC():
    """Drive XSPEC
    
    Attributes:
    -----------
    
    Public Methods:
    ---------------
    do_command:
        Execute any valid XSPEC command
    do_fit:
        Perform a spectral fit
    load_pha:
        Load a data file into XSPEC
    run_script:
        Run an XSPEC script file
    save_script:
        Save the current command stack as a script file
    set_fit_statistic:
        Set the fit statistic
    set_fit_weight:
        Set the fit weight
    set_model:
        Set the model to be fit
    """    
    datafiles = []
    cmd_stack = []
    selnparam = 0
    statistic = 'chi'
    weight = 'standard'
    num_spec = 1
    logger = print
    
    def __init__(self, models=None):
        if models is None:
            models = XspecModelManager()
        self.models = models
        self.my_session = Xi.XspecInterface(self.logger)
        self.my_session.command('query yes')
        self.my_session.response_nowait()
         
    def __del__(self):
        self.my_session.command('cd ..\n')
        self.my_session.command('exit\n')
        self.my_session.capture()
    
    def save_script(self, filename):
        """Save the current command stack as a script file
        
        Parameters:
        --------------
        filename: str
            The filename of the script file
        """
        with open(filename, 'w') as f:
            for command in self.cmd_stack:
                f.write(command+'\n')

    def run_script(self, filename):
        """Run an XSPEC script file
        
        Parameters:
        --------------
        filename: str
            The filename of the script file
        """
        self.my_session.command('@'+filename)
        self.my_session.response_nowait()
    
    def load_pha(self, fullfilename, extension=None, slot=1):
        """Add a data file to XSPEC
        
        Parameters:
        --------------
        fullfilename: str
            The full filename, including path
        extension: int, optional
            The extension number to read if there are multiple extensions
        slot: int, optional
            The slot in which the file is loaded.  Defaults to 1
        """
        # Get the directory name from the selected file
        directory = os.path.dirname(fullfilename)
        self.directory = directory

        # Background and Response filenames saved in PHA files don't have a 
        # path; we must chdir to wherever our data file sits:
        if directory:
            self.my_session.command('cd %s\n' % directory)
            self.my_session.response_nowait()

        # Get the file's basename
        my_file = os.path.basename(fullfilename)

        # Check if we have a PHAII file (num_specs > 1)
        num_specs = parse_fits(fullfilename)

        # print(fullfilename)
        # print(num_specs)
        # print(self.num_spec)

        # If so, does this file have the same number of spectra as any prev. read?
        # if len(self.datafiles) > 0:
        #     if num_specs != self.num_spec:
        #         self.logger(text='**********************************************\n')
        #         self.logger(text='*** Nonmatching number of spectra in file! ***\n')
        #         self.logger(text='***       Found %i, but expected %i!       ***\n' % (num_specs, self.num_spec))
        #         self.logger(text='***       Will reset all data files!       ***\n')
        #         self.logger(text='**********************************************\n')
        #         self.datafiles = []

        self.num_spec = num_specs
        self.datafiles.append(my_file)
        # Loop through each of the data files and 
        if (num_specs > 1):
            command = 'data %i %s{%s}\n' % (slot, my_file, extension)
            self.cmd_stack.extend(command.split('\n'))
            self.my_session.command(command)

        else:
            #command = 'data %i %s\n' % (slot, my_file)
            command = 'data {0}:{1} {2}\n'.format(slot, extension, my_file)
            self.cmd_stack.extend(command.split('\n'))
            self.my_session.command(command)

        self.my_session.response_nowait()
    
    def load_phaii_files(self, filenames):
        # Get the directory name from the selected file
        directory = os.path.dirname(filenames[0])
        self.directory = directory

        # Background and Response filenames saved in PHA files don't have a 
        # path; we must chdir to wherever our data file sits:
        if directory:
            self.my_session.command('cd %s\n' % directory)
            #self.my_session.command('cd %s\nchatter 0\n' % directory)
            self.my_session.response_nowait()
        
        num_specs = parse_fits(filenames[0])
        num_files = len(filenames)
        command = ''
        spec_num = 1
        for ispec in range(num_specs):
            command += 'data'
            for ifile in range(num_files):
                my_file = os.path.basename(filenames[ifile])
                s = '{'+str(ispec+1)+'}'
                command += ' {0}:{1} {2}{3}'.format(ispec+1, spec_num+ifile, my_file, s)
            command += '\n'
            spec_num += (ifile+1)
        self.cmd_stack.extend(command.split('\n'))
        self.my_session.command(command)
        
        self.my_session.response_nowait()
        
    
    def set_fit_statistic(self, statistic):
        """Set the XSPEC fit statistic
        
        Parameters:
        --------------
        statistic: str
            The fit statistic name
        """
        self.statistic = statistic
        command = 'statistic %s\n' % statistic
        self.cmd_stack.extend(command.split('\n'))
        self.my_session.command(command)
        self.my_session.response_nowait()

    def set_fit_weight(self, weight):
        """Set the XSPEC fit weight
        
        Parameters:
        --------------
        weight: str
            The fit weight name
        """
        self.weight = weight
        command = 'weight %s\n' % weight
        self.cmd_stack.extend(command.split('\n'))
        self.my_session.command(command)
        self.my_session.response_nowait()

    def do_command(self, command):
        """Execute an XSPEC command
        
        Parameters:
        --------------
        command: str
            Any valid XSPEC command
        """
        if len(command) != 0:
            self.cmd_stack.append(command)
        self.my_session.command('%s\n' % command)
        output=self.my_session.response_nowait()
        return output

    def do_fit(self, num_iter=100):
        """Perform the spectral fit
        
        Parameters:
        --------------
        num_iter: int, optional
            The number of fit iterations. Default is 100
        """
        command = 'fit {0}\n'.format(num_iter)
        self.cmd_stack.extend(command.split('\n'))        
        self.my_session.command(command)
        self.my_session.response_nowait()
    
    def set_energy_range(self, e_lo, e_hi, num_points=100):
        """Set the energy range for the flux calculation.
        Extends the energy range so that a flux calculation can be 
        performed over the desired range.
        
        Note: XSPEC nullifies any fit results once this is issued, therefore
        the user should set this first and then do a fit.

        Parameters:
        --------------
        e_lo: float
            The lower bound of the energy range
        e_hi: float
            The upper bound of the energy range
        num_points: int, optional
            The number of points used in the extension.  Default is 100
        """
        # Because XSPEC only calculates the flux over the response energy range
        # instead of based on an arbitrary range, we have to set the energy 
        # range before calculating the flux.  Even worse, XSPEC requires you
        # to specify this before doing a fit even though extending the 
        # energy range beyond the response energies has absolutely no bearing
        # on the fit itself.  So, to not be wasteful, we need to do this prior
        # to any fit so that we don't have to re-fit just to calculate the flux.
        cmd = 'energies extend low {0} {1} log\n'.format(e_lo, num_points)
        cmd+= 'energies extend high {0} {1} log\n'.format(e_hi, num_points)
        self.cmd_stack.extend(cmd.split('\n'))
        self.my_session.command(cmd)
        self.my_session.response_nowait()
        
    
    def set_model(self, selected_model, numgroups=1):
        """Set the model to be fitted.
        This reads and sets the default values and frozen/thawed states of the
        parameters for the model. This information is contained in the 
        self.models object and can be set via that object.
        
        Parameters:
        --------------
        selected_model: list
            A list of one or more components that make up the model
        """        
        # get the xspec function names and initialize the Band Epeak and 
        # Comptonized Epeak if either are one of the functions
        num_functions = len(selected_model)
        functions = [self.models.xspec_name(function) for function in selected_model]
        for function in functions:
            if function == 'ngrbep':
                self.set_model_band_epeak()
            elif function == 'cplep':
                self.set_model_comptonized()
        
        # join the functions together into a model
        my_command = 'model {0}\n'.format('+'.join(functions))
        
        # initialize settings for each parameter
        num_params = [len(self.models.parameter_names(function))+1 for function in selected_model]
        params = [self.models.parameter_settings(function) for function in selected_model]
        for igroup in range(numgroups):
            for i in range(num_functions):
                for j in range(num_params[i]-1):
                    my_command += '{0} {1} {2} {3} {4} {5}\n'.format(params[i][j]['default'], 
                                  params[i][j]['fit delta'], params[i][j]['hard min'], 
                                  params[i][j]['soft min'], params[i][j]['soft max'], 
                                  params[i][j]['hard max'])        
                my_command += '{0}\n'.format(igroup/10.0)
        self.cmd_stack.extend(my_command.split('\n'))
        self.my_session.command(my_command, add_terminator=True)
        self.my_session.response_nowait()
        
        # freeze or thaw the parameters
        state_command = ''
        for i in range(num_functions):
            for j in range(num_params[i]-1):
                frozen = params[i][j]['freeze']
                if frozen:
                    state_command += 'freeze {0}\n'.format(j+1)
                else:
                    state_command += 'thaw {0}\n'.format(j+1)

        self.cmd_stack.extend(state_command.split('\n'))
        self.my_session.command(state_command+'\n', add_terminator=True)
        self.my_session.response_nowait()
        
        self.num_params = num_params
    
    def set_model_band_epeak(self):
        """Reparametrize the xspec Band function to be the useful Band function
        parametrized by Epeak
        """
        my_command = 'mdef ngrbep grbm(alpha, beta, tem/(2. + alpha))\n'
        self.cmd_stack.extend(my_command.split('\n'))
        self.my_session.command(my_command, add_terminator=True)
        self.my_session.response_nowait()
    
    def set_model_comptonized(self):
        """Reparametrize the xspec CutOffPL function to be the useful 
        Comptonized function parametrized by Epeak
        """
        my_command = 'mdef cplep cutoffpl(-PhoIndex, foldE/(2.0+PhoIndex))\n'
        self.cmd_stack.extend(my_command.split('\n'))
        self.my_session.command(my_command, add_terminator=True)
        self.my_session.response_nowait()
    
    def get_counts(self, ndets):
        """Retrieve the convolved counts for each detector
        
        Parameters:
        -----------
        ndets: int
            Number of detectors        
        
        Returns:
        -----------
        counts: list
            The convolved counts in each detector
        counts_err: list
            The errors on the convolved counts
        counts_model: list
            The fitted counts model for each detector
        """       
        in_cmd, out_cmd = self._counts_from_xspec(ndets)
        self.my_session.command(in_cmd+'\n')
        self.my_session.capture_nowait()

        self.my_session.command(out_cmd+'\n')
        lines = self.my_session.capture_nowait()
        lines = self._clean_xspec_output(lines, out_cmd)
        return self._format_counts(lines, ndets)
        
    def get_covariance(self, nparams):
        """Retrieve the covariance matrix
        
        Parameters:
        -----------
        nparams: int
            Number of fit parameters      
        
        Returns:
        -----------
        covariance: np.array
            The covariance matrix
        """       
        in_cmd, out_cmd = self._covariance_from_xspec(nparams)
        self.my_session.command(in_cmd+'\n')
        self.my_session.capture_nowait()

        self.my_session.command(out_cmd+'\n')
        lines = self.my_session.capture_nowait()
        return self._format_covariance(lines[1], nparams)

    def get_flux(self, e_lo, e_hi):
        """Retrieve the energy and photon flux over the requested energy range
        
        Parameters:
        -----------
        e_lo: float
            The lower bound of the energy range
        e_hi: float
            The upper bound of the energy range
        
        Returns:
        -----------
        eflux: np.array
            The energy flux (centroid, neg. error, pos. error)
        pflux: np.array
            The photon flux (centroid, neg. error, pos. error)
        """       
        in_cmd, out_cmd = self._flux_from_xspec(e_lo, e_hi)
        self.my_session.command(in_cmd+'\n')
        self.my_session.capture_nowait()

        self.my_session.command(out_cmd+'\n')
        lines = self.my_session.capture_nowait()
        lines = self._clean_xspec_output(lines, out_cmd)
        return self._format_flux(lines)
    
    def get_integrated_counts(self, ndets):
        """Retrieve the cumulative convolved counts for each detector
        
        Parameters:
        -----------
        ndets: int
            Number of detectors        
        
        Returns:
        -----------
        counts_int: list
            The cumulative convolved counts in each detector
        counts_int_model: list
            The fitted cumulative counts model for each detector
        """       
        in_cmd, out_cmd = self._counts_integrated_from_xspec(ndets)
        self.my_session.command(in_cmd+'\n')
        self.my_session.capture_nowait()

        self.my_session.command(out_cmd+'\n')
        lines = self.my_session.capture_nowait()
        lines = self._clean_xspec_output(lines, out_cmd)
        return self._format_counts_integrated(lines, ndets)                
    
    def get_model(self):
        """Retrieve the deconvolved photon model

        Returns:
        -----------
        model_x: np.array
            The energy values of the deconvolved model
        model_y: np.array
            The differential flux values of the deconvolved model
        """       
        in_cmd, out_cmd = self._model_from_xspec()
        self.my_session.command(in_cmd+'\n')
        self.my_session.capture_nowait()

        self.my_session.command(out_cmd+'\n')
        lines = self.my_session.capture_nowait()
        return self._format_model(lines[1:])
    
    def get_parameters(self, nparams):
        """Retrieve the fit info for the model parameters
        
        Parameters:
        -----------
        nparams: int
            Number of parameters        
        
        Returns:
        -----------
        param_name: list
            The names of the parameters
        param_values: list
            The values of the parameters
        param_errs: list
            The 1 sigma errors on the parameters
        """       
        in_cmd, out_cmd = self._parameters_from_xspec(nparams)
        self.my_session.command(in_cmd+'\n')
        self.my_session.capture_nowait()

        self.my_session.command(out_cmd+'\n')
        lines = self.my_session.capture_nowait()
        numstr = len(out_cmd)+8
        badlines = int(numpy.ceil(numstr/80.0))
        lines = lines[badlines:]
        return self._format_parameters(lines, nparams)
    
    def get_residuals(self, ndets):
        """Retrieve the residual sigmas of the fit for each detector
        
        Parameters:
        -----------
        ndets: int
            Number of detectors        
        
        Returns:
        -----------
        residuals: list
            The residual sigmas for each detector
        """       
        in_cmd, out_cmd = self._residuals_from_xspec(ndets)
        self.my_session.command(in_cmd+'\n')
        self.my_session.capture_nowait()

        self.my_session.command(out_cmd+'\n')
        lines = self.my_session.capture_nowait()
        lines = self._clean_xspec_output(lines, out_cmd)
        return self._format_residuals(lines)                

    def get_statistic(self):
        """Retrieve the fit statistic info

        Returns:
        -----------
        stat_name: str
            The name of the fit statistic
        stat_value: float
            The value of the fit statistic
        dof: int
            The fit degrees-of-freedom
        """       
        in_cmd, out_cmd = self._statistic_from_xspec()
        self.my_session.command(in_cmd+'\n')
        self.my_session.capture_nowait()

        self.my_session.command(out_cmd+'\n')
        lines = self.my_session.capture_nowait()

        return self._format_statistic(lines[1:]) 
    
    def get_unfolded_model(self, ndets):
        """Retrieve the deconvolved photon model for each detector
        
        Parameters:
        -----------
        ndets: int
            Number of detectors        
        
        Returns:
        -----------
        model_x: list
            The energy values for each detector
        model_xerr: list
            The energy bin widths for each detector
        model_y: list
            The differential flux values for each detector
        model_yerr: list
            The error on the differential flux
        model: list
            The photon model for each detector
        """       
        in_cmd, out_cmd = self._model_unfolded_from_xspec(ndets)
        self.my_session.command(in_cmd+'\n')
        self.my_session.capture_nowait()

        self.my_session.command(out_cmd+'\n')
        lines = self.my_session.capture_nowait()
        lines = self._clean_xspec_output(lines, out_cmd)
        return self._format_model_unfolded(lines, ndets)                
    
    def _counts_from_xspec(self, ndets):
        in_cmd1 = ''
        in_cmd2 = ''
        in_cmd3 = ''
        out_cmd1 = ''
        out_cmd2 = ''
        out_cmd3 = ''
        for i in range(ndets):
            in_cmd1 += 'set counts{0} [tcloutr plot counts y {0}];'.format(i+1)
            in_cmd2 += 'set counts_err{0} [tcloutr plot counts yerr {0}];'.format(i+1)
            in_cmd3 += 'set counts_model{0} [tcloutr plot counts model {0}];'.format(i+1)
            out_cmd1 += 'echo $counts{0};'.format(i+1)
            out_cmd2 += 'echo $counts_err{0};'.format(i+1)
            out_cmd3 += 'echo $counts_model{0};'.format(i+1)
        return (in_cmd1+in_cmd2+in_cmd3, out_cmd1+out_cmd2+out_cmd3)
    
    def _counts_integrated_from_xspec(self, ndets):
        in_cmd1 = ''
        in_cmd2 = ''
        out_cmd1 = ''
        out_cmd2 = ''
        for i in range(ndets):
            in_cmd1 += 'set counts_integrated{0} [tcloutr plot icounts y {0}];'.format(i+1)
            in_cmd2 += 'set counts_integrated_model{0} [tcloutr plot icounts model {0}];'.format(i+1)
            out_cmd1 += 'echo $counts_integrated{0};'.format(i+1)
            out_cmd2 += 'echo $counts_integrated_model{0};'.format(i+1)
        return (in_cmd1+in_cmd2, out_cmd1+out_cmd2)
    
    def _covariance_from_xspec(self):
        in_cmd = 'set covar [tcloutr covariance];'
        out_cmd = 'echo $covar;'
        return (in_cmd, out_cmd)
    
    def _flux_from_xspec(self, e_lo, e_hi):
        in_cmd = 'flux {0} {1} errorsims;'.format(e_lo, e_hi)
        in_cmd+= 'set fluxes_{0}_{1} [tcloutr flux 1];'.format(int(e_lo), int(e_hi))
        out_cmd = 'echo $fluxes_{0}_{1};'.format(int(e_lo), int(e_hi))
        return (in_cmd, out_cmd)
    
    def _model_from_xspec(self):
        in_cmd = 'set model_x [tcloutr plot model x];'
        in_cmd+= 'set model_y [tcloutr plot model y];'
        
        out_cmd = 'echo $model_x;'
        out_cmd+= 'echo $model_y;'
        return (in_cmd, out_cmd)
    
    def _model_unfolded_from_xspec(self, ndets):
        in_cmd1 = ''
        in_cmd2 = ''
        in_cmd3 = ''
        in_cmd4 = ''
        in_cmd5 = ''
        out_cmd1 = ''
        out_cmd2 = ''
        out_cmd3 = ''
        out_cmd4 = ''
        out_cmd5 = ''
        for i in range(ndets):
            in_cmd1 += 'set model_uf_x{0} [tcloutr plot ufspec x {0}];'.format(i+1)
            in_cmd2 += 'set model_uf_xerr{0} [tcloutr plot ufspec xerr {0}];'.format(i+1)
            in_cmd3 += 'set model_uf_y{0} [tcloutr plot ufspec y {0}];'.format(i+1)
            in_cmd4 += 'set model_uf_yerr{0} [tcloutr plot ufspec yerr {0}];'.format(i+1)
            in_cmd5 += 'set model_uf{0} [tcloutr plot ufspec model {0}];'.format(i+1)
        
            out_cmd1 += 'echo $model_uf_x{0};'.format(i+1)
            out_cmd2 += 'echo $model_uf_xerr{0};'.format(i+1)
            out_cmd3 += 'echo $model_uf_y{0};'.format(i+1)
            out_cmd4 += 'echo $model_uf_yerr{0};'.format(i+1)
            out_cmd5 += 'echo $model_uf{0};'.format(i+1)
        return (in_cmd1+in_cmd2+in_cmd3+in_cmd4+in_cmd5, 
                out_cmd1+out_cmd2+out_cmd3+out_cmd4+out_cmd5)
    
    def _parameters_from_xspec(self, nparams):
        in_cmd1 = ''
        in_cmd2 = ''
        in_cmd3 = ''
        out_cmd1 = ''
        out_cmd2 = ''
        out_cmd3 = ''
        for i in range(nparams):
            in_cmd1 += 'set paramname{0} [tcloutr pinfo {0}];'.format(i+1)
            in_cmd2 += 'set paramvalue{0} [tcloutr param {0}];'.format(i+1)
            in_cmd3 += 'set paramerr{0} [tcloutr sigma {0}];'.format(i+1)
            out_cmd1 += 'echo $paramname{0};'.format(i+1)
            out_cmd2 += 'echo $paramvalue{0};'.format(i+1)
            out_cmd3 += 'echo $paramerr{0};'.format(i+1)
        return (in_cmd1+in_cmd2+in_cmd3, out_cmd1+out_cmd2+out_cmd3)
    
    def _residuals_from_xspec(self, ndets, statistic):
        in_cmd = ''
        out_cmd = ''
        for i in range(ndets):
            # Dont' get the delchi residuals if the fit statistic is not chisqr
            if 'chi' in statistic:
                in_cmd += 'set residuals{0} [tcloutr plot delchi y {0}];'.format(i+1)
                out_cmd += 'echo $residuals{0};'.format(i+1)
            else:
                in_cmd += 'set residuals{0} [tcloutr plot residuals y {0}];'.format(i+1)
                out_cmd += 'echo $residuals{0};'.format(i+1)
        return (in_cmd, out_cmd)
    
    def _statistic_from_xspec(self):
        in_cmd = 'set statname [tcloutr statmethod];'
        in_cmd+= 'set statval [tcloutr stat];'
        in_cmd+= 'set dof [tcloutr dof];'

        out_cmd = 'echo $statname;' # fit statistic name
        out_cmd+= 'echo $statval;' # fit statistic value
        out_cmd+= 'echo $dof;' #dof and number of channels
        return (in_cmd, out_cmd)
        
    def _format_counts(self, xspec_out, ndets):
        counts = xspec_out[:ndets]
        counts = [numpy.array(c.split()).astype(float) for c in counts]
        counts_err = xspec_out[ndets:2*ndets]
        counts_err = [numpy.array(c.split()).astype(float) for c in counts_err]
        counts_model = xspec_out[2*ndets:]
        counts_model = [numpy.array(c.split()).astype(float) for c in counts_model]
        return (counts, counts_err, counts_model)
    
    def _format_counts_integrated(self, xspec_out, ndets):
        counts_int = xspec_out[:ndets]
        counts_int = [numpy.array(c.split()).astype(float) for c in counts_int]
        counts_int_model = xspec_out[ndets:]
        counts_int_model = [numpy.array(c.split()).astype(float) for c in counts_int_model]
        return (counts_int, counts_int_model)
    
    def _format_covariance(self, xspec_out, nparams):
        # because XSPEC only returns the unique covariance elements 
        # (as if it's a huge burden to return the actual matrix :facepalm:)
        # we must reconstruct the matrix from the values we are given
        values = xspec_out.split()
        covar = numpy.zeros((nparams,nparams), dtype=float)
        place = 0
        for i in range(nparams):
            covar[i,0:i+1] = values[place:place+i+1]
            place += (i+1)
        mask = (covar == 0.0)
        covar[mask] = (covar.T)[mask]
        return covar
    
    def _format_flux(self, xspec_out):
        s = xspec_out[0].split()
        eflux = numpy.array(s[:3]).astype(float)  
        eflux[1:] = numpy.abs(eflux[0]-eflux[1:])  
        pflux = numpy.array(s[3:]).astype(float)
        pflux[1:] = numpy.abs(pflux[0]-pflux[1:])  
        return (eflux, pflux)    
    
    def _format_model(self, xspec_out):
        model_x = numpy.array(xspec_out[0].split()).astype(float)
        model_y = numpy.array(xspec_out[1].split()).astype(float)
        return (model_x, model_y)
    
    def _format_model_unfolded(self, xspec_out, ndets):
        model_x = xspec_out[:ndets]
        model_x = [numpy.array(m.split()).astype(float) for m in model_x]
        model_xerr = xspec_out[ndets:2*ndets]
        model_xerr = [numpy.array(m.split()).astype(float) for m in model_xerr]
        model_y = xspec_out[2*ndets:3*ndets]
        model_y = [numpy.array(m.split()).astype(float) for m in model_y]
        model_yerr = xspec_out[3*ndets:4*ndets]
        model_yerr = [numpy.array(m.split()).astype(float) for m in model_yerr]
        model = xspec_out[4*ndets:]
        model = [numpy.array(m.split()).astype(float) for m in model]
        return (model_x, model_xerr, model_y, model_yerr, model)

    def _format_parameters(self, xspec_out, nparams):
        param_names = xspec_out[0:nparams]
        param_values = xspec_out[nparams:2*nparams]
        param_values = [float(val.split()[0]) for val in param_values]
        param_errs = xspec_out[2*nparams:]
        param_errs = [float(err) for err in param_errs]
        return (param_names, param_values, param_errs)
    
    def _format_residuals(self, xspec_out):
        residuals = [numpy.array(r.split()).astype(float) for r in xspec_out]
        return residuals
    
    def _format_statistic(self, xspec_out):
        stat_name, stat_value, dof = tuple(xspec_out)
        dof = dof.split()[0]
        return (stat_name, float(stat_value), int(dof))        

    def _clean_xspec_output(self, xspec_out, out_cmd):
        # Seperate commands (which xspec also returns) from output
        xspec_out = [line for line in xspec_out if ';' not in line]
        return xspec_out
    
    def _clean_xspec_output_old(self, xspec_out, out_cmd):
        # command to retrieve output is also returned, but not in a single 
        # list element.  Each list element contains at most 80 characters,
        # but only for the command
        numstr = len(out_cmd)+8
        badlines = int(numpy.ceil(numstr/80.0))
        xspec_out = xspec_out[badlines:]
        # and sometimes (unpredictably?) XSPEC returns empty strings full
        # of ANSI escape characters...
        xspec_out = [line for line in xspec_out if '\x1b' not in line]
        return xspec_out
     
    def extract_scat_data(self, ndets, nparams, time_range):
        """Extract all of the fit info needed for an SCAT file
        
        Parameters:
        -----------
        ndets: int
            Number of detectors   
        nparams: int
            Number of fit parameters  
        time_range: (float, float)
            Th start and stop time of the fit interval (tstart, tstop)
        
        Returns:
        -----------
        data: dict
            A dictionary containing all of the fit info
        """       

        in_cmd = ''
        out_cmd = ''
        in1, out1 = self._statistic_from_xspec()#3
        in2, out2 = self._covariance_from_xspec()#1
        in3, out3 = self._model_from_xspec()#2
        in4, out4 = self._parameters_from_xspec(nparams)#3*nparams
        in5, out5 = self._counts_from_xspec(ndets)#3*ndets
        in6, out6 = self._counts_integrated_from_xspec(ndets)#2*ndets
        in7, out7 = self._model_unfolded_from_xspec(ndets)#5*ndets
        in8, out8 = self._residuals_from_xspec(ndets, self.statistic)#1*ndets
        in9, out9 = self._flux_from_xspec(10.0, 1000.0)#2
        in10, out10 = self._flux_from_xspec(50.0, 300.0)#2
        in_cmd += in1+in2+in3+in4+in5+in6+in7+in8+in9+in10
        out_cmd += out1+out2+out3+out4+out5+out6+out7+out8+out9+out10
        
        self.my_session.command(in_cmd+'\n')
        self.my_session.capture_nowait()
        
        self.my_session.command(out_cmd+'\n')
        lines = self.my_session.capture_nowait()

        # lines = self._clean_xspec_output(lines, out_cmd)
        lines = self._clean_xspec_output(lines, out_cmd)


        stat_name, stat_val, dof = self._format_statistic(lines[0:3])
        
        covariance = self._format_covariance(lines[3], nparams)
        
        model_x, model_y = self._format_model(lines[4:6])
        row1 = 6+3*nparams
        
        pnames, pvals, perrs = self._format_parameters(lines[6:row1], nparams)
        row2 = row1+3*ndets
        
        counts, counts_err, counts_model = self._format_counts(lines[row1:row2], ndets)
        row1 = row2
        row2 = row1+2*ndets
        
        counts_int, counts_int_model = self._format_counts_integrated(lines[row1:row2],
                                                                      ndets)
        row1 = row2
        row2 = row1+5*ndets
        uf_x, uf_xerr, uf_y, uf_yerr, uf = self._format_model_unfolded(lines[row1:row2], ndets)
        
        row1 = row2
        row2 = row1+ndets

        if 'Chi-Squared' not in stat_name:
            residuals = [numpy.array([ 0.]), numpy.array([ 0.]), numpy.array([ 0.])]
        else:
            residuals = self._format_residuals(lines[row1:row2])
    
    
        eflux, pflux = self._format_flux(lines[row2:row2+1])
        eflux_50_300, pflux_50_300 = self._format_flux(lines[row2+1:])
        
        data = {}
        data['counts'] = numpy.array(counts)
        data['counts_error'] = numpy.array(counts_err)
        data['counts_model'] = numpy.array(counts_model)
        data['counts_integrated'] = numpy.array(counts_int)
        data['counts_integrated_model'] = numpy.array(counts_int_model)
        data['x_ufspec'] = numpy.array(uf_x)
        data['xerr_ufspec'] = numpy.array(uf_xerr)
        data['y_ufspec'] = numpy.array(uf_y)
        data['yerr_ufspec'] = numpy.array(uf_yerr)
        data['y_ufspec_model'] = numpy.array(uf)
        data['energies_model'] = numpy.array(model_x)
        data['ufspec_model'] =  numpy.array(model_y)
        data['residuals_sigma'] = numpy.array(residuals)
        data['stat_name'] = stat_name
        data['stat_value'] = stat_val
        data['dof'] = dof
        data['covar'] = covariance
        data['parameter_names'] = pnames
        data['parameter_values'] = numpy.array(pvals)
        data['parameter_sigmas'] = numpy.array(perrs)
        data['energy_flux'] = eflux
        data['photon_flux'] = pflux
        data['energy_flux_50_300'] = eflux_50_300
        data['photon_flux_50_300'] = pflux_50_300
        data['time_range'] = time_range
        
        return data
    

    def extract_plot_data(self, numberOfDataGroups, numberOfParameters, debug=False):
        return self.extract_scat_data(numberOfDataGroups, numberOfParameters)        
        
        data = {}

        data['counts'] = []
        data['counts_error'] = []
        data['counts_model'] = []
        data['counts_integrated'] = []
        data['counts_integrated_model'] = []
        data['x_ufspec'] = []
        data['xerr_ufspec'] = []
        data['y_ufspec'] = []
        data['yerr_ufspec'] = []
        data['y_ufspec_model'] = []
        data['residuals_sigma'] = []
        data['stat_name'] = []
        data['stat_value'] = []
        data['dof'] = []
        data['covar'] = []
        
        for dataGroup in range(numberOfDataGroups):

            dataGroup = dataGroup + 1

            # Set the tcloutr variables
            extract_plot_data_command = ''
            extract_plot_data_command += 'set counts [tcloutr plot counts y %s];' % dataGroup
            extract_plot_data_command += 'set counts_error [tcloutr plot counts yerr %s];' % dataGroup
            extract_plot_data_command += 'set counts_model [tcloutr plot counts model %s];' % dataGroup
            extract_plot_data_command += 'set counts_integrated [tcloutr plot icounts y %s];' % dataGroup
            extract_plot_data_command += 'set counts_integrated_model [tcloutr plot icounts model %s];' % dataGroup
            extract_plot_data_command += 'set energy [tcloutr plot ufspec x %s];' % dataGroup
            extract_plot_data_command += 'set energy_error [tcloutr plot ufspec xerr %s];' % dataGroup
            extract_plot_data_command += 'set flux [tcloutr plot ufspec y %s];' % dataGroup
            extract_plot_data_command += 'set flux_error [tcloutr plot ufspec yerr %s];' % dataGroup
            extract_plot_data_command += 'set flux_model [tcloutr plot ufspec model %s];' % dataGroup
            extract_plot_data_command += 'set energy_model [tcloutr plot model x];'
            extract_plot_data_command += 'set flux_summed_model [tcloutr plot model y];'
            extract_plot_data_command += 'set residuals_sigma [tcloutr plot delchi y %s];' % dataGroup
            
            for parameterNumber in range(numberOfParameters):
                parameterNumber = parameterNumber + 1
                extract_plot_data_command = extract_plot_data_command + 'set parameterName%s [tcloutr pinfo %s];' % (parameterNumber, parameterNumber)

            for parameterNumber in range(numberOfParameters):
                parameterNumber = parameterNumber + 1  
                extract_plot_data_command = extract_plot_data_command + 'set parameterValue%s [tcloutr param %s];' % (parameterNumber, parameterNumber)

            for parameterNumber in range(numberOfParameters):
                parameterNumber = parameterNumber + 1  
                extract_plot_data_command = extract_plot_data_command + 'set parameterSigma%s [tcloutr sigma %s];' % (parameterNumber, parameterNumber)


            extract_plot_data_command = extract_plot_data_command + '\n'

            if debug == True:
                print("\nIssued Command:")
                print(extract_plot_data_command)
            
            self.my_session.command(extract_plot_data_command)
            # time.sleep(0.25)
            # lines = self.my_session.capture()
            lines = self.my_session.capture_nowait()

            if debug == True:
                print("\nReturned:")
                print(lines)

            # Read the tcloutr variables
            extract_plot_data_command = ''
            extract_plot_data_command = extract_plot_data_command + 'echo $counts;'
            extract_plot_data_command = extract_plot_data_command + 'echo $counts_error;'
            extract_plot_data_command = extract_plot_data_command + 'echo $counts_model;'
            extract_plot_data_command = extract_plot_data_command + 'echo $counts_integrated;'
            extract_plot_data_command = extract_plot_data_command + 'echo $counts_integrated_model;'
            extract_plot_data_command = extract_plot_data_command + 'echo $energy;'
            extract_plot_data_command = extract_plot_data_command + 'echo $energy_error;'
            extract_plot_data_command = extract_plot_data_command + 'echo $flux;'
            extract_plot_data_command = extract_plot_data_command + 'echo $flux_error;'
            extract_plot_data_command = extract_plot_data_command + 'echo $flux_model;'
            extract_plot_data_command = extract_plot_data_command + 'echo $energy_model;'
            extract_plot_data_command = extract_plot_data_command + 'echo $flux_summed_model;'
            extract_plot_data_command = extract_plot_data_command + 'echo $residuals_sigma;'
            
            for parameterNumber in range(numberOfParameters):
                parameterNumber = parameterNumber + 1
                extract_plot_data_command = extract_plot_data_command + 'echo $parameterName%s;' % parameterNumber

            for parameterNumber in range(numberOfParameters):
                parameterNumber = parameterNumber + 1  
                extract_plot_data_command = extract_plot_data_command + 'echo $parameterValue%s;' % parameterNumber

            for parameterNumber in range(numberOfParameters):
                parameterNumber = parameterNumber + 1  
                extract_plot_data_command = extract_plot_data_command + 'echo $parameterSigma%s;' % parameterNumber


            extract_plot_data_command = extract_plot_data_command + '\n'

            if debug == True:
                print("\nIssued Command:")
                print(extract_plot_data_command)

            self.my_session.command(extract_plot_data_command)
            # time.sleep(0.25)
            # lines = self.my_session.capture()
            lines = self.my_session.capture_nowait()

            if debug == True:
                print("\nReturned:")
                print(lines)

            lines = lines[3:]

            # Loop through each line and extract the contained values. Pass on anything that cannot be converted to a float
            values = []
            for line in lines:

                # Skip this line if there's a tcl variable character
                if '$' in line:
                    if debug == True:
                        print("Skipping line:")
                        print(line)
                    continue

                try:
                    line = line.replace('\n', '')
                    value = numpy.array(line.split()).astype(float)
                    values.append(value)

                except:
                    value = line.strip()
                    values.append(value)

                finally:
                    # print("Error: unable to convert tclout value to float")
                    pass

            if debug == True:
                print('\nValues')
                print(values)

            # pickle.dump( values, open( "values.p", "wb" ) )

            # Add the extracted values to the data dictionary
            data['counts'].append(values[0])
            data['counts_error'].append(values[1])
            data['counts_model'].append(values[2])
            data['counts_integrated'].append(values[3])
            data['counts_integrated_model'].append(values[4])
            data['x_ufspec'].append(values[5])
            data['xerr_ufspec'].append(values[6])
            data['y_ufspec'].append(values[7])
            data['yerr_ufspec'].append(values[8])
            data['y_ufspec_model'].append(values[9])  #model_ufspec
            data['energies_model'] = values[10]  # energies
            data['ufspec_model'] = values[11]
            data['residuals_sigma'].append(values[12])

            pararmeterNames = []
            for parameterNumber in range(numberOfParameters):
                parameterNumber = parameterNumber + 1
                pararmeterName = values[12 + parameterNumber]
                pararmeterNames.append(pararmeterName)

            # print("Parameter Names:")
            # print(pararmeterNames)
            data['parameter_names'] = pararmeterNames

            pararmeterValues = []
            for parameterNumber in range(numberOfParameters):
                parameterNumber = parameterNumber + 1
                parameterValue = values[12 + numberOfParameters + parameterNumber][0]
                pararmeterValues.append(parameterValue)

            # print("Parameter Values:")
            # print(pararmeterValues)
            data['parameter_values'] = pararmeterValues

            pararmeterSigmas = []
            for parameterNumber in range(numberOfParameters):
                parameterNumber = parameterNumber + 1
                pararmeterSigma = values[12 + (2*numberOfParameters) + parameterNumber][0]
                pararmeterSigmas.append(pararmeterSigma)

            # print("Parameter Sigmas:")
            # print(pararmeterSigmas)
            data['parameter_sigmas'] = pararmeterSigmas
            

        # Make everything into a numpy array
        data['counts'] = numpy.array(data['counts'])
        data['counts_error'] = numpy.array(data['counts_error'])
        data['counts_model'] = numpy.array(data['counts_model'])
        data['counts_integrated'] = numpy.array(data['counts_integrated'])
        data['counts_integrated_model'] = numpy.array(data['counts_integrated_model'])
        data['x_ufspec'] = numpy.array(data['x_ufspec'])
        data['xerr_ufspec'] = numpy.array(data['xerr_ufspec'])
        data['y_ufspec'] = numpy.array(data['y_ufspec'])
        data['yerr_ufspec'] = numpy.array(data['yerr_ufspec'])
        data['y_ufspec_model'] = numpy.array(data['y_ufspec_model'])
        data['energies_model'] = numpy.array(data['energies_model'])
        data['ufspec_model'] = numpy.array(data['ufspec_model'])
        data['residuals_sigma'] = numpy.array(data['residuals_sigma'])

        # data['parameter_names'] = numpy.array(data['parameter_names'])
        data['parameter_values'] = numpy.array(data['parameter_values'])
        data['parameter_sigmas'] = numpy.array(data['parameter_sigmas'])

        # import pickle
        # pickle.dump( data, open( "data.pickle", "wb" ) )
        #numpy.save("data.npy", [data])
        
        if debug == True:
            print("\n")
            print("Fit Parameters")
            for index in range(len(data['parameter_names'])):
                print(data['parameter_names'][index], data['parameter_values'][index], ' +/- ', data['parameter_sigmas'][index])

        if debug == True:
            print('\nData')
            print(data)
        
        return data
