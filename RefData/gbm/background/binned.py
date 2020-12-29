import numpy as np


class polynomial(object):
    """ A polynomial background fitter
    """
    
    def __init__(self, order):
        """
        Parameters:
        --------------
        order: int (optional)
            The order of the polynomial (0 by default)
        """
        order = int(order)
        assert order >= 0, 'Polynomial order must be non-negative' 
        
        self._order = order
        self.coeff = np.zeros(order+1)
        self.covar = np.zeros((order+1, order+1))
        self.chisq = 0.0
        self.dof = 0
    
    def fit(self, counts, exposure, bounds):
        """Fit a set of binned data with a polynomial.  Model variances are used for 
           chi-squared via two fitting passes.  Adapted from Numerical Recipes LFIT.
        
        Parameters:
        counts: np.array
            The counts in each bin
        exposure: np.array
            The exposure of each bin
        bounds: (np.array, np.array)
            The time bounds for the bins
       
        bin_start: np.array
            The starting edge of each bin
        bin_end: np.array
            The ending edge of each bin
        bin_val: np.array
            The value of each bin

        Returns:
        --------------
        np.array
            The fitted model at each input bin
        np.array
            The model uncertainty at each input bin
        """
        
        assert counts.size == exposure.size, 'Counts array and Exposure array ' \
                                             'must be the same length'
        assert counts.size == bounds[0].size, 'Bounds array must be the size of the ' \
                                              'Counts array'
        assert bounds[0].size == bounds[1].size, 'Bounds array must be of consistent ' \
                                                 'size'
        self._bin_start = bounds[0]
        self._bin_end = bounds[1]
        self._bin_val = counts/exposure
        idx = self._bin_val > 0.0
        self._livetime = exposure#np.sum(exposure[idx])
        numBins = self._bin_val.size
        
        # get basis functions and set up weights array
        basisFunc = self._evalBasis(self._bin_start, self._bin_end)
        weights = np.zeros((numBins, self._order+1))
        
        # Two-pass fitting
        # First pass uses the data variances calculated from data rates
        # Second pass uses the data variances calculated from model rates
        # 1) rate * livetime = counts
        # 2) variance of counts = counts
        # 3) variance of rate = variance of counts/livetime^2 = rate/livetime
        for iPass in range(2):
            if np.max(self._bin_val) <= 0.0:
                continue
            
            if iPass == 0:
                variance = self._bin_val/self._livetime
                idx = variance > 0.0
                for iCoeff in range(self._order+1):
                    weights[idx, iCoeff] = basisFunc[iCoeff, idx]/variance[idx]
            else:
                variance = model/self._livetime
                idx = variance > 0.0
                if np.sum(idx) > 0:
                    for iCoeff in range(self._order+1):
                        weights[idx, iCoeff] = basisFunc[iCoeff, idx]/variance[idx]
                else:
                    raise ValueError('SEVERE ERROR: Background model negative')
            
            # covariance matrix
            covar = np.dot(basisFunc, weights)
            coeff = np.dot(self._bin_val, weights)
            if self._order >= 1:
                self.covar = np.linalg.inv(covar)
            else:
                self.covar = 1.0/covar
            
            # coefficients
            coeff = np.dot(coeff, self.covar)
            
            # evaluate model
            self.coeff = coeff
            model = self._evalModel(self._bin_start, self._bin_end)
         
        # evaluate model uncertainty
        modelUncert = self._evalUncertainty(self._bin_start, self._bin_end)
        # evaluate goodness-of-fit   
        self.chisq, self.dof = self._calcChisq(model)
        
        return model, modelUncert      
        
    def refit(self, order):
        """Refit the data with a different order polynomial
        
        Parameters:
        --------------
        order: int
            The order of the polynomial

        Returns:
        --------------
        np.array
            The fitted model at each input bin
        """
        order = int(order)
        assert order >= 0, 'Polynomial order must be non-negative' 

        self._order = order
        self.coeff = np.zeros(order+1)
        self.covar = np.zeros((order+1, order+1))
        self.chisq = 0.0
        self.dof = 0
        model = self.fit(self._bin_start, self._bin_end, self._bin_val)
        return model

    def interpolate(self, bin_start, bin_end):
        """Interpolation of the fitted polynomial
        
        Parameters:
        --------------
        bin_start: np.array
            The starting edge of each bin
        bin_end: np.array
            The ending edge of each bin

        Returns:
        --------------
        np.array
            The interpolated model at each bin
        np.array
            The model uncertainty at each bin
        """
        interp = self._evalModel(bin_start, bin_end)
        interpUncert = self._evalUncertainty(bin_start, bin_end)
        return interp, interpUncert
    
    def _calcChisq(self, model):
        """Calculate the chi-squared goodness-of-fit for the fitted model.
        
        Parameters:
        --------------
        model: np.array
            The fitted model

        Returns:
        --------------
        float
            The chi-squared goodness-of-fit
        float
            The degrees of freedom
        """
        variance = model/self._livetime
        # do not calculate using bins with value <= 0.0
        idx = self._bin_val > 0.0
        chisq = np.sum((self._bin_val[idx]-model[idx])**2/variance[idx])
        dof = np.sum(idx) - (self._order+1.0)
        return chisq, dof
        
    def _evalBasis(self, bin_start, bin_end):
        """Evaluates basis functions, which are the various polynomials averaged over
           the time bins.
        
        Parameters:
        --------------
        bin_start: np.array
            The starting edge of each bin
        bin_end: np.array
            The ending edge of each bin
        
        Returns:
        --------------
        np.array (nBins, order+1)
            The basis functions for each bin
        """
        numBins = len(bin_start)
        dt = bin_end-bin_start
        zerowidth = (dt == 0.0)
        bin_end[zerowidth] += 2e-6
        dt[zerowidth] = 2e-6
        basisFunc = np.zeros((self._order+1, numBins))
        
        for i in range(self._order+1):
            basisFunc[i,:] = (bin_end**(i+1.0) - bin_start**(i+1.0)) / ((i+1.0)*dt)
        return basisFunc

    def _evalModel(self, bin_start, bin_end):
        """Evaluates the fitted model over the data
        
        Parameters:
        --------------
        bin_start: np.array
            The starting edge of each bin
        bin_end: np.array
            The ending edge of each bin
        
        Returns:
        --------------
        np.array (nBins)
            The model value for each bin
        """
        numBins = len(bin_start)
        dt = bin_end-bin_start
        zerowidth = (dt == 0.0)
        bin_end[zerowidth] += 2e-6
        dt[zerowidth] = 2e-6
        model = np.zeros(numBins)
        for i in range(self._order+1):
            model += (self.coeff[i]*(bin_end**(i+1.0) - \
                                    bin_start**(i+1.0)) / ((i+1.0)*dt)).astype(float)
        return model
    
    def _evalUncertainty(self, bin_start, bin_end):
        """Evaluates the uncertainty in the model-predicted values for the data intervals
           based on the uncertainty in the model coefficients.
        
        Parameters:
        --------------
        bin_start: np.array
            The starting edge of each bin
        bin_end: np.array
            The ending edge of each bin
        
        Returns:
        --------------
        np.array (nBins)
            The model uncertainty for each bin
        """
        numBins = len(bin_start)
        uncertainty = np.zeros(numBins)
        basisFunc = self._evalBasis(bin_start, bin_end)
        basisFuncT = np.transpose(basisFunc)
        
        # formal propagation of uncertainty of fit coefficients to uncertainty of model
        for i in range(numBins):
            uncertainty[i] = np.dot(basisFunc[:,i], np.dot(self.covar, basisFunc[:,i]))
        uncertainty = np.sqrt(uncertainty)
        return uncertainty