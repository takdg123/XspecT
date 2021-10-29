from iminuit import Minuit
import numpy as np
import matplotlib.pyplot as plt

class Linear_fit:
    
    def __init__(self, x, y, x_err = [], y_err = [], logx = False, logy = False, pinit = [-1, 1, 1, 100, 0.8], verbose=False, ext = True):
        self._ext = ext
        self._x = np.asarray(x)
        if len(x) != len(x_err):
            self._x_err = np.zeros(len(x))
            self._ext = False
        else:
            self._x_err = np.asarray(x_err)
        self._y = np.asarray(y)
        if len(y) != len(y_err):
            self._y_err = np.zeros(len(y))
            self._ext = False
        else:
            self._y_err = np.asarray(y_err)
        self._pinit = pinit
        
        if logx: self._x, self._x_err = self.__covLog__(self._x, self._x_err)
        if logy: self._y, self._y_err = self.__covLog__(self._y, self._y_err)
        self._logx=logx
        self._logy=logy
        self._verbose = verbose
        
    @property
    def p(self):
        return self._p

    @property
    def perr(self):
        return self._perr

    @property
    def logx(self):
        return self._logx

    @property
    def logy(self):
        return self._logy
        

    def linear(self, x):
        if hasattr(self, "p"):
            if self.logx and self.logy:
                return 10**(self.__linear__(np.log10(x), *self.p))
            elif self.logx:
                return self.__linear__(np.log10(x), *self.p)
            elif self.logy:
                return 10**self.__linear__(x, *self.p)
            else:
                return self.__linear__(x, *self.p)
        else:
            return None

    def __covLog__(self, t, dt):
        dt = np.log10(t+abs(dt))-np.log10(t)
        t = np.log10(t)
        return t, dt
    
    def __invLog__(self, t, dt):
        t = 10**t
        dt = t*(10**dt-1)
        return t, dt
        
    def __linear_ext_likelihood__(self, m, k, sig_ext):
        likelihood = 1/2.*sum(np.log(m**2.*self._fx_err**2.+self._fy_err**2.+sig_ext**2.))+1/2.*sum((self._fy-m*self._fx-k)**2./(sig_ext**2.+self._fy_err**2.+m**2.*self._fx_err**2.))
        return likelihood

    def __linear_likelihood__(self, m, k):
        likelihood = sum((self._fy-m*self._fx-k)**2./(self._fy_err**2.+m**2.*self._fx_err**2.))
        return likelihood

    def __bknlinear_likelihood__(self, m1, m2, k, bk):
        likelihood = np.zeros(len(self._fx))
        for x, y, xerr, yerr, i in zip(self._fx, self._fy, self._fx_err, self._fy_err, range(len(self._fx))):
            if x <= np.log10(bk):
                likelihood[i] = (y-m1*x-k)**2./(yerr**2.+m1**2.*xerr**2.)
            else:
                likelihood[i] = (y-(m2*x+k+(m1-m2)*np.log10(bk)))**2./(yerr**2.+m2**2.*xerr**2.)
        return sum(likelihood)

    def __sbknlinear_likelihood__(self, m1, m2, k, bk, s):
        diff = -k*(((self._fx/bk)**(s*m1)+(self._fx/bk)**(s*m2))**(-(s+1.)/s)*(m1*(self._fx/bk)**(s*m1)+m2*(self._fx/bk)**(s*m2)))/self._fx
        likelihood = sum((self._fy-self.__sbknlinear__(self._fx, m1, m2, k, bk, s))**2./(self._fy_err**2.+diff**2.*self._fx_err**2.))
        return likelihood

    def __linear__(self, x, m, k):
        return m*x+k
    
    def __bknlinear__(self, x, m1, m2, k, bk):
        temp = np.zeros(len(x))
        for x_c, i in zip(x, range(len(x))):
            if x_c <= np.log10(bk):
                temp[i] = m1*x_c+k
            else:
                temp[i] = m2*x_c+k+(m1-m2)*np.log10(bk)
        return temp

    def __sbknlinear__(self, x, m1, m2, k, bk, s):
        return k*((x/bk)**(s*m1)+(x/bk)**(s*m2))**(-1./s)

    def runFit(self, start = -1, end = -1, returned=False, bkn=False, sbkn=False):
        if start == -1 and end == -1:
            self._fx, self._fy = self._x, self._y
            self._fx_err, self._fy_err = self._x_err, self._y_err
        elif start == -1:
            self._fx, self._fy = self._x[:-end], self._y[:-end]
            self._fx_err, self._fy_err = self._x_err[:-end], self._y_err[:-end]
        elif end == -1:
            self._fx, self._fy = self._x[start:], self._y[start:]
            self._fx_err, self._fy_err = self._x_err[start:], self._y_err[start:]
        else:
            self._fx, self._fy = self._x[start:-end], self._y[start:-end]
            self._fx_err, self._fy_err = self._x_err[start:-end], self._y_err[start:-end]

        if self._ext:
            model=Minuit(self.__linear_ext_likelihood__, 
                  m = self._pinit[0], error_m=0.01, limit_m = [-10, 10],
                  k = self._pinit[1], error_k=0.01, limit_k = [-1e3,1e3],
                  sig_ext = self._pinit[4], error_sig_ext=1e-2, limit_sig_ext = [1e-5, 1e2],
                  print_level=self._verbose, errordef=1)
        elif bkn:
            model=Minuit(self.__bknlinear_likelihood__, 
                  m1 = self._pinit[0], error_m1=0.01, limit_m1 = [-10, 10],
                  m2 = self._pinit[0]+0.1, error_m2=0.01, limit_m2 = [-10, 10],
                  k = self._pinit[1], error_k=0.01, limit_k = [-1e3,1e3],
                  bk = self._pinit[2], error_bk=1, limit_bk = [200, 1e9],
                  #bk = self._pinit[2], error_bk=1, limit_bk = [5, 100],
                  print_level=self._verbose, errordef=1)
        elif sbkn:
            model=Minuit(self.__sbknlinear_likelihood__, 
                  m1 = self._pinit[0], error_m1=0.01, limit_m1 = [-10, 10],
                  m2 = self._pinit[0]+0.1, error_m2=0.01, limit_m2 = [-10, 10],
                  k = self._pinit[1], error_k=0.01, limit_k = [1e-15,1],
                  bk = self._pinit[2], error_bk=1, limit_bk = [1, 1e9],
                  s = self._pinit[3], error_s=0.01, limit_s = [1e-2, 100],
                  print_level=self._verbose, errordef=1)
        else:
            model=Minuit(self.__linear_likelihood__, 
                  m = self._pinit[0],
                  k = self._pinit[1]
                  )

        model.errordef=1
        fit_result = model.migrad()
        model.hesse()


        chisq = fit_result.fval
        dof = len(self._fx) - len(fit_result.parameters)

        self._p = (fit_result.values['m'], fit_result.values['k'])
        self._perr = (fit_result.errors['m'], fit_result.errors['k'])
        self._stat = [chisq, dof]
        if returned: 
            return model

    def plotResult(self, bkn=False, sbkn=False):
        plt.errorbar(self._x, self._y, xerr=self._x_err, yerr=self._y_err, ls='')
        plt.errorbar(self._fx, self._fy, xerr=self._fx_err, yerr=self._fy_err, ls='')
        if bkn:
            x_m = np.linspace(min(self._x), max(self._x), 100)
            plt.plot(x_m, self.__bknlinear__(x_m, self.p[0], self.p[1], self.p[2], self.p[3]))        
            plt.title(r"$\alpha$1 = {:.2f}, $\alpha$2 = {:.2f}, t$_{{bk}}$ = {:.1e}".format(-self.p[0], -self.p[1],self.p[3]))
        elif sbkn:
            x_m = np.logspace(np.log10(min(self._x)), np.log10(max(self._x)), 1000)
            plt.plot(x_m, self.__sbknlinear__(x_m, self.p[0], self.p[1], self.p[2], self.p[3],self.p[4]))        
            plt.title(r"$\alpha$1 = {:.2f}, $\alpha$2 = {:.2f}, t$_{{bk}}$ = {:.1e}, s = {:.2f}".format(self.p[0], self.p[1],self.p[3],self.p[4]))
        else:
            x_m = np.linspace(min(self._x), max(self._x), 100)
            plt.plot(x_m, self.__linear__(x_m, self.p[0], self.p[1]))
            plt.title(r"$\alpha$ = {:.2f}".format(-self.p[0]))
        if sbkn:
            plt.xlabel("time")
            plt.ylabel("flux")
            plt.xscale("log")
            plt.yscale("log")
        else:
            plt.xlabel("log(time)")
            plt.ylabel("log(flux)")
        plt.show(block=False)
