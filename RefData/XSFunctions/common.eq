c===================================================================================
c Global parameters

c Basic array size parameters. NPHENERG is number of energy bins in
c the output spectrum. NLORGAM are the number of Lorentz gamma factor bins.
c MAXITER is the maximum number of iterations when looping.

      INTEGER NPHENERG, NLORGAM
      PARAMETER (NPHENERG=221, NLORGAM=81)

c Geometric factor
c ACHTUNG: there is a geometric factor relating le/s to
c q,s the >dimensionless<  particle injection rates:
      DOUBLE PRECISION GEOZ
c        PARAMETER(GEOZ=16.*atan(1.)/3.)
c              (this ^ is for comparison with big Z,FGBC geoz = 1.)
      PARAMETER(GEOZ=4.1887902047863905)

c Initial Coulomb couling rate
      DOUBLE PRECISION CTURNON
      PARAMETER (CTURNON=1.d0)
 
      DOUBLE PRECISION WMIN, WMAX
c       Non-integer exp initialization not allowed on some compilers
c       PARAMETER (WMIN=10**(-7.025), WMAX=10.**4.025)
      PARAMETER (WMIN=9.44060901e-8, WMAX=10592.537)

      INTEGER IOFF
      DOUBLE PRECISION RIOFF1
      PARAMETER (IOFF=40, RIOFF1=141.5)

c===================================================================================

c For content from ltcz.fits (used in getethd)
c note clav also read from azecomp.fits so overrides

      DOUBLE PRECISION rldisp(61,41,2)
      INTEGER ig(61,41,2)

c For content from comp.fits (used in getethd and mkscat)

      DOUBLE PRECISION ecomp(61,161,3), gcomp(61,161,3)
      INTEGER iecomp(61,161,3), igcomp(61,161,3)

c For content from azntc.fits (used in bbesc and doelec)
c note compr also read from azcomp.fits so overrides

      DOUBLE PRECISION av(NLORGAM,NPHENERG)

c For content from clt2.fits (used in mkscat)

      DOUBLE PRECISION cltr(141,241), clta(141,241), cltd(141,241)

c For content from azcomp.fits (used in dopcomp)

      DOUBLE PRECISION ctab(520483)
      DOUBLE PRECISION compr(NLORGAM,NPHENERG)
      INTEGER ictab(NLORGAM,NPHENERG,2)

c For content from azecomp.fits (used in getethd and doelec)

      DOUBLE PRECISION clav(61,41), ectab(218997)
      INTEGER iectab(NLORGAM,NPHENERG,2)
      
c For content from azpp.fits (used in dopcomp and doelec)

      DOUBLE PRECISION ppr(NLORGAM,NPHENERG), ppf(NLORGAM,NPHENERG)
      INTEGER ipp(NLORGAM,NPHENERG)

c For content from thrm.fits (used in dotherm and bbesc)

      DOUBLE PRECISION thann(3,2), cr(162,2)

c For content from htc.fits (used in dotherm, thcdotherm and bbesc)

      DOUBLE PRECISION htcr(161,4), htcav(161,4), htcd(161,4)
      DOUBLE PRECISION htcavl(NLORGAM,4), htcdl(NLORGAM,4)        

c For content from brem.fits (used in dobrem)

      DOUBLE PRECISION brem(61), brema(7710)
      INTEGER ibrem(61)


c Pair injection function and electron and/or positron injection function 
c (if pair injection is not being used). The latter does nothing at the moment.
c Set in setup and used in solve.

      DOUBLE PRECISION pairinjfn(NLORGAM), epinjfn(NLORGAM)

c (Soft) Photon injection function. Set in setup, used in solve and bbesc.

      DOUBLE PRECISION phinjfn(NPHENERG)

c Energies corresponding to particle Lorentz factor bins.
c Set in setup. 

      DOUBLE PRECISION een(NLORGAM)

c Energies and bin widths corresponding to photon energy bins.
c Set in setup. Used in pairiter, solve, dobrem, thscat, bbesc, doesc, 
c dotherm, mkscat, getethd, thcsolve, thcdotherm.
c Units are m_e c^2 (ie 511 keV).

      DOUBLE PRECISION wgam(NPHENERG), upwgam(NPHENERG)
      DOUBLE PRECISION gwidth(NPHENERG)

c Spectrum of photons which escape (and hence are observed). Calculated 
c and used in pairiter, solve, thcsolve. (Calculated in doesc which is 
c apparently not called from anywhere)

      DOUBLE PRECISION esc(NPHENERG)

c Escape array. Calculated in dotherm. Used in dopcomp, getth, solve, thcdopcomp, thcsolve.

      DOUBLE PRECISION tesc(NPHENERG)

c Temporary escape array. Used in getth. Set in solve, thcsolve.

      DOUBLE PRECISION otesc(NPHENERG)

c Work photon spectrum array. Set in setup. Modified in getth, dopcomp. Used in pairiter, 
c pprod, solve, getethd, doelec, thcdocomp, thcsolve.

      DOUBLE PRECISION gdist(NPHENERG)

c Another work photon spectrum array. Set and used in dopcomp, dobrem, thscat, dotherm,
c thcdopcomp.
 
      DOUBLE PRECISION gder(NPHENERG)

c Thermal electron and positron densities. These replace edist(1) and posdist(1) 
c which were used for this purpose. These are set and used in getth, setup. Set in solve. Used 
c in dobrem, docoul, bbesc, dotherm, mkscat, getethd, thcdotherm.

      DOUBLE PRECISION elecdens, posdens

c Non-thermal electron and positron distributions. Set in setup. Used in dopcomp, dobrem, 
c docoul, bbesc. Used and updated in doelec. Note that edist(1) and posdist(1) 
c are not actually used.

      DOUBLE PRECISION edist(NLORGAM), posdist(NLORGAM)

c Work array. Total cooling rate of pairs ? Set and used in doelec. Updated in dobrem. 
c Updated in docoul. 

      DOUBLE PRECISION gdot(NLORGAM)

c Variables to hold parameter values. 
c       theta  = electron temperature in m_e c^2 (used in setup, getth, dobrem)
c       bbtemp = BB temperature in m_e c^2    (used in setup)
c       rle    = non-thermal particle compactness (used in setup and solve)
c       rls    = soft photon compactness (used in setup and solve)
c       rlth   = thermal particle compactness (used in setup, getth and solve)
c       radius = Source radius in cm    (used in docoul)
c       taup   = Thomson optical depth   (used in solve and getth)

      DOUBLE PRECISION theta, bbtemp, rle, rls, rlth, taup
      DOUBLE PRECISION radius

c Change in thermal energy. Used and set in getth. Updated in dobrem, docoul, getethd.
c Set in doelec.

      DOUBLE PRECISION ethder

c Flag whether or not to calculate thermal particle distribution. Set in solve, thcsolve. 
c Used in dopcomp, getth.

      LOGICAL dothermpairs


      DOUBLE PRECISION diff, dt
      DOUBLE PRECISION eth
      DOUBLE PRECISION prninj, erninj

      INTEGER iter, idump, niter


      LOGICAL qpairinj

      COMMON /eqpbigdata/ av, ppr, ppf, ipp
      COMMON /eqpbrempar/ brem, brema, ibrem
      COMMON /eqpccoul/ radius
      COMMON /eqpchatter/ diff,dt,rle,rls,rlth,iter
      COMMON /eqpchat2/ bbtemp,prninj,erninj,taup,dothermpairs
      COMMON /eqpdist/ edist, gdist, gder, esc, theta, eth, 
     &                 ethder, posdist, elecdens, posdens
      COMMON /eqpeminus/ gdot, idump
      COMMON /eqphthrm/ htcr, htcav, htcd, htcavl, htcdl
      COMMON /eqpiterpar/ niter
      COMMON /eqplotcmp/ cltr, clta, cltd
      COMMON /eqpparams/ pairinjfn, epinjfn, phinjfn,
     &                   een, wgam
      COMMON /eqpscatter/ upwgam, gwidth
      COMMON /eqpthermal/ thann, cr
      COMMON /eqpthstuff/ tesc, otesc
      COMMON /eqplccomp/ clav, rldisp, ig
      COMMON /eqpbigdo/ ecomp, gcomp, iecomp, igcomp
      COMMON /eqpcompton/ compr, ictab, iectab, ctab, ectab
      COMMON /eqpxflags/ qpairinj
