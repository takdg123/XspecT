
      subroutine xsp1tr(ear, ne, param, ifl, photar, photer)

c Power-law with high-energy cutoff absorbed by cold spherical
c matter (isotropic source at centre), taking into account Compton
c scattering as described in Yaqoob 1997 (ApJ, APril 10; or
c see ftp://lheawww.gsfc.nasa.gov/pub/yaqoob/papers/xtrans/PAPER.ps)
c Refer to that paper for the physical limitations of this model.
c
c** parameters **
c par 1       - Radial column density (1e22 cm^-2)
c par 2            - Maximum number of scatterings to consider
c par 3            - Fe abundance
c par 4            - Fe K egde energy (keV)
c par 5       - Photon index
c par 6            - Cutoff energy below which no exponential rollover applied
c par 7            - e-folding energy (keV)
c par 8       - critical albedo for elastic scattering (remmd 0.1)
c par 9       - IFAST=0 for full integration, 1 for simple e-shift
c par 10       - Redshift

      integer ne, ie, ifl, nsmax, namax, kscat
      parameter (nsmax=100, namax=1000)
      real ear(0:ne), param(*), photar(ne), photer(ne)
      real xnh, tausl, taus, pgam, ecut, efold, ans, w1, w2
      real ecen, eobs, wobs, afe, ek, zfac, prb1
      real fexp, abscrs
      real spec, absfac, ewid, taut, sfac, wave
      real wobs1, wobs2, probs(nsmax), critp, spadd, pfchn
      real alb0, acrit, pelas, teff
      integer ifast, nmax, iscat, iwscat
c        integer nact(namax)
      common /col/xnh, taus, tausl
      common /plpars/pgam, ecut, efold
      common /eobs/ecen, wobs
      common /fepars/afe, ek
      common /scat/iscat, nmax, iwscat
      common /wave/wave
      common /wave1/wobs1
      common /wave2/wobs2
      common /probs/probs
c      common /nact/nact

c The above two common blocks are used for passing on the obeserved
c wavelength to intgegrands for Pn and the weighted mean albedo.
c Use wobs1 for p1itgpl and wobs2 for walbit
      external p1itgpl, walbit
c
      data critp/0.01/
c

c suppress a warning message from the compiler
      ie = ifl

c this model does not calculate errors
      DO ie = 1, ne
         photer(ie) = 0.0
      ENDDO



      xnh=10.*param(1)
      taus=0.000798*xnh
      tausl=alog(taus)
      nmax=max(1,int(param(2)))
      afe=param(3)
      ek=param(4)
      pgam=param(5)
      ecut=param(6)
      efold=param(7)
      acrit=param(8)
      ifast=int(param(9))
      zfac=1.+param(10)
c      ecrit=0.01
c      ecrit=5.0
c at the critical albedo ACRIT the switch is made from inelastic
c to elastic scattering. ACRIT=0.185 corresponds to 5 keV for cosmic
c abundances
c if ifast = 1 don't do the integration - simply apply mean compton
c shift
      do 100, ie=1,ne
         ecen = zfac*(ear(ie)+ear(ie-1))/2.
         ewid = ear(ie)-ear(ie-1)
         if (ecen.le.ecut) then
            sfac = 1.0
         else
            sfac = fexp(-(ecen-ecut)/efold)
         endif
         spec = sfac/(ecen**pgam)
         taut = abscrs(ecen,xnh,afe,ek) + taus 
         absfac = fexp(-taut)
c following is the zeroth order single scattering albedo
         alb0=taus/taut
c      write(*,*) 'energy alb0',ecen,alb0
c wavelength at the bin centre
         w2=511./ecen
c minimum wavelength of photons which can make a contribution at ecen.
         eobs=ecen
         wobs1=w2
         photar(ie)=spec*absfac
c       write(*,*) 'Doing energy bin ',ie
         iscat=0
c come back here if doing another scattering
 29      continue
c ignore energy shifts at low energies
         if (alb0.le.acrit) then 
c      write(*,*) 'acrit = ',acrit
            iscat=nmax
c note that PRB1 sets up escape probabilities up to ISCAT scatterings
c in the array probs held in common block /probs/
            pelas = prb1(ecen)
c        write(*,*) ie,pelas
c        write(*,*) (probs(kk),kk=1,nmax)
            kscat=0
 35         continue
            kscat=kscat+1
            spadd = probs(kscat)*spec
c check if the fractional change in the escape spectrum is going to
c be less than CRITP
            if (photar(ie).gt.1e-34) then
               pfchn = spadd/photar(ie)
            else
               pfchn = 1.e34
            endif
c                  write(*,*) 'pfchn = ',pfchn
            if (pfchn.gt.critp) then
               photar(ie)=photar(ie)+spadd
            else
               iscat=kscat
c do next energy bin
               goto 60
            endif
            if (kscat.lt.nmax) goto  35
         else
c don't ignore energy shifts in this part of the IF statement
            iscat=iscat + 1
            if (ifast.lt.1) then  
c do the full integration if IFAST lt 1
c            write(*,*) 'Doing ISCAT = ',is
               w1=w2-2.0*real(iscat)
c        write(*,*) ecen,w1,w2,photar(ie), photer(ie)
c following is the integration. w1 and w2 are the wavelength limits
               call qg32(w1,w2,p1itgpl,ans)
c note the factor 261121 is m_e c^2 in keV squared
               spadd = 261121.*ans/ecen/ecen
c check if the fractional change in the escape spectrum is going to
c be less than CRITP
               if (photar(ie).gt.1.e-34) then
                  pfchn = spadd/photar(ie)
               else
                  pfchn = 1.e34
               endif
c                  write(*,*) 'pfchn = ',pfchn
               if (pfchn.gt.critp) then
                  photar(ie)=photar(ie) + spadd
               else
c do next energy bin
                  goto 60
               endif
c        write(*,*) ecen,w1,w2,photar(ie), photer(ie)
c       write(*,*) 'doing integration'
            else
c this part if the IF statement is for IFAST=1, using the approximation
c exp(-tau_eff) for the attenutation and assumes elastic scattering.
c see Yaqoob 1997 for details.
               teff=taut*sqrt(1.-alb0)
               photar(ie)=spec*exp(-teff)
c this alternative part of the IF statement doesn't do the full integration 
c (IFAST gt 1) but simply shifts photons by the mean energy shift
c for a single compton scattering
c        if (ecen.lt.511.) then
c          edash = ecen/(1.-(ecen/511.))
c         spadd = prb1(edash)*spec
c check if the fractional change in the escape spectrum is going to
c be less than CRITP
c        if (photar(ie).gt.1.e-34) then
c             pfchn = spadd/photar(ie)
c        else
c            pfchn = 1.e34
c        endif
c                  write(*,*) 'pfchn = ',pfchn
c        if (pfchn.gt.critp) then
c               photar(ie)=photar(ie) + spadd
c        else
c do next energy bin
c               goto 60
c        endif
c        endif
            endif
         endif 
c      write(*,*) ecen,ans,photar(ie), photer(ie)
c      write(11,*) ecen,photar(ie), photer(ie)
         if (iscat.lt.nmax) goto 29
 60      continue
       photar(ie)=photar(ie)*ewid/zfac
c actual number of scatterings used for this energy bin
c       nact(ie)=iscat
c      write(11,*) 'NO NO '
 100  continue
      end
c**********************************************************************
      function p1itgpl(wave0)
c T. Yaqoob March 1995
c This is the intgrand for the integral to compute the spectrum
c due to the first compton scattering in spherical cold matter (primary
c source at centre) when the input spectrum is a power law with 
c exponential cutoff at high energies. Be sure to set all the variables
c in the common blocks before calling. The argument to this function
c is the dimensionless wavelength of the primary photons (i.e. 511/e(kev)).
c EOBS is the observed energy and wobs the observed dimensionless
c wavelength. Must multiply the result by (511^2) to get the right norm.
      real xnh, taus, tausl, ecut, efold, fexp
      real wave0, wobs, eobs, prb1, efac, e0
      real p1itgpl, pgam, fdist, wobs1, pll, p1
      real dw, dw1, dw2, dw4, dwn, arg, risc, pi
      integer iscat, nmax, iwscat
      common /col/xnh,taus,tausl
      common /plpars/pgam, ecut, efold
      common /eobs/eobs, wobs
      common /wave1/wobs1
      common /scat/iscat,nmax,iwscat
      data pi/3.14159265/
c      
      e0=511./wave0
c      write(*,*) 'iscat e0 wave0 wobs1 '
c      write(*,*) iscat,e0,wave0,wobs1
c      write(*,*) 'pgam ecut efold e0^pgam'
c      write(*,*) pgam,ecut,efold,e0**pgam
c compute the distribution in wavelength shift.
      dw=wobs1-wave0
      dw1=dw -1.0
      fdist=0.0
      if (iscat.eq.1) then
         fdist = 0.375*(1.+dw1*dw1)
      elseif (iscat.eq.2) then
         if (dw.ge.2.0.and.dw.le.4.) dw=4.-dw
         dw2=dw*dw
         dw4=dw2*dw2
         fdist = 0.140625*(4.*dw - 4.*dw2 +2.*dw2*dw
     &        -(dw4/3.) + (dw4*dw/30.))
      elseif (iscat.ge.3) then
         risc=real(iscat)
         dwn=dw - risc
         arg=-1.25*dwn*dwn/risc
         fdist=sqrt(1.25/pi/risc)*fexp(arg)
      endif
      if (e0.gt.ecut) then 
         efac=exp(-(e0-ecut)/efold)
      else
         efac=1.0      
      endif
c      write(*,*) 'wave0 '
c      write(*,*) wave0
      pll = alog10(fdist)-(pgam*alog10(e0))-(2.*alog10(wave0))+
     &     alog10(efac)
      if (pll.gt.-34.) then
         p1itgpl = 10.**pll
      else
         p1itgpl=0.
      endif
c      write(*,*) 'wave 0 before prb1',wave0
c      write(*,*) wave0
      p1=prb1(e0)
c      write(*,*) 'wave0 after prb1 ',wave0
c      p11=prb1(e0)
c      write(*,*) 'wave0 after second prb1 ',wave0
      p1itgpl=p1itgpl*p1
c      write(*,*) wave0
c      write(*,*) 'p1itgpl= ',p1itgpl
c      write(*,*) 'P1ITGPL: e0 prb1(e0): ',e0,p1
c      write(*,*) 'e0 wave0 wobs1 fdist dw '
c      write(*,*) e0,wave0,wobs1,fdist,dw
c      return
      end
c**************************************************************************
      function prb1(en)
c T. Yaqoob March 1995 ->**
c compute probability that a photon with energy en (keV)
c created at the centre of a sphere with radial column density NH
c (in common block /col/) escapes with exactly ISCAT scatterings.
c Ensure that ISCAT is set before entry.
c MSC is the number of scatterings for which explict expressions 
c for the escape functions have been derived from Monte Carlo 
c calculations. For higher order scatterings all escape functions are
c set equal to the MSC th escape function. See Yaqoob 1995 for details.
      integer msc, nsmax
      parameter (msc=9, nsmax=100)
      real prb1, en, xnh, taus, afe, ek, abscrs
c
c note taus is the radial thomson depth = 0.798 NH_24 (should be
c computed before entry). The iron abundance and Fe K-edge energy 
c should also be set on entry.
c
      real p0fac, u1,u2,u3,v1,v2,v3,vv2,w1,w2,w3
      real z1, z2, sum, pf1, eta(msc+1), pkap(msc+1)
      real fexp, swgtalb, alb(nsmax)
      real albl(nsmax), albll(nsmax), gn(nsmax)
      real tausl, wmatx(nsmax,nsmax), probs(nsmax)
      real tdsh, tdshl, efc, ffc, t0, t1, apwr, sqt, f2, f3
      real t0h, tsfac
      real an(msc), bn(msc), gg(msc), jp(nsmax)
      integer iscat,nmax, i, j, k, l, m, n, n0, iwscat
      common /fepars/afe, ek
      common /col/xnh,taus, tausl      
      common /scat/iscat,nmax,iwscat
      common /probs/probs
c stuff below is from fits to monte carlo runs
      data u1,u2,u3/0.0656774, 1.4967, 0.2635/
      data v1, v2, v3, vv2/2.0956, 1.7885, -2.3559, 0.559128/
      data w1, w2, w3/0.319180, 1.5616, -0.3794/
      data z1, z2/-0.14717, -0.8656/
      data eta/0.,1.2017,1.2850,1.4946,1.8453,1.90,1.90,1.90,1.90,1.90/
      data pkap/0.,7.1695e-2,4.0526e-2,2.0718e-2,9.1339e-3,
     & 5.9456e-3,4.6936e-3,3.7428e-3,3.5835e-3,2.3279e-3/
      data an/0.,-0.50,0.70,1.15,1.40,1.65,1.72,1.80,1.90/
      data bn/0.,-2.60,-3.50,-3.75,-4.00,-4.25,-4.32,-4.45,-4.55/
      data gg/0.,5.7,5.9,5.9,5.9,5.9,5.9,5.9,5.9/
c followig is equal to square-root of (3/4)
      data sqt/0.866025403784/
c      write(10,*) 'PRB1: e0 = ',en
c      write(12,*) 'PRB1: e0 = ',en
      prb1=0.0
      tsfac=fexp(-taus)
c set up ISCAT+1 albedos
c single scattering albedo for zero scatterings
      alb(1) = 1./(1.+(abscrs(en,1.,afe,ek)/0.798e-3))
c      write(12,*) 'i+1 en alb(i+1) '
      do 10, i=1,iscat 
         alb(i+1) = swgtalb(en,i)
c      write(*,*) alb(i+1),swgtalb(en,i)
         albl(i+1) = -alog(alb(i+1))
         if (albl(i+1).gt.0.) albll(i+1)=alog(albl(i+1))
         wmatx(1,i)=alb(i)
c       write(12,*) i+1,en,alb(i+1)
c this bit computes the jumping function
         n0=max(int((511./ek)-(511/en)),0)
         if (i.gt.n0.and.en.gt.ek) then
            apwr=real(i-n0+1)/2.
            jp(i)=(tsfac+(1.-tsfac)*(alb(i)**apwr))
         else
            jp(i)=1.0
         endif
 10   continue
c       write(10,*) (wmatx(1,kk),kk=1,iscat)
c set up remaining albedo matrix
      if (iscat.ge.2) then
         do 30,l=2,iscat
            do 20,m=l,iscat
               wmatx(l,m)=1.0
               do 15,k=1,l
                  wmatx(l,m)=wmatx(l,m)*alb(m-k+1)
 15            continue
 20         continue
c          write(10,*) (wmatx(l,kk),kk=1,iscat)
 30      continue 
      endif
c      write(12,*) 'ALB: ',(alb(kk),kk=1,iscat+1)
c now compute the first scattering function
      t0=taus/alb(1)
      t1=taus/alb(2)
      t0h=0.5*t0
      f2=4.*fexp(-t1)+3.*fexp(-t0)+2.*fexp(-t0h-1.5*t1)+
     &     2.*fexp(-t0h-0.5*t1)+4.*fexp(-t0h-sqt*t1)+fexp(-t0-2.*t1)
      f3=1.+9.8600e-2*t0-2.8717e-4*t0*t0+7.0954e-7*t0*t0*t0
      if (jp(1).le.0.) write(*,*) '**warning** jp zero'
      probs(1)=taus*f2/f3/16./jp(1)       
      p0fac=(1.-fexp(-taus/alb(1)))
      gn(1)=probs(1)/alb(1)/p0fac
c       write(12,*) 'g1 = ',gn(1)
c       write(12,*) 'probs(1) alb(1) p0fac ',probs(1), alb(1),p0fac
      if (iscat.eq.1) then
         prb1=probs(1)
      else
         do 250,n=2,iscat
c compute the Gn
            if (n.le.msc) then
               tdsh = taus*alb(n)/alb(2)/alb(n-1)
               tdshl=alog(tdsh)
               efc = fexp(-tdsh/gg(n))
               ffc = (0.25+0.75*efc)*tdshl + (an(n)+bn(n)*efc)
               if (jp(n).le.0.) write(*,*) '**warning** jp zero'
               gn(n) = gn(1)*(1.+fexp(ffc))/jp(n)
            else
               gn(n) = gn(msc)
            endif
            sum=0.0
            do 100,j=1,n-1
               sum=sum+probs(j)*wmatx(n-j,n)
 100        continue
c     write(*,*) en,gn(n),n,wmatx(n,n),sum
            pf1=p0fac*wmatx(n,n)
            if (pf1.gt.sum) then
               probs(n)=gn(n)*(pf1 - sum)
            else
               probs(n)=0.0
            endif
            prb1=probs(n)
 250     continue
      endif
c      write(11,*) 'PRB1: e0 ISCAT  = ',en,iscat
c      do 300, k=1,iscat 
c      write(11,*) 'n prb1(n) g(n)',k,probs(k),gn(k)
c      300  continue
      end
c*********************************************************************
      function wgtalb(en,isc)
c For photons of initial energy EN, compute the weighted mean albedo
c for ISC Compton scatterings on cold electrons. 
c The maximum wavelength change is 2. 
      real wobs2, wave2, en, wgtalb, ans, walbit
      integer iscat, nmax, iwscat
      common /scat/iscat,nmax,iwscat
      integer isc
      external walbit
      common/wave2/wobs2
      iwscat=isc
      wobs2 = 511./en
      wave2 = wobs2 + 2.*real(iwscat)
      call qg32(wobs2,wave2,walbit,ans)
      wgtalb = ans
      return
      end      
c*********************************************************************
      function swgtalb(en,isc)
c For photons of initial energy EN, compute the weighted mean albedo
c for ISC Compton scatterings on cold electrons. 
c The maximum wavelength change is 2. 
c This routine uses 3-point simpson rule for integration
      real wobs2, wave2, en, swgtalb, walbit, wwid, wnew
      integer iscat, nmax, iwscat, nint, j
      common /scat/iscat,nmax,iwscat
      integer isc
      external walbit
      common/wave2/wobs2
      iwscat=isc
      wobs2 = 511./en
      if (isc.eq.1) then 
         nint=16
      else      
         nint=4      
      endif
      wwid=2./real(nint)
      wave2 = wobs2 + 2.*real(iwscat)
      swgtalb=walbit(wobs2)
      wnew=wobs2
      do 20,j=1,nint*iwscat-1
         wnew=wnew+wwid
         swgtalb=swgtalb+2.*walbit(wnew)
 20   continue
      swgtalb=swgtalb+walbit(wave2)
      swgtalb=0.5*wwid*swgtalb
      return
      end      
c*********************************************************************
      function walbit(y0)
c integrand to compute the weighted mean albedo (by the function wgtalb)
c for the ISCAT'th Compton scattering
c on cold electrons. Y0 is the intial wavelength.
      real walbit
      real wobs2, y0, e0, alb, abscrs, afe, ek, dw, dw1, dw2, dw4
      real risc, dwn, arg, pi, fdist, fexp
      integer iscat,nmax, iwscat
      common /scat/iscat,nmax, iwscat
      common /wave2/wobs2
      common /fepars/afe, ek
      data pi/3.14159265/
      dw=abs(y0-wobs2)
      dw1=dw -1.0
      e0=511./y0
      alb = 1./(1.+(abscrs(e0,1.,afe,ek)/0.798e-3))
c compute the wavelength distribution
      fdist = 0.0
      if (iwscat.eq.1) then
         fdist = 0.375*(1.+dw1*dw1)
      elseif (iwscat.eq.2) then
c       write(13,*) 'iwscat = ',iwscat
c       write(13,*) 'DW before ',dw
         if (dw.ge.2.0.and.dw.le.4.) dw=4.-dw
c       write(13,*) 'DW after ',dw
         dw2=dw*dw
         dw4=dw2*dw2
c       write(13,*) 'DW2 DW4 ',dw2,dw4
         fdist=0.140625*(4.*dw-4.*dw2+2.*dw2*dw-(dw4/3.)+(dw4*dw/30.))
c       write(13,*) 'FDIST ',fdist
      elseif (iwscat.ge.3) then
         risc=real(iwscat)
         dwn=dw - risc
         arg=-1.25*dwn*dwn/risc
         fdist=sqrt(1.25/pi/risc)*fexp(arg)
      endif
      walbit = alb*fdist
      return
      end
C***********************************************************************
      function fndist(y0)
c integrand to compute the non-relativistic distributions of photons
c Compton scattered ISCAT times 
c on cold electrons. Y0 is the intial wavelength.
      real fndist, fexp
      real wave, y0, e0, afe, ek, dw, dw1
      real risc, adw, adw2, adw4, dwn, arg, pi, y01, y02
      integer iscat,nmax,iwscat
      common /scat/iscat,nmax,iwscat
      common /wave/wave
      common /fepars/afe, ek
      data pi/3.14159265/
      dw=wave-y0
      dw1=dw -1.0
      e0=511./y0
      y01=y0
      fndist=0.0
      if (iscat.eq.1) then 
         fndist = 0.375*(1.+dw1*dw1)
      elseif (iscat.eq.2) then 
         adw=abs(dw)
         if (adw.ge.2.0.and.adw.le.4.) adw=4.-adw
         adw2=adw*adw
         adw4=adw2*adw2
         fndist = 0.140625*(4.*adw - 4.*adw2 +2.*adw2*adw 
     &        -(adw4/3.) + (adw4*adw/30.))
      elseif (iscat.ge.3) then
         risc=real(iscat)
         dwn=dw - risc
         arg=-1.25*dwn*dwn/risc
         fndist=sqrt(1.25/pi/risc)*fexp(arg)
      endif
      y02=y0
c      if (y02.ne.y01) write(12,*) 'FNDIST:y0 y01 y02 = ',y0,y01,y02
      return
      end
C*****************************************************************************
  
      SUBROUTINE QG32 (XL, XU, FCT, Y)   
                                       
      EXTERNAL FCT
*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
*                                                                    c 
*    32-point Gaussian quadrature.                                   c 
*    xl  : the lower limit of integration                            c 
*    xu  : the upper limit                                           c 
*    fct : the (external) function                                   c 
*    y   : the value of the integral                                 c 
*                                                                    c 
*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
*                                                                               
      real a, b, c, y, xu, xl, fct
      A = .5 * (XU + XL)                                     
      B = XU - XL                                                               
      C = .498631930924740780       * B  
      Y = .35093050047350483E-2     * (FCT (A + C) + FCT (A - C))  
      C = B * .49280575577263417                                  
      Y = Y + .8137197365452835E-2  * (FCT (A + C) + FCT (A - C))   
      C = B * .48238112779375322                                  
      Y = Y + .1269603265463103E-1  * (FCT (A + C) + FCT (A - C))
      C = B * .46745303796886984                                   
      Y = Y + .17136931456510717E-1 * (FCT (A + C) + FCT (A - C)) 
      C = B * .44816057788302606                                  
      Y = Y + .21417949011113340E-1 * (FCT (A + C) + FCT (A - C)) 
      C = B * .42468380686628499                                  
      Y = Y + .25499029631188088E-1 * (FCT (A + C) + FCT (A - C))  
      C = B * .3972418979839712                                  
      Y = Y + .29342046739267774E-1 * (FCT (A + C) + FCT (A - C)) 
      C = B * .36609105937014484                                  
      Y = Y + .32911111388180923E-1 * (FCT (A + C) + FCT (A - C)) 
      C = B * .3315221334651076                                   
      Y = Y + .36172897054424253E-1 * (FCT (A + C) + FCT (A - C)) 
      C = B * .29385787862038116                                 
      Y = Y + .39096947893535153E-1 * (FCT (A + C) + FCT (A - C)) 
      C = B * .2534499544661147                                  
      Y = Y + .41655962113473378E-1 * (FCT (A + C) + FCT (A - C)) 
      C = B * .21067563806531767                                 
      Y = Y + .43826046502201906E-1 * (FCT (A + C) + FCT (A - C)) 
      C = B * .16593430114106382                                 
      Y = Y + .45586939347881942E-1 * (FCT (A + C) + FCT (A - C)) 
      C = B * .11964368112606854                                 
      Y = Y + .46922199540402283E-1 * (FCT (A + C) + FCT (A - C)) 
      C = B * .7223598079139825E-1                               
      Y = Y + .47819360039637430E-1 * (FCT (A + C) + FCT (A - C))
      C = B * .24153832843869158E-1                               
      Y = B * (Y + .482700442573639E-1 * (FCT (A + C) + FCT (A - C))) 
      RETURN                                                                    
      END                                                                       
C**********************************************************************         
      FUNCTION FEXP(Arg)
      IMPLICIT NONE
C*** Start of declarations inserted by SPAG
      REAL Arg , FEXP
C*** End of declarations inserted by SPAG
      IF ( Arg.GT.80. ) THEN
         FEXP = 1.E36
      ELSEIF ( Arg.LT.-80. ) THEN
         FEXP = 0.0
      ELSE
         FEXP = EXP(Arg)
      ENDIF
      END
