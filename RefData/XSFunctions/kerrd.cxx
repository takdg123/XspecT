// Optically thick extreme-Kerr disk model based on the same tranfer-function used in
// the laor Kerr disk-line model. Local emission is simply assumed to be the diluted
// blackbody. See Laor 1991, ApJ, 376, L90 for explanation of the transfer function.
// See Ebisawa et al. 2003, ApJ, 597, 780 for examples of using this model. 
// Parameters :
//     0 - distance (kpc)
//     1 - spectral hardening factor (Tcol/Teff)
//     2 - mass of the central object (solar mass unit)
//     3 - mass accretion rate in 1e18 erg/s
//     4 - disk inclinatino angle (degree, 0 for face-on)
//     5 - disk inner radius in unit of GM/c^2 (>1.235)
//     6 - disk outer radius 

#include <xsTypes.h>
#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Numerics/Numerics.h>
#include <XSUtil/Numerics/LinearInterp.h>
#include <CCfits/CCfits>
#include <sstream>

// prototype for routine in laor.cxx
void readLaorFile(RealArray& laorRadiusArray, RealArray& laorMidRadiusArray,
		  RealArray& laorAngleArray, RealArray& laorEnergyArray,
		  std::vector<RealArray>& laorTransferFnArray);

// prototypes in the file
Real kerrdBlackbody(const Real kT, const Real en);
Real kerrdTfun(const Real r, const Real m, const Real mdot);

void kerrd(const RealArray& energyArray, const RealArray& params,
	   int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray, 
	   const string& initString)
{
   using namespace Numerics;

  static bool first(true);
  static RealArray laorRadiusArray, laorMidRadiusArray, laorAngleArray;
  static RealArray laorEnergyArray;
  static std::vector<RealArray> laorTransferFnArray;
  
  const size_t nE = energyArray.size();
  const size_t nF = nE - 1;

  fluxArray.resize(nF);
  fluxErrArray.resize(0);
  fluxArray = 0.0;

  // if this is the first time through read in the model file

  if ( first ) {
    readLaorFile(laorRadiusArray, laorMidRadiusArray, laorAngleArray,
		 laorEnergyArray, laorTransferFnArray);
    first = false;
  }

  // Set the physical variables from parameter values.

  Real d = params[0];
  Real fcol = params[1];
  Real m = params[2];
  Real mdot = params[3];
  Real inclin = fmax(0.0,cos(params[4]*M_PI/180.));
  Real Rin = params[5];
  Real Rout = params[6];

  Real cnorm = fcol*fcol*fcol*fcol;

  //  find the tabulated inclinations above and below that requested
  //  and the associated weights

  const size_t nlaorAngles = laorAngleArray.size();

  size_t inclow = 0;
  while ( inclin > laorAngleArray[inclow] && inclow < nlaorAngles-1 ) inclow++;
  if ( inclow > 0 ) inclow--;
  const Real angleBinSize = laorAngleArray[inclow+1] - laorAngleArray[inclow];
  const Real wlow = (laorAngleArray[inclow+1]-inclin)/angleBinSize;
  const Real whgh = (inclin-laorAngleArray[inclow])/angleBinSize;

  // set internal energy range and bin
  const Real ene1 = 0.01;
  const Real ene2 = 20.0;
  const size_t nene = 300;
  const Real dene = log(ene2/ene1)/nene;

  // if the inner radius is within the tabulated radius array range
  // (note that radii are tabulated in decreasing order starting at 400
  // and ending at 1.235)

  if ( Rin < laorRadiusArray[0] ) {

    // set max radius to take into account
    Real Rmax = fmin(Rout, laorRadiusArray[0]);

    // selecting the integer values of the radial bins which are at the inner and
    // outer disk boundary (according to inner and outer radii)

    const size_t nlaorRad = laorRadiusArray.size();
    size_t iRout = nlaorRad-1;
    size_t iRin = 0;
    for (size_t i=0; i<nlaorRad-1; i++) {
      if ( Rmax <= laorRadiusArray[i] ) iRout = i;
      if ( Rin < laorRadiusArray[i] ) iRin = i + 1;
    }

    // outer loop over energy in the emitted spectrum

    for (size_t iemit=0; iemit<nene; iemit++) {

      Real ena = ene1 * exp(dene*iemit);
      Real enb = ene1 * exp(dene*(iemit+1));
      Real enn = sqrt(ena*enb);
  
      // Set up energy array for the line. Note the use of the high energy bin
      // and the low energy bin to suppress error messages from the interpolation
      // program.

      const size_t nlaorE = laorEnergyArray.size();
      RealArray lineEnergyArray(nlaorE+2);

      lineEnergyArray[0] = 0.0;
      lineEnergyArray[1] = (3*laorEnergyArray[0]-laorEnergyArray[1]) * enn / 2.;
      for (size_t i=2; i<nlaorE; i++) {
	lineEnergyArray[i] = (laorEnergyArray[i-1]+laorEnergyArray[i]) * enn / 2.;
      }
      lineEnergyArray[nlaorE] = (3*laorEnergyArray[nlaorE-1]-laorEnergyArray[nlaorE-2]) * enn / 2.;
      lineEnergyArray[nlaorE+1] = 1.0e6;

      // Initialize the flux array
  
      RealArray lineFluxArray(nlaorE+2);
      lineFluxArray = 0.0;

      // Summing the profiles of individual rings according to a given
      // weight function

      for (size_t ir=iRout; ir<=iRin; ir++) {

	Real r = laorRadiusArray[ir];
	Real area = 0.0;
	if ( ir == iRout ) {
	  r = laorMidRadiusArray[iRout];
	  area = Rmax*Rmax - r*r;
	} else if ( ir == iRin ) {
	  r = laorMidRadiusArray[iRin-1];
	  area = r*r - Rin*Rin;
	} else {
	  area = laorMidRadiusArray[ir-1]*laorMidRadiusArray[ir-1] -
	    laorMidRadiusArray[ir]*laorMidRadiusArray[ir];
	}

	area *= M_PI;

	// ems-the radial weight function

	Real kT = fcol * kerrdTfun(r, m, mdot);
	Real ems = 0.5*(kerrdBlackbody(kT,ena)+kerrdBlackbody(kT,enb))*(enb-ena);

	// d is distance in kpc, 1.48*m is the gravitational radius in km.
	// note that the inclination factor is included in Laor's transfer
	// function so it should not appear here.

	ems *= area/cnorm/(d*d)*(1.48*m*1.48*m);

	// sum in the required transfer functions
	// not average of two tabulated inclinations surrounding the input one

	size_t ioff = ir*nlaorAngles + inclow;

	for (size_t i=0; i<nlaorE; i++) {
	  lineFluxArray[i+1] += wlow*laorTransferFnArray[ioff][i]*ems;
	}

	ioff ++;
	for (size_t i=0; i<nlaorE; i++) {
	  lineFluxArray[i+1] += whgh*laorTransferFnArray[ioff][i]*ems;
	}

	// end loop over radii
      }
	
      // divide by d mu. It is necessary to divide by 2d\mu, where \mu is cos(i),
      // where i is the inclination. d\mu is 0.333 in the Laor file. This procedure
      // is required because Laor's transfer function is defined such that adding
      // up all the components (instead of 'integrating' over energy and solid
      // angle will be one. See note attached at end.

      for (size_t i=0; i<nlaorE; i++) {
	lineFluxArray[i+1] *= 1.0/3.3300001e-2/2;
      }

      // Rebin onto passed energies
      RealArray photloc(nF);
      photloc = 0.0;

      size_t inputBin;
      size_t outputBin;
      IntegerVector startBin(nF);
      IntegerVector endBin(nF);
      RealArray startWeight(nF);
      RealArray endWeight(nF);

      const Real FUZZY = 1.0e-6;

      Rebin::findFirstBins(lineEnergyArray, energyArray, FUZZY, inputBin, outputBin);
      Rebin::initializeBins(lineEnergyArray, energyArray, FUZZY, inputBin, outputBin,
			    startBin, endBin, startWeight, endWeight);

      Rebin::rebin(lineFluxArray, startBin, endBin, startWeight, endWeight, photloc);

      // sum into output flux array
      fluxArray += photloc;
      
      // end iemit loop over energies
    }

  }

  // emission from r > 400 Rg if required

  if ( Rout > laorRadiusArray[0] && inclin > 0.0 ) {
    
    size_t nradd = 20;
    Real r1 = fmax(laorRadiusArray[0], Rin);
    Real r2 = Rout;
    Real dr = log(r2/r1)/nradd;
    for (size_t ie=0; ie<nE-1; ie++) {
      Real ena = energyArray[ie];
      Real enb = energyArray[ie+1];
      for (size_t ir=0; ir<nradd; ir++) {
	Real r = r1 * exp(dr*(ir+0.5));
	Real ra = r1 * exp(dr*ir);
	Real rb = r1 * exp(dr*(ir+1));
	Real area = M_PI*(rb*rb-ra*ra);
	Real kT = fcol * kerrdTfun(r, m, mdot);
	Real ems = 0.5*(kerrdBlackbody(kT,ena)+kerrdBlackbody(kT,enb))*(enb-ena);
	// d is distance in kpc, 1.48*m is the gravitational radius in km
	ems *= area/cnorm/(d*d)*inclin*(1.48*m*1.48*m);
	fluxArray[ie] += ems;
      }
    }
  }

}

Real kerrdBlackbody(const Real kT, const Real en)
{

  //     corresponding to photons / sec / km^2 / keV / sr / (1 kpc)^2
  //    
  //    Planck function is 3.14e31* E**2 / ( exp(E/kT)-1 ) [photons/cm^2/s/str],
  //    where E and kT are in keV.
  //
  //    Note, 3.3e-2 = 3.141e31 * (1e5)^2 / (1e3*3.08e18)^2
  //
  //    So, observed flux will be
  //    3.3e-2 * [E**2 / ( exp(E/kT)-1 )] dE dS/d**2 [photons/s/cm^2],
  //    = blackody(kT,En) * dE dS/d**2 [photons/s/cm^2],
  //    where E, kT and dE are in keV, S is the projected emission area in
  //    km**2 and d is the distance in kpc.

  const Real norm = 3.30e-2;

  Real fac = en/kT;
  Real f;
  if ( fac < 1.0e-3 ) {
    f = kT * en * en;
  } else {
    if ( fac < 70.0 ) {
      f = en*en*en/(exp(fac)-1.0);
    } else {
      f = 0.0;
    }
  }

  return norm*f/en;
}

Real kerrdTfun(const Real x, const Real m, const Real mdot)
{
  // x     radius in Rg
  // m     mass in solar mass unit
  // mdot  mass accretion rate in 1e18 g/s

  // this is what Laor's model assumes
  Real a = 0.998;

  // 'special' numbers
  Real y1 = -1.99956;
  Real y2 = 0.963259;
  Real y3 = 1.0363;
      
  // marginally stable orbit
  Real xms = 1.235;
  Real yms = sqrt(xms);

  // Krolik p. 152
  Real y = sqrt(x);
  Real cfac = 1 - yms/y - 1.5*a/y*log(y/yms) -
         3*(y1-a)*(y1-a)/(y*y1*(y1-y2)*(y1-y3))*log((y-y1)/(yms-y1)) - 
         3*(y2-a)*(y2-a)/(y*y2*(y2-y1)*(y2-y3))*log((y-y2)/(yms-y2)) - 
         3*(y3-a)*(y3-a)/(y*y3*(y3-y1)*(y3-y2))*log((y-y3)/(yms-y3));
  Real bfac = 1 - 3/x + 2*a*pow(x,-1.5);

  // radiation temp at this radius in keV
  return 8.24*pow(m,-0.5)*pow(mdot,0.25)*pow(pow(x,-3.0)*cfac/bfac,0.25);

}

//$$$    \documentstyle[psfig]{article}
//$$$    \begin{document}
//$$$    \centerline{Informal note  by K.E. on 2000-4-3}
//$$$
//$$$    The transfer function in the ari.mod file is so normalized that
//$$$    for large enough radius the summation over $g$ and $\mu$ (= $\cos \theta$)
//$$$    will be
//$$$    unity (figure 1).  Namely,
//$$$    \begin{equation}
//$$$    \sum_\mu \sum_g T(g,\mu,r) = 1 \;\;{\rm (for ~ large~enough~ r)}. 
//$$$    \end{equation}
//$$$    This may be rewritten as 
//$$$    \begin{equation}
//$$$    2\pi\sum_\mu \sum_g  \frac{T(g,\mu,r)}{2\Delta \mu} \Delta \mu = \pi \;\;{\rm (for ~ large~enough~ r)}. \label{2}
//$$$    \end{equation}
//$$$
//$$$
//$$$    In figure 1, reduction of the value from 1 as the radius decreases
//$$$    is the effect that the photons are not escaped due to strong gravity.
//$$$    The figure reproduces the values shown in Laor ApJ 376, 90, 1991;
//$$$    0.82 at r=6, 0.5 at r=2.17 and 0.1 at r=1.33.
//$$$
//$$$    On the other hand, integral of a unit specific intensity 
//$$$    ($I  \equiv 1$) over  hemisphere
//$$$    will be 
//$$$    \begin{equation}
//$$$    \int I \cos \theta d\Omega = 2 \pi \int \cos \theta \sin \theta d \theta
//$$$    = 2 \pi \int_0^1 \mu d\mu = \pi
//$$$    \end{equation}
//$$$    \begin{equation}
//$$$    =2 \pi \sum_\mu \mu \Delta\mu.\label{4}
//$$$    \end{equation}
//$$$    This is the total flux from a disk-like emission in the Newtonian case.
//$$$
//$$$    Comparing equations \ref{2} and \ref{4}, we can see 
//$$$    $\sum_g T(g,\mu,r)/2\Delta \mu $ is what tells the polar-angular dependence of 
//$$$    the disk emission, which corresponds to $\mu$ in the Newtonian
//$$$    case.  Figure 2 plots this function for various radii
//$$$    ($r$=400 to $r$=1.235 from top to down at large $\mu$). $\Delta \mu$ is 1/30 in ari.mod
//$$$    file.  We can see that for large radius this function is close to $\mu$,
//$$$    which is expected.  As the  radius gets smaller, the emission will 
//$$$    eventually be isotropic, then reversed (more flux for edge-on disk), 
//$$$    as the Laor 1991 Fig 1 caption says.
//$$$
//$$$    Let's assume $I(r,E)$ is the local specific intensity from the disk.
//$$$    In the Newtonian case, the disk spectrum will be, 
//$$$    \begin{equation}
//$$$    f(E,\theta) = \frac{\mu}{d^2} \times \int_{r_{in}}^{r_{out}}  2 \pi r  I(r,E)  dr 
//$$$    \end{equation}
//$$$    Replacing $\mu$ by $\sum_g T(g,\mu,r)/2\Delta \mu $,
//$$$    we may calculate the observed disk flux as follows.
//$$$    \begin{equation}
//$$$    f(E,\theta) = \frac{1}{d^2} \times \int_{r_{in}}^{r_{out}}  2 \pi r dr \sum_g  I(r,E) T(g,\mu,r)/2\Delta \mu 
//$$$    \end{equation}
//$$$
//$$$    \begin{figure}[b]
//$$$    \centerline{
//$$$    \psfig{file=flux2.ps,angle=270,height=8cm}
//$$$    }
//$$$    \caption{}
//$$$    \end{figure}
//$$$
//$$$    \begin{figure}
//$$$    \centerline{
//$$$    \psfig{file=flux.ps,angle=270,height=8cm}
//$$$    }
//$$$    \caption{}
//$$$    \end{figure}
//$$$    \end{document}
//$$$
//$$$      begin 644 flux2.ps.gz
//$$$      M'XL("*S*Z#@``V9L=7@R+G!S`,5;27,<MQ6^]Z]`#JY*#FUA7W2S:#$7NLPJ
//$$$      MJ2KGL3BFF5`<UG!D)Z7*?\_;@`:&PV4L2BDM'+Q^`+[WO05+#[_[R_F[^8>+
//$$$      MS2_KV7VOU=OS=Z?X8?KNN]/-]K5:_W)UM_IC!<WW5[OK]6MU_O?SLY_?J_/-
//$$$      MW>[=A^W5[4[=7F]V\/QDNU[ML`MK5,G5YN;'U0XZ*C?_<+N=K=9:&?O:!]!X
//$$$      ML_ET<W%U<_EF\^_7ZJ^@=G/Q-Q#_N/GPZ>/Z9G>ZN=G=]0_.5C>7GU:7Z[/U
//$$$      M[^OKU\J`Z.?M%6C2/*\5/+^X^["Z7<.#<]`;.K^]N3C9?,1Q[W#J]>75S?EV
//$$$      M<[VYG%Z=J<\?-[^O=QNUO;ZZP9\?/FVWH'F[N;K9J;O==O.OM6*5_ZI?KFXN
//$$$      MU,7ZU^G5B?I\7(<?VT0:_AS7]]W9/]3GH#Y^NE9WZQUV_>/J8O?;H')RJE[=
//$$$      M;F[5]68EHC?GZC-HWJYVOQT8\^Q<O:HHECYOSSN[KC=W:^J]WOQZ=7W=]_[I
//$$$      MC?I\>;?Z?:UVV]7-W350K7XZ_>%$77RZ5>"'Z[4R`U9E:Q.<!`1(XY]@N!*0
//$$$      MP_AOU:O+[?H.PFK=P3L!VK5R4:O5]H/PU7<['9[O@:8X$+\K"1.()(HEC@D0
//$$$      MO%OO/MU.KR36R4+LJ[_7R2K^G\R;L@E:6?BW$%#`KQL(R/5D!OM,;[I1Z$WX
//$$$      M_^14$+59&=+!S("9-9%VN5W]1\;(KA1@TF@`D;)69]"(UA5UHF9Y=@+_9I%I
//$$$      MY8,RQ?FF#&TP(`WM8ES?=M9F;D_<]CD,SY,>]7.Q?=OKD%L;^/%VG"^8<;R0
//$$$      M.GUH1Q-]/W_TJ0S/DXM#NX3!OF3\@"<Y77H\R8<!3_;C>$7KP?[BS3`^\1E+
//$$$      M\@.A@P`9[05$Z2!`3@<!DLJ":6&UUR!:FZ#R&HLM`[%]%V*V%Q"U@P"Y[:<E
//$$$      M<@<-9+<7$+V#`/GM@1'!/3!BN.]"%/?3$L>B86!XY8T%/GS&=L:_%-;&1HQT
//$$$      M:,Z@8J,!O!94X`%)3U`^.TD"$M1VT^#.)`#Y!$UIB;KTKJ-1<ZZ/9U&?N3=C
//$$$      M3<X\A14H>1CK]$RLZLNP5B@AW8,R'4O;%T(!CTL1,]8F"(*BFR`'C)M%P_IB
//$$$      M0+!H.&/]H.%2B(.&=SD/&D';<8P0@X-I%XUH(7=ZC5A,&Z-X[04I:Q0'('ND
//$$$      MI-$C)8T>*6DPTFG1Z)&21H^4-!AII]$C)0U!&HGEE),R*3J.1>3[!'W!KLM$
//$$$      M/C0C!SX]C35BJ37'Q6,L8'42:&G/+;)KGWF)>@Z>"`!-I2CR".T']L>0TRH;
//$$$      M@(O5DH-Q5*.Q,``8)8\<:]8MH6B;>3(](N1PS-RW,N.TS1TSTT%FU+=@AD-C
//$$$      M.L0,@D;\!-@Z6\LJ6HJ<=?+%,F^+?=+GARV;7M8R">F#ELWLJZQA.0L&<Z..
//$$$      M,W7CM'HS6.(T80'::N$9+*DQDVO=L7W=B7T%%!R#'9/,TJ:O"K7+0G4((?U)
//$$$      MJA\,HNE/42VUX2#5'5Y8GWA9G<>LB2VO(`YKEN6.W)IG><BS.#(%J*AO@[GG
//$$$      M["Z,:H8W]_!\&(!U"'';HJB'D0Y,"L^R@14M%XJGN@)SKIRH<1<]&^,AN[QL
//$$$      M.#(LA6?3G*$`S8&&:_0'$T3F'7PR)'0IP<3L-TA`^`\_V@*N\@3'Q@*]6.J1
//$$$      M!]J96\`R9U:`O=)L-%M5<&"H@.1^##A'.@90S<9SG$3+_^'G@)\3RSWVS32\
//$$$      M<1D'9AWP]6PKBP'GEHQ`N>//:!H>8^@SRF&]0`P:8-I$1AO8CLXVD]5PK("P
//$$$      M)''!?0%S42R3@!^!-.=)"O$VNT!F8!ESL.W!CZ1+'V&_.>,*MN>9B?V%1QA@
//$$$      MP5#-`TN'-C#B0VMW';RG508+A(=GQA=)T"7"NZU)JPP4>X!2>@5:-:$=(PNB
//$$$      M#BR8FR353()_L)CC`F?Z6@!!.18%Z_NB8+VJ1:4F&*6?T;;:DDLYA)U_P,Z;
//$$$      M<.1L&)E`+S%([`/0E!9!EV@=":VC-1".D`!B$;++(LXTJX8-884\%L;[&T+[
//$$$      M".ORN$.`Q702=WB`S9)]R+$11[SHOBX0;$=%;JB5%7:MY]W6H6TG:)\J1;`]
//$$$      MVZMVPX!+F:L:7);1>=;K(JORP)"V<BX7_UD7)<(MN<MZ5]UPK]RULO?5G+&D
//$$$      M@`W!#7'4!(?CJ+.<MYQ/6YY,Z%/5)@SWZ5FIV@5-CJ&/?UMWD$L0B>3I(!J"
//$$$      M'TZR[NL'_Y+MSF`INV_^D.T.C\*\RE,?JR72]Y;[V%/D\'*@.H0%J2RUZ[%C
//$$$      M\#,L>.!H>6]CU7.X0'.QC-YBR71DRKOHRA>F_/),<_X_,^6;+V(*S_!%',NU
//$$$      M2W:?`)(\DX"N5PG/K)2\6-+D1><EI8E'D4@5W!O&#L,M\5ZI77;M=I^R&C(U
//$$$      M=^IZ[;6US(?C=96NMLXFVC"(Q#B\;VO+Z:&)NCVFIZ,TC?K_6`%:.?,N^N?N
//$$$      M/)9X\!Z-G;IVRF-\5,EQ\>&#]\\L@LO<,>S/+9(CYT[VN7.W9=!G5V^EOMTR
//$$$      M."UEN05HT&:Y'N,DJ:*O"V$IU\'0A?/++#A/E^OI4+F67`R6`J#+SN"<[O<O
//$$$      M`=^^?(GCII?93`:?ZAF0]R]-<'#_,I5:%$,RDK?[)]Q'BQH%1EKN,8X/C.EH
//$$$      M$XOU^Z%9[+=%('?'+=>KY-GK$.98Q-/L>/1@T5???>%A7`(Y.G2>7`.B*7#^
//$$$      MEPM"W8VHQ^A9PJF>Y[QMB1L=A1RN7!U$WB#H_?Z\WH@+\#P9O8GR&@$5ZOC+
//$$$      M\5$.D_5LN?`Q,B3O'P*>I;Q!DP+V[=ICK$\]Y\/V;;G-J34C,`JT.>(4N`^B
//$$$      M"H%O2Z%:)%A3/=>'D;;8BM#^K=`<NXWXX1NOMH8W0\>KJ3X"N!`98[!"^7J"
//$$$      MN7>#A/OPD\F3KHX6=#4?RX@[8WRLW?E"V9CBJ]*]7<4(9'D-E-&%>)LK_?@F
//$$$      M"+8W=;YI;K*8%PP`P!-^PYM_X:';P4@T&*L>*EJC(\WRWJE>D%-2<_&'&0*^
//$$$      M!S56CK"QS3?7Z&5>8H@5ID0RBTP>KO'WW*L>=N]TI'LKD)0&_^)?$AT-Y/[E
//$$$      M\W%`LFF.$QPD81SFX"5JFV$?PKB[7#SWZ$6J`'$%')I\P1*2,A7:1>`[@8M+
//$$$      MG'N-K^])9^*X\,;32UQ_9)Q[JV/K%T4$TXJHQKG'HEMEE*G>YUPA/)RIRK/N
//$$$      M`DZF#=8O<'%-8@&7N8-7V?78$&WLD71;H2^H/@_=^(]%OSE3C$BF0IED966)
//$$$      MTP_&4(<SZN8-."K8X_I1W?8E1R%M.JYP/Y+91U(W,1?!Q,$M\%=$%+=?#\?4
//$$$      M^PAGQ3N3>[.F%Y@U/6`]+=K.FIJT[%*1'.'2.A1LS>I0U286'6)R>GF/\D(2
//$$$      M\*9&8!Q\%RQ^C_AMDEHL"&P*;@&;]]R>NT*QG)X61.W%7$>25,N^W]0=AY:A
//$$$      M8Y<`E4V\;ZELRAH:<FB8<4<$Y2YD4VKY$<-*<4.%$@$EV^$*Q1VC<6DO$T3T
//$$$      M8IDP/>X_,3W:8/8"241?"F0Z)I"B)=X>#*1*FRVE6R/`3='YI`XO:=-C2QKL
//$$$      MQ]T8DC%8\^(A.>V%I'HR)&'7HU7.^-4I>96`#SW_<(I>P"5%+R*-X1'PS1S^
//$$$      MS/B1M/`168JO(&M3'B<<#8!RYT0C)AH>]:/*>,..;Q[/2,Q-3?=Z-*WA`U_F
//$$$      M21,/.DZY0)CY^838^!6H7"CE:@FJTXO;7&!GP_.*F8[-]#P3#!V`F0*+H#%T
//$$$      MG`0'JJ)A,\_P<,:!%WP#Z@2-$)$7Y$8@"7OR4`]=>`R_-ZB`2PU4P9>M`@HM
//$$$      MXK:6=H=I:KX2WM"XGKC:3$_YBEYBX\PNUIEA\\4B..ZRJ"E5R;0HP8F'E<0$
//$$$      MBR]@.Q.HS29D,:&Y;C1AS_=N-$%8-G(=V9G@Z8M.'K:6['5\_S[;TMM+7<7_
//$$$      M>TPUYAJ.1USE6@Y-D;]OE?"KB31IIJ@O(5<<`&%/XO!NHP0XQ>7$C&'WZ)PP
//$$$      M1H11\XD\F9[($_5XGHBKDJXI@`"X3<6"M8;IV=P$JU?A[V0Z!2<%++V:^G/0
//$$$      M^19=)O34N^:!-#*?!#.#=WSQ<@*[N27SO(*UT%;?BNY$;JNKA1&W9I[1:/:]
//$$$      MI*6TN$'/?-N`&\;!&_`D[>KOU(,P6N-W_!#%M"#^5BBF!46*>UQ\`Q23<(&^
//$$$      M-1J_EL*O'*E<J5KQ\2?%+J8\Q86_'Q=^+R6?CHNIBPN)9DIEHYVK%:H-Y/FQ
//$$$      MV].B%XH'E7#L@%76XR^BN'831DLN<MT]F"F!09)*E7CYC@\/*BN'%`XCBT!^
//$$$      M8FVKZRWJ=X5G:B%0H5+AD=]"J+\'<??;YH_;U67]%9?WV]75]7I[Z-<&/%YI
//$$$      DXGY4I>3IEQT6Y4<4[_W^S?++-/A+$F]_/IW^!V>$L0HU-```
//$$$      `
//$$$      end
//$$$
//$$$      begin 644 flux.ps.gz
//$$$      M'XL("%?9Z#@``V9L=7@N<',`Q5U);QS'DK[7KZ@Y&)@YE)7[HIM-6W/A@P78
//$$$      MP)QIJ9_,&9I-D)0]@#'_?7*++S*KJ[M)D]*#+;(JUX@OEHR,JBQ^\V_O?UZ^
//$$$      M^[C_=;?H;\7\X_N?W^6+Z9MOWNWOW\Z[7Z\?KOZ\2K>_7#_>[-[.[__S_>5/
//$$$      MO\SO]P^//W^XO[Y[G.]N]H^I_N)^=_68N]065'*]O_WAZC%UG/7RW=W]HH00
//$$$      ML]1O34@MOM]_OOUX??OI^_W_OIW_/36[_?@?J?B'_8?/O^]N']_M;Q\?^HK+
//$$$      MJ]M/GZ\^[2YW?^QNWLXR%?UT?YU:EGG>SJG^X\.'J[M=JGB?V@V=?[S]>+'_
//$$$      M/8_[D*?>?;J^?7^_O]E_FMY<SG_]OO]C][B?[V^N;_/O#Y_O[U/+N_WU[>/\
//$$$      M\'B__Y_=7)O\W_SK]>W'^>/NG].;B_FOYW7X`1.)]-_S^OY\^5_S7W;^_?/-
//$$$      M_+![S%W_O/[X^-O0Y.+=_.9N?S??[*]:T??OY[]2R[NKQ]\VQKQ\/[\A*KC/
//$$$      MC^\[OF[V#[O2>[?_Y_7-3=_[']_/?WUZN/IC-S_>7]T^W"2HYW^\^^YB_OCY
//$$$      M;DYRN-G-<J!U5G2;A)0`:#?_G1B?&Y'#^#_.;S[=[QZ26NTZ\BX2[&+63LQ7
//$$$      M]Q\:7GVW=T/]BNBB!TWN<U.3I$E%EZI.I(*?=X^?[Z8W3=<+A[FO^%9X-=>?
//$$$      MA;TI2"MFE?XQ`#')=9\4<C?)@3_9LR[G+,WT\^)=HPBS5I(V+2/-G(RG_YF&
//$$$      MO/_TZX?$S7T;,^@8$[)2)*)\$/-ENG%*Q_EB7EK=1?JWM#*1Z96J:VSLK/K.
//$$$      MJ5X%NI]*O79CO3%C?ZM,[.N=Z.Y3O0MCO7=C?3!C?<1X$]/KHHH#P2YZ,U#<
//$$$      MMR@D]RT*S7V+0G3?HE!=6TQ,=M^BT-V/40CO6T3%+189$NR)>&-FES0HEZA2
//$$$      M>I'+%]TD4PKHGEI,M7,I*.5A;G>M>>M-HY7;A:J7UGQIO1LIRMNGDS*AQ>N2
//$$$      MX@HE65[6AD2)JX.[VJK\2O"G)E$D8HTJ35)G:I(N\R]E9L*K=J=?I:EHOQ<J
//$$$      MZ/D5Y=]%DEX/CM$$SO0,.;TJ.%,%IZCJ47`J`XEVDYM4:49!9BZ5ZRL:;S9:
//$$$      M?5;PTQ?EK0F^&-EQP1=[F9WTR?JL:0QL:/+%$3[*[W2?M*/QH3H^W)J-JE5E
//$$$      M%('Y`W&A>B[<UO10L<!UWCX!ZZ-Z-+TBUL5=G="C1J]7GNAMUC*1^9#UD)4E
//$$$      M6*M"L&*,#47'+]!Q`W"*K6Y0ZKD-T%21YB4!MND[#\`]5#<`QF4'X9);5DFA
//$$$      M3$8B53=Y%UNYJ+93EL^H4H,IN?5V+[7VJ4=7X(T;"I0L2MH56*?)'FM!+,AR
//$$$      M"ZW+.M<5^!"':8V,86AAG'##H%9(.[2P1HUTV*#':9TR<BS("^3E5`:-VL;*
//$$$      M?6T037*]/?>UH..^%/3<UX*.^S)FSWUI4;F?NH+"?=>BY[X6=-R707ON2XN>
//$$$      M^UK0<5\**O?=M)7[N9JN#SYK2%6.JK9:0-%%O6UWH\N`YI$1UH+6?&F]ES;:
//$$$      MM&!\.(W>TBHE*GCWI4B9B)3Y-"G)8H0DZ:T6G]J_KJ&97BUB;F6;:^CL.;#'
//$$$      MJG?3RI[9UQUZ$Z(H]"2R-U&8+5T&*6<3HVSK1?6GJGCL2;E$D*\AI`PJMU'I
//$$$      M?@G)ZUV42U]<2VH65"L)H5W$V*JD$*U(2FHEE:#:Y%U:F=9T972;75KTL)IZ
//$$$      MV$AEWM!5I/%4XH?*J%8)B5JB12F%*V)"@0(E7:$@UQK4.M32R!I4+<E-YFMY
//$$$      M<M^16T95$4WS"D)4IJ".(-4=I(6$0/-&!4@!E8@;D*)6`PQ#3$H#2`VQ)AW:
//$$$      M>?2-$$)D2-$W^D-(`8N2@6I!2P^X15^B7FD%2"5#:DI$>Q+26B)/`!X(<%^<
//$$$      M7P7<@S``[DEU`ZEWUJ`"003*@NJD!'H*:JC75ZXO,S`)!Y-@%0[H$5`;2`U)
//$$$      MF1-F@I44R(L(1"UJ_0;>Z&&@PL8RWFG9>AK>/>IKO'W#6V=E:G@G@Z()'=$:
//$$$      MY`'>ZA!ONX$W=)215W`/!CV,!][08`?Y>5B3YS(/O$D+9<1XH$H)("H8;XDK
//$$$      MEA5D`*>EX6X2*'GL9^CW-M[6-;Q=6<LKWA9&"[S!9B0MR,YZ6OEH`4@E()4.
//$$$      M,*.=AKO1\1!P]M$.+L.AAPL`/!#@'L:!)42"4!*"&X0`AR(VC,,3_UIT@&<^
//$$$      MS^6&SBIXYK,@;J+#HF@$7)@C8D%#,-U%U7`BE-5+\(H95V6\.F:\><4$W@82
//$$$      M@F%W>%N4>5I4I&.\V=U`5@$C=ZNH/9!&I^!!`V\#O(523\3[E`-/,!:XK2W)
//$$$      MH0JWY2"$.#%Q#7>QK^)K6+]YZ00^')?`:4HXS6X1U7`9&NTL1C8<H7!MA'ZS
//$$$      M]K.[8;FP?K,,]($TND74K9'/@5ZTYQT*(6V/XITYKH!''PAP*PQ63#(HAQ*$
//$$$      M?Q0W,I(NKL`=`!?L6J"074@(+ZP99O;O@!22[P"WB'AXM>T<#T8&[20$AGX$
//$$$      MG.,<.)2T:#]EQ;1G%'S1S:&D?3:63*OA"+%40>BP5&_@PAE3>!MHE>0PD;0Y
//$$$      MLKI"!H(5'$@9J"N[(`.3R-':U!O"X&YX9`O9L])SN\`.!>X&X;9$W)O"MN)`
//$$$      MGQ,1VBV\,QX5\.#AP:W'Y!P3(BP7':8-<"!IQHL!;US%#7S$5EGLL5U+0R,B
//$$$      MU+"7;I/$KBJ@EF,?(.HW8LUN816,=Y;Y$R-"UO(#AY)3(,6A*(T0W/`"CK4(
//$$$      MJR*V0Y'AAAX">"KKW4V_9*I#<6`[)[%<]_LE]%"#UD_K)0&H<3N.0)_JJ]@+
//$$$      MA@[QS-V3/,JI$$4908B7%$-%7",HA`)$^)AABU-C%'>(.'/9[80X'H=.2H[,
//$$$      MQ2&F:AU0NLXQ28F@D&L[^8:^W<'(HF>CR8-F@U91APQXSI$])R;<=BDYKJB`
//$$$      M&\$QH<;&`QXNPH7+7@\JX'+%&D,_+(L,AF3X6"710ZZA&J-ZCCSA5!"#QG`X
//$$$      M&5:'R$7BH-*$`UF8`+1-/!^!BV^U%L//`[0#H:T%AX0:D5H@8@)3"H^B#T(4
//$$$      MWFXJ=AE0ZBU5/J@=R^P*']<%.LVC,4F]3N!"KQG`L@G/#O-4[/^!LRY+V$F<
//$$$      MDQ2<S^7\\]"-4"28\V]PW'"G`:H&7\[[&<F5:^WS=4^Y=NSL.`5'+FZC#++L
//$$$      M,P=M+JKCM!G1X:D$892G<-52KRXOTGIA>[M@A[]T0*L2'YU6Z!KW-)75GD..
//$$$      M%%2:`VW8]JRLF.PD3`_FO-XZ\D8'SI[Q);T,-#_VJMBH8!<#B]8T#:="+2"A
//$$$      M;AXEI!(+`L<%@<R"-.1B2!/2&(7U\SE415M";366-V745F#`ULO`L28AX=/Y
//$$$      M5)CE%I;$RN%2BHT$G#[\GN%=!`U$CE:3`!9D2A3<);1QP>9K0=9V022Q6-1:
//$$$      M]`TA/L$#G-M[JT!(:\W9)<X*=(JIUD'"N$3S4P"^PN8%99P7(7W"_I$5DR8'
//$$$      M*(@[V!42T`N6PX4SL4`KDNHA@;)@0[]@S[+`%!9L>>DJU?IB.R](FT[5UUK*
//$$$      MFUI^,!`[5#F@X="5`RF.$SB;P;X3BP;PA;:2$.`5D<:#M3(1'#>I%4QYXT"5
//$$$      MCLV<<87/A/R*RC?EAY0<7(ICU!W7VB=B?2IE"O^A8M*ARTJ"=!R2\>.2;B/`
//$$$      M:0E>\WFU.XR$D'#%!32*-]`4^[&*L5GS)@Z5>#ZS8/^WN+7&LC/)M68#3]9=
//$$$      M`0FL<2^2.AM%/"%?JBF?I/(.FS3;,'::G_?Q,\"-Q:V+J^`4"`4X7:"!+1>,
//$$$      M9-'000C30`5YT5AXS6>_Z]4AWAVVK*'0_8,KM]&C2:-YH;)BO#1?FM_&JGB;
//$$$      MR$^\NJ2-9Z\\[!#J2FAX=<1#"PZCB4]61#""O/>"87O$P2]R*8MA/:2HN/._
//$$$      M'O+HL+(;\EC["9;,@+B-JZM"U?GM]OF,J1&!'NI&:'C@19(3<:`!X:E$$`H^
//$$$      M$+,@A:\Y6CKDP?!Z!BTT[#'8KN&?$:H4S9QZSS)X8&\V\.3Q]"$M%I)'2+T@
//$$$      M,YCT(KQ"NM3HME*J['4I>]?EOC9"$HZ:!<$=(`HBE%<W=AF&&0($S!HR5YUC
//$$$      M9E@"+VJ^AW3JS6!T&1L.A9VZW:"%5PNS=EKU$?HKY$LM/4*7000.`CFM3KAQ
//$$$      M6IH0TH0VK)JW#)P6@.QZ[4(M^PO/<*M56;X"W-[U,$XK474>(:RN7!]+=X#"
//$$$      MBW%$SFN$YJ#(/R6[<2Y=:IT@N&WD]9*3%'#?R,OAS95%<DC+;*H>CB.@#OK:
//$$$      M@0KF`D.^):2-H(Y!=>R_PHHJ-^P%]2&H6*H7]H((^)>"QHOSI5ZW)3-)D"`G
//$$$      M0J$X3`<C2\%S9\-XCV1!&FI!YG$)7&L/D>5$,./9.6M>$#GD@P/G%V?8G#ID
//$$$      M,3(6\X6U1K.+1BU>KEB0]%WB^<<!YW.E^=6W!C;>R&&'R\$G:RBN`D(E[+\7
//$$$      M;&D6O)34E076;OA,EI#?"MMX66,?S$Z+,@8=HHPR!SF,J&;7Q_K+"21.+VVL
//$$$      M0^XI[OM<MC1(1WCSXQ>%>=DQ(C^T<#YOH>RP:Y=4"D0BK#RR'G,(P2Y[[>3'
//$$$      M<('=+@<QZ*'IL>Z"Q7U!MK=;K1E7+/T]KOS@%KK-^UCD*XOG/(OZN=RI]X%0
//$$$      MMX@)D1U=9$<,VUQ)\K?MB&0?UV4>63#\6#2PKV`OS`O:.DWDL*#Y`6,.U>'!
//$$$      M.9#O7EL2&\ARAM^M\&2,'1EW)AYT8MCS?KSNDVEQ=)93?X+GXBA`\FHBU0&"
//$$$      M(VZ\9G$HN]ZFC!B)#61X$K9J#@ZV=`Y>`F[\$!F.14$&Y#B6Y!$Q77D_^"R>
//$$$      M^7%5Q=/FV*Z)G:U;\C+-3Y$Z2^=M`^L9;SXX2;%IK=V#?9X/WD:O4`(XO2)Q
//$$$      M<`PZ#N#B%Y(.X%('51(&V`_XTI=0<R!>8>Z>6$EV?YS^[W8(W:JC3^+(CX_9
//$$$      M^C;\6AA!ZU7+=<A,1U2+/>^H;#U6ZJ"*>HF3G<X%%>65]'PF4SKCRY'!9+SU
//$$$      M7N0UK1Q(,_EQE9=T@C`IM_=<0,>2QN,RW#&4$QAT.LDD=:8B2G=V)Y!PMF<X
//$$$      M)75X`JD=-.`)N],VM;JC(,KN/7R35DDJ&4=PPTD`)V8ZT)0/P.6C(R)8(CNL
//$$$      MR`XS'TFJAX'HP%`[&$6#TT$AT758#]B=WU(XVY`NLUR4+-BY-4)"M?-V*6C*
//$$$      M3*N<U<NG1I,+SN)2^<7S*@8QC\<5\/O5A3%!.R`,I84?]`@%FWK4<VZE?Q+G
//$$$      MQN732T5^Y=Z69ZPDUB-*T\YH,)VNJ'=W[_RH1%1R0HFF7HF:\JM@_)=7_HI*
//$$$      M@2V49.4!^VK`+'C?SM64/OF%Q\OU`9OVBR&IK^E6PG.!%KH@W3@X<6#X/`?S
//$$$      MF6.6:DMG09K.S\HNITY:M>1))D\H:*W\>12T+DQW]T&NIBXE3YJ:Q:9M>=_E
//$$$      MM-AT?:\11JZMDVSDT]<Q<A99=(2!+X?;^)#KP4GH\41F9QZUZRM1.AVCE*55
//$$$      M#_NQ.T+!MCOBCE&KU;)&13WIOAYB_GL@^RW2.Q6)P9U7D>A"OPX848Y,_JO6
//$$$      M`0XX3#X@-@)(15^6!/971I>G"J^C_.?]U;3EK[2KI)ARL#$5+%1BC>EMVUAO
//$$$      M7B2XE[EC*+YQ(0P6@X+M0!!.P=3CO'_'*;2N?U,QIB>RF%^(7O(C?>WM;/-+
//$$$      M[="-4J+-<"`WGSPO%W2`M1WO;N=;R\9DQK'8J8,G&V=,XP42:).?&]`K>4($
//$$$      MD1Q8@CD*-XD@?$5"U),:^4,QY:&.M>48N<[?.DAXZMF:<J)9UXQG*<&A_#1<
//$$$      M-D_KM'CM,[_3^OCQ?.[,KTV;5>N+5]@\\]O1&^1Y>J?GT3L_A]Y\?KU2$LOI
//$$$      M<!-:VG)8[R='N(UJ^G=/&V-6)]3Q6>=^UNGELU9=SD;A9'E]B8`CI6<;$*.!
//$$$      ML"<:#3<37;Y^8)7*00"M`SJ_:ITB-UV24D09VWMC10ZN+/THQE/H349K\_-1
//$$$      M'K,^N^O5K+>IWE-49ROSXU"O#+:Q;E[)OZ2^3&DKG$IM2^S6ON<@\ULOM?M4
//$$$      M+53F][Q:HX/-8R!"R+[:)"%#FP^XMWXUSR:UP'P+RI+K0EG1;B/S<:_R'2%8
//$$$      M3F=.1*G5#ATK%+5(!MU[*K(*8"&'B'QEQ8/$IB-\*O)ZE9#R=2=&//]?BEY,
//$$$      MR-2KSE,(B8"2Z,@EE0X6TN"/M^RIIVGJR&3*F17BI+8PH9*2-E<)`!/)RM=E
//$$$      M90>A<X(G1=#-6')\D1<T;W*J?*K:U$JJ>]_DH".\K2=+S.]*!!KIA&M8,4)A
//$$$      M5^NLW2OIT6GQ34U\:6OF")^R&NH4WF=\0L$'SE]U*\M74W*=/R?4B".4:M%+
//$$$      M47JJDE=(7'2O#LGT/'D!DKR(-Q5KB-22@LCSS&UB'7V&N55(HGRQEDPO5.%F
//$$$      MON7X,R%2,PJUY+SY3IWYEF]P*:E6^M:*OI153EOZ9E31\2]L@D3)=!+<LE`.
//$$$      MX-:29_A&0M*:8_YM>@4D-RW7V/`T)*<OYVEI2!,D+PX'@T]_<_!C87HW]UI1
//$$$      M\M<#S.!7BQH<"==,*,[/FA:NI1HK`LQD([,OJ*L5SO1=<^PJBQ)$!WUJ)5F?
//$$$      M8`RG]:E(UBIK_^4^N:)I<QYRA>:QX-?FG>0*$JM]Z][\4=IH6FKTG.#7.BO1
//$$$      MS[6BJ#DB:67>QEX1\VE[K8F$LQP$)JY-&UUD<K-(:X$,6*-Z@CM9..'E@4F\
//$$$      M>%EXIDE4)ES^4LL0!+:29SBZUD\I_[Q^1:&=T;*!-GUYIW\\*,Q8.#F()3>J
//$$$      M15\X.)UZ&>4BK^S7#XF==W[T4:WD>?N!#&3T?K",5G#<,IH$?'YK9)1`*_I:
//$$$      M$FA8^'QL?HR(6M'7BHC*[LDK[WG%Z;;EM+EJL.4ODK!ORJD`K>6S7"D)0$?/
//$$$      M?K.P;47'=EBQ38ELUS'9\X9<6Z<O#>:^W]2EQ'EHU_F"M.$7<PC&SK$F4,H3
//$$$      M$%-_Z;GD=/RL2T]91V@O4BRAOGV1;NJ+&:+>TVVMGMI+*ZVS+R/Z,KPH+RV&
//$$$      M*'7YR,-E*:ZWHF3FR[2R/J8+=5(_SB%6)"RM/M-6TU&BQH.!.,FU)E>&_$68
//$$$      M.F]C4U<V39U)E->F0LPOC,ORD"FIPASSLYQ*7IYQP"7]KL!H$!E`>4G`5=Y+
//$$$      MHU8IQBYE#+,:M!'G053,;]$UHC)']5ZT^TU955RFS%P/'-TVV([+2MHVLW8T
//$$$      MLS2M*+_=48K0""7<*)I:-#46\F?$>A;*?64A-!;:@_PPL'!4W5IU0QE:`Q8F
//$$$      M4[X=:H1K4L_O;-5W1\%OZ=KDOT(*R(&.$Z+2;$.N?FW4V_Q)E10-J5"T/MI`
//$$$      M="025B4Z/Z^(-NV+0GU!*+^F'EV*J2IB!;!R^SP[F<[9R3S:21.5%V0"F8!Z
//$$$      M7YQ%;35,7]GU*02/U6GJ.0:376^.QW4S/P/MDK:'7D,"?D3>-YHK\3EQ73+N
//$$$      MMK,\,\?\,+C)%FU-MUKD]^'HGM[%(RGZN=W5FU)G$/BUH6K@Y]L]R=OW1$@A
//$$$      MRF=S-ZF8&A7S5Z#"NQ45TQJ++T=%EJT4^>-J]160XJXJ^^4U=5,_8JR<?IHM
//$$$      M'G7WM`0E$K=MT<P:WG.JKET*K<D_&9JF;]RU<O%8HZR#-GM=DU^7+$Q80_%_
//$$$      MQIXK)@@EK_^-X5''35T-FK?.?K0YL7(?:JU&EV8&;!_D\QL(I@=A&D`@">IF
//$$$      MNZ'.TZ:M`J_NKXJ[#P<(U#K^U/Y0!OVICH??]G_>77VBO\+RR_W5]<WN?NLO
//$$$      I6YC\N'TNFU)ORM_CX,8G&A[\B1C^>R_Y[WC\^-.[Z?\!N@!R(MAF``#N
//$$$      `
//$$$      end

