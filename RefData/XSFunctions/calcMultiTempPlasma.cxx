#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include <XSUtil/Utils/XSutility.h>
#include <Aped.h>
#include <sstream>
#include <cfortran.h>
#include <CCfits/CCfits>
#include <memory>


// clcmek is in clcmek.cxx and is a wrapper to call the Meka(l) code
void clcmek(int plasmaType, const RealArray& energyArray, const Real Temperature,
	    const RealArray& Abundance, const Real Density, const Real z,
            RealArray& fluxArray);

// prototype from tableInterpolate.cxx

void tableInterpolate(const RealArray& energyArray, const RealArray& params, 
		      string fileName, int spectrumNumber, RealArray& fluxArray, 
		      RealArray& fluxErrArray, const string& initString,
		      const string& tableType, const bool readFull);

// wrapper routine for a single temperature
int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType, 
			const IntegerVector& Zarray, const RealArray& abun, 
			const Real dens, const Real z, const Real T, 
			const Real Tb, const Real DEM, const int ifl, 
			const bool qtherm, const Real velocity, RealArray& fluxArray, 
			RealArray& fluxErrArray);

// wrapper routine for a single temperature and no Tb
int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType, 
			const IntegerVector& Zarray, const RealArray& abun, 
			const Real dens, const Real z, const Real T, 
			const Real DEM, const int ifl, const bool qtherm, 
			const Real velocity, RealArray& fluxArray, 
			RealArray& fluxErrArray);

// wrapper routine for no Tbarr
int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType, 
			const IntegerVector& Zarray, const RealArray& abun, 
			const Real dens, const Real z, const RealArray& Tarr, 
			const RealArray& DEMarr, const int ifl, 
			const bool qtherm, const Real velocity, 
			RealArray& fluxArray, RealArray& fluxErrArray);


// routine to return energy array and temperature array from the input file
int plasmaFileInfo(const int plasmaType, RealArray& energyArray, 
		   RealArray& Tarray);

// routine to return bolometric fluxes from the input file
int plasmaBolo(const int plasmaType, const IntegerVector& Zarray, 
	       const RealArray& abun, const Real dens, const Real z,
	       const RealArray& Tarray, const int ifl, const bool qtherm,
	       const Real velocity, RealArray& boloArray);



int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType, 
			const IntegerVector& Zarray, const RealArray& abun, 
			const Real dens, const Real z, const RealArray& Tarr, 
			const RealArray& Tbarr, const RealArray& DEMarr, 
			const int ifl, const bool qtherm, const Real velocity, 
			RealArray& fluxArray, RealArray& fluxErrArray)
{

// Subroutine to calculate the summed emission from plasma with a
// DEM given by the input array.
// Arguments :
//      energyArray   R  i: model energy ranges
//      plasmaType    I  i: type of plasma
//                           1 = R-S  table model (old version)
//                           2 = R-S  APED-style model
//                           3 = Mekal  calculate
//                           4 = Mekal  interpolate
//                           5 = Meka
//                           6 = APEC
//      Zarray        I  i: atomic numbers of elements in abun
//      abun          R  i: abundances
//      dens          R  i: density (cm^-3)
//      z             R  i: redshift
//      Tarr          R  i: temperatures
//      Tbarr         R  i: ion temperatures (assume DEMarr emission measures)
//                          only used for plasmaType=6
//      DEMarr        R  i: emission measures for input temperatures
//      ifl           I  i: dataset number (unused at present)
//      qtherm        B  i: apply thermal broadening (APEC models only)
//      velocity      R  i: gaussian velocity broadening (APEC models only)
//      fluxArray     R  r: spectrum
//      fluxErrArray  R  r: uncertainties on spectrum
// 
// Translated from sumdem.cxx


// At present the max number of elements is 30 (H to Zn)
   const int TOTEL=30;

   int status = 0;

   // Check for plasmaType

   if ( plasmaType < 1 || plasmaType > 6 ) {
     string errMsg = "***calcMultiTempPlasma: Invalid plasmaType";
      FunctionUtility::xsWrite(errMsg,2);
      return -1;
   }

   std::ostringstream oss;
   oss << "calcMultiTempPlasma : plasmaType = " << plasmaType << std::endl;
   FunctionUtility::xsWrite(oss.str(), 25);

   // Set up the inputZ and abunZ arrays which we use because Zarray and abun
   // are input constant.

   IntegerVector inputZ(Zarray.size());
   RealArray abunZ(abun.size());
   inputZ = Zarray;
   abunZ = abun;

   // set abundances - we use the Anders & Grevesse abundance
   // scale internally so shift from the Solar abundance table in use.
   // Note we do not do this for the APEC or new RS cases because this is done
   // internally in the APED class.

   if ( plasmaType != 6 && plasmaType != 2 ) {
     for (size_t i=0; i<abunZ.size(); i++) {
       abunZ[i] *= FunctionUtility::getAbundance(inputZ[i]) / 
	 FunctionUtility::getAbundance("angr", inputZ[i]);
     }
   }

   // Check whether a trace element abundance is being given for the APEC case
   // and if so set the abundances of all the trace elements

   if ( plasmaType == 6 && Zarray.size() != 30 ) {

     Real TraceAbund(1.0);

     string pname("APEC_TRACE_ABUND");
     string pvalue(FunctionUtility::getModelString(pname));

     if (pvalue.length() && pvalue != FunctionUtility::NOT_A_KEY()) {
       int ielt = -1;
       for (size_t i=0; i<Zarray.size(); ++i) {
	 if (pvalue == FunctionUtility::elements(Zarray[i]-1)) ielt = i;
       }

       if ( ielt != -1 ) {
	 TraceAbund = abun[ielt];
       } else {
	 std::istringstream TraceAbundstr(pvalue);
	 if (!(TraceAbundstr >> TraceAbund) || !TraceAbundstr.eof()) {
	   string message = "Failed to read value from " + TraceAbundstr.str() + ": setting trace abundances to Solar";
	   FunctionUtility::xsWrite(message,10);
	   TraceAbund = 1.0;
	 }
       }
     }

     oss.clear();
     oss << "calcMultiTempPlasma : TraceAbund = " << TraceAbund << std::endl;
     FunctionUtility::xsWrite(oss.str(), 25);

     abunZ.resize(TOTEL);
     inputZ.resize(TOTEL);
     for (size_t i=0; i<inputZ.size(); i++) {
       inputZ[i] = i+1;
       bool found = false;
       for (size_t j=0; j<Zarray.size(); j++) {
	 if ( inputZ[i] == Zarray[j] ) {
	   abunZ[i] = abun[j];
	   found = true;
	   break;
	 }
       }
       if ( !found ) abunZ[i] = TraceAbund;
       // H abundance is always 1.0 unless it is zero.
       if ( i == 0 && abunZ[i] != 0.0 ) abunZ[i] = 1.0;
     }

   }

   // Resize the output arrays

   fluxArray.resize(energyArray.size()-1);
   // This will be resized back to 0 before returning if it 
   //   doesn't get filled:
   fluxErrArray.resize(energyArray.size()-1);

   // If calculating meka or mekal...
   if ( plasmaType == 3 || plasmaType == 5 ) {
     for (size_t iTemp=0; iTemp<Tarr.size(); iTemp++) {
       RealArray fluxTemp(energyArray.size()-1);

       clcmek(plasmaType, energyArray, Tarr[iTemp], abunZ, dens, z, fluxTemp);
       fluxArray += fluxTemp * DEMarr[iTemp];
     }

   // Else if the interpolate on table models for RS or Mekal
   } else if ( plasmaType == 1 || plasmaType == 4 ) {

     string filenm(FunctionUtility::modelDataPath());
     if ( plasmaType == 1 ) filenm+= "raysmith.mod";
     if ( plasmaType == 4 ) filenm+= "mekal.mod";

     RealArray tableParams(abunZ.size()+2);
     for (size_t i=1; i<tableParams.size()-1; i++) tableParams[i] = abunZ[i-1];
     tableParams[tableParams.size()-1] = z;

     for (size_t iTemp=0; iTemp<Tarr.size(); iTemp++) {
       tableParams[0] = Tarr[iTemp];
       RealArray fluxTemp(energyArray.size()-1);
       RealArray fluxErrTemp(energyArray.size()-1);

       tableInterpolate(energyArray, tableParams, filenm, ifl, fluxTemp, fluxErrTemp, " ", "add", true);

       fluxArray += fluxTemp * DEMarr[iTemp];
       fluxErrArray += fluxErrTemp*fluxErrTemp * DEMarr[iTemp]*DEMarr[iTemp];
     }
     fluxErrArray = sqrt(fluxErrArray);

   // Else if apec or RS ... call the Aped class method
   } else if ( plasmaType == 2 || plasmaType == 6 ) {

     std::ostringstream message;
     if ( plasmaType == 2 ) {
       status = calcRSSpectrum(energyArray, inputZ, abunZ, z, Tarr, DEMarr, 
			       qtherm, velocity, fluxArray, fluxErrArray);
       message << "***calcMultiTempPlasma : failure in calcRSSpectrum, status = "
     	       << status;
     } else {

       status = calcCEISpectrum(energyArray, inputZ, abunZ, z, Tarr, Tbarr, DEMarr, 
				qtherm, velocity, fluxArray, fluxErrArray);
       message << "***calcMultiTempPlasma : failure in calcCEISpectrum, status = "
	       << status;
     }
     if ( status != 0 ) FunctionUtility::xsWrite(message.str(),5);

   } else {
      string errMsg = "***calcMultiTempPlasma : Unknown value of switch";
      FunctionUtility::xsWrite(errMsg, 5);
   }

   bool anyFluxErrs = false;
   for (size_t iErr=0; !anyFluxErrs && iErr<fluxErrArray.size(); ++iErr)
   {
      if (fluxErrArray[iErr] != 0)
         anyFluxErrs = true;
   }
   if (!anyFluxErrs)
      fluxErrArray.resize(0);
      
   // Save the temperature and DEM arrays for use in plotting

   FunctionUtility::tempsDEM(Tarr);
   FunctionUtility::DEM(DEMarr);

   return status;

}

// wrapper routine for a single temperature
int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType, 
			const IntegerVector& Zarray, const RealArray& abun, 
			const Real dens, const Real z, const Real T, 
			const Real Tb, const Real DEM, const int ifl, 
			const bool qtherm, const Real velocity, RealArray& fluxArray, 
			RealArray& fluxErrArray)
{
  RealArray Tarr(1), DEMarr(1), Tbarr(1);
  Tarr[0] = T;
  Tbarr[0] = Tb;
  DEMarr[0] = DEM;

  return calcMultiTempPlasma(energyArray, plasmaType, Zarray, abun, dens, z, 
			     Tarr, Tbarr, DEMarr, ifl, qtherm, velocity, fluxArray, 
			     fluxErrArray);
}

// wrapper routine for a single temperature and no Tb
int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType, 
			const IntegerVector& Zarray, const RealArray& abun, 
			const Real dens, const Real z, const Real T, 
			const Real DEM, const int ifl, const bool qtherm, 
			const Real velocity, RealArray& fluxArray, 
			RealArray& fluxErrArray)
{
  RealArray Tarr(1), DEMarr(1), Tbarr(1);
  Tarr[0] = T;
  Tbarr[0] = T;
  DEMarr[0] = DEM;

  return calcMultiTempPlasma(energyArray, plasmaType, Zarray, abun, dens, z, 
			     Tarr, Tbarr, DEMarr, ifl, qtherm, velocity, fluxArray, 
			     fluxErrArray);
}

// wrapper routine for no Tbarr
int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType, 
			const IntegerVector& Zarray, const RealArray& abun, 
			const Real dens, const Real z, const RealArray& Tarr, 
			const RealArray& DEMarr, const int ifl, 
			const bool qtherm, const Real velocity, 
			RealArray& fluxArray, RealArray& fluxErrArray)
{
  RealArray Tbarr(Tarr);

  return calcMultiTempPlasma(energyArray, plasmaType, Zarray, abun, dens, z, 
			     Tarr, Tbarr, DEMarr, ifl, qtherm, velocity, fluxArray, 
			     fluxErrArray);
}

// routine to return energy array and temperature array from the input file
int plasmaFileInfo(const int plasmaType, RealArray& energyArray, 
		   RealArray& Tarray)
{
//      plasmaType    I  i: type of plasma
//                           1 = R-S  table model (old version)
//                           2 = R-S  APED-style model
//                           3 = Mekal  calculate
//                           4 = Mekal  interpolate
//                           5 = Meka
//                           6 = APEC

  using namespace CCfits;

  if ( plasmaType == 1 || plasmaType == 3 || 
       plasmaType == 4 || plasmaType == 5 ) {
    // for R-S, Meka or Mekal read the table model file. for the case
    // of mekal calculate or meka this just limits us to using the energies
    // in the file and gives a sensible default temperature array

    string filenm(FunctionUtility::modelDataPath());
    if ( plasmaType == 1 ) {
      filenm += "raysmith.mod";
    } else {
      filenm += "mekal.mod";
    }

    // open the file
    std::unique_ptr<FITS> pInfile((FITS*)0);
    try {
      pInfile.reset(new FITS(filenm, Read, false));

    // Go to the PARAMETERS extension
      ExtHDU& param = pInfile->extension("PARAMETERS");

    // Get the temperature values. The temperature information is in the
    // first row of the PARAMETERS extension
      Column& values = param.column("VALUE");
      values.read(Tarray,1);

    // Go to the ENERGIES extension
      ExtHDU& ener = pInfile->extension("ENERGIES");
      size_t Nenergy (ener.rows());

    // Get the energies
      Column& eLow = ener.column("ENERG_LO");
      Column& eHigh = ener.column("ENERG_HI");

      RealArray enLow(Nenergy), enHigh(Nenergy);
      eLow.read(enLow,1,Nenergy);
      eHigh.read(enHigh,1,Nenergy);

      energyArray.resize(Nenergy+1);
      for (size_t i=0; i<Nenergy; i++) energyArray[i] = enLow[i];
      energyArray[Nenergy] = enHigh[Nenergy-1];
    } catch(...) {
      std::ostringstream msg;
      msg << "plasmaFileInfo: Failed to read " << filenm;
      FunctionUtility::xsWrite(msg.str(),5);
      return 2;
    }

  } else if ( plasmaType == 2 || plasmaType == 6 ) {
    string continuumFile;
    if ( plasmaType == 6 ) {
      // for apec read the Aped continuum file
      string version = FunctionUtility::atomdbVersion();
      string lineFile;
      getApedFileNames(version, continuumFile, lineFile);
    } else {
      continuumFile = FunctionUtility::modelDataPath()+"RS93_coco.fits";
    }

    // open the file
    std::unique_ptr<FITS> pInfile((FITS*)0);
    try {
      pInfile.reset(new FITS(continuumFile, Read, false));

    // Go to the PARAMETERS extension
      ExtHDU& param = pInfile->extension("PARAMETERS");
      int nrows(param.rows());

    // Get the temperature values. The temperature information is in the
    // kT column of the PARAMETERS extension
      Column& values = param.column("kT");
      values.read(Tarray,1,nrows);
    } catch(...) {
      std::ostringstream msg;
      msg << "plasmaFileInfo: Failed to read " << continuumFile;
      FunctionUtility::xsWrite(msg.str(),5);
      return 2;
    }
    // there isn't a single energy array so construct one in 1000 log steps 
    // from 1.05e-2 to 100 keV
    size_t nE(1000);
    Real log10eMin = log10(1.05e-2);
    Real delta = 3.0 - log10eMin;
    energyArray.resize(nE);
    for (size_t i=0; i<nE; i++) {
      energyArray[i] = pow(10.,log10eMin+i*delta/(nE-1.0));
    }
  }

  return 0;
}

// routine to return bolometric luminosities from the input file
int plasmaBolo(const int plasmaType, const IntegerVector& Zarray, 
	       const RealArray& abun, const Real dens, const Real z,
	       const RealArray& Tarray, const int ifl, const bool qtherm,
	       const Real velocity, RealArray& boloArray)
{
  // find the energy ranges in the file. only reread if the plasmaType has changed
  static int savePlasmaType(-1);
  static RealArray energyArray, tabTarray;
  if ( plasmaType != savePlasmaType ) {
    int status = plasmaFileInfo(plasmaType, energyArray, tabTarray);
    if ( status != 0 ) return status;
    savePlasmaType = plasmaType;
  }

  // loop round input temperatures getting the flux for each over the
  // entire energy range
  boloArray.resize(Tarray.size());
  for (size_t it=0; it<Tarray.size(); it++) {
    RealArray fluxArr, fluxErrArr;
    calcMultiTempPlasma(energyArray, plasmaType, Zarray, abun, dens, z, Tarray[it],
			1.0, ifl, qtherm, velocity, fluxArr, fluxErrArr);
    // sum flux multiplied by mean bin energy
    boloArray[it] = 0.0;
    for (size_t ie=0; ie<fluxArr.size(); ie++) boloArray[it] += fluxArr[ie]*0.5*(energyArray[ie]+energyArray[ie+1]);
  }

  return 0;
}

