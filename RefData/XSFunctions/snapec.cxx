#include "Aped.h"
#include <XSFunctions/Utilities/xsFortran.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <fstream>
#include <sstream>


const size_t TOTEL = 30;

extern "C" void 
snapec  (const RealArray& energyArray, 
	 const RealArray& params, 
	 int              spectrumNumber,
	 RealArray&       fluxArray,   
	 RealArray&       fluxErrArray,
	 const string&    initString)
{
  /*
    snapecs version 1.0 11/15/2011 by Esra Bulbul
    XSPEC model subroutine which uses abudances produced by SN models allowing to rescale 
    abundance of elements according to SN type.
    This program calls 'apec'
    Parameters:
    params(0) = kT plasma temperature in keV
    params(1) = Total numb of SNe 
    params(2) = percentage contribution of SN type I Model (SNI/SNII)
    params(3) = Model No of SN Ia in abundance file
    params(4) = Model No of SN Ia in abundance file
    params(5) = redshift 

    modified 12/17/14 by kaa  for new Aped class
  */



  const Real Tinput   = params[0];
  const Real snTot    = params[1]*1e9;   // SN Ia contribution
  const Real snRatio  = params[2];       // SN II contribution
  const Real sn1MoNo  = params[3];       //model no for SN Ia 
  const Real sn2MoNo  = params[4];       //model no for SN II
  const Real Redshift = params[5];

	

  // if first time through read the file and initialize

  static bool qinit(true);
  static vector<vector<Real> > abundance;
  static IntegerVector Zinput(TOTEL);
  static Real AG[TOTEL];
  static Real mass[TOTEL];

  if ( qinit ) {

    // check for existence of file with SN yields

    const string& datadir = FunctionUtility::modelDataPath();
    string fullAbundPath = datadir + string("/") + "snyields.dat";
    ifstream file_exists(fullAbundPath.c_str());
	
    if (!file_exists) {
      string errMsg = "WARNING: Abundance File " + fullAbundPath + " could not be found! ";
      xs_write(const_cast<char*>(errMsg.c_str()),5);
      return;
    }

    // read the file with SN yields
	
    ifstream file;
    file.open(fullAbundPath.c_str());
	
    string line, modelName;
    getline(file,line);  //skip the useless line

    size_t nModels = 0;
    while (true) {	
      file >> modelName;
      if ( modelName.substr(0,1) == "#" ) break;
      nModels += 1;
      vector<Real> tempabund(TOTEL);
      for ( size_t j = 0 ; j < TOTEL; j++) file >> tempabund[j];
      abundance.push_back(tempabund);
      getline(file,line);
      if ( file.eof() ) break;
    }

    file.close();

    /*
      Now we calculate masses of each individual ions in solar units using Anders and Grevesse (1989)
      solar/meteoritic abundances
    
      The total nucleon mass per unit hydrogen atom is m_H = 2.12e-24 g/H-atom
      1 atomic mass unit = 1.66053886e-24 grams
	 
      Then the total ion mass

      mass_i = clusterMass*AG[i]*getAtomicMass(i+1)*amu2gr / m_H [ in the units of Msolar ]
	 
      where AG[i] is the Anders& Grevesse abundance of ion i.
      We assume that the cluster mass is 10^12 Msolar. Clusters with different mass
      needs to be rescaled later see Bulbul et al. (2012)

    */

    const Real amu2gr = 1.66053886e-24;    // conversion factor for atomic mass units to grams
    const Real clusterMass = 1.0e12;       // assumed cluster gas mass in solar units
    const Real mH = 2.27e-24;              // average nucleon mass in grams per H-atom
		
    for (size_t elts = 0; elts < TOTEL; elts++) {
      //AG[elts] = pow(10,(abundance[0][elts]-abundance[0][0])); 
      AG[elts] = pow(10,(abundance[1][elts]-abundance[1][0])); //use Asplund et al. 2009 abundances
      mass[elts] = clusterMass*AG[elts]*getAtomicMass(elts+1)*amu2gr/mH;
      Zinput[elts] = elts+1;
    }

    qinit = false;

  }
	

  // Set the abundance array for the current choice of SNI and SNII models and the SNI/SNII ratio

  RealArray totAbund(TOTEL);

  for (size_t elts = 0; elts < TOTEL; elts++) {

    if (elts < 7 ) {
      totAbund[elts] = AG[elts];
    } else {	
      totAbund[elts] = (snTot/(1.+snRatio))
	*( snRatio * abundance[(int)sn1MoNo+1][elts] + abundance[(int)sn2MoNo+1][elts]) / mass[elts];
    }

  }

  // Calculate the spectrum using the Aped class.
	
  const Real Dem = 1.0;
  bool qtherm = false;      // if true apply thermal broadening
  Real velocity = 0.0;      // gaussian velocity broadening to apply

  int status = calcCEISpectrum(energyArray, Zinput, totAbund, Redshift, Tinput, 
			       Dem, qtherm, velocity, fluxArray, fluxErrArray);

  if ( status != 0 ) {
    ostringstream ostr;
    ostr << "Failure in snapec: could not read APEC input files" << status;
    xs_write(const_cast<char*>(ostr.str().c_str()),5);
  }

  return;
	
}	
