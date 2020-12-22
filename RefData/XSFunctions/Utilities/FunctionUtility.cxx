//   Read the documentation to learn more about C++ code generator
//   versioning.
//	  %X% %Q% %Z% %W%

// FunctionUtility
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSUtil/Utils/XSstream.h>
#include <XSUtil/Utils/IosHolder.h>
#include <XSUtil/Utils/XSutility.h>
#include "xsTypes.h"
#include <fstream>
#include <sstream>
#include <exception>



// Class FunctionUtility::NoInitializer 

FunctionUtility::NoInitializer::NoInitializer()
  : YellowAlert()
{
}

FunctionUtility::NoInitializer::NoInitializer (const string& diag)
  : YellowAlert(" No Initializer for ")
{
  *IosHolder::errHolder() << diag << '\n';
}


// Class FunctionUtility::InvalidAbundanceFile 

FunctionUtility::InvalidAbundanceFile::InvalidAbundanceFile (const string& diag)
  : YellowAlert(" Invalid Abundance File ")
{
  *IosHolder::errHolder() << diag << '\n';
}


// Class FunctionUtility::FunctionException 

FunctionUtility::FunctionException::FunctionException (const string& diag)
  : YellowAlert(" in function ")
{
  *IosHolder::errHolder() << diag << std::endl;
}


// Class FunctionUtility::Cosmology 

FunctionUtility::Cosmology::Cosmology ()
  : H0(.0), q0(.0), lambda0(.0)
{
}


// Class Utility FunctionUtility 
string FunctionUtility::s_XSECT = "bcmc";
string FunctionUtility::s_ABUND = "angr";
const string FunctionUtility::CROSSSECTFILE = "crosssections.dat";
string FunctionUtility::s_abundanceFile = "abundances.dat";
std::vector<string> FunctionUtility::s_elements;
string FunctionUtility::s_managerPath = "";
const size_t FunctionUtility::s_NELEMS = 30;
string FunctionUtility::s_modelDataPath = "";
string FunctionUtility::s_NOT_A_KEY = "$$NOT$$";
FunctionUtility::Cosmology FunctionUtility::s_COSMO = FunctionUtility::Cosmology();
string FunctionUtility::s_abundPath;
string FunctionUtility::s_atomdbVersion;
string FunctionUtility::s_neiVersion;
bool FunctionUtility::s_abundChanged = false;
int FunctionUtility::s_xwriteChatter = 10;
std::vector<double> FunctionUtility::s_tempsDEM;
std::vector<double> FunctionUtility::s_DEM;
std::map<int,std::map<string,Real> > FunctionUtility::s_XFLT;
std::map<string,double> FunctionUtility::s_valueDataBase;

std::map<string,std::vector< float > > FunctionUtility::s_abundanceVectors;

std::map<string,std::string> FunctionUtility::s_crossSections;

std::map<string,std::string> FunctionUtility::s_abundDoc;

std::map<string,std::string> FunctionUtility::s_modelStringDataBase;


bool FunctionUtility::checkXsect (const string& arg)
{
  try
  {
        const string& table = crossSections(arg);       
	if ( table.size() == 0 ) return false;
  }
  catch (NoInitializer&)
  {
        return false;       
  }                

  return true;
}

bool FunctionUtility::checkAbund (const string& arg)
{
  try
  {
        const std::vector<float>& element = abundanceVectors(arg);       
	if ( element.size() == 0 ) return false;
  }
  catch (NoInitializer&)
  {
        return false;       
  }                

  return true;
}

void FunctionUtility::readInitializers (string& xsects, string& abunds)
{

  using namespace std;

  s_abundPath = s_managerPath;

  string fullAbundPath = s_abundPath + string("/") + s_abundanceFile;
  string fullXsectPath = s_managerPath + string("/") + CROSSSECTFILE;

  static const string NULLCHARS(" \t");
  s_elements.reserve(s_NELEMS);

  ifstream abund(fullAbundPath.c_str());    
  abund.exceptions(ios_base::failbit);

  abunds = "[ ";

  try
  {
        string element;
        abund.ignore(256,':');
        string abundLine("");
        getline(abund,abundLine);
        istringstream abundString(abundLine);
        while ( abundString >> element )
        {
                s_elements.push_back(element);     
        }

        bool readingStrings(false);
        while (abund)
        {
                string rawLine("");
                getline(abund, rawLine);       
                vector<float> abundances;
                abundances.reserve(s_NELEMS);
                int first(rawLine.find_first_not_of(NULLCHARS));
                if (first > 0 )rawLine = rawLine.substr(first);
                if ( rawLine.length() != 0 && first != -1)
                {
                        // ignore lines with '#' comment mark as first significant 
                        // character.
                        if ( rawLine[0] == '#') continue;
                        // ignore lines that do not meet syntax specification 
                        // (no delimiter)
                        int delim = rawLine.find_first_of(":");
                        if (delim < 0 ) continue;
                        string abundKey(rawLine.substr(0,rawLine.substr(0,delim).
                                                                find_last_not_of(NULLCHARS)+1));

                        if (abundKey == "References")
                        {
                                readingStrings = true;
                                continue;       
                        }

                        if ( !readingStrings )
                        {
                                abunds += abundKey + " | ";
                                string abundData(rawLine.substr(delim+1));
                                istringstream abVector(abundData);
                                float elementAbundance(0);
                                while (abVector >> elementAbundance)
                                {
                                        abundances.push_back(elementAbundance);      
                                }
                                abundanceVectors(abundKey,abundances);
                        }
                        else
                        {
                        	string descript(rawLine.substr(delim+1));
                                abundDoc(abundKey,descript);
                        }
                }
        }
  }
  catch ( std::exception& )
  {

        if ( !abund.eof())
        {
                string diag = "Cannot read abundance file " + s_abundanceFile
                        + "\n*** in XSPEC data file path (" 
                        + s_managerPath  + ")... exiting\n" ;                            
                throw RedAlert(diag);
        }
  }               


  abunds = abunds.substr(0,abunds.size()-2);
  abunds += "]";

  ifstream xsect(fullXsectPath.c_str());   
  xsect.exceptions(ios_base::failbit);

  xsects = "[ ";
  try
  {
        string rawLine("");
        while (getline(xsect,rawLine))
        {
                int first(rawLine.find_first_not_of(NULLCHARS));
                if (first > 0 )rawLine = rawLine.substr(first);
                if ( rawLine.length() != 0 && first != -1)
                {
                        // ignore lines with '#' comment mark as first significant 
                        // character.
                        if ( rawLine[0] == '#') continue;
                        // ignore lines that do not meet syntax specification 
                        // (no delimiter)
                        int delim = rawLine.find_first_of(":");
                        if (delim < 0 ) continue;
                        string table(rawLine.substr(0,rawLine.substr(0,delim).
                                                                find_last_not_of(NULLCHARS)+1));
                        xsects += table + " | ";
                        string descript(rawLine.substr(delim+1));
                        crossSections(table,descript);
                }
        }
  }
  catch ( std::exception& )
  {

        if ( !xsect.eof())
        {
                string diag = "Cannot read cross section table file " 
                                + CROSSSECTFILE
                                + "\n*** in XSPEC data file path (" 
                                + s_managerPath  + ")... exiting\n" ;                            
                throw RedAlert(diag);
        }
  }     

  xsects = xsects.substr(0,xsects.size()-2);
  xsects += "]";


}

void FunctionUtility::readNewAbundances (const string& file)
{
  using namespace std;
  const string name("file");

  s_abundChanged = true;
  s_abundanceFile = file;
  s_abundPath = XSutility::getRunPath();

  ifstream newAbundanceFile(file.c_str());

  vector<float> newVector(s_NELEMS,0);
  size_t count(0);
  try
  {
        float element;
        newAbundanceFile.exceptions(ios_base::failbit);
        while (count < s_NELEMS && newAbundanceFile >> element)
        {
                newVector[count] = element;
		++count;
        }    
	ABUND(name);
        abundanceVectors(name,newVector);

  }
  catch ( std::exception& )
  {
        if (!newAbundanceFile.eof())
        {
                string diag(file);
                diag += "\nCannot be read from. Either it does not exist or contains invalid data";
                throw InvalidAbundanceFile(diag);       
        }       

  }     

  if (count < s_NELEMS)
  {
        // write a warning to the logfile if file is short.        
     if (XSstream* xscout = dynamic_cast<XSstream*>(IosHolder::outHolder()))
     {
        xscout->setVerbose(0,20);
        *xscout << "\n*** Warning: new element vector has fewer elements than allowed."
              << "\n*** heavier nuclei set with abundance 0" << std::endl;       
        xscout->setVerbose();
     }
  } 
}

float FunctionUtility::getAbundance (const string& element)
{
  string table = s_ABUND;
  return getAbundance(table, element);        
}

float FunctionUtility::getAbundance (const size_t Z)
{
  string table = s_ABUND;
  return getAbundance(table, Z);
}

float FunctionUtility::getAbundance (const string& table, const string& element)
{
  std::vector<string>::iterator f = std::find(s_elements.begin(),s_elements.end(),element);

  if ( f == s_elements.end() ) 
  {
        *IosHolder::errHolder() << "XSPEC::getAbundance: Invalid element: " << element << " entered, returning 0.\n" ;
        return 0;  
  }     
  else
  {
#ifndef STD_COUNT_DEFECT
        size_t index(std::distance(s_elements.begin(),f));
#else
        size_t index(0);
	std::distance(s_elements.begin(),f,index);
#endif

        return getAbundance(table, index+1);        
  }
}

float FunctionUtility::getAbundance (const string& table, const size_t Z)
{
  const std::vector<float>& abVector = abundanceVectors(table);

  if (Z < 1 || Z > abVector.size() ) {
    *IosHolder::errHolder() << "XSPEC::getAbundance: Invalid element atomic number: " << Z << " entered, returning 0.\n" ;
    return 0;  
  } else {
    return abVector[Z-1];
  }

}

const string& FunctionUtility::elements (size_t index)
{

  return s_elements[index];
}

const string& FunctionUtility::getModelString (const string& key)
{
  std::map<string,string>::iterator f(s_modelStringDataBase.find(key));
  return ( f != s_modelStringDataBase.end() ? f->second : s_NOT_A_KEY );
}

void FunctionUtility::setModelString (const string& key, const string& value)
{
  if (!value.length())
  {
     std::map<string,string>::iterator itString =
                s_modelStringDataBase.find(key);
     if (itString != s_modelStringDataBase.end())
        s_modelStringDataBase.erase(itString);
  }
  else
     s_modelStringDataBase[key] = value;
}

void FunctionUtility::eraseModelStringDataBase ()
{
  s_modelStringDataBase.clear();
}

void FunctionUtility::setFunctionCosmoParams (double H0, double q0, double lambda0)
{
  s_COSMO.H0 = static_cast<float>(H0);
  s_COSMO.q0 = static_cast<float>(q0);
  s_COSMO.lambda0 = static_cast<float>(lambda0);
}

void FunctionUtility::XSECT (const string& value)
{
  *IosHolder::outHolder() << " Cross Section Table set to " << value << ": " 
                  << s_crossSections[value] << std::endl;
  s_XSECT = value;
}

void FunctionUtility::ABUND (const string& value)
{
  *IosHolder::outHolder() << " Solar Abundance Vector set to " << value 
        << ": " <<  abundDoc(value) << std::endl;
  s_ABUND = value;
}

const std::vector< float >& FunctionUtility::abundanceVectors (string table)
{
  std::map<string,std::vector<float> >::iterator f(s_abundanceVectors.find(table));

  if ( f == s_abundanceVectors.end())
  {
        // throw quietly to allow this to be caught and the program control
        // sent on to try to read abundance from a file.
        throw NoInitializer();       
  }

  return f->second;
}

const std::string& FunctionUtility::crossSections (string table)
{
  std::map<string,string>::iterator f(s_crossSections.find(table));

  if ( f == s_crossSections.end())
  {
        string diag = " cross section table named "  + table;
        throw NoInitializer(diag);       
  }

  return f->second;
}

const std::string FunctionUtility::abundDoc (string name)
{
  std::map<string,string>::iterator f(s_abundDoc.find(name));

  if ( f == s_abundDoc.end())
  {
        string doc(" User defined abundance vector / no description specified");
        return doc; 
  }
  else return f->second;
}

int FunctionUtility::getNumberXFLT(int ifl)
{
  if ( s_XFLT.count(ifl) > 0 ) {
    return (int)s_XFLT.find(ifl)->second.size();
  } else {
    return 0;
  }
}

bool FunctionUtility::inXFLT(int ifl, int i)
{
   std::ostringstream keystream;
   keystream << "key" << i;
   return inXFLT(ifl, keystream.str());
}


bool FunctionUtility::inXFLT(int ifl, string skey)
{
  if ( getNumberXFLT(ifl) > 0 ) {
    if ( s_XFLT.find(ifl)->second.count(skey) > 0 ) return true;
  }
  return false;
}

double FunctionUtility::getXFLT(int ifl, int i)
{
   std::ostringstream keystream;
   keystream << "key" << i;
   return getXFLT(ifl, keystream.str());
}


double FunctionUtility::getXFLT(int ifl, string skey)
{
  if ( inXFLT(ifl, skey) ) {
    return (double)s_XFLT.find(ifl)->second.find(skey)->second;
  }
  return BADVAL;
}

void FunctionUtility::loadXFLT(int ifl, const std::map<string, Real>& values)
{
  s_XFLT[ifl] = values;
}

void FunctionUtility::clearXFLT()
{
   s_XFLT.clear();
}

double FunctionUtility::getDbValue(const string keyword)
{
   double value = BADVAL;
   std::map<string,double>::const_iterator itMap = s_valueDataBase.find(keyword);
   if (itMap == s_valueDataBase.end()) {
     *IosHolder::errHolder() << "XSPEC::getDbValue: Keyword " << keyword << " not found\n";
   } else {
     value = itMap->second;
   }
   return value;
}

void FunctionUtility::loadDbValue(const string keyword, const double value)
{
  string loKeyword = XSutility::lowerCase(keyword);

  // need to check whether the keyword already exists and if so erase 
  // because the map operator will not overwrite.
  std::map<string,double>::iterator it = s_valueDataBase.find(loKeyword);
  if ( it != s_valueDataBase.end() ) s_valueDataBase.erase(it);

  s_valueDataBase[loKeyword] = value;
}

void FunctionUtility::clearDb()
{
  s_valueDataBase.clear();
}

string FunctionUtility::getDbKeywords()
{
  string keyList;
  for(std::map<string,double>::iterator it=s_valueDataBase.begin(); it!=s_valueDataBase.end(); ++it) {
    keyList += it->first + " ";
  }
  return keyList;
}

void FunctionUtility::xsWrite(const string output, const int chatterLevel)
{
  int absChat = std::abs(chatterLevel);
  XSstream* xsout = dynamic_cast<XSstream*>(IosHolder::outHolder());
  if (xsout) {   
    if (chatterLevel < 0 ) {
      // write only to log file.
      if ( xsout->logChatterLevel() >= absChat ) {
        xsout->setVerbose(0,absChat); 
        *xsout << output << std::endl; 
        xsout->setVerbose();                
      }
    } else {
      if ( xsout->maxChatter() >= chatterLevel ) *xsout << output << std::endl;
    } 
  } else {
    if (FunctionUtility::xwriteChatter() >= absChat)
      *IosHolder::outHolder() << output << std::endl;
  }
  return;
}

// Additional Declarations
