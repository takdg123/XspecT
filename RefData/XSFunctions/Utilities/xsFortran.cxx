
#include <XSFunctions/Utilities/xsFortran.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Utils/IosHolder.h>
#include <XSUtil/Utils/XSstream.h>
#include <XSUtil/Numerics/CosmologyFunction.h>
#include <XSUtil/Numerics/Gamma.h>
#include <cfortran.h>
#include <cstdlib>
#include <string>
#include <unistd.h>

FCALLSCFUN0(STRING,FGDATD,FGDATD,fgdatd)
FCALLSCFUN0(STRING,FGMODF,FGMODF,fgmodf)
FCALLSCFUN0(STRING,FGXSCT,FGXSCT,fgxsct)
FCALLSCFUN0(STRING,FGSOLR,FGSOLR,fgsolr)
FCALLSCFUN1(FLOAT,FGABND,FGABND,fgabnd,STRING)
FCALLSCFUN1(FLOAT,FGABNZ,FGABNZ,fgabnz,INT)
FCALLSCFUN2(FLOAT,FGTABN,FGTABN,ftgabn,STRING,STRING)
FCALLSCFUN2(FLOAT,FGTABZ,FGTABZ,fgtabz,STRING,INT)
FCALLSCFUN1(STRING,FGELTI,FGELTI,fgelti,INT)
FCALLSCFUN0(INT,FGNELT,FGNELT,fgnelt)
FCALLSCFUN0(STRING,FGABFL,FGABFL,fgabfl)
FCALLSCSUB1(FPABFL,FPABFL,fpabfl,STRING)
FCALLSCFUN0(STRING,FGAPTH,FGAPTH,fgapth)
FCALLSCSUB1(FPAPTH,FPAPTH,fpapth,STRING)
FCALLSCFUN1(STRING,FGMSTR,FGMSTR,fgmstr,STRING)
FCALLSCSUB1(FPDATD,FPDATD,fpdatd,STRING)
FCALLSCSUB2(FPSOLR,FPSOLR,fpsolr,STRING,PINT)
FCALLSCSUB2(FPXSCT,FPXSCT,fpxsct,STRING,PINT)
FCALLSCSUB2(FPMSTR,FPMSTR,fpmstr,STRING,STRING)
FCALLSCSUB3(FPSLFL,FPSLFL,fpslfl,FLOATV,INT,PINT)
FCALLSCSUB2(RFLABD,RFLABD,rflabd,STRING,PINT)
FCALLSCSUB0(FNINIT,FNINIT,fninit)
FCALLSCFUN2(INT,xs_write,XWRITE,xwrite,STRING,INT)
FCALLSCFUN3(INT,xs_read,XCREAD,xcread,STRING,PSTRING,PINT)
FCALLSCFUN0(FLOAT,csmgq0,CSMGQ0,csmgq0)
FCALLSCFUN0(FLOAT,csmgh0,CSHMH0,csmgh0)
FCALLSCFUN0(FLOAT,csmgl0,CSMGL0,csmgl0)
FCALLSCSUB1(csmpq0,CSMPQ0,csmpq0,FLOAT)
FCALLSCSUB1(csmph0,CSMPH0,csmph0,FLOAT)
FCALLSCSUB1(csmpl0,CSMPL0,csmpl0,FLOAT)
FCALLSCSUB3(csmpall,CSMPALL,csmpall,FLOAT,FLOAT,FLOAT)
FCALLSCFUN3(FLOAT,fzsq,FZSQ,fzsq,FLOAT,FLOAT,FLOAT)
FCALLSCFUN1(INT,DGNFLT,DGNFLT,dgnflt,INT)
FCALLSCFUN2(FLOAT,DGFILT,DGFILT,dgfilt,INT,STRING)
FCALLSCFUN2(LOGICAL,DGQFLT,DGQFLT,dgqflt,INT,STRING)
FCALLSCSUB4(DPFILT,DPFILT,dpfilt,INT,INT,STRING,FLOAT)
FCALLSCSUB0(DCLFLT,DCLFLT,dclflt)
FCALLSCFUN1(DOUBLE,GDBVAL,GDBVAL,gdbval,STRING)
FCALLSCSUB2(PDBVAL,PDBVAL,pdbval,STRING,DOUBLE)
FCALLSCSUB0(CDBASE,CDBASE,cdbase)
FCALLSCFUN0(STRING,FGATDV,FGATDV,fgatdv)
FCALLSCSUB1(FPATDV,FPATDV,fpatdv,STRING)
FCALLSCFUN0(INT,FGCHAT,FGCHAT,fgchat)
FCALLSCSUB1(FPCHAT,FPCHAT,fpchat,INT)
FCALLSCSUB2(xs_getChat,XGTCHT,xgtcht,PINT,PINT)
FCALLSCFUN2(INT,xs_getVersion,XGVERS,xgvers,PSTRING,INT)
FCALLSCFUN1(FLOAT,xs_erf,ERF,erf,FLOAT)
FCALLSCFUN1(FLOAT,xs_erfc,ERFC,erfc,FLOAT)
FCALLSCFUN2(FLOAT,gammap,GAMMAP,gammap,FLOAT,FLOAT)
FCALLSCFUN2(FLOAT,gammq,GAMMQ,gammq,FLOAT,FLOAT)

char* FGDATD()
{
        return const_cast<char*>(FunctionUtility::managerPath().c_str());

} 

char* FGMODF()
{
        return const_cast<char*>(FunctionUtility::modelDataPath().c_str());       
} 

char* FGXSCT()
{
        return const_cast<char*>(FunctionUtility::XSECT().c_str());

}

char* FGSOLR()
{
        return const_cast<char*>(FunctionUtility::ABUND().c_str());

}

float FGABND(const char* element)
{
  return FunctionUtility::getAbundance(string(element));      
}

float FGABNZ(const int Z)
{
  return FunctionUtility::getAbundance((size_t)Z);      
}

float FGTABN(const char* table, const char* element)
{
  return FunctionUtility::getAbundance(string(table), string(element));      
}

float FGTABZ(const char* table, const int Z)
{
  return FunctionUtility::getAbundance(string(table), (size_t)Z);      
}

char* FGELTI(const int index)
{
  return const_cast<char*>(FunctionUtility::elements((size_t)index).c_str());      
}

int FGNELT()
{
  return (int)FunctionUtility::NELEMS();      
}

char* FGABFL()
{
  return const_cast<char*>(FunctionUtility::abundanceFile().c_str());      
}

void FPABFL(const char* fname)
{
  FunctionUtility::abundanceFile(string(fname));      
  return;
}

char* FGAPTH()
{
  return const_cast<char*>(FunctionUtility::abundPath().c_str());      
}

void FPAPTH(const char* abunDir)
{
  FunctionUtility::abundPath(string(abunDir));      
  return;
}


char* FGMSTR(const char* dname)
{
        static string value;
        value = FunctionUtility::getModelString(string(dname));
        if (value == FunctionUtility::NOT_A_KEY())
        {
           value.erase();
        }
        return const_cast<char*>(value.c_str());      
}

void FPDATD(const char* dataDir)
{
   FunctionUtility::managerPath(string(dataDir));
}

void FPSOLR(const char* table, int* ierr)
{
   string tableName = string(table);
   tableName = XSutility::lowerCase(tableName);
   if (tableName == "file")
   {
      FunctionUtility::ABUND(tableName);
      *ierr = 0;
   }
   else
   {
      if (FunctionUtility::checkAbund(tableName))
      {
         FunctionUtility::ABUND(tableName);
         *ierr = 0;
      }
      else
      {
         *ierr = 1;
      }
   }
}

void FPXSCT(const char* csection, int* ierr)
{
   string tableName = string(csection);
   tableName = XSutility::lowerCase(tableName);
   if (FunctionUtility::checkXsect(tableName))
   {
      FunctionUtility::XSECT(tableName);
      *ierr = 0;
   }
   else
   {
      *ierr = 1;
   }
}

void FPMSTR(const char* value1, const char* value2)
{
   string svalue1 = XSutility::upperCase(string(value1));
   if (svalue1 == "INITIALIZE")
   {
      FunctionUtility::eraseModelStringDataBase();
   }
   else
   {
      FunctionUtility::setModelString(svalue1, string(value2));
   }
}

void FPSLFL(float rvalue[], int nvalue, int *ierr)
{
   // Load the values of the file solar abundance table.
   *ierr = 0;
   size_t nelem = static_cast<size_t>(nvalue);
   if (nelem != FunctionUtility::NELEMS())
   {
      *ierr = 1;
   } 
   else
   {
      std::vector<float> rvect(nelem);
      for (size_t i=0; i<nelem; ++i)
      {
         rvect[i] = rvalue[i];
      }
      FunctionUtility::abundanceVectors("file",rvect);
   }  
}

void RFLABD(const char* fname, int *ierr)
{
   string fNameStr(fname);
   *ierr = 0;
   try
   {
      FunctionUtility::readNewAbundances(fNameStr);
   }
   catch (...)
   {
      *ierr = 1;
   }
}

void FNINIT()
{
   // get directory paths
   const char* cHeadasLoc = getenv("HEADAS");
   const char* cHomeLoc = getenv("HOME");

   if (!cHeadasLoc) {
     *IosHolder::outHolder() << "\n***Error: HEADAS environment variable not set."
			     << std::endl;
     return;
   }
   string HOME;
   if (cHomeLoc) {
      HOME = string(cHomeLoc);
   }
   const string HEADAS(cHeadasLoc);

   string spectralLoc = HEADAS + "/../spectral/";
   if (access(spectralLoc.c_str(),R_OK)) {
     *IosHolder::outHolder() << "\n***Error: Unable to find spectral directory."
			     << std::endl;
     return;          
   }

   // Initialize the directory for the mode data files. Use the
   // XSPEC_MDATA_DIR environment variable if it is set, otherwise
   // use the standard directory.
   string modelDataLoc;
   const char* cMdataLoc = getenv("XSPEC_MDATA_DIR");
   if (cMdataLoc)
   {
      modelDataLoc = string(cMdataLoc);
   }
   else
   {
      modelDataLoc = spectralLoc + "modelData/";
   }

   string managerLoc = spectralLoc + "manager";

   FunctionUtility::modelDataPath(modelDataLoc);
   FunctionUtility::managerPath(managerLoc);

   // read the default settings
   string defaultSettingsFileName(managerLoc+"/Xspec.init");
   std::map<string,string> defaultSettings = XSutility::readSettingsFile(defaultSettingsFileName);

   // read the user's settings
   string userSettingsFileName(HOME+"/.xspec/Xspec.init");
   std::map<string,string> userSettings = XSutility::readSettingsFile(userSettingsFileName);
   
   // Read in the abundance and crosssections.dat files.  
   string dummy1, dummy2;
   FunctionUtility::readInitializers(dummy1, dummy2);

   int ierr=1;

   std::map<string,string>::iterator it;
   string value;

   // Initialize the solar abundance table in use.
   it = userSettings.find("ABUND");
   if ( it != userSettings.end() ) {
     value = it->second;
   } else {
     it = defaultSettings.find("ABUND");
     if ( it != defaultSettings.end() ) {
       value = it->second;
     } else {
       value = "angr";
     }
   }
   
   FPSOLR(value.c_str(), &ierr);
   if (ierr) {
     *IosHolder::outHolder() << "\n***Error: Failed to set Solar abundance table."<<std::endl;
     return;
   }

   // Initialize the photoelectric cross-sections in use.
   it = userSettings.find("XSECT");
   if ( it != userSettings.end() ) {
     value = it->second;
   } else {
     it = defaultSettings.find("XSECT");
     if ( it != defaultSettings.end() ) {
       value = it->second;
     } else {
       value = "bcmc";
     }
   }

   FPXSCT(value.c_str(), &ierr);
   if (ierr) {
     *IosHolder::outHolder() << "\n***Error: Failed to set photoelectric cross-section"<<std::endl;
     return;
   }

   // Initialize the list of model string parameters
   FPMSTR("INITIALIZE"," ");

   // Initialize the ATOMDB_VERSION and NEI version
   it = userSettings.find("ATOMDB_VERSION");
   if ( it != userSettings.end() ) {
     value = it->second;
   } else {
     it = defaultSettings.find("ATOMDB_VERSION");
     if ( it != defaultSettings.end() ) {
       value = it->second;
     } else {
       value = "3.0.9";
     }
   }
   FunctionUtility::atomdbVersion(value);

   it = userSettings.find("NEI_VERSION");
   if ( it != userSettings.end() ) {
     value = it->second;
   } else {
     it = defaultSettings.find("NEI_VERSION");
     if ( it != defaultSettings.end() ) {
       value = it->second;
     } else {
       value = "3.0.4";
     }
   }
   FunctionUtility::neiVersion(value);
}

// Emulates the xanlib routine FXWRITE.
int xs_write(char *wrtstr,int  idest)
{
  FunctionUtility::xsWrite(string(wrtstr), idest);
  return 0;
}

// Emulates XCREAD, etc, but assumes interactive input.
int xs_read(const char *prompt,char * buffer,int* ierr)
{
  XSstream* xsin = dynamic_cast<XSstream*>(IosHolder::inHolder());
  char  bk[] = "/*";
  if (xsin)
  {      
     char* temp_prompt =  new char [strlen(IosHolder::xsPrompt())+1]; 
     XSstream::setPrompter(*xsin,string(prompt));
     *xsin >> buffer;
     XSstream::setPrompter(*xsin,string(temp_prompt));
     delete [] temp_prompt;
  }
  else
  {
     *IosHolder::inHolder() >> buffer;
  }

  if( !strncmp(buffer, bk, strlen(bk)) )
  {
        *ierr = -1;
        return 0;
  }

  *ierr = 0;

  return 0;
}

float csmgq0()
{
   return FunctionUtility::getq0();
}

float csmgh0()
{
   return FunctionUtility::getH0();
}

float csmgl0()
{
   return FunctionUtility::getlambda0();
}

void csmpq0(const float q0)
{
  FunctionUtility::setq0(q0);
}

void csmph0(const float H0)
{
  FunctionUtility::setH0(H0);
}

void csmpl0(const float lambda0)
{
  FunctionUtility::setlambda0(lambda0);
}

void csmpall(const float H0, const float q0, const float lambda0)
{
  FunctionUtility::setFunctionCosmoParams((Real)H0, (Real)q0, (Real)lambda0);
}

float fzsq(const float z, const float q0, const float lambda)
{
   Numerics::FZSQ val;
   return (float)val(z,q0,lambda);
}

int DGNFLT(int ifl)
{
  return FunctionUtility::getNumberXFLT(ifl);
}

float DGFILT(int ifl, const char* key)
{
  return FunctionUtility::getXFLT(ifl, string(key));
}

bool DGQFLT(int ifl, const char* key)
{
  return FunctionUtility::inXFLT(ifl, string(key));
}

void DPFILT(int ifl, int nfilt, const char* key, const float keyval)
{
  std::map<string, Real> filtmap;
  filtmap.insert(std::pair<string,Real>(string(key),keyval));

  FunctionUtility::loadXFLT(ifl, filtmap);
  return;
}

void DCLFLT()
{
  FunctionUtility::clearXFLT();
  return;
}

double GDBVAL(const char* keyword)
{
  return FunctionUtility::getDbValue(string(keyword));
}

void PDBVAL(const char* keyword, double value)
{
  FunctionUtility::loadDbValue(string(keyword), value);
  return;
}

void CDBASE()
{
  FunctionUtility::clearDb();
  return;
}

char* FGATDV()
{
  return const_cast<char*>(FunctionUtility::atomdbVersion().c_str());
}

void FPATDV(const char* version)
{
  FunctionUtility::atomdbVersion(string(version));
  return;
}

int FGCHAT()
{
   return FunctionUtility::xwriteChatter();
}

void FPCHAT(int chat)
{
   FunctionUtility::xwriteChatter(chat);
}

// Emulates xanlib/xparse's xgtcht, but only when this library is
// linked into xspec's executable.  For external programs this
// simply returns levels = 0.
void xs_getChat(int* cons, int* log)
{
  XSstream* xsout = dynamic_cast<XSstream*>(IosHolder::outHolder());
  if (xsout)
  {
     *cons = xsout->consoleChatterLevel();
     *log = xsout->logChatterLevel();
  }
  else
  {
     *cons = 0;
     *log = 0;
  }   
}

int xs_getVersion(char* buffer, int buffSize)
{
   const string& versStr(XSutility::xs_version());
   const int len = static_cast<int>(versStr.length());
   int isOK = -1;
   if (buffSize > len)
   {
      versStr.copy(buffer, len);
      buffer[len] = 0;
      isOK = 0;
   }
   else if (buffSize > 0)
   {
      versStr.copy(buffer, buffSize-1);
      buffer[buffSize-1] = 0;
   }
   return isOK;
}


float xs_erf(float x)
{
        Numerics::Erf E;
        Real y(x);
        return float(E(y));
}

float xs_erfc(float x)
{
        Numerics::Erfc E;
        Real y(x);
        return float(E(y));
}

float gammap(float a, float x)
{
        Numerics::GammaP G;
        Real y(x);
        Real b(a);
        return float(G(b,y));
}



float gammq(float a, float x)
{
        Numerics::GammaQ G;
        Real y(x);
        Real b(a);
        return float(G(b,y));
}

