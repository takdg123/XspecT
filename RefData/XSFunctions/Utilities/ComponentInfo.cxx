//   Read the documentation to learn more about C++ code generator
//   versioning.
//	  %X% %Q% %Z% %W%

// ComponentInfo
#include <XSFunctions/Utilities/ComponentInfo.h>


// Class ComponentInfo 


ComponentInfo::ComponentInfo()
      : m_name(""),
        m_type("null"),
        m_error(false),
        m_infoString(""),
        m_isStandardXspec(false),
        m_isPythonModel(false),
        m_isSpecDependent(false)
{
}

ComponentInfo::ComponentInfo (const string& componentName, const string& componentType, bool errorFlag, string addString, bool isStandard, Real eMin, Real eMax)
      : m_name(componentName),
        m_type(componentType),
        m_error(errorFlag),
        m_infoString(addString),
        m_isStandardXspec(isStandard),
        m_isPythonModel(false),
        m_isSpecDependent(false),
	m_minEnergy(eMin),
	m_maxEnergy(eMax)
{
}


void ComponentInfo::reset ()
{
  m_name = "";
  m_type = "nul";    
  m_error = false;
  m_infoString = "";
  m_isStandardXspec = false;
  m_isPythonModel = false;
  m_isSpecDependent = false;
  m_minEnergy = 0.0;
  m_maxEnergy = 1.e10;
}
