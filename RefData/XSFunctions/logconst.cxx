#include <cmath>
#include <xsTypes.h>
#include <functionMap.h>

void logconst(const RealArray& energy, const RealArray& parameter, 
              int spectrum, RealArray& flux, RealArray& fluxError, 
              const string& init)
{
  Real factor = exp(parameter[0]);

  flux.resize(energy.size()-1);
  for(size_t i=0; i != flux.size(); ++i)
    flux[i] = factor;
}
