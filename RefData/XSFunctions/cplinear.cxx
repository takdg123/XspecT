//
// Piecewise linear model for XSPEC12  
// Patrick Broos, Penn State University, 2007
// $Id$

// WARNING!  When all the used "rate" parameters are thawed, the "norm" 
// parameter supplied by XSPEC is degenerate, and thus MUST be frozen.
// We cannot simply freeze one of the "rate" parameters because if that vertex
// needs to go to zero in the fit then norm will be driven towards zero and the
// other "rate" parameters will be driven towards infinity.

// This code was derived by using Xspec/src/XSFunctions/powerLaw.cxx as an example.

#include <XSFunctions/Utilities/funcType.h>
#include <XSFunctions/functionMap.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include "xsTypes.h"
#include <cmath>
#include <iomanip>

extern "C" void 
cplinear  (const RealArray&       energyArray, 
                 const RealArray& params, 
                 int              spectrumNumber,
                 RealArray&       fluxArray, 
                 RealArray&       fluxErrArray,
                 const string&    initString)
{
  using namespace std;

  size_t N(energyArray.size());                        
  fluxArray.resize(N-1);
  fluxErrArray.resize(0);

// Extract the unignored (energy_j > 0) vertices of the line segments (energy_j,rate_j) from the 
// parameter list into vectors vertexEnergy and vertexRate.
  size_t          MaxVertices(params.size()/2);
  size_t          Nvertices(0);
  RealArray       vertexEnergy;
  RealArray       vertexRate;
  Real            thisEnergy;
  Real            thisRate  ;
  size_t          jj;

  vertexEnergy.resize(MaxVertices);
    vertexRate.resize(MaxVertices);

  for (jj = 0; jj < MaxVertices; jj++)
  {
    thisEnergy = params[jj];
    thisRate   = params[jj+MaxVertices];

    if (thisEnergy >= 0)
    {
      // printf("vertex at (%f,%f)\n", thisEnergy, thisRate);
      vertexEnergy[Nvertices] = thisEnergy;
      vertexRate  [Nvertices] = pow(10.0, thisRate)  ;
      Nvertices++;
      // printf("vertexEnergy[0]=%f)\n", vertexEnergy[0]);
    }
  }

  if (Nvertices == 0)
  {
    string msg("ERROR: cplinear model has no positive energy parameters.");
    xs_write(const_cast<char*>(msg.c_str()),10);
    return;
  }

  Real   mid_energy;
  Real   bin_width;
  Real   slope;
  Real   rate;

  jj = 0;
  for (size_t ii = 0; ii < fluxArray.size(); ii++)
  {
    mid_energy = 0.5 * (energyArray[ii+1] + energyArray[ii]);
    bin_width  =       (energyArray[ii+1] - energyArray[ii]);


    if (mid_energy < vertexEnergy[0])
    {
      // Flat extrapolation from the first vertex.
      rate = vertexRate[0];
//      printf("(E,R) = (%f,%f) extrapolated\n",mid_energy, rate);
    }
    else if (mid_energy > vertexEnergy[Nvertices-1])
    {
      // Flat extrapolation from the last vertex.
      rate = vertexRate[Nvertices-1];
//      printf("(E,R) = (%f,%f) extrapolated\n",mid_energy, rate);
    }
    else
    {
      // Advance jj until vertexEnergy[jj] < mid_energy < vertexEnergy[jj+1]
      while ( vertexEnergy[jj+1] < mid_energy         && // channel energy not in vertexEnergy[jj]:vertexEnergy[jj+1] 
             (jj+2)              < Nvertices             // respect logical bounds of vertex* arrays
            ) jj++; 

      // Linearly interpolate between two vertices.
      slope = (vertexRate[jj+1]   - vertexRate[jj]) / 
              (vertexEnergy[jj+1] - vertexEnergy[jj]);

      rate  = vertexRate[jj] + (mid_energy - vertexEnergy[jj]) * slope;
    }

    // Force rate to be non-negative.
    if (rate < 0) rate = 0;
    fluxArray[ii] = bin_width * rate;
    // printf("(E,R) = (%f,%f); vertecies %f:%f\n",mid_energy, rate, vertexEnergy[jj], vertexEnergy[jj+1]);
  }               
} 

