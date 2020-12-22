#include	<stdio.h>
#include	<math.h>

void xszbabs(const double* energy, int Nflux, const double* parameter, int spectrum, double* flux, double* fluxError, const char* init)

{

  double wav;
  double hcol;
  double heicol;
  double heiicol;
  double z;
  int ie;

  extern double atten();

hcol    = parameter[0]*1e22;
heicol  = parameter[1]*1e22;
heiicol = parameter[2]*1e22;
z       = parameter[3];

for (ie = 0; ie < Nflux; ie = ie+1) {
        wav = 12.39852/energy[ie];
	wav=wav-z*wav;
        flux[ie]=atten(wav,hcol,heicol,heiicol);
    }

 return;

}
