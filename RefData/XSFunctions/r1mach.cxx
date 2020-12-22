/* ======================================================================
* NIST Guide to Available Math Software.
* Fullsource for module R1MACH from package BLAS.
* Retrieved from NETLIB on Tue Dec 13 08:57:38 1994.
* ====================================================================== */

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "cfortran.h"

// Some slight modifications required for C++ compilation: 
//   extern "C" statement added
//   #include <stdlib.h> added to pick up exit function
//   explicit cast to double required for log10 to compile on Solaris
//   - C.G. 8/06 

extern "C" float r1mach(long i)
{
       switch(i){
         case 1: return FLT_MIN;
         case 2: return FLT_MAX;
         case 3: return FLT_EPSILON/FLT_RADIX;
         case 4: return FLT_EPSILON;
         case 5: return log10((double)FLT_RADIX);
         }

       fprintf(stderr, "invalid argument: r1mach(%ld)\n", i);
       exit(1);
       return 0; /* for compilers that complain of missing return values */
}

// Additional wrapper function to convert ints from xspec model functions,
// to the long required by r1mach.  This was first required to prevent
// runtime errors on 64-bit Linux.  5/07

extern "C" float ir1mach(int i)
{
   return r1mach((long)i);
} 

FCALLSCFUN1(FLOAT,ir1mach,R1MACH,r1mach,INT)
