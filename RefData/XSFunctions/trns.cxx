#include <cmath>
#include <XSFunctions/trns.h>

// trns.spg  processed by SPAG 4.50J  at 18:37 on  9 Feb 1995
double TRNS(double E,int Idev)
{

//  E : energy (keV)   0.150-10.000
//  idev : sensor ID ( 0 : SIS0, 1: SIS1, 2: GIS2, 3:GIS3 )

//  orginaly coded by Y. Tawara
//          moddified by H.Awaki  ( Apr.19.'93 )

   double work = .0 , tstrns = .0;
//  thickness of myler and alminuum (unit : cm)
   double mylth[] = {0.22e-4 , 0.22e-4 , 0.54e-4 , 0.54e-4}; 
   double alth[] = {0.28e-5 , 0.28e-5 , 0.37e-5 , 0.37e-5};
   double mesh[] = {0.942 , 0.942 , 0.942 , 0.942};


   work = E*1.0e3;
   if ( work > 0.0 ) 
     tstrns = mesh[Idev]
             *std::exp(-mylth[Idev]*1.41e0*MYLER(work)-alth[Idev]
             *2.69e0*ALMNM(work));
   else
     tstrns = 1.0e-20;
   return tstrns;
}


//  myler.spg  processed by SPAG 4.50J  at 18:37 on  9 Feb 1995

double MYLER(double E)
{
   double norm=.0 , alpha=.0;

   if ( E <= 283.84e0 ) 
   {
      norm = 1.8871e+9;
      alpha = -2.3434;
   }
   else if ( E <= 531.7e0 ) 
   {
      norm = 3.5631e+10;
      alpha = -2.4406;
   }
   else if ( E <= 1560.0e0 ) 
   {
      norm = 2.6550e+11;
      alpha = -2.6547;
   }
   else
   {
      norm = 5.2462e+12;
      alpha = -3.05;
   }
   return norm*std::pow(E,alpha);
}


// almnm.spg  processed by SPAG 4.50J  at 18:37 on  9 Feb 1995

double ALMNM(double E)
{
   double norm=.0 , alpha=.0 , e1=243.534;

   // . . > e1 is the cross point of two functions,
   //    norm1*E**alpha1 and norm2*E**alpha2

   if ( E <= e1 )
   {
      norm = 2.2377e+7;
      alpha = -1.1203;
   }
   else if ( E <= 1559.9e0 ) 
   {
      norm = 7.6486e10;
      alpha = -2.601;
   }
   else
   {
      norm = 2.6897e12;
      alpha = -2.7469;
   }
   return norm*std::pow(E,alpha);
}

