#include <gsl/gsl_math.h>
#include <gsl/gsl_pow_int.h>



#define Pi M_PI
#define Power gsl_pow_int




double born(double *inv)
{
    double s=inv[0], s1=inv[1], s2=inv[2], t1=inv[3], t2=inv[4];
    
    double a=(16*(8 - 22*s2 + 20*Power(s2,2) - 4*Power(s2,3) - 4*Power(s2,4) + 
      2*Power(s2,5) - 10*t1 + 32*s2*t1 - 40*Power(s2,2)*t1 + 
      20*Power(s2,3)*t1 + 2*Power(s2,4)*t1 - 4*Power(s2,5)*t1 - 
      4*Power(t1,2) + 8*s2*Power(t1,2) + 8*Power(s2,2)*Power(t1,2) - 
      22*Power(s2,3)*Power(t1,2) + 8*Power(s2,4)*Power(t1,2) + 
      2*Power(s2,5)*Power(t1,2) + 4*Power(t1,3) - 20*s2*Power(t1,3) + 
      22*Power(s2,2)*Power(t1,3) - 6*Power(s2,4)*Power(t1,3) + 
      4*Power(t1,4) - 2*s2*Power(t1,4) - 8*Power(s2,2)*Power(t1,4) + 
      6*Power(s2,3)*Power(t1,4) - 2*Power(t1,5) + 4*s2*Power(t1,5) - 
      2*Power(s2,2)*Power(t1,5) - 
      Power(s1,6)*(-1 + s2)*(-1 + t1)*(-2*s2 + Power(s2,2) - (-2 + t1)*t1)*
       (-1 + t2) - 10*t2 + 46*s2*t2 - 52*Power(s2,2)*t2 + 
      13*Power(s2,3)*t2 + 8*Power(s2,4)*t2 - 7*Power(s2,5)*t2 + 
      2*Power(s2,6)*t2 + 12*t1*t2 - 86*s2*t1*t2 + 110*Power(s2,2)*t1*t2 - 
      38*Power(s2,3)*t1*t2 + 4*Power(s2,5)*t1*t2 - 2*Power(s2,6)*t1*t2 + 
      22*Power(t1,2)*t2 + 15*s2*Power(t1,2)*t2 - 
      67*Power(s2,2)*Power(t1,2)*t2 + 37*Power(s2,3)*Power(t1,2)*t2 - 
      10*Power(s2,4)*Power(t1,2)*t2 + 3*Power(s2,5)*Power(t1,2)*t2 - 
      38*Power(t1,3)*t2 + 51*s2*Power(t1,3)*t2 - 
      9*Power(s2,2)*Power(t1,3)*t2 - 6*Power(s2,3)*Power(t1,3)*t2 + 
      2*Power(s2,4)*Power(t1,3)*t2 + 12*Power(t1,4)*t2 - 
      21*s2*Power(t1,4)*t2 + 15*Power(s2,2)*Power(t1,4)*t2 - 
      6*Power(s2,3)*Power(t1,4)*t2 + 2*Power(t1,5)*t2 - 
      5*s2*Power(t1,5)*t2 + 3*Power(s2,2)*Power(t1,5)*t2 - 4*Power(t2,2) - 
      34*s2*Power(t2,2) + 85*Power(s2,2)*Power(t2,2) - 
      43*Power(s2,3)*Power(t2,2) - 11*Power(s2,4)*Power(t2,2) + 
      14*Power(s2,5)*Power(t2,2) - 3*Power(s2,6)*Power(t2,2) + 
      22*t1*Power(t2,2) - 13*s2*t1*Power(t2,2) - 
      83*Power(s2,2)*t1*Power(t2,2) + 87*Power(s2,3)*t1*Power(t2,2) - 
      29*Power(s2,4)*t1*Power(t2,2) - 3*Power(s2,5)*t1*Power(t2,2) + 
      3*Power(s2,6)*t1*Power(t2,2) + 26*s2*Power(t1,2)*Power(t2,2) - 
      2*Power(s2,2)*Power(t1,2)*Power(t2,2) - 
      9*Power(s2,3)*Power(t1,2)*Power(t2,2) + 
      18*Power(s2,4)*Power(t1,2)*Power(t2,2) - 
      9*Power(s2,5)*Power(t1,2)*Power(t2,2) + 4*Power(t1,3)*Power(t2,2) - 
      14*s2*Power(t1,3)*Power(t2,2) - 
      6*Power(s2,2)*Power(t1,3)*Power(t2,2) - 
      12*Power(s2,3)*Power(t1,3)*Power(t2,2) + 
      12*Power(s2,4)*Power(t1,3)*Power(t2,2) - 16*Power(t1,4)*Power(t2,2) + 
      20*s2*Power(t1,4)*Power(t2,2) + 
      9*Power(s2,2)*Power(t1,4)*Power(t2,2) - 
      9*Power(s2,3)*Power(t1,4)*Power(t2,2) + 6*Power(t1,5)*Power(t2,2) - 
      9*s2*Power(t1,5)*Power(t2,2) + 
      3*Power(s2,2)*Power(t1,5)*Power(t2,2) + 4*Power(t2,3) + 
      22*s2*Power(t2,3) - 57*Power(s2,2)*Power(t2,3) + 
      9*Power(s2,3)*Power(t2,3) + 22*Power(s2,4)*Power(t2,3) - 
      9*Power(s2,5)*Power(t2,3) + Power(s2,6)*Power(t2,3) - 
      38*t1*Power(t2,3) + 59*s2*t1*Power(t2,3) + 
      103*Power(s2,2)*t1*Power(t2,3) - 122*Power(s2,3)*t1*Power(t2,3) + 
      28*Power(s2,4)*t1*Power(t2,3) + 3*Power(s2,5)*t1*Power(t2,3) - 
      Power(s2,6)*t1*Power(t2,3) + 4*Power(t1,2)*Power(t2,3) - 
      146*s2*Power(t1,2)*Power(t2,3) + 
      100*Power(s2,2)*Power(t1,2)*Power(t2,3) + 
      16*Power(s2,3)*Power(t1,2)*Power(t2,3) - 
      26*Power(s2,4)*Power(t1,2)*Power(t2,3) + 
      4*Power(s2,5)*Power(t1,2)*Power(t2,3) + 56*Power(t1,3)*Power(t2,3) + 
      4*s2*Power(t1,3)*Power(t2,3) - 
      56*Power(s2,2)*Power(t1,3)*Power(t2,3) + 
      34*Power(s2,3)*Power(t1,3)*Power(t2,3) - 
      6*Power(s2,4)*Power(t1,3)*Power(t2,3) - 22*Power(t1,4)*Power(t2,3) + 
      27*s2*Power(t1,4)*Power(t2,3) - 
      18*Power(s2,2)*Power(t1,4)*Power(t2,3) + 
      5*Power(s2,3)*Power(t1,4)*Power(t2,3) - 4*Power(t1,5)*Power(t2,3) + 
      6*s2*Power(t1,5)*Power(t2,3) - 
      2*Power(s2,2)*Power(t1,5)*Power(t2,3) + 4*Power(t2,4) - 
      16*s2*Power(t2,4) + 15*Power(s2,2)*Power(t2,4) + 
      5*Power(s2,3)*Power(t2,4) - 4*Power(s2,4)*Power(t2,4) + 
      12*t1*Power(t2,4) - 3*s2*t1*Power(t2,4) - 
      55*Power(s2,2)*t1*Power(t2,4) + 30*Power(s2,3)*t1*Power(t2,4) - 
      16*Power(t1,2)*Power(t2,4) + 76*s2*Power(t1,2)*Power(t2,4) - 
      27*Power(s2,2)*Power(t1,2)*Power(t2,4) - 
      13*Power(s2,3)*Power(t1,2)*Power(t2,4) + 
      4*Power(s2,4)*Power(t1,2)*Power(t2,4) - 22*Power(t1,3)*Power(t2,4) - 
      15*s2*Power(t1,3)*Power(t2,4) + 
      29*Power(s2,2)*Power(t1,3)*Power(t2,4) - 
      8*Power(s2,3)*Power(t1,3)*Power(t2,4) + 16*Power(t1,4)*Power(t2,4) - 
      16*s2*Power(t1,4)*Power(t2,4) + 
      4*Power(s2,2)*Power(t1,4)*Power(t2,4) - 2*Power(t2,5) + 
      4*s2*Power(t2,5) + Power(s2,2)*Power(t2,5) - 
      4*Power(s2,3)*Power(t2,5) + Power(s2,4)*Power(t2,5) + 
      2*t1*Power(t2,5) - 13*s2*t1*Power(t2,5) + 
      13*Power(s2,2)*t1*Power(t2,5) - Power(s2,3)*t1*Power(t2,5) - 
      Power(s2,4)*t1*Power(t2,5) + 6*Power(t1,2)*Power(t2,5) - 
      3*s2*Power(t1,2)*Power(t2,5) - 
      6*Power(s2,2)*Power(t1,2)*Power(t2,5) + 
      3*Power(s2,3)*Power(t1,2)*Power(t2,5) - 4*Power(t1,3)*Power(t2,5) + 
      6*s2*Power(t1,3)*Power(t2,5) - 
      2*Power(s2,2)*Power(t1,3)*Power(t2,5) - 
      Power(s,6)*(-1 + s1)*(-1 + s2)*(-1 + t1)*(-1 + t2)*(-2 + t1 + t2) - 
      Power(s1,5)*(-1 + s2)*(-1 + t1)*
       (-2*Power(-1 + t2,2) + Power(s2,3)*(1 + t2) - 
         5*t1*(-1 + Power(t2,2)) - 
         Power(s2,2)*(-6 + 4*t1 + 7*t2 + Power(t2,2)) + 
         Power(t1,2)*(-9 + 3*t2 + 4*Power(t2,2)) + 
         s2*(-3 - Power(t1,2)*(-3 + t2) - 4*t2 + 7*Power(t2,2) + 
            t1*(3 + 4*t2 - 3*Power(t2,2)))) + 
      Power(s,5)*(-1 + t1)*(-1 + t2)*
       (-8 + 3*t1 - Power(t1,2) + Power(s2,2)*(3 - 5*t2) + 3*t2 - 
         8*t1*t2 - Power(t2,2) + 
         s2*(3 - 3*t1 + Power(t1,2) + 4*t2 + 8*t1*t2 + Power(t2,2)) + 
         Power(s1,2)*(3 - 5*t1 + s2*(-5 + 5*t1 + 2*t2)) + 
         s1*(3 + 4*t1 + Power(t1,2) - 3*t2 + 8*t1*t2 + Power(t2,2) + 
            Power(s2,2)*(-5 + 2*t1 + 5*t2) - 
            s2*(-6 + 6*t1 + Power(t1,2) + 6*t2 + 8*t1*t2 + Power(t2,2)))) + 
      Power(s1,4)*(-(Power(t1,5)*(-1 + t2)) - 
         Power(s2,5)*(1 + t1)*(-1 + t2) - 2*Power(-1 + t2,2)*(2 + 3*t2) + 
         4*Power(t1,4)*(-1 + Power(t2,2)) + 
         Power(t1,3)*(22 + 28*t2 - 26*Power(t2,2) - 6*Power(t2,3)) + 
         2*t1*(4 - 5*Power(t2,2) + Power(t2,3)) + 
         Power(t1,2)*(-11 - 29*t2 + 18*Power(t2,2) + 12*Power(t2,3)) + 
         Power(s2,4)*(4 + Power(t1,2)*(-1 + t2) - 13*t2 - Power(t2,2) + 
            t1*(-13 + 30*t2 + Power(t2,2))) + 
         Power(s2,3)*(-25 + 33*t2 + 14*Power(t2,2) - 2*Power(t2,3) + 
            t1*(38 - 51*t2 - 3*Power(t2,2)) + 
            Power(t1,2)*(1 - 40*t2 + 3*Power(t2,2))) - 
         Power(s2,2)*(-16 + Power(t1,4)*(-1 + t2) + 7*t2 + 
            11*Power(t2,2) + 8*Power(t2,3) + 
            Power(t1,3)*(-3 - 19*t2 + 4*Power(t2,2)) + 
            t1*(17 + 20*t2 + 5*Power(t2,2) - 20*Power(t2,3)) + 
            Power(t1,2)*(-3 - 79*t2 + 14*Power(t2,2) + 6*Power(t2,3))) + 
         s2*(-7*Power(t1,4)*(-1 + t2) + Power(t1,5)*(-1 + t2) + 
            t1*(-33 + 74*t2 + Power(t2,2) - 22*Power(t2,3)) + 
            Power(t1,2)*(32 - 59*t2 + 17*Power(t2,2) - 6*Power(t2,3)) + 
            Power(t1,3)*(-41 - 15*t2 + 14*Power(t2,2) + 6*Power(t2,3)) + 
            2*(6 - 11*t2 - 3*Power(t2,2) + 8*Power(t2,3)))) - 
      Power(s1,3)*(4 + Power(s2,6)*(-1 + t1)*(-1 + t2) - 20*t2 + 
         22*Power(t2,2) - 6*Power(t2,4) - 
         Power(s2,5)*(-1 + t2)*(-5 + 8*t1 + Power(t1,2) + 4*t2) + 
         Power(t1,5)*(4 + t2 - 3*Power(t2,2)) + 
         Power(t1,4)*(-5 - 30*t2 + 13*Power(t2,2) + 8*Power(t2,3)) - 
         Power(t1,3)*(9 - 122*t2 + 16*Power(t2,2) + 34*Power(t2,3) + 
            5*Power(t2,4)) + t1*
          (-13 + 38*t2 - 37*Power(t2,2) + 6*Power(t2,3) + 6*Power(t2,4)) + 
         Power(t1,2)*(43 - 87*t2 + 9*Power(t2,2) + 12*Power(t2,3) + 
            9*Power(t2,4)) + Power(s2,4)*
          (25 + 2*Power(t1,3) - 38*t2 - Power(t2,2) + 
            Power(t1,2)*(-14 + 3*t2 - 3*Power(t2,2)) + 
            t1*(-33 + 51*t2 + 40*Power(t2,2))) + 
         Power(s2,3)*(-46 + Power(t1,4)*(-5 + t2) + 37*t2 + 
            28*Power(t2,2) + 14*Power(t2,3) - 5*Power(t2,4) + 
            2*Power(t1,3)*(7 + 2*t2 + 3*Power(t2,2)) + 
            Power(t1,2)*(28 - 74*t2 - 48*Power(t2,2) + 6*Power(t2,3)) + 
            t1*(37 - 56*t2 - 74*Power(t2,2) + 4*Power(t2,3) + Power(t2,4))\
) + s2*(29 - 18*t2 - 33*Power(t2,2) + 10*Power(t2,3) + 12*Power(t2,4) + 
            Power(t1,4)*(20 + 19*t2 - 15*Power(t2,2)) + 
            Power(t1,5)*(-7 + 3*Power(t2,2)) + 
            t1*(-101 + 112*t2 + 71*Power(t2,2) - 28*Power(t2,3) - 
               26*Power(t2,4)) + 
            Power(t1,2)*(38 - 133*t2 - 40*Power(t2,2) + 50*Power(t2,3) - 
               3*Power(t2,4)) + 
            Power(t1,3)*(-27 - 84*t2 + 10*Power(t2,2) + 8*Power(t2,3) + 
               5*Power(t2,4))) - 
         Power(s2,2)*(Power(t1,5)*(-3 + t2) + 
            Power(t1,4)*(2 - 2*t2 + 6*Power(t2,2)) + 
            t2*(-23 + 20*t2 + 16*Power(t2,2) + Power(t2,3)) + 
            2*Power(t1,3)*(6 + 5*t2 - 16*Power(t2,2) + 3*Power(t2,3)) + 
            t1*(-71 + 106*t2 - 32*Power(t2,2) + 14*Power(t2,3) - 
               19*Power(t2,4)) + 
            Power(t1,2)*(48 - 244*t2 - 34*Power(t2,2) + 20*Power(t2,3) + 
               6*Power(t2,4)))) + 
      Power(s1,2)*(3*Power(s2,6)*(-1 + t1)*(-1 + t2) + 
         Power(t1,5)*(1 + 13*t2 - 6*Power(t2,2) - 2*Power(t2,3)) - 
         2*Power(-1 + t2,2)*(-10 + 6*Power(t2,2) + Power(t2,3)) + 
         Power(t1,4)*(15 - 55*t2 - 27*Power(t2,2) + 29*Power(t2,3) + 
            4*Power(t2,4)) - Power(t1,3)*
          (57 - 103*t2 - 100*Power(t2,2) + 56*Power(t2,3) + 
            18*Power(t2,4) + 2*Power(t2,5)) + 
         Power(t1,2)*(85 - 83*t2 - 2*Power(t2,2) - 6*Power(t2,3) + 
            9*Power(t2,4) + 3*Power(t2,5)) + 
         t1*(-52 + 110*t2 - 67*Power(t2,2) - 9*Power(t2,3) + 
            15*Power(t2,4) + 3*Power(t2,5)) + 
         Power(s2,5)*(-1 + t2)*
          (9 - 7*t2 - 3*Power(t2,2) + Power(t1,2)*(-8 + 3*t2) + 
            t1*(-3 - 4*t2 + Power(t2,2))) - 
         Power(s2,4)*(-16 + 17*t2 - 3*Power(t2,2) - 3*Power(t2,3) - 
            Power(t2,4) + Power(t1,3)*(8 - 20*t2 + 6*Power(t2,2)) + 
            Power(t1,2)*(11 + 5*t2 + 14*Power(t2,2) + 4*Power(t2,3)) + 
            t1*(7 + 20*t2 - 79*Power(t2,2) - 19*Power(t2,3) + Power(t2,4))\
) + Power(s2,3)*(Power(t1,4)*(1 - 19*t2 + 6*Power(t2,2)) + 
            2*Power(t1,3)*(8 + 7*t2 + 10*Power(t2,2) + 3*Power(t2,3)) + 
            t2*(-71 + 48*t2 + 12*Power(t2,2) + 2*Power(t2,3) - 
               3*Power(t2,4)) + 
            Power(t1,2)*(20 - 32*t2 - 34*Power(t2,2) - 32*Power(t2,3) + 
               6*Power(t2,4)) + 
            t1*(-23 + 106*t2 - 244*Power(t2,2) + 10*Power(t2,3) - 
               2*Power(t2,4) + Power(t2,5))) - 
         Power(s2,2)*(42 - 117*t2 + 66*Power(t2,2) - 3*Power(t2,3) + 
            8*Power(t2,4) - 2*Power(t2,5) + 
            Power(t1,5)*(-2 - 7*t2 + 3*Power(t2,2)) + 
            Power(t1,4)*(8 - 13*t2 + 10*Power(t2,2) + 5*Power(t2,3)) + 
            Power(t1,3)*(-3 + 52*t2 + 42*Power(t2,2) - 36*Power(t2,3) + 
               5*Power(t2,4)) + 
            Power(t1,2)*(66 - 121*t2 - 252*Power(t2,2) + 
               42*Power(t2,3) + 10*Power(t2,4) + 3*Power(t2,5)) - 
            t1*(117 - 112*t2 + 121*Power(t2,2) - 52*Power(t2,3) + 
               13*Power(t2,4) + 7*Power(t2,5))) + 
         s2*(16 + 6*t2 - 13*Power(t2,2) - 29*Power(t2,3) + 
            17*Power(t2,4) + 3*Power(t2,5) + 
            Power(t1,4)*(-4 + 69*t2 + 7*Power(t2,2) - 16*Power(t2,3)) + 
            Power(t1,5)*(-3 - 20*t2 + 9*Power(t2,2) + 2*Power(t2,3)) + 
            Power(t1,2)*(-12 + 58*t2 - 349*Power(t2,2) + 
               132*Power(t2,3) + 19*Power(t2,4)) + 
            Power(t1,3)*(30 - 117*t2 + 24*Power(t2,2) - 18*Power(t2,3) + 
               7*Power(t2,4) + 2*Power(t2,5)) - 
            t1*(51 + 120*t2 - 212*Power(t2,2) + Power(t2,3) + 
               41*Power(t2,4) + 11*Power(t2,5)))) + 
      s1*(Power(s2,6)*(-1 + t1)*(2 - 3*Power(t2,2) + Power(t2,3)) + 
         2*Power(-1 + t2,2)*(-11 - 6*t2 + 3*Power(t2,2) + 2*Power(t2,3)) + 
         Power(t1,5)*(4 - 13*t2 - 3*Power(t2,2) + 6*Power(t2,3)) - 
         Power(t1,4)*(16 + 3*t2 - 76*Power(t2,2) + 15*Power(t2,3) + 
            16*Power(t2,4)) + t1*
          (46 - 86*t2 + 15*Power(t2,2) + 51*Power(t2,3) - 
            21*Power(t2,4) - 5*Power(t2,5)) + 
         Power(t1,3)*(22 + 59*t2 - 146*Power(t2,2) + 4*Power(t2,3) + 
            27*Power(t2,4) + 6*Power(t2,5)) - 
         Power(t1,2)*(34 + 13*t2 - 26*Power(t2,2) + 14*Power(t2,3) - 
            20*Power(t2,4) + 9*Power(t2,5)) - 
         Power(s2,5)*(-1 + t2)*
          (1 + 2*t2 - 12*Power(t2,2) + 4*t1*(2 - t2 + Power(t2,2)) + 
            Power(t1,2)*(-9 - 2*t2 + 4*Power(t2,2))) + 
         Power(s2,3)*(-29 + 101*t2 - 38*Power(t2,2) + 27*Power(t2,3) - 
            20*Power(t2,4) + 7*Power(t2,5) + 
            Power(t1,4)*(-12 + 26*t2 + 3*Power(t2,2) - 5*Power(t2,3)) - 
            2*Power(t1,3)*(5 - 14*t2 + 25*Power(t2,2) + 4*Power(t2,3)) + 
            t1*(18 - 112*t2 + 133*Power(t2,2) + 84*Power(t2,3) - 
               19*Power(t2,4)) + 
            Power(t1,2)*(33 - 71*t2 + 40*Power(t2,2) - 10*Power(t2,3) + 
               15*Power(t2,4) - 3*Power(t2,5))) + 
         Power(s2,4)*(12 - 33*t2 + 32*Power(t2,2) - 41*Power(t2,3) + 
            7*Power(t2,4) - Power(t2,5) + 
            2*Power(t1,3)*(8 - 11*t2 - 3*Power(t2,2) + 3*Power(t2,3)) + 
            Power(t1,2)*(-6 + t2 + 17*Power(t2,2) + 14*Power(t2,3)) + 
            t1*(-22 + 74*t2 - 59*Power(t2,2) - 15*Power(t2,3) - 
               7*Power(t2,4) + Power(t2,5))) + 
         Power(s2,2)*(16 - 51*t2 - 12*Power(t2,2) + 30*Power(t2,3) - 
            4*Power(t2,4) - 3*Power(t2,5) + 
            Power(t1,5)*(3 - 11*t2 + 2*Power(t2,3)) + 
            Power(t1,4)*(17 - 41*t2 + 19*Power(t2,2) + 7*Power(t2,3)) + 
            t1*(6 - 120*t2 + 58*Power(t2,2) - 117*Power(t2,3) + 
               69*Power(t2,4) - 20*Power(t2,5)) + 
            Power(t1,3)*(-29 - t2 + 132*Power(t2,2) - 18*Power(t2,3) - 
               16*Power(t2,4) + 2*Power(t2,5)) + 
            Power(t1,2)*(-13 + 212*t2 - 349*Power(t2,2) + 
               24*Power(t2,3) + 7*Power(t2,4) + 9*Power(t2,5))) + 
         s2*(24 - 58*t2 + 29*Power(t2,2) + Power(t2,3) + 11*Power(t2,4) - 
            7*Power(t2,5) + Power(t1,5)*
             (-7 + 24*t2 + 3*Power(t2,2) - 8*Power(t2,3)) + 
            Power(t1,4)*(11 + 10*t2 - 90*Power(t2,2) + 21*Power(t2,3) + 
               8*Power(t2,4)) + 
            Power(t1,3)*(1 - 32*t2 + 38*Power(t2,2) - 16*Power(t2,3) + 
               21*Power(t2,4) - 8*Power(t2,5)) + 
            Power(t1,2)*(29 - 184*t2 + 308*Power(t2,2) + 38*Power(t2,3) - 
               90*Power(t2,4) + 3*Power(t2,5)) + 
            2*t1*(-29 + 144*t2 - 92*Power(t2,2) - 16*Power(t2,3) + 
               5*Power(t2,4) + 12*Power(t2,5)))) + 
      Power(s,4)*(10 - 17*t1 + 24*Power(t1,2) - 6*Power(t1,3) + 
         Power(t1,4) - 17*t2 + 48*t1*t2 - 48*Power(t1,2)*t2 + 
         18*Power(t1,3)*t2 - Power(t1,4)*t2 + 24*Power(t2,2) - 
         48*t1*Power(t2,2) + 28*Power(t1,2)*Power(t2,2) - 
         10*Power(t1,3)*Power(t2,2) - 6*Power(t2,3) + 18*t1*Power(t2,3) - 
         10*Power(t1,2)*Power(t2,3) + Power(t2,4) - t1*Power(t2,4) + 
         Power(s2,3)*(-1 + t2)*
          (-1 - 3*Power(t1,2) - 10*t2 + 2*t1*(1 + 5*t2)) - 
         Power(s1,3)*(-1 + t1)*
          (1 + 2*Power(s2,2) - 10*t1*(-1 + t2) - 2*t2 + 3*Power(t2,2) + 
            s2*(1 + 10*t1*(-1 + t2) - 6*t2 + Power(t2,2))) + 
         Power(s2,2)*(-7 + 2*t2 + 12*Power(t2,2) + 5*Power(t2,3) + 
            Power(t1,3)*(-1 + 3*t2) - 
            3*Power(t1,2)*(1 - 5*t2 + 6*Power(t2,2)) + 
            t1*(5 - 6*t2 + 6*Power(t2,2) - 5*Power(t2,3))) + 
         s2*(2 + Power(t1,4)*(-1 + t2) - 10*t2 - 12*Power(t2,2) - 
            3*Power(t2,3) - Power(t2,4) + 
            Power(t1,3)*(7 - 21*t2 + 10*Power(t2,2)) + 
            2*Power(t1,2)*(-9 + 12*t2 - 2*Power(t2,2) + 5*Power(t2,3)) + 
            t1*(10 - 22*t2 + 20*Power(t2,2) - 9*Power(t2,3) + Power(t2,4))\
) - Power(s1,2)*(7 + 2*Power(s2,3)*(-1 + t2) + 5*Power(t1,3)*(-1 + t2) - 
            5*t2 + 3*Power(t2,2) + Power(t2,3) + 
            6*Power(t1,2)*(-2 - t2 + 3*Power(t2,2)) - 
            t1*(2 - 6*t2 + 15*Power(t2,2) + 3*Power(t2,3)) + 
            Power(s2,2)*(-6 + 9*Power(t1,2)*(-1 + t2) + 21*t2 - 
               9*Power(t2,2) + t1*(21 - 44*t2 + 9*Power(t2,2))) + 
            s2*(-5*Power(t1,3)*(-1 + t2) + 
               Power(t1,2)*(15 + 9*t2 - 24*Power(t2,2)) + 
               3*(1 - 6*t2 + Power(t2,3)) - 
               t1*(23 - 42*t2 - 10*Power(t2,2) + Power(t2,3)))) + 
         s1*(2 + Power(t1,4)*(-1 + t2) + 10*t2 - 18*Power(t2,2) + 
            7*Power(t2,3) - Power(t2,4) + 
            Power(t1,3)*(-3 - 9*t2 + 10*Power(t2,2)) + 
            2*Power(t1,2)*(-6 + 10*t2 - 2*Power(t2,2) + 5*Power(t2,3)) + 
            t1*(-10 - 22*t2 + 24*Power(t2,2) - 21*Power(t2,3) + 
               Power(t2,4)) - Power(s2,3)*(-1 + t2)*
             (1 + Power(t1,2) - 10*t2 + 2*t1*(-3 + 5*t2)) + 
            Power(s2,2)*(-3 + Power(t1,3)*(-3 + t2) + 23*t2 - 
               15*Power(t2,2) - 5*Power(t2,3) + 
               2*Power(t1,2)*t2*(-5 + 12*t2) + 
               t1*(18 - 42*t2 - 9*Power(t2,2) + 5*Power(t2,3))) + 
            s2*(-8 - Power(t1,4)*(-1 + t2) + 2*t2 - Power(t2,2) + 
               6*Power(t2,3) + Power(t2,4) + 
               Power(t1,3)*(6 + 8*t2 - 10*Power(t2,2)) - 
               Power(t1,2)*(1 - 15*t2 + 32*Power(t2,2) + 10*Power(t2,3)) + 
               t1*(2 + 32*t2 + 15*Power(t2,2) + 8*Power(t2,3) - 
                  Power(t2,4))))) + 
      Power(s,3)*(14 + 19*t1 - 30*Power(t1,2) + 22*Power(t1,3) - 
         2*Power(t1,4) + Power(t1,5) + 19*t2 - 52*t1*t2 + 
         106*Power(t1,2)*t2 - 60*Power(t1,3)*t2 + 12*Power(t1,4)*t2 - 
         Power(t1,5)*t2 - 30*Power(t2,2) + 106*t1*Power(t2,2) - 
         124*Power(t1,2)*Power(t2,2) + 42*Power(t1,3)*Power(t2,2) - 
         6*Power(t1,4)*Power(t2,2) + 22*Power(t2,3) - 60*t1*Power(t2,3) + 
         42*Power(t1,2)*Power(t2,3) - 12*Power(t1,3)*Power(t2,3) - 
         2*Power(t2,4) + 12*t1*Power(t2,4) - 6*Power(t1,2)*Power(t2,4) + 
         Power(t2,5) - t1*Power(t2,5) + 
         Power(s2,4)*(-1 + t2)*
          (1 + t1 + 2*Power(t1,2) + 10*t2 - 10*t1*t2) + 
         Power(s1,4)*(-1 + t1)*
          (1 + 2*s2*(-1 + 5*t1*(-1 + t2) - 3*t2) - 10*t1*(-1 + t2) + t2 + 
            2*Power(t2,2) + Power(s2,2)*(3 + t2)) - 
         2*Power(s2,3)*(3 - 9*t2 + 2*Power(t1,3)*t2 + 13*Power(t2,2) + 
            5*Power(t2,3) + Power(t1,2)*(-2 + 4*t2 - 8*Power(t2,2)) - 
            t1*(9 - 17*t2 + 3*Power(t2,2) + 5*Power(t2,3))) + 
         Power(s2,2)*(32 + 4*Power(t1,4) - 
            10*Power(t1,3)*Power(-1 + t2,2) - 19*t2 + 46*Power(t2,2) + 
            10*Power(t2,3) + 3*Power(t2,4) - 
            2*Power(t1,2)*(-5 + 5*t2 + 8*Power(t2,2) + 10*Power(t2,3)) + 
            t1*(-48 + 81*t2 - 20*Power(t2,2) + 14*Power(t2,3) - 
               3*Power(t2,4))) + 
         s2*(-39 + Power(t1,5)*(-1 + t2) + 3*t2 - 26*Power(t2,2) - 
            6*Power(t2,3) - 3*Power(t2,4) - Power(t2,5) + 
            2*Power(t1,4)*(-1 - 6*t2 + 3*Power(t2,2)) + 
            4*Power(t1,3)*(-1 + 7*t2 - 6*Power(t2,2) + 3*Power(t2,3)) + 
            2*Power(t1,2)*(1 - 23*t2 + 42*Power(t2,2) - 5*Power(t2,3) + 
               3*Power(t2,4)) + 
            t1*(20 - 30*t2 - 56*Power(t2,2) + 24*Power(t2,3) - 
               7*Power(t2,4) + Power(t2,5))) + 
         2*Power(s1,3)*(-3 + 5*Power(t1,3)*(-1 + t2) + 9*t2 + 
            2*Power(t2,2) + 2*Power(s2,3)*(-2 + t1 + t2) + 
            Power(t1,2)*(-13 + 3*t2 + 8*Power(t2,2)) - 
            t1*(-9 + 17*t2 + 4*Power(t2,2) + 2*Power(t2,3)) + 
            2*Power(s2,2)*(1 + 6*t2 - Power(t2,2) + 
               Power(t1,2)*(-5 + 4*t2) + t1*(7 - 17*t2 + Power(t2,2))) + 
            s2*(1 - 5*Power(t1,3)*(-1 + t2) - 13*t2 - 8*Power(t2,2) + 
               2*Power(t2,3) + Power(t1,2)*(17 + t2 - 14*Power(t2,2)) + 
               t1*(-23 + 45*t2 + 8*Power(t2,2)))) + 
         Power(s1,2)*(-3*Power(t1,4)*(-1 + t2) + 
            Power(s2,4)*(3 + t1)*(-1 + t2) + 
            Power(t1,3)*(10 + 14*t2 - 20*Power(t2,2)) - 
            2*Power(t1,2)*(-23 + 10*t2 + 8*Power(t2,2) + 5*Power(t2,3)) + 
            t1*(-19 + 81*t2 - 10*Power(t2,2) + 20*Power(t2,3)) + 
            2*(16 - 24*t2 + 5*Power(t2,2) - 5*Power(t2,3) + 
               2*Power(t2,4)) + 
            4*Power(s2,3)*(1 + Power(t1,2)*(-1 + t2) + 7*t2 - 
               5*Power(t2,2) + t1*(6 - 17*t2 + 4*Power(t2,2))) + 
            2*Power(s2,2)*(8 - 2*Power(t1,3)*(-2 + t2) - 37*t2 + 
               7*Power(t2,2) + 4*Power(t2,3) + 
               Power(t1,2)*(7 + 21*t2 - 28*Power(t2,2)) + 
               t1*(-37 + 82*t2 + 21*Power(t2,2) - 2*Power(t2,3))) + 
            s2*(-27 + 3*Power(t1,4)*(-1 + t2) + 51*t2 + 8*Power(t2,2) + 
               10*Power(t2,3) - 6*Power(t2,4) + 
               2*Power(t1,3)*(-5 - 13*t2 + 14*Power(t2,2)) + 
               Power(t1,2)*(-48 - 30*t2 + 56*Power(t2,2) + 
                  22*Power(t2,3)) + 
               2*t1*(32 - 83*t2 - 26*Power(t2,2) - 10*Power(t2,3) + 
                  Power(t2,4)))) + 
         s1*(-39 + Power(t1,5)*(-1 + t2) + 20*t2 + 2*Power(t2,2) - 
            4*Power(t2,3) - 2*Power(t2,4) - Power(t2,5) + 
            Power(t1,4)*(-3 - 7*t2 + 6*Power(t2,2)) + 
            2*Power(t1,3)*(-3 + 12*t2 - 5*Power(t2,2) + 6*Power(t2,3)) + 
            Power(t1,2)*(-26 - 56*t2 + 84*Power(t2,2) - 24*Power(t2,3) + 
               6*Power(t2,4)) + 
            t1*(3 - 30*t2 - 46*Power(t2,2) + 28*Power(t2,3) - 
               12*Power(t2,4) + Power(t2,5)) + 
            2*Power(s2,4)*(-1 + t2)*(-1 - 5*t2 + t1*(-3 + 5*t2)) + 
            2*Power(s2,3)*(1 + 2*Power(t1,3) - 23*t2 + 17*Power(t2,2) + 
               5*Power(t2,3) - 2*Power(t1,2)*(4 - 4*t2 + 7*Power(t2,2)) + 
               t1*(-13 + 45*t2 + Power(t2,2) - 5*Power(t2,3))) + 
            Power(s2,2)*(-27 + 2*Power(t1,4)*(-3 + t2) + 64*t2 - 
               48*Power(t2,2) - 10*Power(t2,3) - 3*Power(t2,4) + 
               2*Power(t1,3)*(5 - 10*t2 + 11*Power(t2,2)) + 
               4*Power(t1,2)*
                (2 - 13*t2 + 14*Power(t2,2) + 7*Power(t2,3)) + 
               t1*(51 - 166*t2 - 30*Power(t2,2) - 26*Power(t2,3) + 
                  3*Power(t2,4))) + 
            s2*(50 - Power(t1,5)*(-1 + t2) - 42*t2 + 54*Power(t2,2) - 
               24*Power(t2,3) + 9*Power(t2,4) + Power(t2,5) + 
               Power(t1,4)*(9 + 5*t2 - 6*Power(t2,2)) - 
               4*Power(t1,3)*(6 - 7*t2 + 7*Power(t2,2) + 3*Power(t2,3)) + 
               Power(t1,2)*(54 + 28*t2 - 44*Power(t2,2) - 
                  28*Power(t2,3) - 6*Power(t2,4)) + 
               t1*(-42 + 150*t2 + 28*Power(t2,2) + 28*Power(t2,3) + 
                  5*Power(t2,4) - Power(t2,5))))) - 
      Power(s,2)*(-4 - 20*t1 + 17*Power(t1,2) + 3*Power(t1,3) - 
         5*Power(t1,4) - 3*Power(t1,5) - 
         Power(s2,5)*(3 + 5*t1*(-1 + t2) - 5*t2)*(-1 + t2) - 20*t2 - 
         52*t1*t2 + 66*Power(t1,2)*t2 - 60*Power(t1,3)*t2 + 
         20*Power(t1,4)*t2 - 2*Power(t1,5)*t2 + 17*Power(t2,2) + 
         66*t1*Power(t2,2) - 200*Power(t1,2)*Power(t2,2) + 
         130*Power(t1,3)*Power(t2,2) - 22*Power(t1,4)*Power(t2,2) + 
         3*Power(t1,5)*Power(t2,2) + 3*Power(t2,3) - 60*t1*Power(t2,3) + 
         130*Power(t1,2)*Power(t2,3) - 56*Power(t1,3)*Power(t2,3) + 
         5*Power(t1,4)*Power(t2,3) - 5*Power(t2,4) + 20*t1*Power(t2,4) - 
         22*Power(t1,2)*Power(t2,4) + 5*Power(t1,3)*Power(t2,4) - 
         3*Power(t2,5) - 2*t1*Power(t2,5) + 3*Power(t1,2)*Power(t2,5) + 
         Power(s1,5)*(-1 + s2)*(-1 + t1)*
          (3 + 5*t1*(-1 + t2) + (-5 + 2*s2)*t2) + 
         Power(s2,4)*(-11 - 2*Power(t1,3) + 29*t2 - 20*Power(t2,2) - 
            10*Power(t2,3) + Power(t1,2)*(13 - 17*t2 + 10*Power(t2,2)) + 
            2*t1*(5 - 11*t2 + Power(t2,2) + 5*Power(t2,3))) - 
         Power(s2,3)*(-22 + Power(t1,4)*(-5 + t2) + 9*t2 - 
            32*Power(t2,2) - 24*Power(t2,3) - 3*Power(t2,4) + 
            2*Power(t1,3)*(7 - 10*t2 + 3*Power(t2,2)) + 
            2*Power(t1,2)*(6 - 20*t2 + 25*Power(t2,2) + 7*Power(t2,3)) + 
            3*t1*(5 + 8*t2 - 24*Power(t2,2) + 2*Power(t2,3) + Power(t2,4))\
) + s2*(6 + 105*t2 - 75*Power(t2,2) + 35*Power(t2,3) - 5*Power(t2,4) + 
            6*Power(t2,5) + Power(t1,5)*(6 + t2 - 3*Power(t2,2)) + 
            Power(t1,4)*(-3 + 4*t2 + 12*Power(t2,2) - 5*Power(t2,3)) + 
            Power(t1,3)*(25 - 16*t2 - 84*Power(t2,2) + 32*Power(t2,3) - 
               5*Power(t2,4)) + 
            Power(t1,2)*(-67 + 70*t2 + 48*Power(t2,2) - 68*Power(t2,3) + 
               8*Power(t2,4) - 3*Power(t2,5)) - 
            t1*(-57 + 88*t2 - 188*Power(t2,2) + 40*Power(t2,3) - 
               4*Power(t2,4) + Power(t2,5))) + 
         Power(s2,2)*(-20 + Power(t1,5)*(-3 + t2) - 89*t2 + 
            31*Power(t2,2) - 40*Power(t2,3) + Power(t2,4) - 
            3*Power(t2,5) + Power(t1,4)*(-1 - 15*t2 + 6*Power(t2,2)) + 
            Power(t1,3)*(4 + 8*t2 + 8*Power(t2,2) + 8*Power(t2,3)) + 
            Power(t1,2)*(25 - 79*t2 + 98*Power(t2,2) - 4*Power(t2,3) + 
               8*Power(t2,4)) + 
            t1*(-11 + 128*t2 - 271*Power(t2,2) + 72*Power(t2,3) - 
               17*Power(t2,4) + 3*Power(t2,5))) + 
         Power(s1,4)*(-11 + 10*Power(t1,3)*(-1 + t2) + 10*t2 + 
            13*Power(t2,2) - 2*Power(t2,3) + 
            t1*(29 - 22*t2 - 17*Power(t2,2)) + 
            2*Power(t1,2)*(-10 + t2 + 5*Power(t2,2)) + 
            2*Power(s2,3)*(-4 + t1*(3 + t2)) + 
            Power(s2,2)*(t2*(25 + t2) + Power(t1,2)*(-21 + 13*t2) - 
               t1*(-27 + 52*t2 + Power(t2,2))) + 
            s2*(13 - 10*Power(t1,3)*(-1 + t2) - 23*t2 - 20*Power(t2,2) + 
               2*Power(t2,3) + Power(t1,2)*(35 - 3*t2 - 16*Power(t2,2)) + 
               t1*(-58 + 64*t2 + 22*Power(t2,2)))) + 
         Power(s1,3)*(22 - 3*Power(t1,4)*(-1 + t2) - 15*t2 - 
            12*Power(t2,2) - 14*Power(t2,3) + 5*Power(t2,4) + 
            2*Power(s2,4)*(-4 + (3 + t1)*t2) - 
            2*Power(t1,3)*(-12 + 3*t2 + 7*Power(t2,2)) + 
            Power(t1,2)*(32 + 72*t2 - 50*Power(t2,2) - 6*Power(t2,3)) - 
            t1*(9 + 24*t2 - 40*Power(t2,2) - 20*Power(t2,3) + 
               Power(t2,4)) + 
            2*Power(s2,3)*(2 + 21*t2 - 5*Power(t2,2) + 
               Power(t1,2)*(-5 + 3*t2) + t1*(21 - 56*t2 + 3*Power(t2,2))) \
+ 2*Power(s2,2)*(23 + Power(t1,3)*(5 - 3*t2) - 57*t2 - 10*Power(t2,2) + 
               4*Power(t2,3) + 
               Power(t1,2)*(19 + 33*t2 - 24*Power(t2,2)) + 
               t1*(-65 + 119*t2 + 20*Power(t2,2))) + 
            s2*(-52 + 3*Power(t1,4)*(-1 + t2) + 69*t2 + 30*Power(t2,2) + 
               18*Power(t2,3) - 5*Power(t2,4) + 
               2*Power(t1,3)*(-9 - 10*t2 + 15*Power(t2,2)) + 
               2*Power(t1,2)*(-40 - 46*t2 + 27*Power(t2,2) + 
                  9*Power(t2,3)) + 
               t1*(105 - 128*t2 - 62*Power(t2,2) - 28*Power(t2,3) + 
                  Power(t2,4)))) + 
         Power(s1,2)*(-20 + 2*Power(s2,5)*t1*(-1 + t2) + 
            3*Power(t1,5)*(-1 + t2) - 11*t2 + 25*Power(t2,2) + 
            4*Power(t2,3) - Power(t2,4) - 3*Power(t2,5) + 
            Power(t1,4)*(1 - 17*t2 + 8*Power(t2,2)) + 
            Power(t1,3)*(-40 + 72*t2 - 4*Power(t2,2) + 8*Power(t2,3)) + 
            Power(t1,2)*(31 - 271*t2 + 98*Power(t2,2) + 8*Power(t2,3) + 
               6*Power(t2,4)) + 
            t1*(-89 + 128*t2 - 79*Power(t2,2) + 8*Power(t2,3) - 
               15*Power(t2,4) + Power(t2,5)) + 
            Power(s2,4)*(-(Power(t1,2)*(-1 + t2)) + 3*(9 - 7*t2)*t2 + 
               t1*(25 - 52*t2 + 13*Power(t2,2))) + 
            2*Power(s2,3)*(23 + 4*Power(t1,3) - 65*t2 + 19*Power(t2,2) + 
               5*Power(t2,3) - 
               2*Power(t1,2)*(5 - 10*t2 + 12*Power(t2,2)) + 
               t1*(-57 + 119*t2 + 33*Power(t2,2) - 3*Power(t2,3))) + 
            Power(s2,2)*(-122 + 173*t2 - 10*Power(t2,2) + 
               20*Power(t2,3) - 13*Power(t2,4) + 
               Power(t1,4)*(-13 + 5*t2) + 
               Power(t1,3)*(20 - 42*t2 + 42*Power(t2,2)) + 
               2*Power(t1,2)*
                (-5 - 64*t2 + 12*Power(t2,2) + 21*Power(t2,3)) + 
               t1*(173 - 380*t2 - 128*Power(t2,2) - 42*Power(t2,3) + 
                  5*Power(t2,4))) - 
            s2*(-86 + 3*Power(t1,5)*(-1 + t2) + 63*t2 + 46*Power(t2,3) - 
               8*Power(t2,4) - 3*Power(t2,5) + 
               4*Power(t1,4)*(-2 - 5*t2 + 3*Power(t2,2)) + 
               2*Power(t1,3)*
                (2 + 7*t2 + 11*Power(t2,2) + 12*Power(t2,3)) + 
               2*Power(t1,2)*(-16 - 150*t2 + 31*Power(t2,2) + 
                  15*Power(t2,3) + 6*Power(t2,4)) + 
               t1*(5 - 88*t2 - 120*Power(t2,2) - 32*Power(t2,3) - 
                  14*Power(t2,4) + Power(t2,5)))) + 
         s1*(6 + 57*t2 - 67*Power(t2,2) + 25*Power(t2,3) - 3*Power(t2,4) + 
            6*Power(t2,5) - Power(t1,5)*(-6 + t2 + 3*Power(t2,2)) + 
            Power(t1,4)*(-5 + 4*t2 + 8*Power(t2,2) - 5*Power(t2,3)) + 
            Power(t1,3)*(35 - 40*t2 - 68*Power(t2,2) + 32*Power(t2,3) - 
               5*Power(t2,4)) + 
            Power(t1,2)*(-75 + 188*t2 + 48*Power(t2,2) - 84*Power(t2,3) + 
               12*Power(t2,4) - 3*Power(t2,5)) + 
            t1*(105 - 88*t2 + 70*Power(t2,2) - 16*Power(t2,3) + 
               4*Power(t2,4) + Power(t2,5)) + 
            Power(s2,5)*(-1 + t2)*(3 - 5*t2 + t1*(-7 + 5*t2)) + 
            Power(s2,4)*(13 + 2*Power(t1,3) - 58*t2 + 35*Power(t2,2) + 
               10*Power(t2,3) - 
               2*Power(t1,2)*(10 - 11*t2 + 8*Power(t2,2)) - 
               t1*(23 - 64*t2 + 3*Power(t2,2) + 10*Power(t2,3))) + 
            Power(s2,3)*(-52 + Power(t1,4)*(-5 + t2) + 105*t2 - 
               80*Power(t2,2) - 18*Power(t2,3) - 3*Power(t2,4) + 
               2*Power(t1,3)*(9 - 14*t2 + 9*Power(t2,2)) + 
               Power(t1,2)*(30 - 62*t2 + 54*Power(t2,2) + 
                  30*Power(t2,3)) + 
               t1*(69 - 128*t2 - 92*Power(t2,2) - 20*Power(t2,3) + 
                  3*Power(t2,4))) - 
            Power(s2,2)*(-86 + Power(t1,5)*(-3 + t2) + 5*t2 - 
               32*Power(t2,2) + 4*Power(t2,3) - 8*Power(t2,4) - 
               3*Power(t2,5) + 
               2*Power(t1,4)*(-4 - 7*t2 + 6*Power(t2,2)) + 
               2*Power(t1,2)*t2*
                (-60 + 31*t2 + 11*Power(t2,2) + 6*Power(t2,3)) + 
               Power(t1,3)*(46 - 32*t2 + 30*Power(t2,2) + 
                  24*Power(t2,3)) + 
               t1*(63 - 88*t2 - 300*Power(t2,2) + 14*Power(t2,3) - 
                  20*Power(t2,4) + 3*Power(t2,5))) + 
            s2*(-42 - 111*t2 + 81*Power(t2,2) - 25*Power(t2,3) + 
               10*Power(t2,4) - 9*Power(t2,5) + 
               Power(t1,5)*(-9 + 2*t2 + 3*Power(t2,2)) + 
               Power(t1,4)*(10 - 35*t2 + 12*Power(t2,2) + 5*Power(t2,3)) + 
               Power(t1,3)*(-25 + 100*t2 + 24*Power(t2,3) + 
                  5*Power(t2,4)) + 
               t1*(-111 + 132*t2 - 352*Power(t2,2) + 100*Power(t2,3) - 
                  35*Power(t2,4) + 2*Power(t2,5)) + 
               Power(t1,2)*(81 - 352*t2 + 108*Power(t2,2) + 
                  12*Power(t2,4) + 3*Power(t2,5))))) + 
      s*(12 - 20*t1 + 14*Power(t1,2) - 18*Power(t1,3) + 14*Power(t1,4) - 
         2*Power(t1,5) + Power(s1,6)*(-1 + s2)*(-1 + t1)*(-2 + s2 + t1)*
          (-1 + t2) - 20*t2 + 82*t1*t2 - 25*Power(t1,2)*t2 - 
         21*Power(t1,3)*t2 - 3*Power(t1,4)*t2 + 11*Power(t1,5)*t2 + 
         14*Power(t2,2) - 25*t1*Power(t2,2) - 10*Power(t1,2)*Power(t2,2) + 
         94*Power(t1,3)*Power(t2,2) - 46*Power(t1,4)*Power(t2,2) - 
         3*Power(t1,5)*Power(t2,2) - 18*Power(t2,3) - 21*t1*Power(t2,3) + 
         94*Power(t1,2)*Power(t2,3) - 92*Power(t1,3)*Power(t2,3) + 
         27*Power(t1,4)*Power(t2,3) - 2*Power(t1,5)*Power(t2,3) + 
         14*Power(t2,4) - 3*t1*Power(t2,4) - 46*Power(t1,2)*Power(t2,4) + 
         27*Power(t1,3)*Power(t2,4) - 2*Power(t2,5) + 11*t1*Power(t2,5) - 
         3*Power(t1,2)*Power(t2,5) - 2*Power(t1,3)*Power(t2,5) - 
         Power(s2,6)*(-1 + t1)*(2 - 3*t2 + Power(t2,2)) + 
         Power(s1,5)*(-1 + s2)*(-1 + t1)*
          (-3 - 5*Power(t1,2)*(-1 + t2) - 4*t2 + 7*Power(t2,2) + 
            Power(s2,2)*(1 + 3*t2) + t1*(6 + 2*t2 - 4*Power(t2,2)) + 
            s2*(9 + 4*t1*(-2 + t2) - 12*t2 - Power(t2,2))) + 
         Power(s2,5)*(-1 + t2)*
          (3 - 6*t2 - 5*Power(t2,2) + Power(t1,2)*(-7 + 4*t2) + 
            t1*(4 - 2*t2 + 5*Power(t2,2))) + 
         Power(s2,4)*(-6 + 16*t2 - 12*Power(t2,2) + 25*Power(t2,3) + 
            Power(t2,4) + Power(t1,2)*t2*(9 - 13*t2 - 8*Power(t2,2)) - 
            2*Power(t1,3)*(5 - 10*t2 + 3*Power(t2,2)) - 
            t1*(-16 + 65*t2 - 63*Power(t2,2) + 13*Power(t2,3) + 
               Power(t2,4))) + 
         Power(s2,3)*(17 - 91*t2 + 48*Power(t2,2) - 50*Power(t2,3) + 
            7*Power(t2,4) - 3*Power(t2,5) + 
            Power(t1,4)*(6 - 20*t2 + 6*Power(t2,2)) + 
            2*Power(t1,3)*(3 - 11*t2 + 13*Power(t2,2) + Power(t2,3)) + 
            Power(t1,2)*(-19 + 43*t2 - 28*Power(t2,2) + 42*Power(t2,3) - 
               2*Power(t2,4)) + 
            t1*(-10 + 118*t2 - 152*Power(t2,2) + 2*Power(t2,3) - 
               9*Power(t2,4) + 3*Power(t2,5))) - 
         Power(s2,2)*(12 - 97*t2 + 26*Power(t2,2) - 12*Power(t2,3) + 
            8*Power(t2,4) - 9*Power(t2,5) + 
            Power(t1,5)*(1 - 8*t2 + 3*Power(t2,2)) + 
            Power(t1,4)*(7 - 30*t2 + 13*Power(t2,2) + 2*Power(t2,3)) + 
            Power(t1,3)*(-23 + 2*t2 + 78*Power(t2,2) - 6*Power(t2,3) - 
               3*Power(t2,4)) - 
            t1*(6 - 2*t2 + 39*Power(t2,2) + 114*Power(t2,3) - 
               38*Power(t2,4) + Power(t2,5)) + 
            Power(t1,2)*(9 + 119*t2 - 209*Power(t2,2) + 114*Power(t2,3) - 
               27*Power(t2,4) + 6*Power(t2,5))) + 
         s2*(Power(t1,4)*(-13 + t2 + 37*Power(t2,2) - 17*Power(t2,3)) + 
            Power(t1,5)*(3 - 19*t2 + 6*Power(t2,2) + 2*Power(t2,3)) + 
            t1*(14 - 174*t2 + 155*Power(t2,2) - 135*Power(t2,3) + 
               59*Power(t2,4) - 15*Power(t2,5)) + 
            Power(t1,3)*(-1 - 7*t2 + 36*Power(t2,2) + 36*Power(t2,3) - 
               22*Power(t2,4) + 2*Power(t2,5)) - 
            2*(5 + 20*Power(t2,2) - 22*Power(t2,3) + 7*Power(t2,4) + 
               2*Power(t2,5)) + 
            Power(t1,2)*(7 + 151*t2 - 274*Power(t2,2) + 66*Power(t2,3) + 
               5*Power(t2,4) + 9*Power(t2,5))) + 
         Power(s1,4)*(-(Power(t1,4)*(-1 + t2)) + 
            2*Power(s2,4)*(-3 + t1 + t2 + t1*t2) + 
            Power(t1,3)*(25 - 13*t2 - 8*Power(t2,2)) - 
            2*(3 - 8*t2 + 5*Power(t2,3)) - 
            Power(t1,2)*(12 - 63*t2 + 13*Power(t2,2) + 6*Power(t2,3)) + 
            t1*(16 - 65*t2 + 9*Power(t2,2) + 20*Power(t2,3)) + 
            Power(s2,3)*(Power(t1,2)*(-11 + 3*t2) + 
               t1*(33 - 67*t2 - 2*Power(t2,2)) + 
               2*(-3 + 16*t2 + Power(t2,2))) + 
            Power(s2,2)*(39 - 4*Power(t1,3)*(-2 + t2) - 55*t2 - 
               32*Power(t2,2) + 4*Power(t2,3) + 
               Power(t1,2)*(28 + 39*t2 - 19*Power(t2,2)) + 
               t1*(-87 + 104*t2 + 23*Power(t2,2))) + 
            s2*(-21 + Power(t1,4)*(-1 + t2) + 5*t2 + 30*Power(t2,2) + 
               6*Power(t2,3) + Power(t1,3)*(-25 + t2 + 16*Power(t2,2)) + 
               t1*(44 + 10*t2 - 22*Power(t2,2) - 20*Power(t2,3)) + 
               Power(t1,2)*(-21 - 73*t2 + 16*Power(t2,2) + 6*Power(t2,3)))) \
+ Power(s1,3)*(17 + 3*Power(t1,5)*(-1 + t2) + 
            Power(s2,5)*(1 + 3*t1)*(-1 + t2) - 10*t2 - 19*Power(t2,2) + 
            6*Power(t2,3) + 6*Power(t2,4) + 
            Power(t1,4)*(7 - 9*t2 - 2*Power(t2,2)) + 
            2*Power(t1,3)*(-25 + t2 + 21*Power(t2,2) + Power(t2,3)) + 
            t1*(-91 + 118*t2 + 43*Power(t2,2) - 22*Power(t2,3) - 
               20*Power(t2,4)) + 
            2*Power(t1,2)*(24 - 76*t2 - 14*Power(t2,2) + 13*Power(t2,3) + 
               3*Power(t2,4)) + 
            Power(s2,4)*(-6 - 2*Power(t1,2)*(-1 + t2) + 33*t2 - 
               11*Power(t2,2) + t1*(32 - 67*t2 + 3*Power(t2,2))) + 
            4*Power(s2,3)*(17 + Power(t1,3) - 31*t2 - 2*Power(t2,2) + 
               Power(t2,3) + Power(t1,2)*(-2 + 17*t2 - 6*Power(t2,2)) + 
               t1*(-31 + 51*t2 + 17*Power(t2,2))) + 
            2*Power(s2,2)*(-50 + 2*Power(t1,4)*(-2 + t2) + 52*t2 + 
               23*Power(t2,2) + 16*Power(t2,3) - 5*Power(t2,4) + 
               2*Power(t1,3)*(2 - 11*t2 + 6*Power(t2,2)) + 
               Power(t1,2)*(-7 - 92*t2 + Power(t2,2) + 12*Power(t2,3)) + 
               t1*(63 - 82*t2 - 58*Power(t2,2) - 10*Power(t2,3) + 
                  Power(t2,4))) - 
            s2*(3*Power(t1,5)*(-1 + t2) + 
               Power(t1,4)*(7 - 21*t2 + 6*Power(t2,2)) + 
               6*Power(t1,3)*(-9 + t2 + 3*Power(t2,2) + 3*Power(t2,3)) - 
               2*(7 + 6*t2 - 8*Power(t2,2) - 21*Power(t2,3) + 
                  2*Power(t2,4)) + 
               2*Power(t1,2)*(22 - 167*t2 + 15*Power(t2,2) + 
                  9*Power(t2,3) + 3*Power(t2,4)) - 
               2*t1*(38 - 71*t2 + 25*Power(t2,2) + 13*Power(t2,3) + 
                  9*Power(t2,4)))) + 
         Power(s1,2)*(-12 + Power(s2,6)*(-1 + t1)*(-1 + t2) + 6*t2 - 
            9*Power(t2,2) + 23*Power(t2,3) - 7*Power(t2,4) - Power(t2,5) - 
            Power(s2,5)*(-1 + t2)*
             (-8 + 15*t1 + Power(t1,2) + 8*t2 - 4*t1*t2) + 
            Power(t1,5)*(9 + t2 - 6*Power(t2,2)) + 
            Power(t1,4)*(-8 - 38*t2 + 27*Power(t2,2) + 3*Power(t2,3)) - 
            2*Power(t1,3)*(-6 - 57*t2 + 57*Power(t2,2) - 3*Power(t2,3) + 
               Power(t2,4)) - Power(t1,2)*
             (26 - 39*t2 - 209*Power(t2,2) + 78*Power(t2,3) + 
               13*Power(t2,4) + 3*Power(t2,5)) + 
            t1*(97 - 2*t2 - 119*Power(t2,2) - 2*Power(t2,3) + 
               30*Power(t2,4) + 8*Power(t2,5)) + 
            Power(s2,4)*(39 + 4*Power(t1,3) - 87*t2 + 28*Power(t2,2) + 
               8*Power(t2,3) + 
               Power(t1,2)*(-32 + 23*t2 - 19*Power(t2,2)) + 
               t1*(-55 + 104*t2 + 39*Power(t2,2) - 4*Power(t2,3))) + 
            2*Power(s2,3)*(-50 + Power(t1,4)*(-5 + t2) + 63*t2 - 
               7*Power(t2,2) + 4*Power(t2,3) - 4*Power(t2,4) + 
               2*Power(t1,3)*(8 - 5*t2 + 6*Power(t2,2)) + 
               Power(t1,2)*(23 - 58*t2 + Power(t2,2) + 12*Power(t2,3)) + 
               2*t1*(26 - 41*t2 - 46*Power(t2,2) - 11*Power(t2,3) + 
                  Power(t2,4))) - 
            2*Power(s2,2)*(-31 + Power(t1,5)*(-3 + t2) - 24*t2 + 
               33*Power(t2,2) + 20*Power(t2,3) - Power(t2,4) - 
               3*Power(t2,5) + Power(t1,4)*(-1 - 8*t2 + 9*Power(t2,2)) + 
               2*Power(t1,3)*
                (10 + 5*t2 - 4*Power(t2,2) + 6*Power(t2,3)) + 
               Power(t1,2)*(33 - 200*t2 + 16*Power(t2,2) - 
                  8*Power(t2,3) + 9*Power(t2,4)) + 
               t1*(-24 + 57*t2 - 200*Power(t2,2) + 10*Power(t2,3) - 
                  8*Power(t2,4) + Power(t2,5))) + 
            s2*(34 - 132*t2 + 69*Power(t2,2) + 9*Power(t2,3) + 
               13*Power(t2,4) - 5*Power(t2,5) + 
               Power(t1,5)*(-15 + t2 + 6*Power(t2,2)) + 
               Power(t1,4)*(32 - 4*t2 - 9*Power(t2,2) + 5*Power(t2,3)) + 
               2*Power(t1,3)*(-32 + 19*t2 + 13*Power(t2,2) + 
                  Power(t2,3) + 5*Power(t2,4)) + 
               Power(t1,2)*(157 - 521*t2 - 64*Power(t2,2) + 
                  54*Power(t2,3) + 15*Power(t2,4) + 3*Power(t2,5)) - 
               2*t1*(132 - 153*t2 + 94*Power(t2,2) - 27*Power(t2,3) + 
                  21*Power(t2,4) + 3*Power(t2,5)))) + 
         s1*(-10 + 14*t2 + 7*Power(t2,2) - Power(t2,3) - 13*Power(t2,4) + 
            3*Power(t2,5) + Power(s2,6)*(-1 + t1)*
             (3 - 4*t2 + Power(t2,2)) + 
            Power(t1,4)*(-14 + 59*t2 + 5*Power(t2,2) - 22*Power(t2,3)) + 
            Power(t1,5)*(-4 - 15*t2 + 9*Power(t2,2) + 2*Power(t2,3)) + 
            t1*t2*(-174 + 151*t2 - 7*Power(t2,2) + Power(t2,3) - 
               19*Power(t2,4)) + 
            Power(t1,3)*(44 - 135*t2 + 66*Power(t2,2) + 36*Power(t2,3) - 
               17*Power(t2,4) + 2*Power(t2,5)) + 
            Power(t1,2)*(-40 + 155*t2 - 274*Power(t2,2) + 36*Power(t2,3) + 
               37*Power(t2,4) + 6*Power(t2,5)) - 
            Power(s2,5)*(-1 + t2)*
             (12 + 4*Power(t1,2)*(-2 + t2) - 14*t2 - 5*Power(t2,2) + 
               t1*(-8 + 2*t2 + 5*Power(t2,2))) + 
            Power(s2,4)*(-21 + 44*t2 - 21*Power(t2,2) - 25*Power(t2,3) - 
               Power(t2,4) + Power(t1,3)*(6 - 20*t2 + 6*Power(t2,2)) + 
               2*Power(t1,2)*
                (15 - 11*t2 + 8*Power(t2,2) + 8*Power(t2,3)) + 
               t1*(5 + 10*t2 - 73*Power(t2,2) + Power(t2,3) + Power(t2,4))) \
- Power(s2,3)*(-14 - 76*t2 + 44*Power(t2,2) - 54*Power(t2,3) + 
               7*Power(t2,4) - 3*Power(t2,5) + 
               2*Power(t1,4)*(-2 - 9*t2 + 3*Power(t2,2)) + 
               2*Power(t1,3)*
                (21 - 13*t2 + 9*Power(t2,2) + 9*Power(t2,3)) + 
               2*Power(t1,2)*(8 - 25*t2 + 15*Power(t2,2) + 
                  9*Power(t2,3) + 3*Power(t2,4)) + 
               t1*(-12 + 142*t2 - 334*Power(t2,2) + 6*Power(t2,3) - 
                  21*Power(t2,4) + 3*Power(t2,5))) + 
            Power(s2,2)*(34 - 264*t2 + 157*Power(t2,2) - 64*Power(t2,3) + 
               32*Power(t2,4) - 15*Power(t2,5) + 
               Power(t1,5)*(-5 - 6*t2 + 3*Power(t2,2)) + 
               Power(t1,4)*(13 - 42*t2 + 15*Power(t2,2) + 
                  10*Power(t2,3)) + 
               Power(t1,3)*(9 + 54*t2 + 54*Power(t2,2) + 2*Power(t2,3) + 
                  5*Power(t2,4)) + 
               t1*(-132 + 306*t2 - 521*Power(t2,2) + 38*Power(t2,3) - 
                  4*Power(t2,4) + Power(t2,5)) + 
               Power(t1,2)*(69 - 188*t2 - 64*Power(t2,2) + 
                  26*Power(t2,3) - 9*Power(t2,4) + 6*Power(t2,5))) + 
            s2*(-34 + 152*t2 - 83*Power(t2,2) + 15*Power(t2,3) - 
               11*Power(t2,4) + 9*Power(t2,5) + 
               Power(t1,5)*(9 + 21*t2 - 12*Power(t2,2) - 2*Power(t2,3)) - 
               Power(t1,4)*(11 + 35*t2 - 10*Power(t2,2) + 4*Power(t2,3)) + 
               Power(t1,2)*(-83 + 25*t2 + 468*Power(t2,2) - 
                  188*Power(t2,3) + 10*Power(t2,4) - 12*Power(t2,5)) + 
               Power(t1,3)*(15 + 59*t2 - 188*Power(t2,2) + 
                  60*Power(t2,3) - 4*Power(t2,4) - 2*Power(t2,5)) + 
               t1*(152 - 22*t2 + 25*Power(t2,2) + 59*Power(t2,3) - 
                  35*Power(t2,4) + 21*Power(t2,5)))))))/
  (Power(-1 + s1,2)*Power(-1 + s2,2)*Power(-1 + t1,2)*Power(s - s2 + t1,2)*
    Power(-1 + t2,2)*Power(s - s1 + t2,2));
  return a/4.; // 1/4 from polarization of 2 photons
};
