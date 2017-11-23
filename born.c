#include <gsl/gsl_math.h>
#include <gsl/gsl_pow_int.h>



#define Pi M_PI
#define Power gsl_pow_int




double born(double *inv)
{
    double s=inv[0], s1=inv[1], s2=inv[2], t1=inv[3], t2=inv[4];
    
    double a=(-8*(-20 + 64*s2 - 76*Power(s2,2) + 36*Power(s2,3) - 4*Power(s2,5) + 
      16*t1 - 64*s2*t1 + 108*Power(s2,2)*t1 - 80*Power(s2,3)*t1 + 
      12*Power(s2,4)*t1 + 8*Power(s2,5)*t1 + 20*Power(t1,2) - 
      60*s2*Power(t1,2) + 28*Power(s2,2)*Power(t1,2) + 
      40*Power(s2,3)*Power(t1,2) - 24*Power(s2,4)*Power(t1,2) - 
      4*Power(s2,5)*Power(t1,2) - 4*Power(t1,3) + 48*s2*Power(t1,3) - 
      72*Power(s2,2)*Power(t1,3) + 16*Power(s2,3)*Power(t1,3) + 
      12*Power(s2,4)*Power(t1,3) - 16*Power(t1,4) + 20*s2*Power(t1,4) + 
      8*Power(s2,2)*Power(t1,4) - 12*Power(s2,3)*Power(t1,4) + 
      4*Power(t1,5) - 8*s2*Power(t1,5) + 4*Power(s2,2)*Power(t1,5) + 
      2*Power(s1,6)*(-1 + s2)*(-1 + t1)*
       (-2*s2 + Power(s2,2) - (-2 + t1)*t1)*(-1 + t2) + 22*t2 - 86*s2*t2 + 
      89*Power(s2,2)*t2 - 31*Power(s2,3)*t2 + 3*Power(s2,4)*t2 + 
      5*Power(s2,5)*t2 - 2*Power(s2,6)*t2 - 40*t1*t2 + 201*s2*t1*t2 - 
      207*Power(s2,2)*t1*t2 + 57*Power(s2,3)*t1*t2 - 9*Power(s2,4)*t1*t2 - 
      4*Power(s2,5)*t1*t2 + 2*Power(s2,6)*t1*t2 - 44*Power(t1,2)*t2 - 
      58*s2*Power(t1,2)*t2 + 120*Power(s2,2)*Power(t1,2)*t2 - 
      30*Power(s2,3)*Power(t1,2)*t2 + 13*Power(s2,4)*Power(t1,2)*t2 - 
      Power(s2,5)*Power(t1,2)*t2 + 94*Power(t1,3)*t2 - 
      82*s2*Power(t1,3)*t2 - 5*Power(s2,3)*Power(t1,3)*t2 - 
      7*Power(s2,4)*Power(t1,3)*t2 - 38*Power(t1,4)*t2 + 
      28*s2*Power(t1,4)*t2 + Power(s2,2)*Power(t1,4)*t2 + 
      9*Power(s2,3)*Power(t1,4)*t2 + 6*Power(t1,5)*t2 - 
      3*s2*Power(t1,5)*t2 - 3*Power(s2,2)*Power(t1,5)*t2 + 22*Power(t2,2) + 
      5*s2*Power(t2,2) - 103*Power(s2,2)*Power(t2,2) + 
      71*Power(s2,3)*Power(t2,2) + 9*Power(s2,4)*Power(t2,2) - 
      17*Power(s2,5)*Power(t2,2) + 3*Power(s2,6)*Power(t2,2) - 
      27*t1*Power(t2,2) + 82*s2*t1*Power(t2,2) + 
      65*Power(s2,2)*t1*Power(t2,2) - 137*Power(s2,3)*t1*Power(t2,2) + 
      53*Power(s2,4)*t1*Power(t2,2) + 5*Power(s2,5)*t1*Power(t2,2) - 
      3*Power(s2,6)*t1*Power(t2,2) - 63*Power(t1,2)*Power(t2,2) - 
      45*s2*Power(t1,2)*Power(t2,2) + 
      108*Power(s2,2)*Power(t1,2)*Power(t2,2) - 
      30*Power(s2,3)*Power(t1,2)*Power(t2,2) - 
      34*Power(s2,4)*Power(t1,2)*Power(t2,2) + 
      10*Power(s2,5)*Power(t1,2)*Power(t2,2) + 55*Power(t1,3)*Power(t2,2) - 
      67*s2*Power(t1,3)*Power(t2,2) + 
      4*Power(s2,2)*Power(t1,3)*Power(t2,2) + 
      56*Power(s2,3)*Power(t1,3)*Power(t2,2) - 
      14*Power(s2,4)*Power(t1,3)*Power(t2,2) + 25*Power(t1,4)*Power(t2,2) + 
      10*s2*Power(t1,4)*Power(t2,2) - 
      53*Power(s2,2)*Power(t1,4)*Power(t2,2) + 
      10*Power(s2,3)*Power(t1,4)*Power(t2,2) - 20*Power(t1,5)*Power(t2,2) + 
      23*s2*Power(t1,5)*Power(t2,2) - 
      3*Power(s2,2)*Power(t1,5)*Power(t2,2) - 30*Power(t2,3) + 
      19*s2*Power(t2,3) + 93*Power(s2,2)*Power(t2,3) - 
      71*Power(s2,3)*Power(t2,3) + 2*Power(s2,4)*Power(t2,3) + 
      4*Power(s2,5)*Power(t2,3) - Power(s2,6)*Power(t2,3) + 
      87*t1*Power(t2,3) - 238*s2*t1*Power(t2,3) - 
      43*Power(s2,2)*t1*Power(t2,3) + 170*Power(s2,3)*t1*Power(t2,3) - 
      40*Power(s2,4)*t1*Power(t2,3) - Power(s2,5)*t1*Power(t2,3) + 
      Power(s2,6)*t1*Power(t2,3) + 47*Power(t1,2)*Power(t2,3) + 
      309*s2*Power(t1,2)*Power(t2,3) - 
      283*Power(s2,2)*Power(t1,2)*Power(t2,3) - 
      7*Power(s2,3)*Power(t1,2)*Power(t2,3) + 
      37*Power(s2,4)*Power(t1,2)*Power(t2,3) - 
      7*Power(s2,5)*Power(t1,2)*Power(t2,3) - 197*Power(t1,3)*Power(t2,3) + 
      68*s2*Power(t1,3)*Power(t2,3) + 
      129*Power(s2,2)*Power(t1,3)*Power(t2,3) - 
      79*Power(s2,3)*Power(t1,3)*Power(t2,3) + 
      15*Power(s2,4)*Power(t1,3)*Power(t2,3) + 67*Power(t1,4)*Power(t2,3) - 
      98*s2*Power(t1,4)*Power(t2,3) + 
      60*Power(s2,2)*Power(t1,4)*Power(t2,3) - 
      13*Power(s2,3)*Power(t1,4)*Power(t2,3) + 12*Power(t1,5)*Power(t2,3) - 
      16*s2*Power(t1,5)*Power(t2,3) + 
      4*Power(s2,2)*Power(t1,5)*Power(t2,3) + 6*Power(t2,4) - 
      13*s2*Power(t2,4) - 19*Power(s2,2)*Power(t2,4) + 
      52*Power(s2,3)*Power(t2,4) - 34*Power(s2,4)*Power(t2,4) + 
      6*Power(s2,5)*Power(t2,4) - 29*t1*Power(t2,4) + 
      62*s2*t1*Power(t2,4) - 16*Power(s2,2)*t1*Power(t2,4) - 
      11*Power(s2,3)*t1*Power(t2,4) + 10*Power(s2,4)*t1*Power(t2,4) - 
      2*Power(s2,5)*t1*Power(t2,4) + 21*Power(t1,2)*Power(t2,4) - 
      126*s2*Power(t1,2)*Power(t2,4) + 
      74*Power(s2,2)*Power(t1,2)*Power(t2,4) + 
      9*Power(s2,3)*Power(t1,2)*Power(t2,4) - 
      8*Power(s2,4)*Power(t1,2)*Power(t2,4) + 64*Power(t1,3)*Power(t2,4) + 
      15*s2*Power(t1,3)*Power(t2,4) - 
      75*Power(s2,2)*Power(t1,3)*Power(t2,4) + 
      22*Power(s2,3)*Power(t1,3)*Power(t2,4) - 46*Power(t1,4)*Power(t2,4) + 
      50*s2*Power(t1,4)*Power(t2,4) - 
      12*Power(s2,2)*Power(t1,4)*Power(t2,4) + 11*s2*Power(t2,5) - 
      16*Power(s2,2)*Power(t2,5) + Power(s2,3)*Power(t2,5) - 
      7*t1*Power(t2,5) + 21*s2*t1*Power(t2,5) + 
      3*Power(s2,2)*t1*Power(t2,5) - 5*Power(s2,3)*t1*Power(t2,5) - 
      13*Power(t1,2)*Power(t2,5) - 14*s2*Power(t1,2)*Power(t2,5) + 
      17*Power(s2,2)*Power(t1,2)*Power(t2,5) - 
      2*Power(s2,3)*Power(t1,2)*Power(t2,5) + 14*Power(t1,3)*Power(t2,5) - 
      12*s2*Power(t1,3)*Power(t2,5) + 
      2*Power(s2,2)*Power(t1,3)*Power(t2,5) + 
      2*Power(s,6)*(-1 + s1)*(-1 + s2)*(-1 + t1)*(-1 + t2)*(-1 + t1 + t2) + 
      Power(s1,5)*(-2*Power(-1 + t2,2) - 8*t1*(-1 + Power(t2,2)) + 
         Power(s2,4)*(3 + 9*t1*(-1 + t2) - 9*t2 + 4*Power(t2,2)) + 
         Power(t1,4)*(1 - 5*t2 + 4*Power(t2,2)) - 
         3*Power(t1,3)*(2 - 9*t2 + 7*Power(t2,2)) + 
         Power(t1,2)*(-3 - 32*t2 + 33*Power(t2,2)) - 
         Power(s2,3)*(22 + 13*Power(t1,2)*(-1 + t2) - 33*t2 + 
            7*Power(t2,2) + t1*(-33 + 22*t2 + 7*Power(t2,2))) + 
         Power(s2,2)*(35 + 7*Power(t1,3)*(-1 + t2) - 40*t2 + 
            3*Power(t2,2) + Power(t1,2)*(-28 + 26*t2) + 
            t1*(-32 + 5*t2 + 19*Power(t2,2))) - 
         s2*(4 + 3*Power(t1,4)*(-1 + t2) + 8*t2 - 12*Power(t2,2) + 
            Power(t1,2)*(-60 + 65*t2 - 9*Power(t2,2)) + 
            Power(t1,3)*(9 - 10*t2 + Power(t2,2)) + 
            t1*(34 - 76*t2 + 38*Power(t2,2)))) + 
      Power(s1,4)*(3*Power(t1,5)*(-1 + t2) + 
         2*Power(-1 + t2,2)*(3 + 2*t2) - 
         t1*Power(-1 + t2,2)*(-23 + 5*t2) - 
         4*Power(t1,4)*t2*(-5 + 3*t2 + 2*Power(t2,2)) + 
         Power(s2,5)*(1 + 3*t1*(-1 + t2) - 7*t2 + 4*Power(t2,2)) + 
         Power(t1,3)*(68 - 163*t2 + 43*Power(t2,2) + 38*Power(t2,3)) - 
         Power(t1,2)*(98 - 201*t2 + 46*Power(t2,2) + 47*Power(t2,3)) + 
         Power(s2,4)*(3 - 6*Power(t1,2)*(-1 + t2) + 38*t2 - 
            19*Power(t2,2) - 8*Power(t2,3) + 
            t1*(17 - 6*t2 - 21*Power(t2,2))) + 
         Power(s2,2)*(-1 + 2*Power(t1,4)*(-1 + t2) - 54*t2 + 
            97*Power(t2,2) - 32*Power(t2,3) + 
            Power(t1,3)*(18 - 17*t2 - 15*Power(t2,2)) - 
            t1*(5 - 194*t2 + 158*Power(t2,2) + Power(t2,3)) - 
            Power(t1,2)*(-72 + 97*t2 + 12*Power(t2,2) + 5*Power(t2,3))) + 
         Power(s2,3)*(16 + 4*Power(t1,3)*(-1 + t2) - 44*t2 - 
            5*Power(t2,2) + 11*Power(t2,3) + 
            Power(t1,2)*(-19 + 25*t2 + 20*Power(t2,2)) + 
            t1*(-71 - 29*t2 + 89*Power(t2,2) + 11*Power(t2,3))) + 
         s2*(-3*Power(t1,5)*(-1 + t2) + (-55 + t2)*Power(-1 + t2,2) + 
            2*Power(t1,4)*(-7 + t2 + 6*Power(t2,2)) + 
            Power(t1,2)*(-99 + 57*t2 + 92*Power(t2,2) - 44*Power(t2,3)) + 
            2*Power(t1,3)*(-2 + 34*t2 - 23*Power(t2,2) + 5*Power(t2,3)) + 
            t1*(145 - 243*t2 + 3*Power(t2,2) + 75*Power(t2,3)))) + 
      Power(s1,3)*(Power(s2,6)*(-1 + t1)*(-1 + t2) + 
         Power(t1,5)*(5 - 7*t2)*t2 - 
         2*Power(-1 + t2,2)*(-8 + 8*t2 + Power(t2,2)) + 
         2*Power(s2,5)*(2 + (8 - 3*t1 + Power(t1,2))*t2 - 4*Power(t2,3)) + 
         Power(t1,3)*(-189 + 228*t2 + 180*Power(t2,2) - 148*Power(t2,3) - 
            23*Power(t2,4)) + Power(t1,4)*
          (38 - 92*t2 + 3*Power(t2,2) + 39*Power(t2,3) + 4*Power(t2,4)) + 
         t1*(-83 + 139*t2 - 75*Power(t2,2) + Power(t2,3) + 
            18*Power(t2,4)) + Power(t1,2)*
          (232 - 202*t2 - 171*Power(t2,2) + 118*Power(t2,3) + 
            21*Power(t2,4)) + Power(s2,4)*
          (-44 + Power(t1,3)*(2 - 8*t2) + 
            3*Power(t1,2)*Power(-1 + t2,2) - 52*t2 + 7*Power(t2,2) + 
            59*Power(t2,3) + 4*Power(t2,4) + 
            t1*(5 - 14*t2 + 24*Power(t2,2) + 21*Power(t2,3))) + 
         Power(s2,3)*(53 - 25*t2 + 6*Power(t1,4)*t2 + 50*Power(t2,2) - 
            65*Power(t2,3) + 3*Power(t2,4) + 
            Power(t1,3)*(7 + 18*t2 - 3*Power(t2,2)) - 
            Power(t1,2)*(62 - 31*t2 + 48*Power(t2,2) + 11*Power(t2,3)) + 
            t1*(86 + 124*t2 - 141*Power(t2,2) - 114*Power(t2,3) - 
               5*Power(t2,4))) - 
         Power(s2,2)*(85 - 195*t2 + 132*Power(t2,2) + 15*Power(t2,3) - 
            35*Power(t2,4) + Power(t1,5)*(1 + t2) + 
            Power(t1,4)*(12 + t2 + 7*Power(t2,2)) - 
            Power(t1,3)*(59 - 117*t2 + 67*Power(t2,2) + 13*Power(t2,3)) + 
            Power(t1,2)*(127 - 159*t2 - 147*Power(t2,2) + 
               13*Power(t2,3) - 4*Power(t2,4)) + 
            t1*(-98 + 321*t2 - 83*Power(t2,2) - 185*Power(t2,3) + 
               29*Power(t2,4))) + 
         s2*(85 - 77*t2 - 59*Power(t2,2) + 73*Power(t2,3) - 
            22*Power(t2,4) + Power(t1,5)*(1 - 4*t2 + 7*Power(t2,2)) - 
            Power(t1,4)*(6 - 79*t2 + 36*Power(t2,2) + 15*Power(t2,3)) + 
            t1*(-215 + 61*t2 + 413*Power(t2,2) - 213*Power(t2,3) - 
               42*Power(t2,4)) + 
            Power(t1,3)*(31 - 105*t2 - 36*Power(t2,2) + 31*Power(t2,3) - 
               7*Power(t2,4)) + 
            Power(t1,2)*(104 + 16*t2 - 315*Power(t2,2) + 74*Power(t2,3) + 
               41*Power(t2,4)))) + 
      Power(s1,2)*(3*Power(s2,6)*(-1 + t1 + t2 - t1*t2) + 
         2*Power(-1 + t2,2)*(-31 + t2 + 5*Power(t2,2)) + 
         Power(t1,5)*(17 - 29*t2 + 14*Power(t2,2) + 4*Power(t2,3)) + 
         Power(t1,4)*(-97 + 79*t2 + 128*Power(t2,2) - 80*Power(t2,3) - 
            22*Power(t2,4)) + t1*
          (112 - 223*t2 + 114*Power(t2,2) + 52*Power(t2,3) - 
            46*Power(t2,4) - 9*Power(t2,5)) - 
         Power(t1,2)*(185 + 49*t2 - 249*Power(t2,2) + 4*Power(t2,3) + 
            60*Power(t2,4) + Power(t2,5)) + 
         Power(t1,3)*(207 + 54*t2 - 481*Power(t2,2) + 90*Power(t2,3) + 
            98*Power(t2,4) + 4*Power(t2,5)) + 
         Power(s2,5)*(-8 - 27*t2 - 10*Power(t2,2) + 17*Power(t2,3) + 
            4*Power(t2,4) - 3*Power(t1,2)*(3 - 4*t2 + 3*Power(t2,2)) + 
            t1*(15 + 3*t2 + 5*Power(t2,2) - 5*Power(t2,3))) + 
         Power(s2,4)*(53 + 70*t2 + 7*Power(t2,2) - 101*Power(t2,3) - 
            31*Power(t2,4) + Power(t1,3)*(15 - 20*t2 + 23*Power(t2,2)) + 
            Power(t1,2)*(-32 - 7*t2 + 20*Power(t2,2) + 3*Power(t2,3)) - 
            t1*(22 - 39*t2 + 24*Power(t2,3) + 9*Power(t2,4))) + 
         Power(s2,3)*(-48 + 64*t2 - 59*Power(t2,2) + 93*Power(t2,3) + 
            31*Power(t2,4) - 5*Power(t2,5) + 
            Power(t1,4)*(-15 + 16*t2 - 19*Power(t2,2)) + 
            Power(t1,3)*t2*(15 - 54*t2 + 5*Power(t2,2)) + 
            Power(t1,2)*(126 - 5*t2 - 88*Power(t2,2) + 47*Power(t2,3) + 
               2*Power(t2,4)) - 
            t1*(93 + 258*t2 - 218*Power(t2,2) - 169*Power(t2,3) - 
               43*Power(t2,4) + Power(t2,5))) + 
         s2*(44 - 175*t2 + 188*Power(t2,2) - 82*Power(t2,3) + 
            12*Power(t2,4) + 13*Power(t2,5) - 
            Power(t1,5)*(23 - 34*t2 + 19*Power(t2,2) + 4*Power(t2,3)) + 
            Power(t1,3)*(-14 - 26*t2 + 87*Power(t2,2) + 20*Power(t2,3) - 
               9*Power(t2,4)) + 
            Power(t1,4)*(61 - 71*t2 - 101*Power(t2,2) + 71*Power(t2,3) + 
               6*Power(t2,4)) + 
            t1*(57 + 564*t2 - 709*Power(t2,2) + 23*Power(t2,3) + 
               160*Power(t2,4) + 5*Power(t2,5)) - 
            Power(t1,2)*(117 + 282*t2 - 674*Power(t2,2) + 
               10*Power(t2,3) + 121*Power(t2,4) + 12*Power(t2,5))) + 
         Power(s2,2)*(14 - 105*t2 + 2*Power(t2,2) + 115*Power(t2,3) - 
            64*Power(t2,4) - 12*Power(t2,5) + 
            Power(t1,5)*(6 - 5*t2 + 5*Power(t2,2)) + 
            Power(t1,4)*(43 - 48*t2 + 48*Power(t2,2) + Power(t2,3)) + 
            Power(t1,3)*(-174 + 93*t2 + 185*Power(t2,2) - 
               115*Power(t2,3) - 3*Power(t2,4)) + 
            Power(t1,2)*(163 + 127*t2 - 462*Power(t2,2) + 
               12*Power(t2,3) + 17*Power(t2,4) + Power(t2,5)) + 
            t1*(-34 + 34*t2 + 100*Power(t2,2) - 279*Power(t2,3) - 
               18*Power(t2,4) + 17*Power(t2,5)))) - 
      s1*(2*Power(-1 + t2,2)*(-31 - 14*t2 + 7*Power(t2,2)) + 
         Power(s2,6)*(-1 + t1)*(2 - 3*Power(t2,2) + Power(t2,3)) + 
         Power(t1,5)*(18 - 15*t2 - 13*Power(t2,2) + 16*Power(t2,3)) - 
         t1*Power(-1 + t2,2)*(-72 + 27*t2 + 89*Power(t2,2) + 
            16*Power(t2,3)) - 2*Power(t1,4)*
          (33 + 34*t2 - 98*Power(t2,2) + 7*Power(t2,3) + 28*Power(t2,4)) - 
         2*Power(t1,2)*(10 + 82*t2 - 25*Power(t2,2) - 41*Power(t2,3) + 
            5*Power(t2,4) + 7*Power(t2,5)) + 
         Power(t1,3)*(58 + 306*t2 - 320*Power(t2,2) - 153*Power(t2,3) + 
            123*Power(t2,4) + 18*Power(t2,5)) + 
         Power(s2,4)*(23 + 18*t2 + 56*Power(t2,2) - 80*Power(t2,3) - 
            53*Power(t2,4) + Power(t1,2)*
             (-39 - 38*t2 + 37*Power(t2,2) + 8*Power(t2,3)) + 
            Power(t1,3)*(29 - 35*t2 + 9*Power(t2,2) + 15*Power(t2,3)) + 
            t1*(-13 + 83*t2 - 40*Power(t2,2) + 21*Power(t2,3) - 
               15*Power(t2,4))) + 
         Power(s2,5)*(-7 - 13*t2 - 23*Power(t2,2) + 13*Power(t2,3) + 
            10*Power(t2,4) + Power(t1,2)*
             (-13 + 13*t2 + Power(t2,2) - 7*Power(t2,3)) - 
            2*t1*(-10 + 2*t2 - 5*Power(t2,2) + 3*Power(t2,3) + 
               Power(t2,4))) + 
         Power(s2,3)*(21 + 59*t2 - 46*Power(t2,2) + 32*Power(t2,3) + 
            70*Power(t2,4) - 4*Power(t2,5) - 
            Power(t1,4)*(27 - 31*t2 + 9*Power(t2,2) + 13*Power(t2,3)) + 
            Power(t1,3)*(3 + 96*t2 - 97*Power(t2,2) - 10*Power(t2,3) + 
               6*Power(t2,4)) + 
            t1*(-111 - 190*t2 + 118*Power(t2,2) + 172*Power(t2,3) + 
               43*Power(t2,4) - 6*Power(t2,5)) - 
            Power(t1,2)*(-114 + 56*t2 + 50*Power(t2,2) + 35*Power(t2,3) - 
               27*Power(t2,4) + 2*Power(t2,5))) + 
         Power(s2,2)*(-111 + 59*t2 - 85*Power(t2,2) + 129*Power(t2,3) - 
            40*Power(t2,4) - 28*Power(t2,5) + 
            Power(t1,5)*(9 - 9*t2 + 2*Power(t2,2) + 4*Power(t2,3)) + 
            Power(t1,4)*(45 - 78*t2 + 36*Power(t2,2) + 29*Power(t2,3) - 
               4*Power(t2,4)) + 
            Power(t1,3)*(-160 - 98*t2 + 337*Power(t2,2) - 
               37*Power(t2,3) - 62*Power(t2,4) + 2*Power(t2,5)) + 
            Power(t1,2)*(60 + 527*t2 - 507*Power(t2,2) - 
               97*Power(t2,3) + 47*Power(t2,4) + 18*Power(t2,5)) + 
            t1*(157 - 365*t2 + 205*Power(t2,2) - 202*Power(t2,3) - 
               47*Power(t2,4) + 20*Power(t2,5))) + 
         s2*(Power(t1,5)*(-27 + 24*t2 + 11*Power(t2,2) - 20*Power(t2,3)) + 
            Power(-1 + t2,2)*(138 + 37*t2 + 25*Power(t2,2) + 
               24*Power(t2,3)) + 
            Power(t1,4)*(48 + 99*t2 - 211*Power(t2,2) + 22*Power(t2,3) + 
               40*Power(t2,4)) + 
            Power(t1,3)*(70 - 201*t2 + 33*Power(t2,2) + 65*Power(t2,3) + 
               15*Power(t2,4) - 12*Power(t2,5)) + 
            t1*(-127 + 723*t2 - 345*Power(t2,2) - 289*Power(t2,3) + 
               164*Power(t2,4) + 26*Power(t2,5)) - 
            Power(t1,2)*(102 + 390*t2 - 511*Power(t2,2) - 
               265*Power(t2,3) + 190*Power(t2,4) + 26*Power(t2,5)))) - 
      Power(s,5)*(Power(s1,2)*
          (6 - 29*t2 + 25*Power(t2,2) - 4*Power(t2,3) + 
            t1*(-23 + 46*t2 - 19*Power(t2,2)) + 
            Power(t1,2)*(11 - 15*t2 + 4*Power(t2,2)) + 
            Power(s2,2)*(2 + 4*(-2 + t1)*t2 + 4*Power(t2,2)) + 
            s2*(-8 + 7*Power(t1,2)*(-1 + t2) + 33*t2 - 21*Power(t2,2) + 
               t1*(15 - 34*t2 + 11*Power(t2,2)))) + 
         s1*(6 + 4*Power(t1,3)*(-1 + t2) + 39*t2 - 51*Power(t2,2) + 
            6*Power(t2,3) + Power(t1,2)*(5 - 7*t2 + 4*Power(t2,2)) + 
            t1*(5 - 40*t2 + 27*Power(t2,2) + 2*Power(t2,3)) + 
            2*Power(s2,2)*(-4 + 13*t2 - 9*Power(t2,2) + 
               Power(t1,2)*(-1 + 2*t2) + t1*(3 - 11*t2 + 5*Power(t2,2))) \
+ s2*(2 - 4*Power(t1,3)*(-1 + t2) - 57*t2 + 53*Power(t2,2) + 
               2*Power(t2,3) + 
               Power(t1,2)*(-11 + 19*t2 - 12*Power(t2,2)) + 
               t1*(5 + 30*t2 - 21*Power(t2,2) - 2*Power(t2,3)))) + 
         2*(-8 - 2*Power(t1,3)*(-1 + t2) - t2 + 13*Power(t2,2) - 
            Power(t2,3) + Power(t1,2)*(-10 + 15*t2 - 4*Power(t2,2)) - 
            t1*(-13 + 11*t2 + 4*Power(t2,2) + Power(t2,3)) + 
            Power(s2,2)*(1 - 5*t2 + 7*Power(t2,2) + 
               Power(t1,2)*(-1 + 2*t2) + t1*(1 + t2 - 5*Power(t2,2))) + 
            s2*(7 + 2*Power(t1,3)*(-1 + t2) + 4*t2 - 16*Power(t2,2) - 
               Power(t2,3) + Power(t1,2)*(13 - 21*t2 + 6*Power(t2,2)) + 
               t1*(-18 + 18*t2 + 5*Power(t2,2) + Power(t2,3))))) + 
      Power(s,4)*(-14 + 62*t1 - 79*Power(t1,2) + 40*Power(t1,3) - 
         5*Power(t1,4) - 12*t2 - 36*t1*t2 + 86*Power(t1,2)*t2 - 
         49*Power(t1,3)*t2 + 5*Power(t1,4)*t2 + 35*Power(t2,2) - 
         37*t1*Power(t2,2) - 13*Power(t1,2)*Power(t2,2) + 
         7*Power(t1,3)*Power(t2,2) - 61*Power(t2,3) + 37*t1*Power(t2,3) + 
         6*Power(t1,2)*Power(t2,3) + 8*Power(t2,4) + 
         Power(s2,3)*(-8 - 17*t2 + 31*Power(t2,2) + 
            Power(t1,2)*(-9 + 11*t2) + t1*(19 - 6*t2 - 19*Power(t2,2))) + 
         Power(s2,2)*(52 + Power(t1,3)*(8 - 10*t2) + 15*t2 - 
            103*Power(t2,2) - 20*Power(t2,3) + 
            Power(t1,2)*(21 - 39*t2 + 14*Power(t2,2)) + 
            t1*(-75 + 56*t2 + 45*Power(t2,2) + 12*Power(t2,3))) - 
         s2*(26 + 5*Power(t1,4)*(-1 + t2) + 14*t2 - 89*Power(t2,2) - 
            45*Power(t2,3) + Power(t1,3)*(56 - 75*t2 + 15*Power(t2,2)) + 
            Power(t1,2)*(-91 + 114*t2 - 39*Power(t2,2) + 
               14*Power(t2,3)) + 
            t1*(34 - 66*t2 + 65*Power(t2,2) + 25*Power(t2,3))) + 
         Power(s1,3)*(22 - 78*t2 + 58*Power(t2,2) - 8*Power(t2,3) + 
            t1*(-59 + 111*t2 - 44*Power(t2,2)) + 
            Power(t1,2)*(23 - 35*t2 + 12*Power(t2,2)) + 
            Power(s2,2)*(5 - 23*t2 + 12*Power(t2,2) + t1*(-7 + 15*t2)) + 
            s2*(-21 + 11*Power(t1,2)*(-1 + t2) + 81*t2 - 48*Power(t2,2) + 
               t1*(44 - 82*t2 + 22*Power(t2,2)))) + 
         Power(s1,2)*(-17 + 131*t2 - 81*Power(t2,2) - 31*Power(t2,3) + 
            8*Power(t2,4) - 2*Power(t1,3)*(9 - 13*t2 + 4*Power(t2,2)) + 
            Power(t1,2)*(21 - 62*t2 + 45*Power(t2,2) - 8*Power(t2,3)) + 
            t1*(46 - 97*t2 - 16*Power(t2,2) + 37*Power(t2,3)) + 
            2*Power(s2,3)*(2 - 11*t2 + 6*Power(t2,2) + t1*(-1 + 5*t2)) + 
            Power(s2,2)*(-13 + 4*Power(t1,2)*(-2 + t2) + 134*t2 - 
               91*Power(t2,2) - 8*Power(t2,3) + 
               t1*(31 - 100*t2 + 23*Power(t2,2))) + 
            s2*(18 - 10*Power(t1,3)*(-1 + t2) - 231*t2 + 
               168*Power(t2,2) + 19*Power(t2,3) + 
               Power(t1,2)*(-13 + 50*t2 - 29*Power(t2,2)) - 
               t1*(59 - 179*t2 + 39*Power(t2,2) + 13*Power(t2,3)))) + 
         s1*(-15 - 5*Power(t1,4)*(-1 + t2) - 13*t2 - 32*Power(t2,2) + 
            100*Power(t2,3) - 16*Power(t2,4) + 
            Power(t1,3)*(-14 + 7*t2 + Power(t2,2)) + 
            t1*(7 - 50*t2 + 137*Power(t2,2) - 74*Power(t2,3)) + 
            Power(t1,2)*(-5 + 71*t2 - 64*Power(t2,2) + 2*Power(t2,3)) + 
            Power(s2,3)*(-4 + 55*t2 - 43*Power(t2,2) + 
               Power(t1,2)*(1 + 5*t2) + t1*(-1 - 36*t2 + 19*Power(t2,2))) \
- Power(s2,2)*(52 + 130*t2 + 6*Power(t1,3)*t2 - 162*Power(t2,2) - 
               28*Power(t2,3) + 
               Power(t1,2)*(37 - 63*t2 + 34*Power(t2,2)) + 
               t1*(-75 - 21*t2 + 28*Power(t2,2) + 12*Power(t2,3))) + 
            s2*(69 + 5*Power(t1,4)*(-1 + t2) + 124*t2 - 169*Power(t2,2) - 
               64*Power(t2,3) + 
               3*Power(t1,3)*(10 - 11*t2 + 5*Power(t2,2)) + 
               Power(t1,2)*(5 - 51*t2 + 30*Power(t2,2) + 14*Power(t2,3)) + 
               t1*(-47 - 51*t2 + 2*Power(t2,2) + 38*Power(t2,3))))) + 
      Power(s,3)*(-38 - 7*t1 + 101*Power(t1,2) - 90*Power(t1,3) + 
         27*Power(t1,4) - 3*Power(t1,5) - 47*t2 + 99*t1*t2 - 
         199*Power(t1,2)*t2 + 155*Power(t1,3)*t2 - 33*Power(t1,4)*t2 + 
         3*Power(t1,5)*t2 + 29*Power(t2,2) - 87*t1*Power(t2,2) + 
         118*Power(t1,2)*Power(t2,2) - 76*Power(t1,3)*Power(t2,2) + 
         8*Power(t1,4)*Power(t2,2) + 26*Power(t2,3) - 71*t1*Power(t2,3) + 
         30*Power(t1,2)*Power(t2,3) - 3*Power(t1,3)*Power(t2,3) - 
         44*Power(t2,4) + 36*t1*Power(t2,4) + 2*Power(t1,2)*Power(t2,4) + 
         4*Power(t2,5) + Power(s2,4)*
          (12 - 2*Power(t1,2)*(-2 + t2) + 23*t2 - 29*Power(t2,2) + 
            t1*(-14 - 9*t2 + 17*Power(t2,2))) + 
         Power(s2,3)*(-57 - 41*t2 + 113*Power(t2,2) + 49*Power(t2,3) + 
            Power(t1,3)*(-11 + 7*t2) - 2*Power(t1,2)*(-9 + Power(t2,2)) + 
            t1*(30 + 4*t2 - 61*Power(t2,2) - 25*Power(t2,3))) + 
         Power(s2,2)*(11 - 2*Power(t1,4)*(-2 + t2) + 37*t2 - 
            81*Power(t2,2) - 177*Power(t2,3) - 6*Power(t2,4) + 
            Power(t1,3)*(40 - 37*t2 - 9*Power(t2,2)) + 
            2*Power(t1,2)*(-79 + 71*t2 + 6*Power(t2,2) + 9*Power(t2,3)) + 
            t1*(143 - 134*t2 + 6*Power(t2,2) + 75*Power(t2,3) + 
               2*Power(t2,4))) + 
         s2*(-3*Power(t1,5)*(-1 + t2) + 
            Power(t1,4)*(-35 + 43*t2 - 12*Power(t2,2)) + 
            Power(t1,3)*(69 - 157*t2 + 125*Power(t2,2) - 
               13*Power(t2,3)) + 
            Power(t1,2)*(35 + 115*t2 - 244*Power(t2,2) + 
               16*Power(t2,3) - 6*Power(t2,4)) + 
            2*(33 + 36*t2 - 71*Power(t2,2) + 109*Power(t2,3) + 
               Power(t2,4)) - 
            t1*(150 + 52*t2 - 327*Power(t2,2) + 115*Power(t2,3) + 
               14*Power(t2,4))) + 
         Power(s1,4)*(-38 + 89*t2 - 49*Power(t2,2) + 4*Power(t2,3) + 
            Power(t1,2)*(-23 + 35*t2 - 12*Power(t2,2)) + 
            t1*(71 - 114*t2 + 39*Power(t2,2)) - 
            2*Power(s2,2)*(3 - 12*t2 + 6*Power(t2,2) + 
               2*t1*(-4 + 5*t2)) + 
            s2*(32 - 11*Power(t1,2)*(-1 + t2) - 85*t2 + 41*Power(t2,2) + 
               t1*(-67 + 94*t2 - 19*Power(t2,2)))) - 
         Power(s1,3)*(-68 + 139*t2 + 11*Power(t2,2) - 70*Power(t2,3) + 
            8*Power(t2,4) - 8*Power(t1,3)*(4 - 7*t2 + 3*Power(t2,2)) + 
            Power(t1,2)*(50 - 159*t2 + 111*Power(t2,2) - 
               16*Power(t2,3)) + 
            2*t1*(41 - 12*t2 - 77*Power(t2,2) + 33*Power(t2,3)) + 
            2*Power(s2,3)*(5 - 29*t2 + 16*Power(t2,2) + 
               t1*(-9 + 17*t2)) - 
            Power(s2,2)*(27 - 246*t2 + 151*Power(t2,2) + 
               16*Power(t2,3) + 2*Power(t1,2)*(1 + 6*t2) + 
               t1*(-97 + 170*t2 - 11*Power(t2,2))) + 
            s2*(87 - 8*Power(t1,3)*(-1 + t2) - 375*t2 + 206*Power(t2,2) + 
               26*Power(t2,3) + 
               Power(t1,2)*(4 + 51*t2 - 27*Power(t2,2)) + 
               t1*(-215 + 312*t2 + Power(t2,2) - 22*Power(t2,3)))) - 
         Power(s1,2)*(25 + 57*t2 - 216*Power(t2,2) + 139*Power(t2,3) + 
            15*Power(t2,4) - 4*Power(t2,5) + 
            4*Power(t1,4)*(4 - 5*t2 + Power(t2,2)) + 
            Power(s2,4)*(3 - 5*t1 - 21*t2 + 9*t1*t2 + 12*Power(t2,2)) + 
            2*Power(t1,3)*(-15 + 7*t2 - 4*Power(t2,2) + 8*Power(t2,3)) + 
            t1*(23 - 316*t2 + 417*Power(t2,2) - 37*Power(t2,3) - 
               29*Power(t2,4)) + 
            Power(t1,2)*(-58 + 293*t2 - 143*Power(t2,2) - 
               48*Power(t2,3) + 4*Power(t2,4)) + 
            Power(s2,3)*(9 - 2*Power(t1,2)*(-3 + t2) + 200*t2 - 
               133*Power(t2,2) - 24*Power(t2,3) + 
               t1*(21 - 120*t2 + 19*Power(t2,2))) - 
            Power(s2,2)*(85 + 445*t2 - 322*Power(t2,2) - 
               134*Power(t2,3) - 4*Power(t2,4) + 
               Power(t1,3)*(1 + 7*t2) + 
               Power(t1,2)*(66 - 167*t2 + 61*Power(t2,2)) + 
               t1*(-6 - 163*t2 - 55*Power(t2,2) + 18*Power(t2,3))) + 
            s2*(14 + 12*Power(t1,4)*(-1 + t2) + 345*t2 - 
               131*Power(t2,2) - 233*Power(t2,3) + 13*Power(t2,4) + 
               Power(t1,3)*(71 - 71*t2 + 16*Power(t2,2)) + 
               2*Power(t1,2)*(-11 - 101*t2 + 58*Power(t2,2) + 
                  8*Power(t2,3)) + 
               t1*(121 - 68*t2 - 349*Power(t2,2) + 103*Power(t2,3) + 
                  5*Power(t2,4)))) + 
         s1*(45 - 3*Power(t1,5)*(-1 + t2) + 110*t2 - 133*Power(t2,2) + 
            19*Power(t2,3) + 67*Power(t2,4) - 8*Power(t2,5) + 
            Power(t1,4)*(-7 + 5*t2 - 4*Power(t2,2)) + 
            Power(s2,4)*(-5 + t1 - 52*t2 + 34*t1*t2 - 6*Power(t1,2)*t2 + 
               41*Power(t2,2) - 17*t1*Power(t2,2)) + 
            Power(t1,3)*(-28 - 13*t2 + 4*Power(t2,2) + 19*Power(t2,3)) + 
            t1*(-31 - 181*t2 + 167*Power(t2,2) + 140*Power(t2,3) - 
               65*Power(t2,4)) + 
            2*Power(t1,2)*(13 + 67*t2 - 3*Power(t2,2) - 57*Power(t2,3) + 
               Power(t2,4)) + Power(s2,3)*
             (116 + 143*t2 - 174*Power(t2,2) - 73*Power(t2,3) + 
               Power(t1,3)*(3 + 9*t2) + 
               Power(t1,2)*(44 - 74*t2 + 42*Power(t2,2)) + 
               t1*(-115 + 6*t2 + 25*Power(t2,3))) + 
            Power(s2,2)*(-197 - 200*t2 - 6*Power(t1,4)*t2 + 
               236*Power(t2,2) + 275*Power(t2,3) + 10*Power(t2,4) + 
               Power(t1,3)*(-81 + 70*t2 - 31*Power(t2,2)) + 
               Power(t1,2)*(78 + 17*t2 - 21*Power(t2,2) - 
                  38*Power(t2,3)) + 
               t1*(72 + 51*t2 + 76*Power(t2,2) - 53*Power(t2,3) - 
                  2*Power(t2,4))) + 
            s2*(27 + 3*Power(t1,5)*(-1 + t2) + 15*t2 + 112*Power(t2,2) - 
               385*Power(t2,3) + 11*Power(t2,4) + 
               3*Power(t1,4)*(5 - 5*t2 + 4*Power(t2,2)) + 
               Power(t1,3)*(114 - 50*t2 - 29*Power(t2,2) + 
                  13*Power(t2,3)) + 
               Power(t1,2)*(-224 - 15*t2 + 109*Power(t2,2) + 
                  40*Power(t2,3) + 6*Power(t2,4)) + 
               t1*(163 + 42*t2 - 448*Power(t2,2) + 116*Power(t2,3) + 
                  19*Power(t2,4))))) + 
      Power(s,2)*(-42 - 22*t1 + 83*Power(t1,2) - 9*Power(t1,3) - 
         21*Power(t1,4) + 3*Power(t1,5) - 3*t2 - 165*t1*t2 + 
         243*Power(t1,2)*t2 - 167*Power(t1,3)*t2 + 66*Power(t1,4)*t2 - 
         8*Power(t1,5)*t2 + 165*t1*Power(t2,2) - 
         431*Power(t1,2)*Power(t2,2) + 265*Power(t1,3)*Power(t2,2) - 
         46*Power(t1,4)*Power(t2,2) + 7*Power(t1,5)*Power(t2,2) + 
         38*Power(t2,3) - 63*t1*Power(t2,3) + 
         129*Power(t1,2)*Power(t2,3) - 81*Power(t1,3)*Power(t2,3) - 
         3*Power(t1,4)*Power(t2,3) - 12*Power(t2,4) - 42*t1*Power(t2,4) + 
         38*Power(t1,2)*Power(t2,4) - 4*Power(t1,3)*Power(t2,4) - 
         13*Power(t2,5) + 13*t1*Power(t2,5) + 2*Power(t1,2)*Power(t2,5) - 
         Power(s2,5)*(2 + 15*t2 - 11*Power(t2,2) + 
            Power(t1,2)*(-1 + 3*t2) + t1*(1 - 14*t2 + 7*Power(t2,2))) + 
         Power(s2,4)*(24 + 47*t2 - 32*Power(t2,2) - 47*Power(t2,3) + 
            Power(t1,3)*(-1 + 7*t2) + 
            2*Power(t1,2)*(-1 - 10*t2 + 5*Power(t2,2)) + 
            t1*(-7 - 16*t2 + 14*Power(t2,2) + 23*Power(t2,3))) - 
         Power(s2,3)*(12 + 39*t2 + 14*Power(t2,2) - 215*Power(t2,3) - 
            18*Power(t2,4) + Power(t1,4)*(1 + 5*t2) + 
            Power(t1,3)*(5 - 5*t2 - 2*Power(t2,2)) + 
            Power(t1,2)*(-70 - 31*t2 + 94*Power(t2,2) + Power(t2,3)) + 
            2*t1*(41 + 16*t2 - 56*Power(t2,2) + 42*Power(t2,3) + 
               3*Power(t2,4))) + 
         Power(s2,2)*(-81 - 122*t2 + 187*Power(t2,2) - 242*Power(t2,3) - 
            62*Power(t2,4) + Power(t1,5)*(1 + t2) + 
            2*Power(t1,4)*(6 - 3*t2 + Power(t2,2)) + 
            Power(t1,3)*(-70 + 69*t2 + 2*Power(t2,2) - 19*Power(t2,3)) + 
            Power(t1,2)*(-23 - 158*t2 + 250*Power(t2,2) + 
               23*Power(t2,3) + 4*Power(t2,4)) + 
            t1*(179 + 256*t2 - 495*Power(t2,2) + 26*Power(t2,3) + 
               38*Power(t2,4))) - 
         s2*(-103 - 136*t2 + 78*Power(t2,2) + 104*Power(t2,3) - 
            148*Power(t2,4) + 7*Power(t2,5) + 
            Power(t1,5)*(4 - 7*t2 + 7*Power(t2,2)) + 
            Power(t1,4)*(-2 + 47*t2 - 52*Power(t2,2) + 5*Power(t2,3)) + 
            Power(t1,3)*(-119 - 38*t2 + 283*Power(t2,2) - 
               136*Power(t2,3) + 4*Power(t2,4)) + 
            Power(t1,2)*(183 - 3*t2 - 275*Power(t2,2) + 
               243*Power(t2,3) + 2*Power(t2,4) + 2*Power(t2,5)) + 
            t1*(29 + 117*t2 - 125*Power(t2,2) - 306*Power(t2,3) + 
               98*Power(t2,4) + 5*Power(t2,5))) + 
         Power(s1,5)*(24 - 40*t2 + 14*Power(t2,2) + 
            t1*(-37 + 49*t2 - 12*Power(t2,2)) + 
            Power(t1,2)*(11 - 15*t2 + 4*Power(t2,2)) + 
            Power(s2,2)*(5 + 11*t1*(-1 + t2) - 11*t2 + 4*Power(t2,2)) + 
            s2*(-23 + 7*Power(t1,2)*(-1 + t2) + 39*t2 - 12*Power(t2,2) + 
               6*t1*(7 - 8*t2 + Power(t2,2)))) + 
         Power(s1,4)*(-47 + 35*t2 + 55*Power(t2,2) - 33*Power(t2,3) - 
            4*Power(t1,3)*(7 - 13*t2 + 6*Power(t2,2)) + 
            Power(t1,2)*(73 - 176*t2 + 103*Power(t2,2) - 8*Power(t2,3)) + 
            t1*(18 + 111*t2 - 162*Power(t2,2) + 31*Power(t2,3)) + 
            2*Power(s2,3)*(6 - 27*t2 + 14*Power(t2,2) + 
               t1*(-17 + 21*t2)) - 
            Power(s2,2)*(61 - 204*t2 + 97*Power(t2,2) + 8*Power(t2,3) + 
               4*Power(t1,2)*(-3 + 5*t2) + 
               t1*(-139 + 154*t2 + 3*Power(t2,2))) + 
            s2*(124 - 4*Power(t1,3)*(-1 + t2) - 273*t2 + 
               106*Power(t2,2) + 9*Power(t2,3) + 
               Power(t1,2)*(-13 + 44*t2 - 15*Power(t2,2)) + 
               t1*(-215 + 205*t2 + 33*Power(t2,2) - 11*Power(t2,3)))) + 
         Power(s1,3)*(-15 + 147*t2 - 189*Power(t2,2) + 29*Power(t2,3) + 
            26*Power(t2,4) + 6*Power(t1,4)*(3 - 5*t2 + 2*Power(t2,2)) + 
            Power(s2,4)*(7 - 17*t1 - 49*t2 + 25*t1*t2 + 28*Power(t2,2)) + 
            2*Power(t1,3)*(-4 + 9*t2 - 20*Power(t2,2) + 16*Power(t2,3)) + 
            t1*(197 - 560*t2 + 283*Power(t2,2) + 116*Power(t2,3) - 
               28*Power(t2,4)) + 
            Power(t1,2)*(-200 + 437*t2 - 56*Power(t2,2) - 
               115*Power(t2,3) + 4*Power(t2,4)) + 
            Power(s2,3)*(9 + Power(t1,2)*(20 - 30*t2) + 271*t2 - 
               160*Power(t2,2) - 40*Power(t2,3) + 
               t1*(79 - 131*t2 - 34*Power(t2,2))) + 
            Power(s2,2)*(36 - 521*t2 + 230*Power(t2,2) + 
               131*Power(t2,3) + 4*Power(t2,4) + 
               Power(t1,3)*(-9 + 11*t2) + 
               Power(t1,2)*(-51 + 163*t2 - 22*Power(t2,2)) + 
               t1*(-250 + 213*t2 + 196*Power(t2,2) - 3*Power(t2,3))) + 
            s2*(-89 + 6*Power(t1,4)*(-1 + t2) + 228*t2 + 
               141*Power(t2,2) - 240*Power(t2,3) + 16*Power(t2,4) + 
               Power(t1,3)*(27 - 17*t2 - 14*Power(t2,2)) + 
               Power(t1,2)*(117 - 446*t2 + 176*Power(t2,2) + 
                  3*Power(t2,3)) + 
               t1*(155 + 281*t2 - 615*Power(t2,2) + 87*Power(t2,3) + 
                  6*Power(t2,4)))) + 
         Power(s1,2)*(30 + 9*Power(t1,5)*(-1 + t2) - 201*t2 + 
            120*Power(t2,2) + 72*Power(t2,3) - 64*Power(t2,4) - 
            7*Power(t2,5) + Power(s2,5)*
             (2 + 4*t1*(-1 + t2) - 8*t2 + 4*Power(t2,2)) - 
            Power(t1,4)*(-29 + 5*t2 + 4*Power(t2,2) + 8*Power(t2,3)) - 
            Power(t1,3)*(28 + 6*t2 - 89*Power(t2,2) + 39*Power(t2,3) + 
               8*Power(t2,4)) + 
            Power(t1,2)*(239 - 292*t2 - 390*Power(t2,2) + 
               272*Power(t2,3) + 33*Power(t2,4)) + 
            t1*(-293 + 395*t2 + 185*Power(t2,2) - 371*Power(t2,3) + 
               15*Power(t2,4) + 9*Power(t2,5)) + 
            Power(s2,4)*(21 + 134*t2 - 77*Power(t2,2) - 24*Power(t2,3) + 
               6*Power(t1,2)*(1 + t2) + t1*(7 - 74*t2 + 13*Power(t2,2))) \
+ Power(s2,3)*(-197 - 357*t2 + 236*Power(t2,2) + 212*Power(t2,3) + 
               12*Power(t2,4) - 3*Power(t1,3)*(3 + 5*t2) - 
               2*Power(t1,2)*(28 - 63*t2 + 23*Power(t2,2)) + 
               t1*(82 + 22*t2 + 120*Power(t2,2) - 2*Power(t2,3))) + 
            Power(s2,2)*(178 + 310*t2 - 21*Power(t2,2) - 
               454*Power(t2,3) - 29*Power(t2,4) + 
               2*Power(t1,4)*(-1 + 7*t2) + 
               Power(t1,3)*(145 - 120*t2 + 31*Power(t2,2)) + 
               Power(t1,2)*(-196 - 54*t2 - 9*Power(t2,2) + 
                  37*Power(t2,3)) + 
               t1*(195 + 68*t2 - 579*Power(t2,2) - 27*Power(t2,3) + 
                  Power(t2,4))) - 
            s2*(48 + 9*Power(t1,5)*(-1 + t2) - 286*t2 + 588*Power(t2,2) - 
               390*Power(t2,3) - 81*Power(t2,4) + 13*Power(t2,5) + 
               Power(t1,4)*(59 - 47*t2 + 12*Power(t2,2)) + 
               Power(t1,3)*(22 + 75*t2 - 46*Power(t2,2) - 
                  11*Power(t2,3)) + 
               Power(t1,2)*(51 - 550*t2 + 21*Power(t2,2) + 
                  153*Power(t2,3) + Power(t2,4)) + 
               t1*(-23 + 791*t2 - 935*Power(t2,2) - 148*Power(t2,3) + 
                  80*Power(t2,4) + Power(t2,5)))) + 
         s1*(42 + 94*t2 - 48*Power(t2,2) - 74*Power(t2,3) + 
            42*Power(t2,4) + 20*Power(t2,5) - 
            Power(t1,5)*(-6 + t2 + 7*Power(t2,2)) + 
            Power(t1,4)*(-50 - 3*t2 + 18*Power(t2,2) + 11*Power(t2,3)) + 
            Power(t1,3)*(145 - 41*t2 - 146*Power(t2,2) + 48*Power(t2,3) + 
               12*Power(t2,4)) + 
            t1*(177 + 18*t2 - 259*Power(t2,2) + 183*Power(t2,3) + 
               71*Power(t2,4) - 22*Power(t2,5)) - 
            Power(t1,2)*(286 - 39*t2 - 494*Power(t2,2) + 
               166*Power(t2,3) + 83*Power(t2,4) + 2*Power(t2,5)) + 
            Power(s2,5)*((23 - 15*t2)*t2 + Power(t1,2)*(-1 + 3*t2) + 
               t1*(5 - 18*t2 + 7*Power(t2,2))) + 
            Power(s2,4)*(-76 + Power(t1,3)*(1 - 7*t2) - 104*t2 + 
               61*Power(t2,2) + 71*Power(t2,3) + 
               Power(t1,2)*(-28 + 42*t2 - 30*Power(t2,2)) + 
               t1*(65 + 9*t2 + 13*Power(t2,2) - 23*Power(t2,3))) + 
            Power(s2,3)*(212 + 211*t2 - 154*Power(t2,2) - 
               347*Power(t2,3) - 30*Power(t2,4) + 
               Power(t1,4)*(1 + 5*t2) + 
               Power(t1,3)*(62 - 46*t2 + 38*Power(t2,2)) + 
               Power(t1,2)*(-106 + 17*t2 - 4*Power(t2,2) + 
                  41*Power(t2,3)) + 
               3*t1*(-15 - 7*t2 - 10*Power(t2,2) + 2*Power(t2,3) + 
                  2*Power(t2,4))) - 
            Power(s2,2)*(61 - 24*t2 + 163*Power(t2,2) - 525*Power(t2,3) - 
               79*Power(t2,4) + Power(t1,5)*(1 + t2) + 
               Power(t1,4)*(34 - 20*t2 + 22*Power(t2,2)) + 
               Power(t1,3)*(90 - 8*t2 - 31*Power(t2,2) + 
                  21*Power(t2,3)) + 
               Power(t1,2)*(-394 + 67*t2 + 147*Power(t2,2) + 
                  28*Power(t2,3) + 12*Power(t2,4)) + 
               t1*(356 + 138*t2 - 625*Power(t2,2) - 60*Power(t2,3) + 
                  23*Power(t2,4))) + 
            s2*(-75 - 392*t2 + 423*Power(t2,2) - 79*Power(t2,3) - 
               229*Power(t2,4) + 20*Power(t2,5) + 
               Power(t1,5)*(-5 + 2*t2 + 7*Power(t2,2)) + 
               Power(t1,4)*(111 - 62*t2 + 5*Power(t2,3)) + 
               Power(t1,3)*(-224 + 290*t2 + 3*Power(t2,2) - 
                  67*Power(t2,3) + 4*Power(t2,4)) + 
               Power(t1,2)*(177 - 430*t2 - 47*Power(t2,2) + 
                  209*Power(t2,3) + 19*Power(t2,4) + 2*Power(t2,5)) + 
               2*t1*(20 + 271*t2 - 318*Power(t2,2) - 201*Power(t2,3) + 
                  70*Power(t2,4) + 3*Power(t2,5))))) - 
      s*(42 - 24*t1 - 16*Power(t1,2) - 46*Power(t1,3) + 58*Power(t1,4) - 
         14*Power(t1,5) + 2*Power(s1,6)*(-1 + s2)*(-1 + t1)*(-2 + s2 + t1)*
          (-1 + t2) - 56*t2 + 157*t1*t2 + 8*Power(t1,2)*t2 - 
         118*Power(t1,3)*t2 + 4*Power(t1,4)*t2 + 21*Power(t1,5)*t2 + 
         31*Power(t2,2) - 9*t1*Power(t2,2) - 149*Power(t1,2)*Power(t2,2) + 
         290*Power(t1,3)*Power(t2,2) - 118*Power(t1,4)*Power(t2,2) - 
         7*Power(t1,5)*Power(t2,2) - 39*Power(t2,3) - 87*t1*Power(t2,3) + 
         266*Power(t1,2)*Power(t2,3) - 202*Power(t1,3)*Power(t2,3) + 
         64*Power(t1,4)*Power(t2,3) - 4*Power(t1,5)*Power(t2,3) + 
         15*Power(t2,4) + t1*Power(t2,4) - 68*Power(t1,2)*Power(t2,4) + 
         36*Power(t1,3)*Power(t2,4) + 6*Power(t1,4)*Power(t2,4) + 
         7*Power(t2,5) + 26*t1*Power(t2,5) - 23*Power(t1,2)*Power(t2,5) - 
         2*Power(t1,3)*Power(t2,5) - 
         Power(s2,6)*(-1 + t1)*(2 - 3*t2 + Power(t2,2)) + 
         Power(s2,5)*(3 + 15*t2 + 11*Power(t2,2) - 17*Power(t2,3) + 
            Power(t1,2)*(9 - 15*t2 + 10*Power(t2,2)) + 
            t1*(-12 + 4*t2 - 13*Power(t2,2) + 9*Power(t2,3))) - 
         Power(s2,4)*(15 + 9*t2 + 49*Power(t2,2) - 89*Power(t2,3) - 
            18*Power(t2,4) + Power(t1,3)*(17 - 27*t2 + 22*Power(t2,2)) + 
            Power(t1,2)*(-23 - 36*t2 + 25*Power(t2,2) + 4*Power(t2,3)) + 
            t1*(-9 + 82*t2 - 66*Power(t2,2) + 33*Power(t2,3) + 
               6*Power(t2,4))) + 
         Power(s2,3)*(-33 - 51*t2 + 85*Power(t2,2) - 83*Power(t2,3) - 
            86*Power(t2,4) + 3*Power(t1,4)*(5 - 7*t2 + 6*Power(t2,2)) - 
            Power(t1,3)*(3 + 83*t2 - 87*Power(t2,2) + 17*Power(t2,3)) - 
            2*Power(t1,2)*(45 - 34*t2 + 11*Power(t2,2) - 41*Power(t2,3) + 
               5*Power(t2,4)) + 
            t1*(111 + 147*t2 - 148*Power(t2,2) - 92*Power(t2,3) + 
               34*Power(t2,4))) + 
         Power(s2,2)*(107 + 36*t2 + 48*Power(t2,2) - 159*Power(t2,3) + 
            156*Power(t2,4) - 2*Power(t2,5) + 
            Power(t1,5)*(-5 + 6*t2 - 5*Power(t2,2)) + 
            Power(t1,4)*(-29 + 73*t2 - 62*Power(t2,2) + 8*Power(t2,3)) + 
            2*Power(t1,3)*(68 + 11*t2 - 95*Power(t2,2) + 20*Power(t2,3) + 
               9*Power(t2,4)) - 
            Power(t1,2)*(80 + 257*t2 - 391*Power(t2,2) + 
               203*Power(t2,3) - 9*Power(t2,4) + 4*Power(t2,5)) - 
            t1*(129 - 84*t2 + 142*Power(t2,2) - 416*Power(t2,3) + 
               73*Power(t2,4) + 10*Power(t2,5))) + 
         s2*(-106 + 88*t2 - 153*Power(t2,2) + 177*Power(t2,3) - 
            41*Power(t2,4) - 29*Power(t2,5) + 
            Power(t1,5)*(19 - 27*t2 + 12*Power(t2,2) + 4*Power(t2,3)) - 
            2*Power(t1,4)*(22 + 20*t2 - 67*Power(t2,2) + 32*Power(t2,3) + 
               Power(t2,4)) + Power(t1,3)*
             (-70 + 84*t2 - 45*Power(t2,2) + 143*Power(t2,3) - 
               70*Power(t2,4) + 2*Power(t2,5)) + 
            t1*(47 - 389*t2 + 373*Power(t2,2) - 205*Power(t2,3) - 
               46*Power(t2,4) + 16*Power(t2,5)) + 
            Power(t1,2)*(154 + 268*t2 - 397*Power(t2,2) - 89*Power(t2,3) + 
               109*Power(t2,4) + 19*Power(t2,5))) + 
         Power(s1,5)*(-4 - 8*t2 + 12*Power(t2,2) + 
            t1*(-23 + 76*t2 - 49*Power(t2,2)) - 
            4*Power(t1,3)*(3 - 5*t2 + 2*Power(t2,2)) + 
            4*Power(s2,3)*(2 + 5*t1*(-1 + t2) - 5*t2 + 2*Power(t2,2)) + 
            Power(t1,2)*(43 - 76*t2 + 33*Power(t2,2)) - 
            Power(s2,2)*(45 + 8*Power(t1,2)*(-1 + t2) - 72*t2 + 
               19*Power(t2,2) + t1*(-73 + 68*t2 + Power(t2,2))) + 
            s2*(59 - 4*Power(t1,3)*(-1 + t2) - 80*t2 + 17*Power(t2,2) + 
               Power(t1,2)*(-23 + 28*t2 - 5*Power(t2,2)) + 
               4*t1*(-17 + 12*t2 + 3*Power(t2,2)))) + 
         Power(s1,4)*(-(Power(-1 + t2,2)*(41 + 3*t2)) + 
            4*Power(t1,4)*(2 - 5*t2 + 3*Power(t2,2)) + 
            Power(t1,2)*(-130 + 153*t2 + 64*Power(t2,2) - 
               65*Power(t2,3)) + 
            2*Power(t1,3)*(-1 + 21*t2 - 26*Power(t2,2) + 8*Power(t2,3)) + 
            t1*(163 - 268*t2 + Power(t2,2) + 84*Power(t2,3)) + 
            Power(s2,4)*(7 - 37*t2 + 20*Power(t2,2) + t1*(-21 + 25*t2)) - 
            Power(s2,3)*(26 - 157*t2 + 75*Power(t2,2) + 16*Power(t2,3) + 
               4*Power(t1,2)*(-7 + 9*t2) + 
               t1*(-90 + 67*t2 + 43*Power(t2,2))) + 
            Power(s2,2)*(102 - 228*t2 + 52*Power(t2,2) + 24*Power(t2,3) + 
               Power(t1,3)*(-15 + 19*t2) + 
               Power(t1,2)*(-48 + 81*t2 + 5*Power(t2,2)) + 
               t1*(-215 + 66*t2 + 157*Power(t2,2))) + 
            s2*(-52 - 4*Power(t1,4)*(-1 + t2) - 11*t2 + 148*Power(t2,2) - 
               65*Power(t2,3) + 
               Power(t1,3)*(-19 + 27*t2 - 16*Power(t2,2)) + 
               Power(t1,2)*(178 - 318*t2 + 87*Power(t2,2) + Power(t2,3)) + 
               t1*(-11 + 340*t2 - 325*Power(t2,2) + 24*Power(t2,3)))) + 
         Power(s1,3)*(75 + 9*Power(t1,5)*(-1 + t2) - 83*t2 - 
            9*Power(t2,2) + 31*Power(t2,3) - 14*Power(t2,4) + 
            Power(s2,5)*(3 + 7*t1*(-1 + t2) - 15*t2 + 8*Power(t2,2)) + 
            Power(t1,4)*(17 + 25*t2 - 20*Power(t2,2) - 16*Power(t2,3)) + 
            t1*(-255 + 99*t2 + 384*Power(t2,2) - 175*Power(t2,3) - 
               49*Power(t2,4)) + 
            Power(t1,3)*(52 - 231*t2 + 143*Power(t2,2) + 18*Power(t2,3) - 
               8*Power(t2,4)) + 
            Power(t1,2)*(96 + 129*t2 - 454*Power(t2,2) + 
               110*Power(t2,3) + 43*Power(t2,4)) + 
            Power(s2,4)*(19 - 6*Power(t1,2)*(-2 + t2) + 122*t2 - 
               57*Power(t2,2) - 32*Power(t2,3) + 
               t1*(25 - 48*t2 - 23*Power(t2,2))) + 
            Power(s2,3)*(-73 - 276*t2 + 89*Power(t2,2) + 
               156*Power(t2,3) + 8*Power(t2,4) - 2*Power(t1,3)*(5 + t2) + 
               2*Power(t1,2)*(-20 + 45*t2 + 7*Power(t2,2)) + 
               t1*(-85 - 8*t2 + 177*Power(t2,2) + 40*Power(t2,3))) + 
            s2*(-118 - 9*Power(t1,5)*(-1 + t2) + 392*t2 - 
               367*Power(t2,2) + 28*Power(t2,3) + 61*Power(t2,4) + 
               Power(t1,4)*(-53 + 29*t2 + 12*Power(t2,2)) + 
               2*Power(t1,3)*(16 + 9*t2 - 22*Power(t2,2) + 
                  17*Power(t2,3)) + 
               t1*(329 - 940*t2 + 381*Power(t2,2) + 308*Power(t2,3) - 
                  54*Power(t2,4)) + 
               Power(t1,2)*(-299 + 498*t2 + 202*Power(t2,2) - 
                  174*Power(t2,3) + Power(t2,4))) + 
            Power(s2,2)*(40 + 48*t2 + 188*Power(t2,2) - 235*Power(t2,3) + 
               11*Power(t2,4) + 2*Power(t1,4)*(-2 + 5*t2) + 
               Power(t1,3)*(90 - 77*t2 - 15*Power(t2,2)) - 
               Power(t1,2)*(17 + 195*t2 + 6*Power(t2,2) + 4*Power(t2,3)) + 
               t1*(171 + 414*t2 - 603*Power(t2,2) - 89*Power(t2,3) + 
                  Power(t2,4)))) + 
         s1*(Power(s2,6)*(-1 + t1)*(3 - 4*t2 + Power(t2,2)) + 
            Power(t1,5)*(20 - 37*t2 + 21*Power(t2,2) + 4*Power(t2,3)) - 
            Power(-1 + t2,2)*(76 + 43*t2 + 57*Power(t2,2) + 
               16*Power(t2,3)) + 
            Power(t1,4)*(-110 + 163*t2 + 45*Power(t2,2) - 
               78*Power(t2,3) - 10*Power(t2,4)) + 
            t1*(22 - 486*t2 + 280*Power(t2,2) + 140*Power(t2,3) - 
               74*Power(t2,4) - 34*Power(t2,5)) + 
            Power(t1,3)*(170 - 215*t2 - 117*Power(t2,2) + 
               184*Power(t2,3) - 22*Power(t2,4) + 2*Power(t2,5)) + 
            Power(t1,2)*(-42 + 390*t2 - 270*Power(t2,2) - 
               269*Power(t2,3) + 123*Power(t2,4) + 32*Power(t2,5)) + 
            Power(s2,5)*(-10 - 38*t2 - 5*Power(t2,2) + 25*Power(t2,3) - 
               2*Power(t1,2)*(4 - 5*t2 + 5*Power(t2,2)) + 
               t1*(14 + 12*t2 + 7*Power(t2,2) - 9*Power(t2,3))) - 
            Power(s2,3)*(64 + 8*t2 + 131*Power(t2,2) - 320*Power(t2,3) - 
               127*Power(t2,4) + 
               2*Power(t1,4)*(8 - 5*t2 + 9*Power(t2,2)) + 
               Power(t1,3)*(21 - 19*t2 + 19*Power(t2,2) + 
                  23*Power(t2,3)) + 
               t1*(183 + 211*t2 - 346*Power(t2,2) - 81*Power(t2,3) - 
                  11*Power(t2,4)) + 
               2*Power(t1,2)*(-112 + 9*t2 + 75*Power(t2,2) - 
                  3*Power(t2,3) + 3*Power(t2,4))) + 
            Power(s2,4)*(85 + 109*t2 - 13*Power(t2,2) - 155*Power(t2,3) - 
               30*Power(t2,4) + 
               2*Power(t1,3)*(7 - 6*t2 + 11*Power(t2,2)) + 
               Power(t1,2)*(-26 - 20*t2 + 24*Power(t2,3)) + 
               t1*(-45 + 23*t2 + 11*Power(t2,2) - 5*Power(t2,3) + 
                  6*Power(t2,4))) + 
            Power(s2,2)*(-55 - 313*t2 + 292*Power(t2,2) - 
               72*Power(t2,3) - 232*Power(t2,4) + 16*Power(t2,5) + 
               Power(t1,5)*(7 - 4*t2 + 5*Power(t2,2)) + 
               Power(t1,4)*(63 - 56*t2 + 39*Power(t2,2) + 
                  12*Power(t2,3)) + 
               Power(t1,3)*(-240 + 219*t2 + 109*Power(t2,2) - 
                  68*Power(t2,3) - 2*Power(t2,4)) + 
               Power(t1,2)*(120 - 204*t2 - 57*Power(t2,2) + 
                  93*Power(t2,3) - 6*Power(t2,4) + 4*Power(t2,5)) + 
               t1*(141 + 494*t2 - 536*Power(t2,2) - 449*Power(t2,3) + 
                  32*Power(t2,4) + 12*Power(t2,5))) + 
            s2*(103 + 97*t2 + 83*Power(t2,2) - 281*Power(t2,3) + 
               102*Power(t2,4) + 48*Power(t2,5) - 
               Power(t1,5)*(27 - 41*t2 + 26*Power(t2,2) + 4*Power(t2,3)) + 
               Power(t1,4)*(47 - 133*t2 + 6*Power(t2,2) + 
                  34*Power(t2,3) + 2*Power(t2,4)) + 
               t1*(124 + 268*t2 - 579*Power(t2,2) + 526*Power(t2,3) + 
                  103*Power(t2,4) - 42*Power(t2,5)) + 
               Power(t1,3)*(145 + 57*t2 - 311*Power(t2,2) + 
                  63*Power(t2,3) + 48*Power(t2,4) - 2*Power(t2,5)) - 
               Power(t1,2)*(376 + 266*t2 - 1023*Power(t2,2) + 
                  154*Power(t2,3) + 147*Power(t2,4) + 20*Power(t2,5)))) + 
         Power(s1,2)*(Power(s2,6)*(-1 + t1)*(-1 + t2) + 
            Power(t1,5)*(3 + 7*t2 - 14*Power(t2,2)) + 
            Power(t1,4)*(15 - 128*t2 + 29*Power(t2,2) + 50*Power(t2,3) + 
               4*Power(t2,4)) + 
            2*Power(t1,3)*(-61 + 175*t2 - 28*Power(t2,2) - 
               60*Power(t2,3) + 5*Power(t2,4)) + 
            t2*(-37 + 48*t2 - 44*Power(t2,2) + 24*Power(t2,3) + 
               9*Power(t2,4)) + 
            Power(t1,2)*(3 - 430*t2 + 532*Power(t2,2) + 106*Power(t2,3) - 
               130*Power(t2,4) - 9*Power(t2,5)) + 
            t1*(139 + 352*t2 - 511*Power(t2,2) - 26*Power(t2,3) + 
               138*Power(t2,4) + 8*Power(t2,5)) + 
            Power(s2,5)*(4 + 38*t2 - 14*Power(t2,2) - 8*Power(t2,3) + 
               Power(t1,2)*(-1 + 5*t2) + t1*(5 - 23*t2 + 6*Power(t2,2))) + 
            Power(s2,4)*(-108 + Power(t1,3)*(3 - 15*t2) - 141*t2 + 
               47*Power(t2,2) + 118*Power(t2,3) + 12*Power(t2,4) + 
               Power(t1,2)*(-21 + 34*t2 - 27*Power(t2,2)) - 
               2*t1*(-28 + 3*t2 - 25*Power(t2,2) + Power(t2,3))) + 
            Power(s2,3)*(196 + 174*t2 + 32*Power(t2,2) - 353*Power(t2,3) - 
               65*Power(t2,4) + Power(t1,4)*(1 + 11*t2) + 
               Power(t1,3)*(58 - 22*t2 + 36*Power(t2,2)) + 
               2*Power(t1,2)*(-81 + 24*t2 - 21*Power(t2,2) + 
                  8*Power(t2,3)) + 
               t1*(95 + 79*t2 - 244*Power(t2,2) - 117*Power(t2,3) - 
                  13*Power(t2,4))) - 
            Power(s2,2)*(131 - 303*t2 + 421*Power(t2,2) - 
               334*Power(t2,3) - 97*Power(t2,4) + 14*Power(t2,5) + 
               2*Power(t1,5)*(1 + t2) + 
               Power(t1,4)*(42 - 17*t2 + 29*Power(t2,2)) - 
               Power(t1,3)*(21 - 159*t2 + 88*Power(t2,2) + 
                  4*Power(t2,3)) + 
               Power(t1,2)*(-89 - 319*t2 + 21*Power(t2,2) + 
                  6*Power(t2,3) + 3*Power(t2,4)) + 
               t1*(111 + 712*t2 - 733*Power(t2,2) - 354*Power(t2,3) + 
                  8*Power(t2,4) + 2*Power(t2,5))) + 
            s2*(104 - 428*t2 + 176*Power(t2,2) + 205*Power(t2,3) - 
               138*Power(t2,4) - 19*Power(t2,5) + 
               Power(t1,5)*(-1 - 5*t2 + 14*Power(t2,2)) + 
               Power(t1,4)*(70 + 60*t2 - 48*Power(t2,2) - 10*Power(t2,3)) - 
               2*Power(t1,3)*(74 - 17*t2 - 60*Power(t2,2) + 
                  36*Power(t2,3) + 5*Power(t2,4)) + 
               Power(t1,2)*(392 - 316*t2 - 726*Power(t2,2) + 
                  264*Power(t2,3) + 85*Power(t2,4) + Power(t2,5)) + 
               t1*(-405 + 593*t2 + 242*Power(t2,2) - 693*Power(t2,3) - 
                  3*Power(t2,4) + 26*Power(t2,5)))))))/
  (Power(-1 + s1,2)*Power(-1 + s2,2)*Power(-1 + t1,2)*Power(s - s2 + t1,2)*
    Power(-1 + t2,2)*Power(s - s1 + t2,2));
  
  return a;
};