
#include <gsl/gsl_math.h>
#include <gsl/gsl_pow_int.h>
#include "nlo_functions.h"

#define Pi M_PI
#define Sqrt sqrt
#define Power gsl_pow_int
#define ln gsl_sf_log_abs
#define Log gsl_sf_log_abs

double Box1(double s,double s1,double s2,double t1,double t2){
    double a=(((-64*(-1 + s2 - t1 + t2)*(-2*Power(s,2)*(-1 + s1) + 
            s*(-4 + t1 + s1*(2*s1 + 2*s2 - t1 - t2) + t2) - 
            (1 + s1)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2))*
          ((4 - Power(s1,4) + Power(s,3)*(3 + s1) + 15*s2 + 
               Power(s2,2) - 4*t1 - s2*t1 - 
               Power(s,2)*(2 + 3*Power(s1,2) + 4*s2 + 
                  s1*(3 + 2*s2 - t2) - 5*t2) - 4*t2 - s2*t2 + 
               2*Power(t2,2) - 3*s2*Power(t2,2) + 
               Power(s1,3)*(5 - 2*s2 + t2) + 
               Power(s1,2)*(9 - Power(s2,2) + t1 - 3*t2 + 
                  2*s2*(2 + t2)) + 
               s*(-19 + 3*Power(s1,3) + 3*s2 + t1 + 
                  Power(s1,2)*(-5 + 4*s2 - 2*t2) + 5*t2 - 7*s2*t2 + 
                  2*Power(t2,2) + 
                  s1*(-7 + 2*s2 + Power(s2,2) - t1 - 2*t2 - 2*s2*t2)) \
+ s1*(s2*(4 + t1 - 3*t2) + Power(s2,2)*(-1 + t2) + 
                  3*(5 + t1 - 3*t2 - Power(t2,2))))/
             ((-1 + s1)*(-1 + t2)) + 
            (-7 + Power(s,3)*s1 - Power(s1,4) + 11*s2 - 
               6*Power(s2,2) - t1 + 7*s2*t1 + 9*t2 - 5*s2*t2 + 
               Power(s2,2)*t2 - 2*s2*t1*t2 + 2*s2*Power(t2,2) + 
               Power(s1,3)*(5 - 2*s2 + t2) - 
               Power(s,2)*(5 + 3*Power(s1,2) + s2 - 2*t1 + 
                  s1*(-8 + 2*s2 - t2) + 2*t2) + 
               s1*(8 + s2*(8 + 3*t1 - 11*t2) + t1*(9 - 2*t2) + 
                  Power(s2,2)*(-3 + t2) - 7*t2 + 2*Power(t2,2)) + 
               Power(s1,2)*(17 - Power(s2,2) + 2*t1 - 11*t2 + 
                  s2*(3 + 2*t2)) + 
               s*(-2 + 3*Power(s1,3) + Power(s2,2) - 7*t1 + 
                  Power(s1,2)*(-13 + 4*s2 - 2*t2) + t2 + 2*t1*t2 - 
                  2*Power(t2,2) + s2*(7 - 2*t1 + t2) + 
                  s1*(-11 + Power(s2,2) - 5*t1 + 14*t2 - 
                     2*s2*(2 + t2))))/((s - s2 + t1)*(s - s1 + t2)) + 
            (-1 + Power(s,3)*(-1 + s1) - Power(s1,4) + 12*s2 + 
               2*Power(s2,2) + t1 - 2*s2*t1 + 3*t2 - 2*s2*t2 - 
               Power(s2,2)*t2 - 2*Power(t2,2) + 3*s2*Power(t2,2) + 
               Power(s1,3)*(7 - 2*s2 + t2) + 
               s1*(11 - 4*t1 + 2*s2*(4 + t1 - 5*t2) + 
                  Power(s2,2)*(-1 + t2) + 3*Power(t2,2)) + 
               Power(s,2)*(2 - 3*Power(s1,2) + 3*s2 - 3*t2 + 
                  s1*(6 - 2*s2 + t2)) + 
               Power(s1,2)*(5 - Power(s2,2) + 3*t1 - 10*t2 + 
                  s2*(5 + 2*t2)) + 
               s*(-9 + 3*Power(s1,3) - Power(s2,2) + 2*t1 + 
                  2*Power(s1,2)*(-6 + 2*s2 - t2) + 6*s2*(-1 + t2) - 
                  2*t2 - 2*Power(t2,2) + 
                  s1*(Power(s2,2) - 2*(4 + t1 - 6*t2) - 2*s2*(3 + t2))\
))/((-1 + s1)*(-s + s1 - t2)) + 
            (-9 - Power(s1,4) + Power(s,3)*(1 + s1) + 8*s2 - 
               2*Power(s2,2) - 2*Power(s2,3) + 3*t1 - s2*t1 + 
               6*Power(s2,2)*t1 - 4*s2*Power(t1,2) + 
               Power(s1,3)*(2 - 3*s2 + t1) + 
               Power(s,2)*(1 - 3*Power(s1,2) - 3*s2 + 4*t1 + 
                  s1*(2 - 3*s2 + t1) - t2) + 3*t2 - 3*s2*t2 - 
               2*Power(s2,2)*t2 + 3*s2*t1*t2 - s2*Power(t2,2) + 
               Power(s1,2)*(16 - 2*Power(s2,2) - 5*t1 + 
                  s2*(3 + t1 + t2)) + 
               s1*(12 - 4*Power(t1,2) + s2*(15 + t1 - 2*t2) + 
                  Power(s2,2)*(-1 + t2) - Power(t2,2) + 
                  t1*(-4 + 3*t2)) + 
               s*(-3 + 3*Power(s1,3) + 5*Power(s2,2) + 
                  Power(s1,2)*(-5 + 6*s2 - 2*t1) + t1 + 
                  4*Power(t1,2) - t2 - 3*t1*t2 + Power(t2,2) + 
                  s2*(-3 - 10*t1 + 4*t2) + 
                  s1*(-17 + 2*Power(s2,2) + t1 + t2 - 
                     s2*(2 + t1 + t2))))/((-1 + s2)*(-s + s2 - t1)) + 
            (-6 - 2*Power(s1,4) + Power(s,3)*(-3 + 2*s1) + 13*s2 + 
               7*Power(s2,2) - 2*t1 - 17*s2*t1 - 2*Power(s2,2)*t1 + 
               4*s2*Power(t1,2) + 2*t2 + 6*s2*t2 - Power(s2,2)*t2 - 
               3*s2*t1*t2 + s2*Power(t2,2) + 
               Power(s1,3)*(2 - 4*s2 + 2*t2) + 
               s1*(17 + 4*Power(t1,2) + s2*(17 - 4*t1 - 3*t2) + 
                  2*Power(s2,2)*t2 + Power(t2,2) - 3*t1*(5 + t2)) - 
               2*Power(s1,2)*
                (-8 + Power(s2,2) + t1 + t2 - s2*(1 + 2*t2)) + 
               Power(s,2)*(16 - 6*Power(s1,2) + 4*s2 - 5*t1 + 
                  s1*(7 - 4*s2 + 2*t2)) + 
               s*(-11 + 6*Power(s1,3) - Power(s2,2) + 17*t1 - 
                  4*Power(t1,2) + Power(s1,2)*(-6 + 8*s2 - 4*t2) - 
                  8*t2 + 3*t1*t2 - Power(t2,2) + 
                  s2*(-25 + 7*t1 + t2) + 
                  s1*(-32 + 2*Power(s2,2) + 7*t1 + 2*t2 - 
                     s2*(5 + 4*t2))))/((-1 + s2)*(-1 + t1)) + 
            (2 - Power(s,4) - Power(s1,4) + 19*s2 - 2*Power(s2,2) - 
               3*t1 - 5*s2*t1 + Power(t1,2) + 
               Power(s,3)*(4*s1 + 2*s2 - t1 - t2) - 8*t2 - 8*s2*t2 + 
               3*Power(s2,2)*t2 + t1*t2 + s2*t1*t2 - 
               2*s2*Power(t2,2) + Power(s1,3)*(1 - 2*s2 + t2) - 
               Power(s1,2)*(-17 + Power(s2,2) + t1 - 2*s2*(-1 + t2) + 
                  t2) + s1*(25 - Power(t1,2) + 
                  Power(s2,2)*(-3 + t2) - 12*t2 - 2*Power(t2,2) + 
                  t1*(-3 + 2*t2) + s2*(11 + 2*t2)) + 
               Power(s,2)*(16 - 6*Power(s1,2) - Power(s2,2) + t1 - 
                  2*t2 - t1*t2 + s2*(-3 + t1 + 2*t2) + 
                  s1*(2 - 6*s2 + 2*t1 + 3*t2)) + 
               s*(-27 + 4*Power(s1,3) + 6*t1 + 
                  s2*(-11 + t1*(-1 + t2) - t2) - 
                  Power(s2,2)*(-3 + t2) + 11*t2 - t1*t2 + 
                  2*Power(t2,2) + 
                  Power(s1,2)*(6*s2 - t1 - 3*(1 + t2)) + 
                  s1*(-33 + 2*Power(s2,2) + (3 + t1)*t2 - 
                     s2*(-4 + t1 + 4*t2))))/((-1 + t1)*(-1 + t2))))/
        (Power(s,3)*(-1 + s1)*(-1 + s2 - t1 + t2) - 
          Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
          Power(s,2)*(3 + s2*(-3 + t1) + 2*t1 - (3 + t1)*t2 + 
             Power(t2,2) + 2*Power(s1,2)*(-1 + s2 - t1 + t2) + 
             s1*(s2 - t2)*(-1 + s2 - t1 + t2)) + 
          s*(2*(-1 + s2)*(-1 + t1) - 
             (2*(2 + t1) + s2*(-5 + s2 + t1))*t2 + 
             (2 + s2)*Power(t2,2) + Power(s1,3)*(-1 + s2 - t1 + t2) + 
             Power(s1,2)*(1 + s2 - t2)*(-1 + s2 - t1 + t2) + 
             s1*(2 + s2*(-3 + s2 - t1) + 4*t1 - t2 - 
                (1 + s2)*(s2 - t1)*t2 - (1 + s2)*Power(t2,2)))) + 
       8*((-8*(-11 + 2*Power(s,3)*(-1 + s1)*s1 + 11*s2 - 
               2*Power(s2,2) - 2*Power(s1,4)*(s2 - t1) - 6*t1 - 
               s2*t1 + 3*Power(t1,2) + 6*t2 - 10*s2*t2 + 
               3*Power(s2,2)*t2 + 6*t1*t2 - s2*t1*t2 - 
               2*Power(t1,2)*t2 + 5*Power(t2,2) - 3*s2*Power(t2,2) + 
               2*t1*Power(t2,2) - 
               2*Power(s,2)*(-9 + 2*Power(s1,3) + s2 + 
                  Power(s1,2)*(-7 + 3*s2 - 2*t1) + 2*t1 - 2*t2 + 
                  s1*(13 - 3*s2 + t1 + t2)) - 
               Power(s1,3)*(3*Power(s2,2) + Power(t1,2) + 
                  2*(-1 + t2) + t1*(5 + t2) - s2*(5 + 4*t1 + t2)) + 
               Power(s1,2)*(-5 - 2*Power(s2,3) + 2*Power(t1,3) + 
                  Power(t1,2)*(5 - 2*t2) + 4*t2 + Power(t2,2) + 
                  Power(s2,2)*(6*t1 + t2) + 6*t1*(-3 + 2*t2) + 
                  s2*(19 - 6*Power(t1,2) + t1*(-5 + t2) - 14*t2 + 
                     Power(t2,2))) + 
               s1*(2*Power(s2,3) + 2*Power(t1,3) - 
                  Power(s2,2)*(9 + 2*t1 + 2*t2) - 
                  Power(t1,2)*(3 + 4*t2) + 
                  t1*(-7 - 11*t2 + 2*Power(t2,2)) - 
                  4*(4 - 7*t2 + 3*Power(t2,2)) + 
                  s2*(11 - 2*Power(t1,2) + 5*t2 + 6*t1*(2 + t2))) + 
               s*(-4 + 2*Power(s1,4) + 2*Power(s2,2) + 11*t1 - 
                  3*Power(t1,2) + Power(s1,3)*(8*s2 - 6*(1 + t1)) + 
                  s2*(-10 + t1 - 7*t2) + 21*t2 + 3*Power(t2,2) + 
                  Power(s1,2)*
                   (-16 + 6*Power(s2,2) + 17*t1 + 3*Power(t1,2) + 
                     15*t2 - Power(t2,2) - s2*(14 + 9*t1 + t2)) - 
                  2*s1*(-8 + 3*Power(s2,2) + 15*t2 + t1*(11 + t2) - 
                     s2*(11 + 3*t1 + 2*t2)))))/
           ((s - s2 + t1)*(1 + s*(-1 + s1) + s1*(-s2 + t1) - t2)*
             (s - s1 + t2)) - 
          (8*(-19 + 2*Power(s,3)*Power(-1 + s1,2) - s2 - 
               2*Power(s2,2) - 2*Power(s1,4)*(s2 - t1) + 10*t1 + 
               3*s2*t1 - Power(t1,2) + 
               Power(s1,3)*(2 - 5*Power(s2,2) - 3*Power(t1,2) + 
                  s2*(11 + 8*t1 - t2) + t1*(-11 + t2) - 2*t2) + 
               10*t2 + 2*s2*t2 + 3*Power(s2,2)*t2 - 8*t1*t2 - 
               5*s2*t1*t2 + 2*Power(t1,2)*t2 + 7*Power(t2,2) - 
               3*s2*Power(t2,2) + 2*Power(t2,3) - 
               2*Power(s,2)*(-2 + 2*Power(s1,3) + 3*s2 - t1 + 
                  Power(s1,2)*(3*s2 - 2*(2 + t1)) - 3*t2 + 
                  s1*(4 - 6*s2 + 5*t1 + t2)) - 
               Power(s1,2)*(11 + 2*Power(s2,3) - 2*Power(t1,3) - 
                  12*t2 - 12*t1*t2 + Power(t2,2) + 
                  Power(s2,2)*(-6 - 6*t1 + t2) + 
                  Power(t1,2)*(-3 + 2*t2) + 
                  s2*(1 + 6*Power(t1,2) - 3*t1*(-3 + t2) + 10*t2 + 
                     Power(t2,2))) + 
               s*(16 + 2*Power(s1,4) + 2*Power(s2,2) + 
                  2*Power(s1,3)*(-5 + 4*s2 - 3*t1) - 11*t1 + 
                  Power(t1,2) + s2*(6 - 3*t1 - 11*t2) + 13*t2 + 
                  4*t1*t2 + 5*Power(t2,2) + 
                  Power(s1,2)*
                   (12 + 6*Power(s2,2) + 15*t1 + 5*Power(t1,2) + 
                     7*t2 - 4*t1*t2 + Power(t2,2) + 
                     s2*(-16 - 11*t1 + t2)) - 
                  2*s1*(11 + 5*Power(s2,2) - 3*t1 + 4*Power(t1,2) + 
                     12*t2 + 2*t1*t2 - 3*s2*(1 + 3*t1 + t2))) + 
               s1*(6 + 2*Power(s2,3) - 2*Power(t1,3) + 
                  Power(s2,2)*(3 - 6*t1 - 4*t2) + 
                  Power(t1,2)*(11 - 2*t2) + 2*t2 - 8*Power(t2,2) - 
                  t1*(17 + 13*t2) + 
                  s2*(13 + 6*Power(t1,2) + 17*t2 + 2*t1*(-7 + 3*t2)))))/
           ((-1 + s1)*(-s + s1 - t2)*
             (1 + s*(-1 + s1) + s1*(-s2 + t1) - t2)) - 
          (8*(-27 + 2*Power(s,3)*(-3 + 2*s1 + Power(s1,2)) - s2 - 
               3*Power(s1,4)*(s2 - t1) - 4*t1 + s2*t1 - Power(t1,2) + 
               45*t2 - 4*s2*t2 - 3*Power(s2,2)*t2 + t1*t2 + 
               5*s2*t1*t2 - 2*Power(t1,2)*t2 - 16*Power(t2,2) + 
               8*s2*Power(t2,2) - 2*Power(t2,3) + 
               Power(s,2)*(4 - 4*Power(s1,3) + 8*s2 + t1 + 
                  s1*(2 - 10*s2 + 12*t1) - 17*t2 + 
                  Power(s1,2)*(6 - 6*s2 + 3*t1 + t2)) - 
               Power(s1,3)*(5*Power(s2,2) + Power(t1,2) + 
                  3*(-1 + t2) + 2*t1*(6 + t2) - s2*(13 + 6*t1 + t2)) + 
               s1*(2*Power(t1,3) + 
                  s2*(23 - 4*Power(t1,2) - 4*t1*(-1 + t2) - 11*t2) + 
                  9*(-1 + t2) + Power(t1,2)*(-7 + 2*t2) + 
                  Power(s2,2)*(3 + 2*t1 + 2*t2) + t1*(-34 + 22*t2)) + 
               s*(23 + 2*Power(s1,4) + Power(t1,2) - 
                  s2*(4 + t1 - 22*t2) + 
                  Power(s1,3)*(-14 + 9*s2 - 6*t1 - t2) - 12*t2 - 
                  5*t1*t2 - 10*Power(t2,2) + 
                  s1*(-22 + s2 + 6*Power(s2,2) - 16*s2*t1 + 
                     10*Power(t1,2) + 17*t2 - 2*s2*t2 + 2*t1*t2) + 
                  Power(s1,2)*
                   (3 + 6*Power(s2,2) - Power(t1,2) + 6*t2 + 
                     3*t1*t2 - s2*(10 + 5*t1 + 2*t2))) + 
               Power(s1,2)*(-13 - 2*Power(s2,3) - 2*Power(t1,3) + 
                  t1*(-11 + t2) + 11*t2 + 2*Power(t2,2) + 
                  Power(s2,2)*(4 + 2*t1 + t2) + 
                  Power(t1,2)*(-3 + 2*t2) + 
                  s2*(14 + 2*Power(t1,2) - 4*t2 - t1*(1 + 3*t2)))))/
           ((-1 + s1)*(1 + s*(-1 + s1) + s1*(-s2 + t1) - t2)*(-1 + t2)) \
- (8*(-2*Power(s,3)*(-1 + Power(s1,2)) + 3*Power(s1,4)*(s2 - t1) + 
               Power(s1,3)*(-3 + 8*Power(s2,2) + 4*t1 + 
                  4*Power(t1,2) + 3*t2 - t1*t2 + 
                  s2*(-5 - 12*t1 + 2*t2)) + 
               Power(s,2)*(-2 + 4*Power(s1,3) - 6*s2 + 8*t1 + 
                  Power(s1,2)*(-12 + 8*s2 - 5*t1 + t2) + 
                  s1*(6 + 2*s2 - 9*t1 + 5*t2)) + 
               Power(s1,2)*(5 + 4*Power(s2,3) + 33*t1 - 
                  12*Power(t1,2) - 2*Power(s2,2)*(5 + 4*t1 - t2) - 
                  6*t2 + Power(t2,2) + 
                  s2*(-38 + 4*Power(t1,2) - 2*t1*(-11 + t2) + 4*t2 + 
                     Power(t2,2))) + 
               s1*(6*Power(s2,3) - 6*Power(t1,3) - 
                  2*Power(s2,2)*(6 + 9*t1 - 5*t2) + 
                  2*Power(t1,2)*t2 + t1*(-6 + 19*t2) + 
                  s2*(8 + 18*Power(t1,2) - 12*t1*(-1 + t2) - 22*t2 + 
                     Power(t2,2)) - 5*(-6 + 5*t2 + Power(t2,2))) + 
               2*(1 + 3*Power(s2,2)*(-1 + t2) + 
                  3*Power(t1,2)*(-1 + t2) + 2*t2 - 3*Power(t2,2) + 
                  t1*(2 + 2*t2 - Power(t2,2)) + 
                  s2*(4 - 6*t1*(-1 + t2) - 9*t2 + 2*Power(t2,2))) - 
               s*(2*Power(s1,4) + 
                  Power(s1,3)*(-8 + 11*s2 - 8*t1 + t2) - 
                  2*(2 + 3*Power(s2,2) + 3*Power(t1,2) + 
                     3*t1*(-2 + t2) + t2 - Power(t2,2) - 
                     s2*(1 + 6*t1 + t2)) + 
                  Power(s1,2)*
                   (-23 + 10*Power(s2,2) + 4*Power(t1,2) + 
                     t1*(23 - 3*t2) + Power(t2,2) + 
                     s2*(-25 - 14*t1 + 3*t2)) + 
                  s1*(35 + 8*Power(s2,2) + 16*Power(t1,2) - 13*t2 + 
                     Power(t2,2) - t1*(7 + 9*t2) + 
                     s2*(-2 - 24*t1 + 15*t2)))))/
           ((-1 + s2)*(-s + s2 - t1)*
             (-1 + s - s*s1 + s1*(s2 - t1) + t2)) + 
          (8*(32*(-1 + s1)*(-1 + s + t1)*(-1 + s2 - t1 + t2) + 
               ((13 - 6*Power(s1,2) + 2*s*(-5 + 3*s1) + 3*s2 - t1 - 
                    t2 + 2*s1*(5 - 2*s2 + t2))*
                  (2*Power(s,2)*(-1 + s1) - 
                    s*(-4 + 2*Power(s1,2) + t1 + 
                       s1*(2*s2 - t1 - t2) + t2) + 
                    (1 + s1)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2)))/
                (1 + s*(-1 + s1) + s1*(-s2 + t1) - t2) + 
               16*(6 - 2*Power(s1,3) - 6*s2 - Power(s2,2) + t1 + 
                  9*s2*t1 + Power(s2,2)*t1 - 9*Power(t1,2) - 
                  s2*Power(t1,2) - 8*t2 - s2*t2 - Power(s2,2)*t2 + 
                  7*t1*t2 + 2*s2*t1*t2 + 2*Power(t2,2) - 
                  s2*Power(t2,2) + Power(s1,2)*(3 - 2*t1 + 3*t2) + 
                  s1*(Power(s2,2) + 10*Power(t1,2) + t1*(2 - 9*t2) - 
                     Power(-3 + t2,2) + s2*(8 - 11*t1 + t2)) + 
                  s*(-9 + 2*Power(s1,2) - 13*t1 + Power(t1,2) + 
                     s1*(9 - 13*s2 + 14*t1 - 14*t2) + 13*t2 - 
                     2*t1*t2 + Power(t2,2) + s2*(11 - t1 + t2))) + 
               8*(-18 + 17*Power(s1,3) + 19*s2 + 4*Power(s2,2) - 
                  11*t1 - 35*s2*t1 - 4*Power(s2,2)*t1 + 
                  36*Power(t1,2) + 4*s2*Power(t1,2) + 
                  Power(s1,2)*(-17 + s2 + 16*t1 - 25*t2) + 
                  Power(s,2)*(s1 - 2*s2 + 4*t1 - 4*t2) + 38*t2 + 
                  10*s2*t2 + 5*Power(s2,2)*t2 - 25*t1*t2 - 
                  10*s2*t1*t2 - 14*Power(t2,2) + 6*s2*Power(t2,2) - 
                  s*(-33 + 18*Power(s1,2) - 2*Power(s2,2) - 60*t1 + 
                     2*Power(t1,2) + 
                     s1*(31 - 60*s2 + 71*t1 - 72*t2) + 63*t2 - 
                     8*t1*t2 + 6*Power(t2,2) + s2*(42 + t1 + t2)) - 
                  s1*(-28 + 7*Power(s2,2) + 41*Power(t1,2) + 
                     t1*(10 - 33*t2) + 17*t2 - 8*Power(t2,2) + 
                     s2*(30 - 48*t1 + 9*t2))) + 
               4*(40 + 2*Power(s,3) - 48*Power(s1,3) - 23*s2 - 
                  9*Power(s2,2) + 21*t1 + 64*s2*t1 + 
                  5*Power(s2,2)*t1 - 63*Power(t1,2) - 
                  5*s2*Power(t1,2) - 84*t2 - 29*s2*t2 - 
                  8*Power(s2,2)*t2 + 43*t1*t2 + 16*s2*t1*t2 + 
                  34*Power(t2,2) - 12*s2*Power(t2,2) + 
                  Power(s,2)*(-8*s1 + 4*s2 - 11*t1 + 15*t2) + 
                  Power(s1,2)*(23 - 6*s2 - 44*t1 + 69*t2) + 
                  s*(-71 + 54*Power(s1,2) - 6*Power(s2,2) - 116*t1 - 
                     Power(t1,2) + 
                     s1*(55 - 113*s2 + 148*t1 - 155*t2) + 
                     s2*(73 + 10*t1 - 7*t2) + 128*t2 - 9*t1*t2 + 
                     12*Power(t2,2)) + 
                  s1*(-34 + 15*Power(s2,2) + 74*Power(t1,2) + 
                     t1*(15 - 53*t2) + 36*t2 - 22*Power(t2,2) + 
                     s2*(52 - 91*t1 + 24*t2))) - 
               2*(28 + 4*Power(s,3) - 38*Power(s1,3) - 19*s2 - 
                  11*Power(s2,2) + 23*t1 + 48*s2*t1 + 
                  2*Power(s2,2)*t1 - 41*Power(t1,2) - 
                  2*s2*Power(t1,2) - 
                  Power(s,2)*(9 + 6*s1 + 6*t1 - 14*t2) - 64*t2 - 
                  27*s2*t2 - 4*Power(s2,2)*t2 + 27*t1*t2 + 
                  8*s2*t1*t2 + 28*Power(t2,2) - 8*s2*Power(t2,2) + 
                  Power(s1,2)*(-6 - 40*t1 + 58*t2) + 
                  2*s1*(-21 + 6*Power(s2,2) + 24*Power(t1,2) + 
                     t1*(4 - 15*t2) + 20*t2 - 10*Power(t2,2) + 
                     s2*(14 - 31*t1 + 8*t2)) + 
                  s*(-36 + 40*Power(s1,2) - 4*Power(s2,2) - 87*t1 - 
                     2*Power(t1,2) + s2*(57 + 8*t1 - 10*t2) + 92*t2 - 
                     2*t1*t2 + 8*Power(t2,2) - 
                     2*s1*(-29 + 41*s2 - 54*t1 + 57*t2))) - 
               ((-2 + s1)*(4*(-1 + s2 - t1 + t2)*
                     Power(2*Power(s,2)*(-1 + s1) - 
                       s*(-4 + 2*Power(s1,2) + t1 + 
                       s1*(2*s2 - t1 - t2) + t2) + 
                       (1 + s1)*
                       (-1 + s1*(s2 - t1) + t1 + t2 - s2*t2),2) - 
                    2*(6*Power(s,4)*Power(-1 + s1,2)*
                        (-1 + s2 - t1 + t2) + 
                       Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2)*
                        (-3 + s2 - t1 + 3*t2 + 
                        2*s1*(-1 + 2*s2 - 2*t1 + t2) + 
                        Power(s1,2)*(-1 + s2 - t1 + t2)) - 
                       2*Power(s,3)*(-1 + s1)*
                        (12 + 9*t1 - 2*Power(t1,2) - 15*t2 - 
                        2*t1*t2 + 4*Power(t2,2) + 
                        6*Power(s1,2)*(-1 + s2 - t1 + t2) + 
                        3*s2*(-4 + t1 + t2) + 
                        3*s1*
                        (2*Power(s2,2) + t1 + Power(t1,2) + t2 - 
                        Power(t2,2) + s2*(-2 - 3*t1 + t2))) + 
                       Power(s,2)*
                        (-30 + 30*s2 - 8*t1 - 18*s2*t1 + 
                        11*Power(t1,2) + s2*Power(t1,2) - 
                        Power(t1,3) + 52*t2 - 32*s2*t2 + 
                        6*Power(s2,2)*t2 + 8*t1*t2 + 2*s2*t1*t2 - 
                        Power(t1,2)*t2 - 25*Power(t2,2) + 
                        3*s2*Power(t2,2) - t1*Power(t2,2) + 
                        3*Power(t2,3) + 
                        6*Power(s1,4)*(-1 + s2 - t1 + t2) + 
                        6*Power(s1,3)*
                       (3*Power(s2,2) + 2*t1 + 2*Power(t1,2) + t2 - 
                       t1*t2 - Power(t2,2) + s2*(-3 - 5*t1 + 2*t2)) \
+ 2*s1*(Power(t1,3) + 3*Power(s2,2)*(-5 + t1 + t2) + 
                       Power(t1,2)*(-7 + 2*t2) - 
                       2*t2*(4 - 5*t2 + Power(t2,2)) - 
                       t1*(15 - 3*t2 + Power(t2,2)) + 
                       s2*(15 - 4*Power(t1,2) + t1*(23 - 6*t2) - 
                       6*t2 + 4*Power(t2,2))) + 
                        Power(s1,2)*
                        (30 + 6*Power(s2,3) - Power(t1,3) - 42*t2 + 
                        11*Power(t2,2) + Power(t2,3) - 
                        3*Power(t1,2)*(3 + t2) - 
                        6*Power(s2,2)*(1 + 2*t1 + t2) + 
                        t1*(26 - 8*t2 + 3*Power(t2,2)) + 
                        s2*(-30 + 7*Power(t1,2) + 20*t2 - 
                        11*Power(t2,2) + 2*t1*(7 + 5*t2)))) - 
                       2*s*(-7 + 6*s2 + 5*t1 - 7*s2*t1 + 
                        3*Power(t1,2) + s2*Power(t1,2) - 
                        Power(t1,3) + 17*t2 - 18*s2*t2 + 
                        5*Power(s2,2)*t2 - 4*t1*t2 + 4*s2*t1*t2 - 
                        Power(s2,2)*t1*t2 - Power(t1,2)*t2 + 
                        s2*Power(t1,2)*t2 - 13*Power(t2,2) + 
                        12*s2*Power(t2,2) - 
                        3*Power(s2,2)*Power(t2,2) - t1*Power(t2,2) - 
                        s2*t1*Power(t2,2) + 3*Power(t2,3) + 
                        3*Power(s1,4)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                        Power(s1,3)*
                        (3*Power(s2,3) - Power(t1,3) - 
                        Power(t1,2)*(-1 + t2) + 3*Power(-1 + t2,2) - 
                        Power(s2,2)*(-1 + 7*t1 + t2) + 
                        s2*(-6 + 5*Power(t1,2) + 2*t1*(-1 + t2) + 
                       10*t2 - 4*Power(t2,2)) + 
                        2*t1*(2 - 3*t2 + Power(t2,2))) + 
                        s1*(Power(t1,3) - 3*Power(s2,3)*t2 + 
                        2*Power(t1,2)*(3 + t2) - 
                        Power(-1 + t2,2)*(5 + t2) - 
                        2*t1*(7 - 8*t2 + Power(t2,2)) + 
                        s2*(17 + t1*(2 - 18*t2) + 
                       Power(t1,2)*(-6 + t2) - 23*t2 + 
                       7*Power(t2,2) - Power(t2,3)) + 
                        Power(s2,2)*
                        (-9 + 18*t2 - Power(t2,2) + t1*(5 + 2*t2))) + 
                        Power(s1,2)*
                        (-13*Power(t1,2) + Power(t1,3) - 
                        3*Power(s2,3)*(-1 + t2) - 
                        Power(-1 + t2,2)*(-3 + 2*t2) + 
                        t1*(-4 + 3*t2 + Power(t2,2)) + 
                        Power(s2,2)*
                        (-13 + 2*t2 - 2*Power(t2,2) + t1*(-4 + 5*t2)) \
+ s2*(4 - 2*t2 - 2*Power(t1,2)*t2 - 3*Power(t2,2) + Power(t2,3) + 
                        t1*(25 + Power(t2,2))))))))/
                (Power(-1 + s - s*s1 + s1*(s2 - t1) + t2,2)*
                  (-1 + s2 - t1 + t2))))/((-1 + s2)*(-1 + t1)) + 
          (4*(-64*s*(-1 + s1)*(-1 + s2 - t1 + t2) - 
               (2*(-15 + 3*Power(s,2) + 3*Power(s1,2) + s2 + 
                    Power(s2,2) - 3*t1 - s2*t1 + 
                    s1*(-7 + 4*s2 - 2*t1 - t2) + 4*t2 - s2*t2 + 
                    t1*t2 + s*(6 - 6*s1 - 4*s2 + 3*t1 + t2))*
                  (2*Power(s,2)*(-1 + s1) - 
                    s*(-4 + 2*Power(s1,2) + t1 + 
                       s1*(2*s2 - t1 - t2) + t2) + 
                    (1 + s1)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2)))/
                (1 + s*(-1 + s1) + s1*(-s2 + t1) - t2) + 
               4*(32 + Power(s,3) + 3*Power(s1,3) + 6*s2 - 
                  3*Power(s2,2) - 5*t1 + 13*s2*t1 - 9*Power(t1,2) - 
                  63*t2 - 2*s2*t2 - 3*Power(s2,2)*t2 + 18*t1*t2 + 
                  3*s2*t1*t2 + 22*Power(t2,2) - 6*s2*Power(t2,2) + 
                  Power(s,2)*(-31 + 41*s1 + s2 - t1 + 8*t2) + 
                  Power(s1,2)*(7 + 37*s2 - 33*t1 + 12*t2) + 
                  s*(-59 - 45*Power(s1,2) - 2*Power(s2,2) - 62*t1 - 
                     2*Power(t1,2) + 
                     s1*(67 - 120*s2 + 96*t1 - 62*t2) + 
                     s2*(66 + 4*t1 - 5*t2) + 56*t2 + 2*t1*t2 + 
                     6*Power(t2,2)) + 
                  s1*(-7 + 10*Power(s2,2) + 12*Power(t1,2) + 
                     s2*(14 - 21*t1 - 25*t2) + 36*t2 - 
                     16*Power(t2,2) + t1*(7 + t2))) - 
               32*(4 + 2*Power(s,2)*(-1 + s1) - s2 - t1 + s2*t1 - 
                  Power(t1,2) - 6*t2 + 2*s2*t2 - Power(s2,2)*t2 + 
                  t1*t2 + s2*t1*t2 + 2*Power(t2,2) - s2*Power(t2,2) + 
                  Power(s1,2)*(-1 + 2*s2 - 2*t1 + t2) + 
                  s*(-12 - 2*Power(s1,2) - 11*t1 + 
                     s1*(15 - 15*s2 + 14*t1 - 12*t2) + 10*t2 - 
                     t1*t2 + Power(t2,2) + s2*(12 + t2)) + 
                  s1*(-2 + Power(s2,2) + 2*t1 + Power(t1,2) + 3*t2 - 
                     Power(t2,2) - s2*(1 + 2*t1 + t2))) + 
               16*(21 + Power(s1,3) - 4*s2 - 2*Power(s2,2) - 5*t1 + 
                  7*s2*t1 - 5*Power(t1,2) + 
                  Power(s,2)*(17*s1 + 2*(-9 + t2)) - 34*t2 + 7*s2*t2 - 
                  4*Power(s2,2)*t2 + 8*t1*t2 + 4*s2*t1*t2 + 
                  12*Power(t2,2) - 5*s2*Power(t2,2) + 
                  Power(s1,2)*(-5 + 16*s2 - 15*t1 + 6*t2) + 
                  s*(-45 - 18*Power(s1,2) - 49*t1 - Power(t1,2) + 
                     s1*(71 - 75*s2 + 66*t1 - 51*t2) + 41*t2 - 
                     2*t1*t2 + 5*Power(t2,2) + s2*(55 + t1 + 2*t2)) + 
                  s1*(-17 + 6*Power(s2,2) + 7*Power(t1,2) - 
                     t1*(-11 + t2) + 23*t2 - 7*Power(t2,2) - 
                     s2*(3 + 13*t1 + 9*t2))) + 
               8*(-32 + Power(s,3) - 5*Power(s1,3) + 6*s2 + 
                  4*Power(s2,2) + 7*t1 - 16*s2*t1 + 11*Power(t1,2) + 
                  Power(s1,2)*(5 - 45*s2 + 39*t1 - 13*t2) + 
                  Power(s,2)*(45 - 51*s1 - 3*s2 + 2*t1 - 7*t2) + 
                  68*t2 - 5*s2*t2 + 6*Power(s2,2)*t2 - 21*t1*t2 - 
                  6*s2*t1*t2 - 27*Power(t2,2) + 9*s2*Power(t2,2) + 
                  s1*(46 - 14*Power(s2,2) - 16*Power(t1,2) + 
                     t1*(-17 + t2) - 54*t2 + 18*Power(t2,2) + 
                     s2*(-6 + 29*t1 + 28*t2)) + 
                  s*(64 + 55*Power(s1,2) + 2*Power(s2,2) + 91*t1 + 
                     3*Power(t1,2) - 75*t2 - 9*Power(t2,2) + 
                     s2*(-101 - 5*t1 + t2) + 
                     s1*(-128 + 163*s2 - 134*t1 + 91*t2))) + 
               ((3 + s - s1 - s2 + t1)*
                  (4*(-1 + s2 - t1 + t2)*
                     Power(2*Power(s,2)*(-1 + s1) - 
                       s*(-4 + 2*Power(s1,2) + t1 + 
                       s1*(2*s2 - t1 - t2) + t2) + 
                       (1 + s1)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2),
                      2) - 2*
                     (6*Power(s,4)*Power(-1 + s1,2)*
                        (-1 + s2 - t1 + t2) + 
                       Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2)*
                        (-3 + s2 - t1 + 3*t2 + 
                        2*s1*(-1 + 2*s2 - 2*t1 + t2) + 
                        Power(s1,2)*(-1 + s2 - t1 + t2)) - 
                       2*Power(s,3)*(-1 + s1)*
                        (12 + 9*t1 - 2*Power(t1,2) - 15*t2 - 
                        2*t1*t2 + 4*Power(t2,2) + 
                        6*Power(s1,2)*(-1 + s2 - t1 + t2) + 
                        3*s2*(-4 + t1 + t2) + 
                        3*s1*
                        (2*Power(s2,2) + t1 + Power(t1,2) + t2 - 
                        Power(t2,2) + s2*(-2 - 3*t1 + t2))) + 
                       Power(s,2)*
                        (-30 + 30*s2 - 8*t1 - 18*s2*t1 + 
                        11*Power(t1,2) + s2*Power(t1,2) - 
                        Power(t1,3) + 52*t2 - 32*s2*t2 + 
                        6*Power(s2,2)*t2 + 8*t1*t2 + 2*s2*t1*t2 - 
                        Power(t1,2)*t2 - 25*Power(t2,2) + 
                        3*s2*Power(t2,2) - t1*Power(t2,2) + 
                        3*Power(t2,3) + 
                        6*Power(s1,4)*(-1 + s2 - t1 + t2) + 
                        6*Power(s1,3)*
                        (3*Power(s2,2) + 2*t1 + 2*Power(t1,2) + t2 - 
                        t1*t2 - Power(t2,2) + s2*(-3 - 5*t1 + 2*t2)) \
+ 2*s1*(Power(t1,3) + 3*Power(s2,2)*(-5 + t1 + t2) + 
                        Power(t1,2)*(-7 + 2*t2) - 
                        2*t2*(4 - 5*t2 + Power(t2,2)) - 
                        t1*(15 - 3*t2 + Power(t2,2)) + 
                        s2*(15 - 4*Power(t1,2) + t1*(23 - 6*t2) - 
                        6*t2 + 4*Power(t2,2))) + 
                        Power(s1,2)*
                        (30 + 6*Power(s2,3) - Power(t1,3) - 42*t2 + 
                        11*Power(t2,2) + Power(t2,3) - 
                        3*Power(t1,2)*(3 + t2) - 
                        6*Power(s2,2)*(1 + 2*t1 + t2) + 
                        t1*(26 - 8*t2 + 3*Power(t2,2)) + 
                        s2*(-30 + 7*Power(t1,2) + 20*t2 - 
                        11*Power(t2,2) + 2*t1*(7 + 5*t2)))) - 
                       2*s*(-7 + 6*s2 + 5*t1 - 7*s2*t1 + 
                        3*Power(t1,2) + s2*Power(t1,2) - Power(t1,3) + 
                        17*t2 - 18*s2*t2 + 5*Power(s2,2)*t2 - 
                        4*t1*t2 + 4*s2*t1*t2 - Power(s2,2)*t1*t2 - 
                        Power(t1,2)*t2 + s2*Power(t1,2)*t2 - 
                        13*Power(t2,2) + 12*s2*Power(t2,2) - 
                        3*Power(s2,2)*Power(t2,2) - t1*Power(t2,2) - 
                        s2*t1*Power(t2,2) + 3*Power(t2,3) + 
                        3*Power(s1,4)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                        Power(s1,3)*
                        (3*Power(s2,3) - Power(t1,3) - 
                        Power(t1,2)*(-1 + t2) + 3*Power(-1 + t2,2) - 
                        Power(s2,2)*(-1 + 7*t1 + t2) + 
                        s2*(-6 + 5*Power(t1,2) + 2*t1*(-1 + t2) + 
                        10*t2 - 4*Power(t2,2)) + 
                        2*t1*(2 - 3*t2 + Power(t2,2))) + 
                        s1*(Power(t1,3) - 3*Power(s2,3)*t2 + 
                        2*Power(t1,2)*(3 + t2) - 
                        Power(-1 + t2,2)*(5 + t2) - 
                        2*t1*(7 - 8*t2 + Power(t2,2)) + 
                        s2*(17 + t1*(2 - 18*t2) + 
                        Power(t1,2)*(-6 + t2) - 23*t2 + 
                        7*Power(t2,2) - Power(t2,3)) + 
                        Power(s2,2)*
                        (-9 + 18*t2 - Power(t2,2) + t1*(5 + 2*t2))) + 
                        Power(s1,2)*
                        (-13*Power(t1,2) + Power(t1,3) - 
                        3*Power(s2,3)*(-1 + t2) - 
                        Power(-1 + t2,2)*(-3 + 2*t2) + 
                        t1*(-4 + 3*t2 + Power(t2,2)) + 
                        Power(s2,2)*
                        (-13 + 2*t2 - 2*Power(t2,2) + t1*(-4 + 5*t2)) \
+ s2*(4 - 2*t2 - 2*Power(t1,2)*t2 - 3*Power(t2,2) + Power(t2,3) + 
                        t1*(25 + Power(t2,2))))))))/
                (Power(-1 + s - s*s1 + s1*(s2 - t1) + t2,2)*
                  (-1 + s2 - t1 + t2))))/((-1 + t1)*(-1 + t2))) + 
       (((-64*(2*Power(s,2)*(-1 + s1) - 
                  s*(-4 + 2*Power(s1,2) + t1 + 
                     s1*(2*s2 - t1 - t2) + t2) + 
                  (1 + s1)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2))*
                (3 + 2*Power(s,2) - s2 + t1 + 3*t2 - 2*s2*t2 + 
                  2*t1*t2 - Power(s1,2)*(-1 + s2 - t1 + t2) - 
                  s1*(-2 + Power(s2,2) + t1 + Power(t1,2) - t1*t2 + 
                     s2*(-1 - 2*t1 + t2)) + 
                  s*(s1*(-1 + s2 - t1 + t2) + 2*(2 - s2 + t1 + t2))))/
              (1 + s*(-1 + s1) + s1*(-s2 + t1) - t2) + 
             128*(-18 + 2*Power(s,3) - 15*s2 - Power(s2,2) + 19*t1 + 
                2*s2*t1 - Power(t1,2) - 5*t2 + 4*s2*t2 - 7*t1*t2 - 
                s2*t1*t2 + Power(t1,2)*t2 + 7*Power(t2,2) - 
                3*t1*Power(t2,2) + 2*Power(t2,3) + 
                Power(s1,3)*(-1 + s2 - t1 + t2) + 
                Power(s1,2)*(-11 + Power(s2,2) + 4*t1 + Power(t1,2) - 
                   2*s2*(1 + t1) + 2*t2 - Power(t2,2)) + 
                Power(s,2)*(2 - 2*s2 + 4*t2 + 
                   s1*(-3 + s2 - t1 + t2)) - 
                s*(-4 - Power(t1,2) - 4*t2 + 2*s2*t2 - 
                   4*Power(t2,2) + 
                   2*Power(s1,2)*(-1 + s2 - t1 + t2) + 
                   t1*(3 + s2 + 3*t2) + 
                   s1*(-6 + s2 + Power(s2,2) + t1 - 2*s2*t1 + 
                      Power(t1,2) + 9*t2 - Power(t2,2))) - 
                s1*(7 + Power(s2,2)*(-1 + t2) - 4*t2 + 
                   Power(t1,2)*t2 + 5*Power(t2,2) - 
                   t1*(5 - 2*t2 + Power(t2,2)) + 
                   s2*(4 + t1 - t2 - 2*t1*t2 + Power(t2,2)))))/
           ((-1 + s1)*(-1 + t2)) + 
          (64*(((2*Power(s,2)*(-1 + s1) - 
                    s*(-4 + 2*Power(s1,2) + t1 + 
                       s1*(2*s2 - t1 - t2) + t2) + 
                    (1 + s1)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2))*
                  (-6 - s2 + Power(s2,2) + t1 - 2*s2*t1 + 
                    Power(t1,2) + s*(-3 + 2*s1 - s2 + t1 - t2) - 
                    4*t2 + s2*t2 - t1*t2 - 
                    Power(s1,2)*(-1 + s2 - t1 + t2) + 
                    s1*(Power(s2,2) + t1 + Power(t1,2) + t2 - 
                       2*t1*t2 + Power(t2,2) + s2*(-1 - 2*t1 + 2*t2))\
))/(1 + s*(-1 + s1) + s1*(-s2 + t1) - t2) - 
               2*(9*s2 - 6*Power(s2,2) - 2*Power(s2,3) - 8*t1 + 
                  10*s2*t1 + 8*Power(s2,2)*t1 - 4*Power(t1,2) - 
                  10*s2*Power(t1,2) + 4*Power(t1,3) + 
                  Power(s,2)*(-3 + 2*s1 - s2 + t1 - t2) + 10*t2 - 
                  5*s2*t2 - 4*Power(s2,2)*t2 + 4*t1*t2 + 10*s2*t1*t2 - 
                  6*Power(t1,2)*t2 - 2*s2*Power(t2,2) + 
                  2*t1*Power(t2,2) + Power(s1,3)*(-1 + s2 - t1 + t2) - 
                  Power(s1,2)*
                   (-5 + t1 + Power(t1,2) + 2*t2 - 2*t1*t2 + 
                     Power(t2,2) + s2*(4 - t1 + t2)) + 
                  s*(-1 + 3*Power(s2,2) - 3*t1 + 5*Power(t1,2) - 
                     7*t1*t2 + 2*Power(t2,2) - 
                     Power(s1,2)*(1 + s2 - t1 + t2) + 
                     s2*(4 - 8*t1 + 5*t2) + 
                     s1*(-7 + Power(s2,2) + 3*t1 - 2*s2*t1 + 
                        Power(t1,2) + t2 + 2*s2*t2 - 2*t1*t2 + 
                        Power(t2,2))) + 
                  s1*(24 + 2*Power(s2,2) + 5*Power(t1,2) + 7*t2 + 
                     Power(t2,2) - 4*t1*(5 + t2) + 
                     s2*(-7*t1 + 4*(3 + t2))))))/
           ((-1 + s2)*(-s + s2 - t1)) + 
          (64*(((2 - Power(s1,2) + s*(3 + s1) - 3*s2 + 3*t1 + 
                    s1*(-2 + t2))*(s - s1 + t2)*
                  (2*Power(s,2)*(-1 + s1) - 
                    s*(-4 + 2*Power(s1,2) + t1 + 
                       s1*(2*s2 - t1 - t2) + t2) + 
                    (1 + s1)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2)))/
                (1 + s*(-1 + s1) + s1*(-s2 + t1) - t2) - 
               2*(12 - Power(s1,4) + Power(s,3)*(3 + s1) + 14*s2 + 
                  2*Power(s2,2) - 14*t1 - 4*s2*t1 + 2*Power(t1,2) - 
                  Power(s1,3)*(s2 - 2*t2) + 11*t2 + s2*t2 - 4*t1*t2 - 
                  s2*t1*t2 + Power(t1,2)*t2 + 5*Power(t2,2) - 
                  3*t1*Power(t2,2) + 2*Power(t2,3) + 
                  Power(s,2)*
                   (3 - 3*Power(s1,2) - 3*s2 + t1 - 
                     s1*(5 + t1 - 3*t2) + 5*t2) + 
                  Power(s1,2)*
                   (13 + Power(t1,2) + 5*t2 - Power(t2,2) - 
                     t1*(4 + t2) + s2*(3 - t1 + 2*t2)) + 
                  s*(-3 + 3*Power(s1,3) + Power(t1,2) + 
                     Power(s1,2)*(2 + s2 + t1 - 5*t2) + 6*t2 - 
                     2*t1*t2 + 4*Power(t2,2) - s2*(3 + t1 + 3*t2) + 
                     s1*(-14 + t1 - Power(t1,2) + s2*(2 + t1 - t2) - 
                        7*t2 + 2*Power(t2,2))) + 
                  s1*(5 - 2*Power(s2,2) - 14*t2 - 7*Power(t2,2) - 
                     Power(t1,2)*(3 + t2) + t1*t2*(9 + t2) + 
                     s2*(3 - 5*t2 - Power(t2,2) + t1*(5 + t2))))))/
           ((-1 + s1)*(-s + s1 - t2)) + 
          (64*(-8 - 12*s2 + 12*Power(s2,2) + 12*t1 - 24*s2*t1 + 
               12*Power(t1,2) + Power(s1,5)*(-s2 + t1) - 19*t2 + 
               32*s2*t2 - 16*Power(s2,2)*t2 - 25*t1*t2 + 28*s2*t1*t2 - 
               12*Power(t1,2)*t2 + 27*Power(t2,2) - 
               23*s2*Power(t2,2) + 6*Power(s2,2)*Power(t2,2) + 
               18*t1*Power(t2,2) - 8*s2*t1*Power(t2,2) + 
               2*Power(t1,2)*Power(t2,2) - 2*t1*Power(t2,3) - 
               Power(s1,4)*(-1 - 2*Power(t1,2) + 
                  s2*(-4 + 2*t1 - 3*t2) + t2 + t1*(3 + 4*t2)) + 
               Power(s,3)*(8 - 4*s2 - 2*t1 + 
                  Power(s1,2)*(-4 + t1 - t2) + 2*t2 + 
                  s1*(-4 + 4*s2 - 5*t1 + 5*t2)) + 
               Power(s1,3)*(-2*Power(s2,2)*(-1 + t1) - 
                  2*Power(t1,3) + 
                  s2*(24 - 2*t1 + 4*Power(t1,2) - 5*t2 - 
                     3*Power(t2,2)) + 4*(-1 + Power(t2,2)) + 
                  t1*(-20 + t2 + 3*Power(t2,2))) + 
               Power(s,2)*(12 + 4*Power(s2,2) - 7*t1 + 29*t2 - 
                  2*t1*t2 + 2*Power(t2,2) + 
                  Power(s1,3)*(10 + s2 - 3*t1 + 4*t2) - 
                  2*s2*(-3 + 2*t1 + 6*t2) + 
                  Power(s1,2)*
                   (25 + t1 + 2*Power(t1,2) - 10*t2 - 
                     2*Power(t2,2) + s2*(5 - 2*t1 + t2)) - 
                  s1*(53 + 8*Power(s2,2) - 15*t1 + 8*Power(t1,2) + 
                     23*t2 - 2*t1*t2 - 6*Power(t2,2) + 
                     s2*(2 - 16*t1 + 5*t2))) + 
               Power(s1,2)*(-8*Power(s2,3) + 2*Power(t1,3)*(3 + t2) - 
                  2*Power(t1,2)*(-5 + 6*t2 + Power(t2,2)) + 
                  t1*(-29 + 13*t2 + 8*Power(t2,2)) - 
                  3*(8 - 9*t2 + Power(t2,3)) + 
                  2*Power(s2,2)*(10 - 9*t2 + t1*(11 + t2)) + 
                  s2*(27 - 14*t2 - 6*Power(t2,2) + Power(t2,3) - 
                     4*Power(t1,2)*(5 + t2) + 
                     2*t1*(-15 + 15*t2 + Power(t2,2)))) + 
               s1*(-29 - 2*Power(t1,3)*(-6 + t2) + 
                  4*Power(s2,3)*(-3 + t2) - 
                  16*Power(t1,2)*(-1 + t2) + 45*t2 - 7*Power(t2,2) - 
                  9*Power(t2,3) + 
                  Power(s2,2)*
                   (20 + t1*(36 - 10*t2) - 30*t2 + 4*Power(t2,2)) + 
                  t1*(5 - 57*t2 + 13*Power(t2,2) + 2*Power(t2,3)) - 
                  s2*(12 + Power(t1,2)*(36 - 8*t2) - 64*t2 + 
                     14*Power(t2,2) + Power(t2,3) + 
                     t1*(36 - 46*t2 + 4*Power(t2,2)))) + 
               s*(-7 - t1 - 12*Power(t1,2) + 27*t2 + 11*t1*t2 + 
                  2*Power(t1,2)*t2 + 21*Power(t2,2) - 
                  2*t1*Power(t2,2) + 2*Power(s2,2)*(-8 + 5*t2) + 
                  Power(s1,3)*
                   (-18 + 3*t1 - 4*Power(t1,2) + 
                     s2*(-9 + 4*t1 - 2*t2) + 12*t2 + 2*t1*t2 + 
                     4*Power(t2,2)) + Power(s1,4)*(t1 - 3*(2 + t2)) + 
                  s2*(8 - 17*t2 - 8*Power(t2,2) - 4*t1*(-7 + 3*t2)) + 
                  Power(s1,2)*
                   (1 + 8*Power(t1,2) + 2*Power(t1,3) + 
                     2*Power(s2,2)*(3 + t1) + 9*t2 - Power(t2,2) - 
                     Power(t2,3) - t1*(-39 + 13*t2 + Power(t2,2)) + 
                     s2*(-50 - 14*t1 - 4*Power(t1,2) + 23*t2 + 
                        2*Power(t2,2))) + 
                  s1*(52 + 4*Power(s2,3) - 2*Power(t1,3) + 
                     Power(t1,2)*(22 - 8*t2) - 77*t2 - 
                     14*Power(t2,2) + Power(t2,3) - 
                     2*Power(s2,2)*(-4 + 5*t1 + 2*t2) + 
                     3*t1*(-8 - 4*t2 + 3*Power(t2,2)) + 
                     s2*(27 + 8*Power(t1,2) + 24*t2 - 10*Power(t2,2) + 
                        6*t1*(-5 + 2*t2))))))/
           ((s - s2 + t1)*(1 + s*(-1 + s1) + s1*(-s2 + t1) - t2)*
             (s - s1 + t2)) + 
          (32*((-2*(2*Power(s,2)*(-1 + s1) - 
                    s*(-4 + 2*Power(s1,2) + t1 + 
                       s1*(2*s2 - t1 - t2) + t2) + 
                    (1 + s1)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2))*
                  (18 + 2*Power(s,2)*(-1 + s1) + 2*Power(s1,3) + 
                    10*s2 - Power(s2,2) - 12*t1 - s2*t1 + 
                    2*Power(t1,2) + 
                    Power(s1,2)*(-7 + s2 + t1 - 3*t2) + 5*t2 + 
                    2*s2*t2 - 2*t1*t2 + Power(t2,2) - 
                    s*(13 + 4*Power(s1,2) - 3*t1 + 
                       s1*(-9 + s2 + t1 - 3*t2) + 4*t2) + 
                    s1*(11 + 4*t2 + Power(t2,2) - s2*(3 + t2) - 
                       t1*(3 + t2))))/
                (1 + s*(-1 + s1) + s1*(-s2 + t1) - t2) + 
               4*(-12 + Power(s,3)*(-1 + s1) - Power(s1,4) - 25*s2 - 
                  5*Power(s2,2) + 39*t1 + 26*s2*t1 + 
                  2*Power(s2,2)*t1 - 23*Power(t1,2) - 
                  6*s2*Power(t1,2) + 4*Power(t1,3) - 
                  Power(s,2)*
                   (3*Power(s1,2) + s1*(-5 + s2 + t1 - 2*t2) - 
                     2*(-4 + s2 + t1 - t2)) - 
                  Power(s1,3)*(-3 + s2 + t1 - 2*t2) - 15*t2 - 
                  17*s2*t2 - Power(s2,2)*t2 + 26*t1*t2 + 5*s2*t1*t2 - 
                  6*Power(t1,2)*t2 - 7*Power(t2,2) - s2*Power(t2,2) + 
                  2*t1*Power(t2,2) - 
                  Power(s1,2)*
                   (6 + s2*(-3 + t1 - 2*t2) + 3*t2 + Power(t2,2) - 
                     t1*(4 + t2)) + 
                  s1*(-21 - Power(t1,2) - t1*(-4 + t2) + t2 + 
                     s2*(1 + t1 - 3*t2 + t1*t2 - Power(t2,2))) + 
                  s*(29 + 3*Power(s1,3) + 14*s2 - Power(s2,2) - 
                     25*t1 - 5*s2*t1 + 6*Power(t1,2) + 
                     Power(s1,2)*(-7 + 2*s2 + 2*t1 - 4*t2) + 15*t2 + 
                     5*s2*t2 - 8*t1*t2 + 3*Power(t2,2) + 
                     s1*(12 + s2*(-7 + t1 - 2*t2) + 3*t2 + 
                        Power(t2,2) - t1*(4 + t2)))) + 
               ((-5 + s*(-1 + s1) - Power(s1,2) - 2*s2 + 2*t1 - 
                    3*t2 + s1*(4 + t2))*
                  (4*(-1 + s2 - t1 + t2)*
                     Power(2*Power(s,2)*(-1 + s1) - 
                       s*(-4 + 2*Power(s1,2) + t1 + 
                       s1*(2*s2 - t1 - t2) + t2) + 
                       (1 + s1)*
                       (-1 + s1*(s2 - t1) + t1 + t2 - s2*t2),2) - 
                    2*(6*Power(s,4)*Power(-1 + s1,2)*
                        (-1 + s2 - t1 + t2) + 
                       Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2)*
                        (-3 + s2 - t1 + 3*t2 + 
                        2*s1*(-1 + 2*s2 - 2*t1 + t2) + 
                        Power(s1,2)*(-1 + s2 - t1 + t2)) - 
                       2*Power(s,3)*(-1 + s1)*
                        (12 + 9*t1 - 2*Power(t1,2) - 15*t2 - 
                        2*t1*t2 + 4*Power(t2,2) + 
                        6*Power(s1,2)*(-1 + s2 - t1 + t2) + 
                        3*s2*(-4 + t1 + t2) + 
                        3*s1*
                        (2*Power(s2,2) + t1 + Power(t1,2) + t2 - 
                        Power(t2,2) + s2*(-2 - 3*t1 + t2))) + 
                       Power(s,2)*
                        (-30 + 30*s2 - 8*t1 - 18*s2*t1 + 
                        11*Power(t1,2) + s2*Power(t1,2) - 
                        Power(t1,3) + 52*t2 - 32*s2*t2 + 
                        6*Power(s2,2)*t2 + 8*t1*t2 + 2*s2*t1*t2 - 
                        Power(t1,2)*t2 - 25*Power(t2,2) + 
                        3*s2*Power(t2,2) - t1*Power(t2,2) + 
                        3*Power(t2,3) + 
                        6*Power(s1,4)*(-1 + s2 - t1 + t2) + 
                        6*Power(s1,3)*
                       (3*Power(s2,2) + 2*t1 + 2*Power(t1,2) + t2 - 
                       t1*t2 - Power(t2,2) + s2*(-3 - 5*t1 + 2*t2)) \
+ 2*s1*(Power(t1,3) + 3*Power(s2,2)*(-5 + t1 + t2) + 
                       Power(t1,2)*(-7 + 2*t2) - 
                       2*t2*(4 - 5*t2 + Power(t2,2)) - 
                       t1*(15 - 3*t2 + Power(t2,2)) + 
                       s2*(15 - 4*Power(t1,2) + t1*(23 - 6*t2) - 
                       6*t2 + 4*Power(t2,2))) + 
                        Power(s1,2)*
                        (30 + 6*Power(s2,3) - Power(t1,3) - 42*t2 + 
                        11*Power(t2,2) + Power(t2,3) - 
                        3*Power(t1,2)*(3 + t2) - 
                        6*Power(s2,2)*(1 + 2*t1 + t2) + 
                        t1*(26 - 8*t2 + 3*Power(t2,2)) + 
                        s2*(-30 + 7*Power(t1,2) + 20*t2 - 
                        11*Power(t2,2) + 2*t1*(7 + 5*t2)))) - 
                       2*s*(-7 + 6*s2 + 5*t1 - 7*s2*t1 + 
                        3*Power(t1,2) + s2*Power(t1,2) - 
                        Power(t1,3) + 17*t2 - 18*s2*t2 + 
                        5*Power(s2,2)*t2 - 4*t1*t2 + 4*s2*t1*t2 - 
                        Power(s2,2)*t1*t2 - Power(t1,2)*t2 + 
                        s2*Power(t1,2)*t2 - 13*Power(t2,2) + 
                        12*s2*Power(t2,2) - 
                        3*Power(s2,2)*Power(t2,2) - t1*Power(t2,2) - 
                        s2*t1*Power(t2,2) + 3*Power(t2,3) + 
                        3*Power(s1,4)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                        Power(s1,3)*
                        (3*Power(s2,3) - Power(t1,3) - 
                        Power(t1,2)*(-1 + t2) + 3*Power(-1 + t2,2) - 
                        Power(s2,2)*(-1 + 7*t1 + t2) + 
                        s2*(-6 + 5*Power(t1,2) + 2*t1*(-1 + t2) + 
                       10*t2 - 4*Power(t2,2)) + 
                        2*t1*(2 - 3*t2 + Power(t2,2))) + 
                        s1*(Power(t1,3) - 3*Power(s2,3)*t2 + 
                        2*Power(t1,2)*(3 + t2) - 
                        Power(-1 + t2,2)*(5 + t2) - 
                        2*t1*(7 - 8*t2 + Power(t2,2)) + 
                        s2*(17 + t1*(2 - 18*t2) + 
                       Power(t1,2)*(-6 + t2) - 23*t2 + 
                       7*Power(t2,2) - Power(t2,3)) + 
                        Power(s2,2)*
                        (-9 + 18*t2 - Power(t2,2) + t1*(5 + 2*t2))) + 
                        Power(s1,2)*
                        (-13*Power(t1,2) + Power(t1,3) - 
                        3*Power(s2,3)*(-1 + t2) - 
                        Power(-1 + t2,2)*(-3 + 2*t2) + 
                        t1*(-4 + 3*t2 + Power(t2,2)) + 
                        Power(s2,2)*
                        (-13 + 2*t2 - 2*Power(t2,2) + t1*(-4 + 5*t2)) \
+ s2*(4 - 2*t2 - 2*Power(t1,2)*t2 - 3*Power(t2,2) + Power(t2,3) + 
                        t1*(25 + Power(t2,2))))))))/
                (Power(-1 + s - s*s1 + s1*(s2 - t1) + t2,2)*
                  (-1 + s2 - t1 + t2))))/((-1 + s2)*(-1 + t1)) + 
          (32*((2*(2*Power(s,2)*(-1 + s1) - 
                    s*(-4 + 2*Power(s1,2) + t1 + 
                       s1*(2*s2 - t1 - t2) + t2) + 
                    (1 + s1)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2))*
                  (-16 + 2*Power(s,3) - 6*s2 + Power(s2,2) + 9*t1 - 
                    2*s2*t1 + Power(t1,2) - 3*t2 + 4*s2*t2 - 
                    6*t1*t2 + s2*t1*t2 - Power(t1,2)*t2 + 
                    2*Power(t2,2) - s2*Power(t2,2) + t1*Power(t2,2) + 
                    Power(s1,2)*(s2 - t1 + 2*t2) + 
                    Power(s,2)*(2 - 4*s1 - 2*s2 + t1 + 3*t2) + 
                    s*(4 + 2*Power(s1,2) + s2 - 3*t1 + s2*t1 - 
                       Power(t1,2) + s1*(-2 + s2 - 5*t2) + 2*t2 - 
                       3*s2*t2 + 2*t1*t2 + Power(t2,2)) + 
                    s1*(-10 + Power(s2,2) + 3*t1 + Power(t1,2) - 
                       2*t2 - t1*t2 - Power(t2,2) + 
                       s2*(-1 - 2*t1 + 2*t2))))/
                (1 + s*(-1 + s1) + s1*(-s2 + t1) - t2) - 
               4*(18 + Power(s,4) + 14*s2 - Power(s2,2) - 23*t1 - 
                  3*s2*t1 + 5*Power(t1,2) - 
                  Power(s,3)*(3*s1 + s2 - 2*t2) - 2*t2 - 7*s2*t2 + 
                  Power(s2,2)*t2 + 10*t1*t2 - s2*t1*t2 - 
                  Power(t1,2)*t2 - 8*Power(t2,2) + 2*s2*Power(t2,2) + 
                  t1*Power(t2,2) - Power(s1,3)*(s2 - t1 + t2) - 
                  Power(s1,2)*
                   (-7 + s2 + Power(s2,2) + 2*t1 - 2*s2*t1 + 
                     Power(t1,2) + 2*t2 - Power(t2,2)) + 
                  Power(s,2)*
                   (2 + 3*Power(s1,2) - 3*t1 - Power(t1,2) + 
                     s1*(s2 + t1 - 5*t2) + s2*(1 + t1 - 2*t2) - 
                     2*t2 + t1*t2 + Power(t2,2)) + 
                  s1*(23 + Power(s2,2)*(-2 + t2) + 6*t2 + 
                     Power(t1,2)*t2 + 2*Power(t2,2) - 
                     t1*(19 - 3*t2 + Power(t2,2)) + 
                     s2*(12 - 2*t1*(-1 + t2) - 2*t2 + Power(t2,2))) + 
                  s*(-21 - Power(s1,3) + Power(s2,2) + 16*t1 + 
                     Power(t1,2) - 4*t2 - 8*t1*t2 - Power(t1,2)*t2 + 
                     2*Power(t2,2) + t1*Power(t2,2) + 
                     Power(s1,2)*(s2 - 2*t1 + 4*t2) + 
                     s2*(-8 + t1*(-3 + t2) + 7*t2 - Power(t2,2)) + 
                     s1*(-11 + Power(s2,2) + 7*t1 + 2*Power(t1,2) + 
                        2*t2 - t1*t2 - 2*Power(t2,2) + 
                        s2*(-2 - 3*t1 + 2*t2)))) - 
               ((3 + Power(s,2) - (-2 + s1 + s2 - t1)*t2 + 
                    s*(2 - s1 - s2 + t1 + t2))*
                  (4*(-1 + s2 - t1 + t2)*
                     Power(2*Power(s,2)*(-1 + s1) - 
                       s*(-4 + 2*Power(s1,2) + t1 + 
                       s1*(2*s2 - t1 - t2) + t2) + 
                       (1 + s1)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2),
                      2) - 2*
                     (6*Power(s,4)*Power(-1 + s1,2)*
                        (-1 + s2 - t1 + t2) + 
                       Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2)*
                        (-3 + s2 - t1 + 3*t2 + 
                        2*s1*(-1 + 2*s2 - 2*t1 + t2) + 
                        Power(s1,2)*(-1 + s2 - t1 + t2)) - 
                       2*Power(s,3)*(-1 + s1)*
                        (12 + 9*t1 - 2*Power(t1,2) - 15*t2 - 
                        2*t1*t2 + 4*Power(t2,2) + 
                        6*Power(s1,2)*(-1 + s2 - t1 + t2) + 
                        3*s2*(-4 + t1 + t2) + 
                        3*s1*
                        (2*Power(s2,2) + t1 + Power(t1,2) + t2 - 
                        Power(t2,2) + s2*(-2 - 3*t1 + t2))) + 
                       Power(s,2)*
                        (-30 + 30*s2 - 8*t1 - 18*s2*t1 + 
                        11*Power(t1,2) + s2*Power(t1,2) - 
                        Power(t1,3) + 52*t2 - 32*s2*t2 + 
                        6*Power(s2,2)*t2 + 8*t1*t2 + 2*s2*t1*t2 - 
                        Power(t1,2)*t2 - 25*Power(t2,2) + 
                        3*s2*Power(t2,2) - t1*Power(t2,2) + 
                        3*Power(t2,3) + 
                        6*Power(s1,4)*(-1 + s2 - t1 + t2) + 
                        6*Power(s1,3)*
                        (3*Power(s2,2) + 2*t1 + 2*Power(t1,2) + t2 - 
                        t1*t2 - Power(t2,2) + s2*(-3 - 5*t1 + 2*t2)) \
+ 2*s1*(Power(t1,3) + 3*Power(s2,2)*(-5 + t1 + t2) + 
                        Power(t1,2)*(-7 + 2*t2) - 
                        2*t2*(4 - 5*t2 + Power(t2,2)) - 
                        t1*(15 - 3*t2 + Power(t2,2)) + 
                        s2*(15 - 4*Power(t1,2) + t1*(23 - 6*t2) - 
                        6*t2 + 4*Power(t2,2))) + 
                        Power(s1,2)*
                        (30 + 6*Power(s2,3) - Power(t1,3) - 42*t2 + 
                        11*Power(t2,2) + Power(t2,3) - 
                        3*Power(t1,2)*(3 + t2) - 
                        6*Power(s2,2)*(1 + 2*t1 + t2) + 
                        t1*(26 - 8*t2 + 3*Power(t2,2)) + 
                        s2*(-30 + 7*Power(t1,2) + 20*t2 - 
                        11*Power(t2,2) + 2*t1*(7 + 5*t2)))) - 
                       2*s*(-7 + 6*s2 + 5*t1 - 7*s2*t1 + 
                        3*Power(t1,2) + s2*Power(t1,2) - Power(t1,3) + 
                        17*t2 - 18*s2*t2 + 5*Power(s2,2)*t2 - 
                        4*t1*t2 + 4*s2*t1*t2 - Power(s2,2)*t1*t2 - 
                        Power(t1,2)*t2 + s2*Power(t1,2)*t2 - 
                        13*Power(t2,2) + 12*s2*Power(t2,2) - 
                        3*Power(s2,2)*Power(t2,2) - t1*Power(t2,2) - 
                        s2*t1*Power(t2,2) + 3*Power(t2,3) + 
                        3*Power(s1,4)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                        Power(s1,3)*
                        (3*Power(s2,3) - Power(t1,3) - 
                        Power(t1,2)*(-1 + t2) + 3*Power(-1 + t2,2) - 
                        Power(s2,2)*(-1 + 7*t1 + t2) + 
                        s2*(-6 + 5*Power(t1,2) + 2*t1*(-1 + t2) + 
                        10*t2 - 4*Power(t2,2)) + 
                        2*t1*(2 - 3*t2 + Power(t2,2))) + 
                        s1*(Power(t1,3) - 3*Power(s2,3)*t2 + 
                        2*Power(t1,2)*(3 + t2) - 
                        Power(-1 + t2,2)*(5 + t2) - 
                        2*t1*(7 - 8*t2 + Power(t2,2)) + 
                        s2*(17 + t1*(2 - 18*t2) + 
                        Power(t1,2)*(-6 + t2) - 23*t2 + 
                        7*Power(t2,2) - Power(t2,3)) + 
                        Power(s2,2)*
                        (-9 + 18*t2 - Power(t2,2) + t1*(5 + 2*t2))) + 
                        Power(s1,2)*
                        (-13*Power(t1,2) + Power(t1,3) - 
                        3*Power(s2,3)*(-1 + t2) - 
                        Power(-1 + t2,2)*(-3 + 2*t2) + 
                        t1*(-4 + 3*t2 + Power(t2,2)) + 
                        Power(s2,2)*
                        (-13 + 2*t2 - 2*Power(t2,2) + t1*(-4 + 5*t2)) \
+ s2*(4 - 2*t2 - 2*Power(t1,2)*t2 - 3*Power(t2,2) + Power(t2,3) + 
                        t1*(25 + Power(t2,2))))))))/
                (Power(-1 + s - s*s1 + s1*(s2 - t1) + t2,2)*
                  (-1 + s2 - t1 + t2))))/((-1 + t1)*(-1 + t2)))/
        (s - s1 + t2))*B1(s1,1 - s2 + t1 - t2,1 - s + s1 - t2))/
   (128.*Power(Pi,2)) - (((-8*
          (-4 - 9*t2 + 4*s*t2 - 12*s2*t2 - 13*Power(t2,2) - 
            9*s*Power(t2,2) + 5*s2*Power(t2,2) + 
            t1*(4 + s*(-4 + t2) + 11*t2) + 
            Power(s1,3)*(s2*t2 + t1*(4 + t2)) + 
            s1*(-4 - 12*Power(t1,2) + s2*(4 - 4*t1*(-2 + t2) - 5*t2) + 
               6*t2 + 4*s*t2 + 10*Power(t2,2) + 
               t1*(14 - 4*s + 9*t2 + 8*s*t2)) - 
            Power(s1,2)*(-4*Power(t1,2) + t2*(5 + t2 - s*t2) + 
               t1*(6 + (9 + s)*t2) + s2*(-4 + Power(t2,2)))))/
        ((-1 + s1)*(s1*t1 - t2)*(-1 + t2)) - 
       (8*(-3 - 3*t2 + 5*s*t2 - 9*s2*t2 - (1 + s)*Power(t1,2)*t2 - 
            8*Power(t2,2) - 6*s*Power(t2,2) + 10*s2*Power(t2,2) + 
            4*Power(t2,3) + 2*Power(s1,3)*(s2 - t1*(1 + 2*t2)) + 
            t1*(3 + t2 + 2*s2*t2 + (1 - 3*s2)*Power(t2,2) + 
               3*s*(-1 + Power(t2,2))) + 
            Power(s1,2)*(s2*(5 - t1*(-4 + t2) - 4*t2) + 
               Power(t1,2)*(-4 + t2) + 
               t1*(7 + 2*s*(-1 + t2) + 4*t2 + 2*Power(t2,2)) + 
               2*(-1 + t2 + s*t2 + 2*Power(t2,2))) + 
            s1*(-5 + 2*Power(t1,3) + Power(t1,2)*(-4 + s*(2 - 3*t2)) - 
               3*t2 + 5*s*t2 - 4*Power(t2,2) - 4*s*Power(t2,2) - 
               2*Power(t2,3) + 
               t1*(10 + 9*t2 - 5*Power(t2,2) + 
                  s*(-7 + 6*t2 + Power(t2,2))) + 
               s2*(3 + 2*Power(t1,2)*(-1 + t2) - 11*t2 + 
                  4*Power(t2,2) - t1*(-6 + 7*t2 + Power(t2,2))))))/
        ((-1 + t1)*(s1*t1 - t2)*(-1 + t2)) + 
       (8*(-1 + s1 + t1 - t2)*
          (1 + (-1 + s)*t1 + Power(s1,2)*(-s2 + t1) - t2 + s*t2 + 
            s2*t2 - s1*(-1 + s2 + t2 - s2*t2 + s*(t1 + t2)))*
          (-((5 + Power(s1,4) - 2*s2 + 3*Power(s2,2) + 7*t1 - 
                 8*s2*t1 - Power(s2,2)*t1 - 2*Power(t1,2) + 
                 2*s2*Power(t1,2) - 2*t2 + 2*Power(s2,2)*t2 - 
                 3*t1*t2 - s2*t1*t2 + Power(t2,2) + 
                 Power(s1,3)*(-4 + 2*s2 + t1 + t2) + 
                 Power(s,2)*
                  (1 + 2*Power(s1,2) - 6*t1 + Power(t1,2) + 
                    s1*(-5 + 3*t1) + 2*t2) + 
                 s1*(5 + Power(s2,2)*t1 + 4*Power(t1,2) + 
                    s2*(t1*(-6 + t2) + 2*(-2 + t2)) - 8*t2 + 
                    Power(t2,2) - 4*t1*(3 + t2)) + 
                 Power(s1,2)*
                  (Power(s2,2) + t1*(-2 + t2) - 2*(5 + t2) + 
                    s2*(-5 + 2*t1 + t2)) - 
                 s*(1 + 3*Power(s1,3) - 7*t1 + 2*Power(t1,2) - t2 - 
                    2*t1*t2 + Power(s1,2)*(-10 + 3*s2 + 4*t1 + t2) + 
                    s1*(-7 + Power(t1,2) + s2*(-5 + 4*t1) + 
                       t1*(-11 + t2) + t2) + 
                    s2*(4 - 7*t1 + Power(t1,2) + 4*t2)))/
               ((-1 + s2)*(-s + s2 - t1))) + 
            (-11 + 3*Power(s1,3) + 6*s2 + Power(s2,2) + 11*t1 - 
               s2*t1 - 4*Power(s,2)*(-2 + t2) + 8*t2 + s2*t2 - 
               2*Power(s2,2)*t2 + 3*t1*t2 + Power(t2,2) - 
               Power(s1,2)*(7 - 3*t1 + s2*(-3 + t2) + 3*t2) + 
               s1*(8 - 14*t1 + s2*(7 + t1) - Power(s2,2)*(-2 + t2) + 
                  14*t2 - 3*t1*t2 + 3*Power(t2,2)) - 
               s*(6 + Power(s1,2) + 9*s2 - t1 + t2 - 6*s2*t2 + 
                  s1*(14 + 2*s2 + t1 - 2*t2 - s2*t2)))/
             ((-1 + s1)*(-s + s1 - t2)) + 
            (-2 + Power(s1,3) + 11*s2 + 2*t1 - 2*s2*t1 + 3*t2 - 
               4*Power(s,2)*t2 + 4*s2*t2 + Power(s2,2)*t2 - 3*t1*t2 - 
               Power(t2,2) + 
               Power(s1,2)*(-3 + t1 - s2*(-2 + t2) + 4*t2) + 
               s1*(12 + 3*t1*(-1 + t2) + 5*t2 - Power(s2,2)*t2 - 
                  3*Power(t2,2) + s2*(-1 + 2*t1 + t2)) + 
               s*(-11 - 2*Power(s1,2) + 2*t1 - 4*t2 + 3*s2*t2 + 
                  s1*(1 - 2*t1 + 4*t2 + s2*t2)))/((-1 + s1)*(-1 + t2)) \
- (1 - 8*s2 + Power(s2,2) - 4*t1 + s2*t1 + Power(t1,2) + 
               2*Power(s1,3)*(-1 + t2) - 4*t2 + s2*t2 + 2*t1*t2 - 
               s2*t1*t2 + 3*Power(t2,2) + s2*Power(t2,2) + 
               Power(s,2)*(-1 + (3 + s1 - t1)*t2) - 
               Power(s1,2)*
                (6 + 3*t1 - 3*(-1 + s2)*t2 + Power(t2,2)) - 
               s1*(6 + Power(t1,2) + t1*(-6 + t2) + 4*t2 - 
                  Power(s2,2)*t2 - 2*Power(t2,2) + 
                  s2*(8 + 4*t2 + Power(t2,2))) + 
               s*(8 - 3*Power(s1,2)*t2 - 3*s2*t2 - 2*Power(t2,2) + 
                  t1*(-1 + t2 + s2*t2) + 
                  s1*(10 + (1 - 2*s2 + t1)*t2 + Power(t2,2))))/
             ((-1 + t1)*(-1 + t2)) + 
            (-4 + 2*s2 - 3*Power(s2,2) - 2*t1 + 7*s2*t1 + 
               2*Power(t1,2) - Power(s1,3)*(-5 + t2) + 9*t2 + 
               4*s2*t2 - 4*Power(s2,2)*t2 + 3*t1*t2 - 2*s2*t1*t2 + 
               3*Power(t2,2) + s2*Power(t2,2) - 
               Power(s,2)*(3 + t1*(-6 + t2) + 2*s1*(-3 + t2) + 
                  2*t2) - Power(s1,2)*
                (-5 - 2*t1 + 2*s2*(-4 + t2) + 4*t2 + Power(t2,2)) - 
               s1*(11 + 2*Power(t1,2) + t1*(-10 + t2) - 9*t2 + 
                  Power(s2,2)*t2 - 2*Power(t2,2) + 
                  s2*(-2 - 7*t1 + 10*t2 + Power(t2,2))) + 
               s*(1 - 7*t1 - 8*t2 + 2*t1*t2 - 2*Power(t2,2) + 
                  Power(s1,2)*(-14 + 3*t2) + 
                  s2*(t1*(-6 + t2) + 6*(1 + t2)) + 
                  s1*(t1*(-13 + t2) + 3*s2*(-2 + t2) + t2*(8 + t2))))/
             ((s - s2 + t1)*(s - s1 + t2)) + 
            (-11 - 2*Power(s1,4) + 3*s2 + 4*Power(s2,2) + 3*t1 - 
               4*s2*t1 - Power(s2,2)*t1 - 2*Power(t1,2) + 
               2*s2*Power(t1,2) + 8*t2 + 6*s2*t2 - 3*Power(s2,2)*t2 - 
               3*t1*t2 - s2*t1*t2 + Power(t2,2) + 
               Power(s1,3)*(4 - 3*s2 - 2*t1 + t2) + 
               Power(s,2)*(9 - s1 - Power(s1,2) - 5*t1 + 
                  Power(t1,2) + t2) + 
               Power(s1,2)*(11 - Power(s2,2) - t2 + t1*(4 + t2) + 
                  s2*(4 - 3*t1 + t2)) + 
               s1*(1 - Power(s2,2)*(-2 + t1) + 4*Power(t1,2) + 
                  s2*(12 + t1*(-1 + t2) - t2) + 7*t2 + Power(t2,2) - 
                  t1*(5 + 4*t2)) + 
               s*(-3 + 3*Power(s1,3) + 3*t1 - 2*Power(t1,2) - 
                  s2*(13 - 6*t1 + Power(t1,2) - 2*t2) + 
                  Power(s1,2)*(-1 + 2*s2 + 2*t1 - t2) - 6*t2 + 
                  2*t1*t2 + s1*
                   (s2*(-1 + t1) - Power(t1,2) - t1*(-5 + t2) - 
                     2*(9 + t2))))/((-1 + s2)*(-1 + t1))))/
        (Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
          Power(s,2)*(-1 + s1)*(t1*(-1 + t2) + (s1 - t2)*t2) - 
          s*(-1 + s1)*(Power(s1,2)*(-1 + t2) + 2*(1 + t1 - t2)*t2 + 
             s2*(1 + t1*(-1 + t2) - 2*t2 - Power(t2,2)) + 
             s1*(t1*(-3 + t2) - Power(-1 + t2,2) + s2*(1 + t2)))) + 
       (8*(-1 + Power(s1,6)*s2*(-s2 + t1) + 12*t2 + 3*s*t2 - 4*s2*t2 + 
            2*s*s2*t2 + 15*Power(t2,2) - 5*s*Power(t2,2) - 
            4*Power(s,2)*Power(t2,2) + 25*s2*Power(t2,2) + 
            3*s*s2*Power(t2,2) - Power(s2,2)*Power(t2,2) + 
            26*Power(t2,3) - 16*s*Power(t2,3) - 
            2*Power(s,2)*Power(t2,3) + 25*s2*Power(t2,3) + 
            5*s*s2*Power(t2,3) + Power(s2,2)*Power(t2,3) + 
            16*Power(t2,4) + 4*s2*Power(t2,4) + 
            (-1 + s)*Power(t1,3)*(-1 + s + 3*t2) - 
            Power(t1,2)*(3 + 2*(-9 + 2*s2)*t2 - 
               (4 + 9*s2)*Power(t2,2) + Power(s,2)*(1 + 4*t2) + 
               s*(-4 + (9 - 4*s2)*t2 + 14*Power(t2,2))) + 
            t1*(3 + (-27 + 8*s2)*t2 + 
               (-19 - 34*s2 + Power(s2,2))*Power(t2,2) - 
               (17 + 13*s2)*Power(t2,3) + 5*Power(s,2)*t2*(1 + t2) + 
               s*(-2 + (3 - 6*s2)*t2 + (35 - 9*s2)*Power(t2,2) + 
                  11*Power(t2,3))) - 
            Power(s1,5)*(t1*(1 + t1*(-7 + t2) + s*(-2 + t2)) + 
               Power(s2,2)*(t1 - 3*t2) + 
               s2*(-2 - 3*Power(t1,2) + t2 + 2*s*t2 + 
                  t1*(-2 + s + 3*t2))) + 
            Power(s1,4)*(-1 - Power(t1,3)*(-18 + t2) + t2 + 
               s*Power(t2,2) - Power(s,2)*Power(t2,2) + 
               Power(t1,2)*(-12 - 6*s - 18*t2 + Power(t2,2)) + 
               t1*(-2 + s*(-1 + t2) - 11*t2 + Power(s,2)*t2 + 
                  2*Power(t2,2)) + 
               Power(s2,2)*(2 + 2*Power(t1,2) + t2 - 3*Power(t2,2) + 
                  t1*(3 + 2*t2)) + 
               s2*(2*Power(t1,3) - Power(t1,2)*(4 + 3*s + 4*t2) + 
                  t2*(-8 + 3*s + 3*t2 + 4*s*t2) - 
                  t1*(10 + 3*s + 12*t2 + 2*s*t2 - 3*Power(t2,2)))) + 
            Power(s1,3)*(9*Power(t1,4) + 
               Power(t1,3)*(-17 + s*(-21 + t2) - 21*t2) + 
               t2*(5 + 4*t2 - Power(t2,2) + Power(s,2)*t2*(1 + t2) - 
                  s*(1 + 5*t2)) + 
               Power(t1,2)*(-12 - 39*t2 + 14*Power(t2,2) + 
                  Power(s,2)*(3 + t2) + s*(7 + 17*t2 - 2*Power(t2,2))) \
+ Power(s2,2)*(2*Power(t1,3) - 2*Power(t1,2)*(1 + t2) + 
                  t2*(-9 - 2*t2 + Power(t2,2)) - 
                  t1*(1 + 10*t2 + Power(t2,2))) + 
               s2*(-4 + (5 - 2*s)*Power(t1,3) + (10 + 7*s)*t2 + 
                  (15 - 2*s)*Power(t2,2) - (3 + 2*s)*Power(t2,3) + 
                  Power(t1,2)*
                   (-9 + s*(-2 + t2) - 10*t2 + Power(t2,2)) + 
                  t1*(-4 + s + 32*t2 + 18*s*t2 + 14*Power(t2,2) + 
                     3*s*Power(t2,2) - Power(t2,3))) + 
               t1*(11 + 30*t2 + 33*Power(t2,2) - 2*Power(t2,3) - 
                  2*Power(s,2)*t2*(2 + t2) + 
                  s*(-1 + 6*t2 - 4*Power(t2,2) + Power(t2,3)))) + 
            s1*(-3*(-1 + s)*Power(t1,4) + 
               Power(s2,2)*t2*
                (4 + 2*Power(t1,2) + t1*(-6 + t2) - 4*t2 - 
                  7*Power(t2,2)) + 
               Power(t1,3)*(2*Power(s,2) - 6*(3 + t2) + 
                  s*(11 + 23*t2)) - 
               Power(t1,2)*(-27 - 27*t2 - 47*Power(t2,2) + 
                  Power(s,2)*(3 + 11*t2) + 
                  s*(7 + 43*t2 + 24*Power(t2,2))) + 
               t2*(-7 + 14*t2 + 4*Power(t2,2) + 11*Power(t2,3) - 
                  Power(s,2)*t2*(1 + 7*t2) + 
                  s*(1 + 13*t2 + 2*Power(t2,2) - 3*Power(t2,3))) + 
               t1*(-12 - 14*t2 - 88*Power(t2,2) - 55*Power(t2,3) + 
                  4*Power(s,2)*t2*(1 + 4*t2) + 
                  s*(-1 - 13*t2 + 6*Power(t2,2) + 7*Power(t2,3))) + 
               s2*(2 + Power(t1,3)*(2 - 2*s - 15*t2) - (13 + 5*s)*t2 + 
                  (-31 + 10*s)*Power(t2,2) + 2*(6 + 7*s)*Power(t2,3) + 
                  3*Power(t2,4) - 
                  t1*(2 + (21 - 6*s)*t2 + (29 + 23*s)*Power(t2,2) + 
                     12*Power(t2,3)) + 
                  Power(t1,2)*
                   (-2 + 49*t2 + 24*Power(t2,2) + s*(2 + 11*t2)))) + 
            Power(s1,2)*(2 + (2 - 9*s)*Power(t1,4) - 11*t2 - 3*s*t2 - 
               21*Power(t2,2) + 4*s*Power(t2,2) + 
               5*Power(s,2)*Power(t2,2) - 15*Power(t2,3) + 
               6*s*Power(t2,3) + Power(t2,4) - s*Power(t2,4) + 
               Power(t1,3)*(-10 + 5*Power(s,2) - 39*t2 + 
                  s*(12 + 13*t2)) + 
               Power(t1,2)*(5 + Power(s,2)*(1 - 10*t2) + 79*t2 + 
                  60*Power(t2,2) + s*(10 + 27*t2 - 8*Power(t2,2))) + 
               Power(s2,2)*(-1 - 2*Power(t1,3) - 
                  4*Power(t1,2)*(-1 + t2) + t2 + 14*Power(t2,2) + 
                  Power(t2,3) + t1*(-1 + 6*t2 + 7*Power(t2,2))) + 
               t1*(1 - 2*t2 + 18*Power(t2,2) - 24*Power(t2,3) + 
                  Power(s,2)*t2*(-6 + 5*t2) + 
                  s*(3 - 12*t2 - 21*Power(t2,2) + 5*Power(t2,3))) + 
               s2*(6*Power(t1,4) - Power(t1,3)*(15 + 4*s + 11*t2) + 
                  Power(t1,2)*
                   (-4 + 3*t2 + 8*Power(t2,2) + s*(-5 + 16*t2)) + 
                  t1*(13 + 28*t2 - 4*Power(t2,2) - 4*Power(t2,3) + 
                     s*(3 - 11*Power(t2,2))) + 
                  t2*(16 - 28*t2 - 12*Power(t2,2) + Power(t2,3) - 
                     s*(5 + 23*t2 + Power(t2,2)))))))/
        ((-1 + s1)*(-s + s1 - t2)*(-1 + s1 + t1 - t2)*
          Power(-(s1*t1) + t2,2)) + 
       (8*(-1 + 2*Power(s1,6)*t1*(s2 + 3*t1) + 13*t2 + s*t2 - 6*s2*t2 + 
            2*s*s2*t2 - 14*s*Power(t2,2) - 2*Power(s,2)*Power(t2,2) + 
            36*s2*Power(t2,2) - s*s2*Power(t2,2) - 
            7*Power(s2,2)*Power(t2,2) + Power(t2,3) - 
            28*s*Power(t2,3) + 37*s2*Power(t2,3) + s*s2*Power(t2,3) - 
            5*Power(s2,2)*Power(t2,3) + 7*Power(t2,4) - 
            5*s*Power(t2,4) + 3*s2*Power(t2,4) - 
            (-1 + s)*Power(t1,4)*(-1 + s + t2) + 
            Power(t1,3)*(4 + (-7 + 3*s2)*t2 - (4 + 3*s2)*Power(t2,2) + 
               Power(s,2)*(2 + t2) + 
               s*(-6 - 3*(-1 + s2)*t2 + 4*Power(t2,2))) + 
            Power(t1,2)*(-6 - 12*(-2 + s2)*t2 + 
               (-2 + 25*s2 - 2*Power(s2,2))*Power(t2,2) + 
               (7 + 3*s2)*Power(t2,3) + 
               Power(s,2)*(-1 - 4*t2 + Power(t2,2)) + 
               s*(6 + 11*(-1 + s2)*t2 + 2*(-13 + s2)*Power(t2,2) - 
                  3*Power(t2,3))) + 
            t1*(4 + (-31 + 15*s2)*t2 + 
               (6 - 58*s2 + 9*Power(s2,2))*Power(t2,2) - 
               (7 + 20*s2)*Power(t2,3) - 4*Power(t2,4) + 
               Power(s,2)*t2*(3 + 2*t2 - Power(t2,2)) + 
               s*(-2 + (8 - 10*s2)*t2 + (45 - 4*s2)*Power(t2,2) + 
                  (24 + s2)*Power(t2,3))) + 
            Power(s1,5)*(Power(s2,2) + 
               s2*(9*Power(t1,2) - 2*t2 - t1*(1 + s + 5*t2)) + 
               t1*(9*Power(t1,2) - t1*(10 + 5*s + 9*t2) + 
                  2*(-1 + (-6 + s)*t2))) + 
            Power(s1,4)*(2*Power(t1,4) - 
               Power(t1,3)*(3 + 3*s + 8*t2) + 2*t2*(1 + 3*t2 - s*t2) + 
               Power(t1,2)*(-28 - s + Power(s,2) - 14*t2 + 6*s*t2 + 
                  3*Power(t2,2)) - 
               t1*(-1 + s - 25*t2 - 9*s*t2 + Power(s,2)*t2 - 
                  18*Power(t2,2) + 3*s*Power(t2,2)) + 
               Power(s2,2)*(1 + Power(t1,2) - 4*t2 + t1*(4 + t2)) + 
               s2*(-2 + 7*Power(t1,3) + t2 + 3*s*t2 + 5*Power(t2,2) - 
                  2*Power(t1,2)*(1 + 2*s + 6*t2) + 
                  t1*(-11 + s*(-1 + t2) - 15*t2 + 4*Power(t2,2)))) + 
            Power(s1,3)*(1 - Power(t1,5) + 
               Power(t1,4)*(3 + 5*s - t2) - t2 - s*t2 - 
               15*Power(t2,2) - 4*s*Power(t2,2) + 
               2*Power(s,2)*Power(t2,2) - 9*Power(t2,3) + 
               3*s*Power(t2,3) - 
               Power(t1,3)*(37 + 5*t2 - 2*Power(t2,2) + 
                  4*s*(5 + t2)) + 
               Power(s2,2)*(-1 + Power(t1,3) - 8*t2 + 5*Power(t2,2) + 
                  Power(t1,2)*(5 + t2) - t1*t2*(11 + 2*t2)) + 
               Power(t1,2)*(22 + 2*Power(s,2) + 35*t2 + 
                  19*Power(t2,2) + s*(21 + 20*t2)) - 
               s2*(2 + 2*Power(t1,4) - (19 + s)*t2 + 
                  (-6 + 7*s)*Power(t2,2) + 4*Power(t2,3) + 
                  Power(t1,3)*(-15 + 2*s + t2) + 
                  Power(t1,2)*(30 + 2*s + 24*t2 - 3*Power(t2,2)) - 
                  t1*(-1 + s + 27*t2 + 17*s*t2 + 21*Power(t2,2) + 
                     2*s*Power(t2,2) - Power(t2,3))) + 
               t1*(11 + 53*t2 - 4*Power(s,2)*t2 - 3*Power(t2,2) - 
                  6*Power(t2,3) + 
                  s*(1 - 5*t2 - 11*Power(t2,2) + Power(t2,3)))) + 
            s1*(-1 + (-1 + s)*Power(t1,5) - t2 + s*t2 + 
               52*Power(t2,2) + 16*s*Power(t2,2) - 
               2*Power(s,2)*Power(t2,2) + 26*Power(t2,3) - 
               4*Power(s,2)*Power(t2,3) + 2*Power(t2,4) + 
               Power(t1,4)*(6 + 8*t2 - s*(3 + 7*t2)) + 
               Power(s2,2)*t2*
                (6 + Power(t1,3) - 9*Power(t2,2) - 
                  2*Power(t1,2)*(3 + t2) + 
                  t1*(-1 + 17*t2 + 3*Power(t2,2))) + 
               Power(t1,3)*(-20 + 5*t2 - 16*Power(t2,2) + 
                  2*Power(s,2)*(1 + t2) + 
                  s*(11 + 40*t2 + 7*Power(t2,2))) - 
               Power(t1,2)*(-25 + 15*t2 - 24*Power(t2,2) - 
                  7*Power(t2,3) + 
                  Power(s,2)*(2 + 8*t2 + 4*Power(t2,2)) + 
                  s*(8 + 59*t2 + 42*Power(t2,2) + Power(t2,3))) + 
               t1*(-9 + 3*t2 - 63*Power(t2,2) - 14*Power(t2,3) + 
                  2*Power(t2,4) + 
                  2*Power(s,2)*t2*(2 + 5*t2 + Power(t2,2)) + 
                  s*(-1 + 7*t2 + 40*Power(t2,2) + 11*Power(t2,3))) + 
               s2*(2 - (17 + s)*t2 + (-42 + 13*s)*Power(t2,2) + 
                  (12 + 13*s)*Power(t2,3) + 2*Power(t2,4) + 
                  Power(t1,4)*(-1 + s + 5*t2) - 
                  t1*(2 + (30 - 7*s)*t2 + (46 + 27*s)*Power(t2,2) + 
                     (14 + 5*s)*Power(t2,3)) - 
                  Power(t1,3)*
                   (-4 + 40*t2 + 6*Power(t2,2) + s*(5 + 3*t2)) + 
                  Power(t1,2)*
                   (-3 + 82*t2 + 34*Power(t2,2) + Power(t2,3) + 
                     s*(4 + 3*t2 + 7*Power(t2,2))))) + 
            Power(s1,2)*(1 + (-4 + 3*s)*Power(t1,5) - 15*t2 - s*t2 - 
               25*Power(t2,2) + 12*s*Power(t2,2) + 
               2*Power(s,2)*Power(t2,2) + 8*Power(t2,3) + 
               5*s*Power(t2,3) - 2*Power(s,2)*Power(t2,3) + 
               3*Power(t2,4) - s*Power(t2,4) - 
               Power(t1,4)*(2 + Power(s,2) - 10*t2 + 2*s*(9 + 2*t2)) + 
               Power(t1,2)*(3 + 95*t2 + 10*Power(t2,2) - 
                  4*Power(t2,3) - 3*Power(s,2)*t2*(2 + t2) + 
                  s*(-5 + 8*t2)) + 
               t1*(-3 - 68*t2 - 58*Power(t2,2) - 13*Power(t2,3) + 
                  s*(3 - 37*t2 - 23*Power(t2,2)) + 
                  Power(s,2)*t2*(-2 + 6*t2 + Power(t2,2))) + 
               Power(t1,3)*(5 - 19*t2 - 2*Power(t2,2) + 
                  Power(s,2)*(2 + 3*t2) + s*(26 + 13*t2 + Power(t2,2))) \
+ Power(s2,2)*(-1 + 4*t2 + 16*Power(t2,2) - 2*Power(t2,3) + 
                  Power(t1,3)*(1 + t2) + 
                  Power(t1,2)*(2 - 9*t2 - 2*Power(t2,2)) + 
                  t1*(-2 - 11*t2 + 4*Power(t2,2) + Power(t2,3))) + 
               s2*(2 - 2*Power(t1,5) + (9 - 5*s)*t2 - 
                  (37 + 13*s)*Power(t2,2) + (-9 + 4*s)*Power(t2,3) + 
                  Power(t2,4) + Power(t1,4)*(17 + s + 3*t2) - 
                  Power(t1,3)*
                   (32 + 14*t2 + Power(t2,2) + s*(-1 + 4*t2)) + 
                  Power(t1,2)*
                   (2*(3 + t2 + 5*Power(t2,2)) + 
                     s*(-6 + 16*t2 + 5*Power(t2,2))) + 
                  t1*(9 + 60*t2 + 13*Power(t2,2) - 5*Power(t2,3) + 
                     s*(1 + t2 - 13*Power(t2,2) - 2*Power(t2,3)))))))/
        ((-1 + s2)*(-1 + t1)*(-1 + s1 + t1 - t2)*Power(-(s1*t1) + t2,2)) \
+ (8*(5 + Power(s1,6)*(2*Power(s2,2) + s2*t1 + 5*Power(t1,2)) - 5*t2 + 
            2*s*t2 + 12*s2*t2 - 10*s*s2*t2 - Power(t2,2) - 
            22*s*Power(t2,2) + 3*Power(s,2)*Power(t2,2) + 
            4*s2*Power(t2,2) + 24*s*s2*Power(t2,2) + 
            Power(s2,2)*Power(t2,2) - 21*Power(t2,3) + 
            18*s*Power(t2,3) - 7*Power(s,2)*Power(t2,3) + 
            32*s2*Power(t2,3) + 14*s*s2*Power(t2,3) - 
            9*Power(s2,2)*Power(t2,3) + 10*Power(t2,4) + 
            2*s*Power(t2,4) + (-1 + s)*Power(t1,4)*(-1 + s + t2) - 
            Power(t1,3)*(8 + 3*(-3 + s2)*t2 - (4 + 3*s2)*Power(t2,2) + 
               Power(s,2)*(6 + t2) + 
               s*(-14 + (5 - 3*s2)*t2 + 4*Power(t2,2))) - 
            Power(t1,2)*(-18 + (20 - 18*s2)*t2 + 
               (17 + 21*s2 - 2*Power(s2,2))*Power(t2,2) + 
               (7 + 3*s2)*Power(t2,3) + 
               Power(s,2)*(-5 - 5*t2 + Power(t2,2)) + 
               s*(22 + (4 + 17*s2)*t2 + (-27 + 2*s2)*Power(t2,2) - 
                  3*Power(t2,3))) + 
            t1*(-16 + (17 - 27*s2)*t2 + 
               (14 + 14*s2 - 3*Power(s2,2))*Power(t2,2) + 
               (7 + 13*s2)*Power(t2,3) + 4*Power(t2,4) + 
               Power(s,2)*t2*(-8 + 8*t2 + Power(t2,2)) - 
               s*(-10 - 6*(1 + 4*s2)*t2 + (28 + 5*s2)*Power(t2,2) + 
                  (20 + s2)*Power(t2,3))) + 
            Power(s1,5)*(Power(s2,2)*(-3 + 6*t1 - 6*t2) + 
               s2*(-4 + 10*Power(t1,2) - t2 + 4*s*t2 - t1*(3 + t2)) + 
               t1*(-1 - 23*t1 + 10*Power(t1,2) - 10*t2 - 2*t1*t2 + 
                  s*(-4 - 9*t1 + 5*t2))) + 
            Power(s1,4)*(2 - (39 + 19*s)*Power(t1,3) + 6*Power(t1,4) + 
               t2 + 5*Power(t2,2) - 5*s*Power(t2,2) + 
               2*Power(s,2)*Power(t2,2) + 
               Power(t1,2)*(-16 + 2*Power(s,2) - 7*t2 - 
                  3*Power(t2,2) + 17*s*(1 + t2)) + 
               Power(s2,2)*(-7 + 7*Power(t1,2) + 4*t2 + 
                  6*Power(t2,2) - t1*(6 + 13*t2)) + 
               t1*(3 + 47*t2 - 4*Power(s,2)*t2 + 4*Power(t2,2) + 
                  s*(10 + 8*t2 - 4*Power(t2,2))) + 
               s2*(6 + 17*Power(t1,3) + 15*t2 - 10*s*t2 + 
                  Power(t2,2) - 8*s*Power(t2,2) - 
                  3*Power(t1,2)*(5 + s + 3*t2) - 
                  t1*(16 + s*(2 - 13*t2) + 14*t2 + Power(t2,2)))) + 
            Power(s1,3)*(-3 + Power(t1,5) - 9*t2 - 24*Power(t2,2) + 
               9*s*Power(t2,2) - 3*Power(s,2)*Power(t2,2) - 
               2*Power(t2,3) + 4*s*Power(t2,3) - 
               2*Power(s,2)*Power(t2,3) + 
               Power(t1,4)*(-11 - 13*s + t2) + 
               Power(t1,3)*(3*Power(s,2) + 2*s*(20 + 7*t2) - 
                  2*(10 + 4*t2 + Power(t2,2))) + 
               Power(t1,2)*(-(Power(s,2)*(1 + 8*t2)) + 
                  s*(25 + t2 - 2*Power(t2,2)) + 
                  2*(22 + 54*t2 + Power(t2,2))) + 
               Power(s2,2)*(3 + 3*Power(t1,3) + 
                  Power(t1,2)*(5 - 9*t2) + 28*t2 + Power(t2,2) - 
                  2*Power(t2,3) + t1*(-2 - t2 + 8*Power(t2,2))) + 
               t1*(10 + 26*t2 - 15*Power(t2,2) + 6*Power(t2,3) + 
                  Power(s,2)*t2*(4 + 7*t2) - 
                  s*(-10 + 68*t2 + 19*Power(t2,2) + Power(t2,3))) + 
               s2*(14 + 10*Power(t1,4) - 2*(2 + 3*s)*t2 + 
                  (-8 + 9*s)*Power(t2,2) + (1 + 4*s)*Power(t2,3) - 
                  Power(t1,3)*(25 + 4*s + 9*t2) + 
                  t1*(3 + (59 - 15*s)*t2 - 15*(-1 + s)*Power(t2,2) + 
                     Power(t2,3)) - 
                  Power(t1,2)*
                   (44 + 22*t2 + Power(t2,2) - 3*s*(-4 + 5*t2)))) + 
            s1*(3 - (-1 + s)*Power(t1,5) + 19*t2 + 48*Power(t2,2) - 
               19*s*Power(t2,2) + 3*Power(s,2)*Power(t2,2) + 
               25*Power(t2,3) - 44*s*Power(t2,3) + 
               10*Power(s,2)*Power(t2,3) + 3*Power(t2,4) - 
               s*Power(t2,4) + 
               Power(t1,4)*(-8*(1 + t2) + s*(5 + 7*t2)) + 
               Power(s2,2)*t2*
                (-Power(t1,3) + 2*Power(t1,2)*(-6 + t2) + 
                  t1*(25 + 13*t2 - 3*Power(t2,2)) + 
                  3*(-4 + 4*t2 + 3*Power(t2,2))) + 
               Power(t1,3)*(Power(s,2)*(3 - 2*t2) - 
                  s*t2*(39 + 7*t2) + 2*(7 + 16*t2 + 8*Power(t2,2))) + 
               Power(t1,2)*(-5 - 5*t2 - 34*Power(t2,2) - 
                  7*Power(t2,3) + Power(s + 2*s*t2,2) + 
                  s*(2 - t2 + 27*Power(t2,2) + Power(t2,3))) - 
               t1*(5 + 38*t2 - 33*Power(t2,2) + 16*Power(t2,3) + 
                  2*Power(t2,4) + 
                  Power(s,2)*t2*(4 + 17*t2 + 2*Power(t2,2)) - 
                  s*(-6 + 87*t2 + 59*Power(t2,2) + 2*Power(t2,3))) - 
               s2*(10 - (13 + 2*s)*t2 + (-9 + 51*s)*Power(t2,2) + 
                  3*(-3 + 7*s)*Power(t2,3) + Power(t2,4) + 
                  Power(t1,4)*(-1 + s + 5*t2) + 
                  t1*(-20 + (44 + 57*s)*t2 - 
                     2*(-65 + 8*s)*Power(t2,2) - 5*s*Power(t2,3)) - 
                  Power(t1,3)*
                   (-2 + 30*t2 + 6*Power(t2,2) + 3*s*(1 + t2)) + 
                  Power(t1,2)*
                   (9 - 6*t2 + 11*Power(t2,2) + Power(t2,3) + 
                     s*(2 - 18*t2 + 7*Power(t2,2))))) + 
            Power(s1,2)*(-7 + (4 - 3*s)*Power(t1,5) - 2*s*t2 - 
               4*Power(t2,2) + 29*s*Power(t2,2) - 
               5*Power(s,2)*Power(t2,2) + 12*Power(t2,3) - 
               2*s*Power(t2,3) + Power(s,2)*Power(t2,3) - 
               3*Power(t2,4) + s*Power(t2,4) + 
               Power(t1,4)*(Power(s,2) + 4*s*(4 + t2) - 2*(8 + 5*t2)) + 
               Power(t1,3)*(1 + 37*t2 + 2*Power(t2,2) - 
                  Power(s,2)*(2 + 3*t2) + s*(1 + 6*t2 - Power(t2,2))) + 
               Power(t1,2)*(15 + 14*t2 + 8*Power(t2,2) + 
                  4*Power(t2,3) + 
                  Power(s,2)*(-7 + 5*t2 + 3*Power(t2,2)) - 
                  s*(21 + 109*t2 + 20*Power(t2,2))) - 
               t1*(-3 + 104*t2 + 96*Power(t2,2) + 5*Power(t2,3) + 
                  Power(s,2)*t2*(-12 + 4*t2 + Power(t2,2)) + 
                  s*(20 + 22*t2 - 76*Power(t2,2) - 3*Power(t2,3))) - 
               Power(s2,2)*(-5 + Power(t1,3)*(-7 + t2) + 8*t2 + 
                  30*Power(t2,2) + 2*Power(t2,3) + 
                  Power(t1,2)*(8 + 13*t2 - 2*Power(t2,2)) + 
                  t1*(4 + 5*t2 - 10*Power(t2,2) + Power(t2,3))) + 
               s2*(-6 + 2*Power(t1,5) + (-47 + 20*s)*t2 + 
                  (-20 + 34*s)*Power(t2,2) + (-2 + s)*Power(t2,3) - 
                  Power(t2,4) - Power(t1,4)*(11 + s + 3*t2) + 
                  Power(t1,3)*
                   (-8 - 10*t2 + Power(t2,2) + s*(-9 + 4*t2)) + 
                  Power(t1,2)*
                   (16 + 103*t2 + 11*Power(t2,2) + 
                     s*(25 + 2*t2 - 5*Power(t2,2))) + 
                  t1*(7 + 83*t2 - 16*Power(t2,2) + 2*Power(t2,3) + 
                     s*(2 + 19*t2 - 2*Power(t2,2) + 2*Power(t2,3)))))))/
        ((-1 + s2)*(-s + s2 - t1)*(-1 + s1 + t1 - t2)*
          Power(-(s1*t1) + t2,2)) + 
       (8*(5 + Power(-1 + s,2)*Power(t1,4) + 
            Power(s1,6)*(2*Power(s2,2) - s2*t1 + Power(t1,2)) - 3*t2 + 
            2*s*t2 + 12*s2*t2 - 10*s*s2*t2 + 11*Power(t2,2) - 
            22*s*Power(t2,2) + 3*Power(s,2)*Power(t2,2) + 
            10*s2*Power(t2,2) + 24*s*s2*Power(t2,2) + 
            Power(s2,2)*Power(t2,2) + 9*Power(t2,3) + 
            14*s*Power(t2,3) - 7*Power(s,2)*Power(t2,3) + 
            40*s2*Power(t2,3) + 14*s*s2*Power(t2,3) - 
            9*Power(s2,2)*Power(t2,3) + 34*Power(t2,4) - 
            2*s*Power(t2,4) + 2*s2*Power(t2,4) + 4*Power(t2,5) - 
            Power(t1,3)*(8 + (-4 + 3*s2)*t2 + 5*Power(t2,2) + 
               Power(s,2)*(6 + t2) - 
               s*(14 + (-1 + 3*s2)*t2 + Power(t2,2))) - 
            Power(t1,2)*(-18 + (11 - 18*s2)*t2 + 
               (1 + 7*s2 - 2*Power(s2,2))*Power(t2,2) - 
               3*(4 + s2)*Power(t2,3) + 
               Power(s,2)*(-5 - 5*t2 + Power(t2,2)) + 
               s*(22 + (9 + 17*s2)*t2 + (-13 + 2*s2)*Power(t2,2) + 
                  4*Power(t2,3))) + 
            t1*(-16 + (10 - 27*s2)*t2 - 
               (5 + 3*s2 + 3*Power(s2,2))*Power(t2,2) - 
               (28 + 3*s2)*Power(t2,3) - (11 + 3*s2)*Power(t2,4) + 
               Power(s,2)*t2*(-8 + 8*t2 + Power(t2,2)) + 
               s*(10 + 8*(1 + 3*s2)*t2 - (19 + 5*s2)*Power(t2,2) - 
                  (6 + s2)*Power(t2,3) + 3*Power(t2,4))) + 
            Power(s1,5)*(Power(s2,2)*(-3 + 6*t1 - 6*t2) + 
               s2*(-4 + t1 - Power(t1,2) + t2 + 4*s*t2 + 5*t1*t2) + 
               t1*(1 - 15*t1 + Power(t1,2) - 2*t2 + 4*t1*t2 + 
                  s*(-4 - 3*t1 + 3*t2))) + 
            Power(s1,4)*(2 - t2 + Power(t2,2) - 3*s*Power(t2,2) + 
               2*Power(s,2)*Power(t2,2) + 
               Power(t1,3)*(-29 - 4*s + 6*t2) + 
               Power(t1,2)*(3 + 9*s + 2*Power(s,2) + 11*t2 - 
                  3*Power(t2,2)) + 
               t1*(-1 + 25*t2 - 4*Power(s,2)*t2 - 8*Power(t2,2) + 
                  2*s*(5 + t2)) + 
               Power(s2,2)*(-7 + 7*Power(t1,2) + 4*t2 + 
                  6*Power(t2,2) - t1*(6 + 13*t2)) - 
               s2*(-6 + Power(t1,3) + Power(t1,2)*(-4 + 3*s - 14*t2) + 
                  (-11 + 10*s)*t2 + (5 + 8*s)*Power(t2,2) + 
                  t1*(16 + s*(2 - 13*t2) + 5*t2 + 7*Power(t2,2)))) + 
            Power(s1,3)*(-3 - 5*t2 - 10*Power(t2,2) + 
               9*s*Power(t2,2) - 3*Power(s,2)*Power(t2,2) + 
               4*Power(t2,3) - 2*Power(s,2)*Power(t2,3) + 
               Power(t1,4)*(-17 - s + t2) + 
               Power(t1,3)*(2 + 3*Power(s,2) + s*(27 - 8*t2) + 
                  22*t2 + Power(t2,2)) - 
               Power(t1,2)*(-15 - 33*t2 + 22*Power(t2,2) + 
                  2*Power(t2,3) + Power(s,2)*(1 + 8*t2) + 
                  s*(-23 + 9*t2 - 10*Power(t2,2))) + 
               Power(s2,2)*(3 + 3*Power(t1,3) + 
                  Power(t1,2)*(5 - 9*t2) + 28*t2 + Power(t2,2) - 
                  2*Power(t2,3) + t1*(-2 - t2 + 8*Power(t2,2))) + 
               t1*(10 + t2 - 18*Power(t2,2) + 6*Power(t2,3) + 
                  Power(s,2)*t2*(4 + 7*t2) + 
                  s*(10 - 56*t2 + 2*Power(t2,2) - 3*Power(t2,3))) - 
               s2*(-14 + Power(t1,4) + Power(t1,3)*(1 + 4*s - 13*t2) + 
                  (4 + 6*s)*t2 + (6 - 9*s)*Power(t2,2) - 
                  (7 + 4*s)*Power(t2,3) + 
                  t1*(1 + (-22 + 15*s)*t2 + (16 + 15*s)*Power(t2,2) - 
                     3*Power(t2,3)) + 
                  Power(t1,2)*
                   (51 + 9*t2 + 14*Power(t2,2) - 3*s*(-4 + 5*t2)))) + 
            s1*(3 + 15*t2 + 20*Power(t2,2) - 17*s*Power(t2,2) + 
               3*Power(s,2)*Power(t2,2) - 15*Power(t2,3) - 
               34*s*Power(t2,3) + 10*Power(s,2)*Power(t2,3) - 
               7*Power(t2,4) + 3*s*Power(t2,4) - 2*Power(t2,5) + 
               Power(t1,4)*(-3 + s + 9*t2 - s*t2) + 
               Power(t1,2)*(2 + 18*t2 + 33*Power(t2,2) + 
                  24*Power(t2,3) + Power(s + 2*s*t2,2) + 
                  s*t2*(-5 + 15*t2 - 7*Power(t2,2))) + 
               Power(s2,2)*t2*
                (-Power(t1,3) + 2*Power(t1,2)*(-6 + t2) + 
                  t1*(25 + 13*t2 - 3*Power(t2,2)) + 
                  3*(-4 + 4*t2 + 3*Power(t2,2))) + 
               Power(t1,3)*(5 + Power(s,2)*(3 - 2*t2) + 7*t2 - 
                  24*Power(t2,2) + s*(5 - 21*t2 + 7*Power(t2,2))) - 
               t1*(7 + 49*t2 + 8*Power(t2,2) + 53*Power(t2,3) + 
                  7*Power(t2,4) + 
                  Power(s,2)*t2*(4 + 17*t2 + 2*Power(t2,2)) - 
                  s*(-6 + 81*t2 + 46*Power(t2,2) - 4*Power(t2,3) + 
                     Power(t2,4))) - 
               s2*(10 + (-1 + s)*Power(t1,4) - (11 + 2*s)*t2 + 
                  (11 + 51*s)*Power(t2,2) + (11 + 21*s)*Power(t2,3) + 
                  7*Power(t2,4) + 
                  t1*(-20 + (49 + 57*s)*t2 + 
                     (107 - 16*s)*Power(t2,2) - 
                     5*(3 + s)*Power(t2,3) + Power(t2,4)) + 
                  Power(t1,3)*
                   (2 - 7*t2 + 5*Power(t2,2) - 3*s*(1 + t2)) + 
                  Power(t1,2)*
                   (9 - 31*t2 - 3*Power(t2,2) - 6*Power(t2,3) + 
                     s*(2 - 18*t2 + 7*Power(t2,2))))) - 
            Power(s1,2)*(7 + 4*Power(t1,5) + 2*s*t2 - 2*Power(t2,2) - 
               25*s*Power(t2,2) + 5*Power(s,2)*Power(t2,2) - 
               6*Power(t2,3) + 6*s*Power(t2,3) - 
               Power(s,2)*Power(t2,3) + 3*Power(t2,4) - 
               3*s*Power(t2,4) - 
               Power(t1,4)*(-7 + Power(s,2) - 3*s*(-4 + t2) + 12*t2) + 
               Power(t1,3)*(3 - 11*t2 + 14*Power(t2,2) + 
                  Power(s,2)*(2 + 3*t2) + s*(4 + 8*t2 - 4*Power(t2,2))) \
+ Power(t1,2)*(-14 - 3*t2 + 3*Power(t2,2) - 2*Power(t2,3) + 
                  Power(s,2)*(7 - 5*t2 - 3*Power(t2,2)) + 
                  s*(15 + 79*t2 - 12*Power(t2,2) + Power(t2,3))) + 
               Power(s2,2)*(-5 + Power(t1,3)*(-7 + t2) + 8*t2 + 
                  30*Power(t2,2) + 2*Power(t2,3) + 
                  Power(t1,2)*(8 + 13*t2 - 2*Power(t2,2)) + 
                  t1*(4 + 5*t2 - 10*Power(t2,2) + Power(t2,3))) + 
               t1*(-7 + 47*t2 - 9*Power(t2,2) - 23*Power(t2,3) - 
                  4*Power(t2,4) + 
                  Power(s,2)*t2*(-12 + 4*t2 + Power(t2,2)) + 
                  s*(20 + 22*t2 - 61*Power(t2,2) + 13*Power(t2,3))) + 
               s2*(6 + Power(t1,4)*(2 + s - 2*t2) + (43 - 20*s)*t2 + 
                  (2 - 34*s)*Power(t2,2) - (6 + s)*Power(t2,3) + 
                  3*Power(t2,4) + 
                  Power(t1,3)*
                   (16 + s*(9 - 4*t2) - 3*t2 + 3*Power(t2,2)) - 
                  Power(t1,2)*
                   (15 + 48*t2 - 28*Power(t2,2) + Power(t2,3) + 
                     s*(25 + 2*t2 - 5*Power(t2,2))) - 
                  t1*(9 + 110*t2 + 9*Power(t2,2) + 21*Power(t2,3) + 
                     s*(2 + 19*t2 - 2*Power(t2,2) + 2*Power(t2,3)))))))/
        ((s - s2 + t1)*(-1 + s1 + t1 - t2)*(s - s1 + t2)*
          Power(-(s1*t1) + t2,2)) - 
       (4*((2*(-6 + 7*Power(t1,2) - Power(t1,3) - 2*t2 - 10*t1*t2 - 
                 6*s2*t1*t2 - 6*Power(t1,2)*t2 + s2*Power(t1,2)*t2 + 
                 Power(t1,3)*t2 - 13*t1*Power(t2,2) + 
                 4*s2*t1*Power(t2,2) + 2*Power(t1,2)*Power(t2,2) - 
                 s2*Power(t1,2)*Power(t2,2) + 10*Power(t2,3) - 
                 2*s2*Power(t2,3) + 3*t1*Power(t2,3) + 
                 s2*t1*Power(t2,3) - 2*Power(t2,4) + 
                 Power(s1,4)*(s2 + t1 - 2*t1*t2) + 
                 s*(-(Power(t1,3)*(-1 + t2)) + 2*t2*(-4 + 3*t2) + 
                    Power(t1,2)*(-6 + 5*t2 + 2*Power(t2,2)) - 
                    t1*(6 - 4*t2 + 2*Power(t2,2) + Power(t2,3))) - 
                 Power(s1,3)*
                  (1 + Power(t1,2)*(-2 + t2) + t2 - s*t2 - 
                    2*Power(t2,2) - 
                    t1*(4 + 3*s + 4*t2 + 2*s*t2 + 3*Power(t2,2)) + 
                    s2*(-9 + 5*t2 + Power(t2,2) + t1*(2 + t2))) - 
                 s1*(14 - 10*t2 + 19*Power(t2,2) + 6*s*Power(t2,2) + 
                    Power(t2,3) - 3*s*Power(t2,3) - 2*Power(t2,4) + 
                    Power(t1,3)*(s + 2*t2 + s*t2) + 
                    t1*(-16 + (13 - 2*s)*t2 + 
                       (13 + 4*s)*Power(t2,2) + (-1 + s)*Power(t2,3)) \
- Power(t1,2)*(15 + 11*t2 - 5*Power(t2,2) + 
                       2*s*(1 - t2 + Power(t2,2))) + 
                    s2*(-6 + Power(t1,2)*Power(-1 + t2,2) + 18*t2 - 
                       10*Power(t2,2) + Power(t2,3) - 
                       t1*t2*(-2 + 4*t2 + Power(t2,2)))) + 
                 Power(s1,2)*
                  (-9 + t2 + 5*s*t2 - 3*Power(t2,2) - 
                    6*s*Power(t2,2) - 3*Power(t2,3) - s*Power(t2,3) + 
                    Power(t1,3)*(1 + t2) + 
                    t1*(9 + 5*s + 11*t2 + Power(t2,2) - 
                       2*Power(t2,3)) + 
                    Power(t1,2)*
                     (4 + 6*t2 + Power(t2,2) + s*(2 + t2)) + 
                    s2*(14 - 17*t2 + 7*Power(t2,2) + Power(t2,3) - 
                       Power(t1,2)*(3 + t2) - 
                       t1*(-4 + t2 + 2*Power(t2,2))))))/
             ((-1 + t1)*(s1*t1 - t2)) + 
            (2*(-7 + Power(s1,4)*(s2 - t1) - 3*t2 - 15*s*t2 + s2*t2 - 
                 Power(t2,2) + 4*s*Power(t2,2) + 2*s2*Power(t2,2) + 
                 11*Power(t2,3) + 5*s*Power(t2,3) - s2*Power(t2,3) + 
                 Power(t1,2)*(3 + 2*t2 + s*(-3 + 2*t2)) + 
                 Power(s1,3)*
                  (-1 + 3*Power(t1,2) + t2 + s*t2 + 
                    t1*(8 + s + 4*t2 - Power(t2,2)) + 
                    s2*(9 + t1 - 5*t2 - 2*t1*t2 - Power(t2,2))) - 
                 t1*(-4 + (5 + 3*s2)*t2 + 19*Power(t2,2) + 
                    s*(7 - 3*t2 + 7*Power(t2,2))) + 
                 Power(s1,2)*
                  (-9 + (5 + s)*Power(t1,2) + 4*Power(t1,3) - 3*t2 + 
                    7*s*t2 - 3*Power(t2,2) - 4*s*Power(t2,2) + 
                    Power(t2,3) - s*Power(t2,3) + 
                    s2*(15 + t1*(12 - 5*t2) - 19*t2 - 
                       2*Power(t1,2)*t2 + 6*Power(t2,2) + Power(t2,3)\
) + t1*(1 + 9*t2 + 3*Power(t2,2) + s*(-1 + 7*t2 + Power(t2,2)))) + 
                 s1*(-15 - 4*Power(t1,3) + 17*t2 + 7*s*t2 - 
                    16*Power(t2,2) - 12*s*Power(t2,2) - 
                    4*Power(t2,3) + 
                    Power(t1,2)*(5 + 14*t2 + s*(2 + 6*t2)) + 
                    t1*(20 - 12*t2 - 11*Power(t2,2) + 
                       s*(7 + 2*t2 - 6*Power(t2,2))) + 
                    s2*(7 - 29*t2 + 2*Power(t1,2)*t2 + 
                       13*Power(t2,2) + t1*(-5 - 2*t2 + 4*Power(t2,2)))\
)))/((-1 + s1)*(s1*t1 - t2)) + 
            ((-1 + t2)*((-2*
                    (-14 + 2*Power(s1,3) - 8*t1 + 2*s2*t1 - 
                      3*Power(t1,2) + Power(t1,3) + 
                      Power(s1,2)*(-1 + s2 + 5*t1 - 3*t2) + 8*t2 + 
                      2*s2*t2 - s2*t1*t2 - Power(t1,2)*t2 + 
                      2*Power(t2,2) - 
                      s*(-9 + Power(s1,2) + s1*(5 + t1 - 3*t2) - 
                       2*t1*(-1 + t2) + 7*t2) + 
                      s1*(-18 + 4*Power(t1,2) + 
                       s2*(5 + t1 - 2*t2) + 7*t2 - Power(t2,2) - 
                       4*t1*(2 + t2)))*
                    (-1 + Power(s1,2)*(s2 - t1) + t1 - s*t1 + t2 - 
                      s*t2 - s2*t2 + 
                      s1*(-1 + s2 + t2 - s2*t2 + s*(t1 + t2))))/
                  (s1*t1 - t2) + 
                 4*(-2 - Power(s1,4) + s2 + 5*t1 - s2*t1 + 
                    6*Power(t1,2) - Power(t1,3) + 
                    Power(s,2)*
                     (5 + s1*(-3 + t2) + t1*(-1 + t2) - 3*t2) - 
                    s2*t2 + 2*Power(s2,2)*t2 - 4*t1*t2 - 4*s2*t1*t2 + 
                    4*Power(t1,2)*t2 - 2*Power(t2,2) + 
                    s2*Power(t2,2) - t1*Power(t2,2) + 
                    Power(s1,3)*(2 - 3*t1 + t2) + 
                    Power(s1,2)*
                     (12 - 6*s2 - 3*Power(t1,2) - 4*t2 + 
                       Power(t2,2) + t1*(9 + t2)) + 
                    s1*(15 - 2*Power(s2,2) + 3*Power(t1,2) - 
                       Power(t1,3) - 9*t2 - Power(t2,2) + 
                       s2*(1 - t1 - 3*t2 + Power(t2,2)) + 
                       t1*(15 + t2 + Power(t2,2))) + 
                    s*(-6 + Power(s1,3) - t1 - 4*Power(t1,2) + 
                       Power(t1,3) + Power(s1,2)*(2 + 3*t1 - 2*t2) + 
                       7*t2 + 5*t1*t2 - Power(t1,2)*t2 + 
                       s2*(-5 + t1 + t2 - t1*t2) - 
                       s1*(10 + 7*t1 - 3*Power(t1,2) + s2*(-5 + t2) - 
                        9*t2 + 3*t1*t2 + Power(t2,2)))) - 
                 (2*(-9 + Power(s1,2) + t1 + s1*(1 + t1 - 2*t2) + 
                      4*t2 - t1*t2)*
                    (-1 + Power(s1,5)*Power(s2 - t1,2) + 
                      Power(-1 + s,2)*Power(t1,3) + 3*t2 - 2*s*t2 - 
                      2*s2*t2 + 2*s*s2*t2 - 3*Power(t2,2) + 
                      4*s*Power(t2,2) - Power(s,2)*Power(t2,2) + 
                      4*s2*Power(t2,2) - 6*s*s2*Power(t2,2) - 
                      Power(s2,2)*Power(t2,2) + Power(t2,3) - 
                      2*s*Power(t2,3) + Power(s,2)*Power(t2,3) - 
                      2*s2*Power(t2,3) - 4*s*s2*Power(t2,3) + 
                      Power(s2,2)*Power(t2,3) + 
                      Power(t1,2)*
                       (-3 + Power(s,2)*(-1 + t2) + (3 - 2*s2)*t2 + 
                        2*s*(2 + (-1 + s2)*t2)) + 
                      t1*(3 + (-6 + 4*s2)*t2 - 
                        3*Power(s,2)*Power(t2,2) + 
                        (3 - 4*s2 + Power(s2,2))*Power(t2,2) + 
                        2*s*
                        (-1 - 2*(-1 + s2)*t2 + (3 + s2)*Power(t2,2))) \
+ Power(s1,4)*(Power(s2,2)*(1 + t1 - 3*t2) + 
                        t1*(2 + Power(t1,2) - 2*s*(1 + t1) - 2*t2 - 
                       t1*(1 + t2)) + 
                        2*s2*
                        (-1 - Power(t1,2) + t2 + 2*t1*t2 + 
                        s*(t1 + t2))) + 
                      Power(s1,3)*
                       (1 - 2*Power(t1,3) - 2*t2 - 2*t1*t2 - 
                        2*Power(t1,2)*t2 + Power(t2,2) + 
                        2*t1*Power(t2,2) - 
                        2*s*t1*
                       (1 + t1 + Power(t1,2) - 4*t2 - t1*t2) + 
                        Power(s,2)*(Power(t1,2) + Power(t2,2)) + 
                        Power(s2,2)*
                       (-1 - 2*(2 + t1)*t2 + 3*Power(t2,2)) + 
                        2*s2*
                        (-1 - (-3 + s)*t2 - 2*(1 + s)*Power(t2,2) + 
                        Power(t1,2)*(1 + s + t2) + 
                        t1*(1 + 2*t2 - Power(t2,2)))) + 
                      s1*(-(Power(-1 + t1 + t2,2)*
                       (1 + 2*t1 + 2*t2)) - 
                        Power(s2,2)*t2*
                       (-2 + 2*t1 + 3*t2 + 2*Power(t2,2)) + 
                        2*s*
                       (t1 + Power(t1,3) - 8*t1*t2 + 
                       (-1 + t2)*Power(t2,2) - 
                       2*Power(t1,2)*(1 + t2)) + 
                        Power(s,2)*
                       (Power(t1,2) - 2*Power(t1,3) + Power(t2,2) + 
                       2*t1*t2*(2 + t2)) + 
                        2*s2*
                        (1 + (-3 + s)*t2 + 5*s*Power(t2,2) + 
                        (2 + s)*Power(t2,3) + 
                        Power(t1,2)*(1 + t2 - s*t2) + 
                        t1*(-2 + (2 + 4*s)*t2 + 3*Power(t2,2)))) + 
                      Power(s1,2)*
                       (1 + (2 + 2*s + Power(s,2))*Power(t1,3) - 
                        3*t2 + 2*s*t2 + 3*Power(t2,2) - 
                        2*s*Power(t2,2) - Power(s,2)*Power(t2,2) - 
                        Power(t2,3) - Power(s,2)*Power(t2,3) + 
                        Power(s2,2)*
                        (-1 + t1 + 3*t2 + 5*Power(t2,2) + 
                        t1*Power(t2,2) - Power(t2,3)) - 
                        Power(t1,2)*
                        (-6*t2 - 4*s*(1 + t2) + Power(s,2)*(1 + t2)) \
+ t1*(Power(s,2)*(-4 + t2)*t2 + s*(4 + 4*t2 - 6*Power(t2,2)) + 
                        3*(-1 + Power(t2,2))) - 
                        2*s2*
                        (-1 - 2*t2 + 4*Power(t2,2) - Power(t2,3) + 
                        Power(t1,2)*(1 + t2) + 2*t1*t2*(2 + t2) + 
                        s*(Power(t1,2) + t1*Power(1 + t2,2) - 
                        t2*(-2 + Power(t2,2)))))))/
                  ((-1 + s1 + t1 - t2)*Power(-(s1*t1) + t2,2))))/
             ((-1 + s2)*(-s + s2 - t1)) + 
            ((-1 + t2)*((-2*
                    (-1 + Power(s1,2)*(s2 - t1) + t1 - s*t1 + t2 - 
                      s*t2 - s2*t2 + 
                      s1*(-1 + s2 + t2 - s2*t2 + s*(t1 + t2)))*
                    (-17 + 2*Power(s1,3) - 5*t1 + 2*s2*t1 - 
                      2*Power(t1,2) + 
                      Power(s1,2)*(-4 + s2 + 2*t1 - t2) + 4*t2 + 
                      2*s2*t2 - s2*t1*t2 + Power(t1,2)*t2 + 
                      Power(t2,2) - t1*Power(t2,2) - 
                      s*(-9 + Power(s1,2) + s1*(5 + t1 - 3*t2) - 
                       2*t1*(-1 + t2) + 7*t2) + 
                      s1*(s2*(5 + t1 - 2*t2) + 2*t1*(-4 + t2) - 
                        3*(4 - 3*t2 + Power(t2,2)))))/(s1*t1 - t2) + 
                 4*(-1 + 2*s2 + 8*t1 - 4*s2*t1 + 9*Power(t1,2) + 
                    Power(s,2)*
                     (5 + s1*(-3 + t2) + t1*(-1 + t2) - 3*t2) + 
                    Power(s1,3)*(5 - 2*t2) + 13*t2 - 5*s2*t2 + 
                    2*Power(s2,2)*t2 - 10*t1*t2 - 4*s2*t1*t2 - 
                    3*Power(t1,2)*t2 - 3*Power(t2,2) + 
                    2*s2*Power(t2,2) + 6*t1*Power(t2,2) - 
                    Power(t2,3) + 
                    Power(s1,2)*
                     (13 - 3*s2 + t1*(13 - 3*t2) - 9*t2 + 
                       2*Power(t2,2)) + 
                    s1*(2 - 2*Power(s2,2) + t1*(19 - 9*t2) - 
                       Power(t1,2)*(-7 + t2) - 12*t2 + 
                       s2*(-5 + t2)*t2 + 4*Power(t2,2) + Power(t2,3)) \
- s*(3 + Power(s1,2) + 2*t1 + 6*Power(t1,2) + 
                       s2*(5 + t1*(-1 + t2) - t2) - 3*t2 - 8*t1*t2 - 
                       Power(t1,2)*t2 + Power(t2,2) + 
                       t1*Power(t2,2) + 
                       s1*(9 + 10*t1 + s2*(-5 + t2) - 14*t2 - t1*t2 + 
                        2*Power(t2,2)))) - 
                 (2*(-9 + Power(s1,2) + t1 + s1*(1 + t1 - 2*t2) + 
                      4*t2 - t1*t2)*
                    (-1 + Power(s1,5)*Power(s2 - t1,2) + 
                      Power(-1 + s,2)*Power(t1,3) + 3*t2 - 2*s*t2 - 
                      2*s2*t2 + 2*s*s2*t2 - 3*Power(t2,2) + 
                      4*s*Power(t2,2) - Power(s,2)*Power(t2,2) + 
                      4*s2*Power(t2,2) - 6*s*s2*Power(t2,2) - 
                      Power(s2,2)*Power(t2,2) + Power(t2,3) - 
                      2*s*Power(t2,3) + Power(s,2)*Power(t2,3) - 
                      2*s2*Power(t2,3) - 4*s*s2*Power(t2,3) + 
                      Power(s2,2)*Power(t2,3) + 
                      Power(t1,2)*
                       (-3 + Power(s,2)*(-1 + t2) + (3 - 2*s2)*t2 + 
                        2*s*(2 + (-1 + s2)*t2)) + 
                      t1*(3 + (-6 + 4*s2)*t2 - 
                        3*Power(s,2)*Power(t2,2) + 
                        (3 - 4*s2 + Power(s2,2))*Power(t2,2) + 
                        2*s*
                        (-1 - 2*(-1 + s2)*t2 + (3 + s2)*Power(t2,2))) \
+ Power(s1,4)*(Power(s2,2)*(1 + t1 - 3*t2) + 
                        t1*(2 + Power(t1,2) - 2*s*(1 + t1) - 2*t2 - 
                       t1*(1 + t2)) + 
                        2*s2*
                        (-1 - Power(t1,2) + t2 + 2*t1*t2 + 
                        s*(t1 + t2))) + 
                      Power(s1,3)*
                       (1 - 2*Power(t1,3) - 2*t2 - 2*t1*t2 - 
                        2*Power(t1,2)*t2 + Power(t2,2) + 
                        2*t1*Power(t2,2) - 
                        2*s*t1*
                       (1 + t1 + Power(t1,2) - 4*t2 - t1*t2) + 
                        Power(s,2)*(Power(t1,2) + Power(t2,2)) + 
                        Power(s2,2)*
                       (-1 - 2*(2 + t1)*t2 + 3*Power(t2,2)) + 
                        2*s2*
                        (-1 - (-3 + s)*t2 - 2*(1 + s)*Power(t2,2) + 
                        Power(t1,2)*(1 + s + t2) + 
                        t1*(1 + 2*t2 - Power(t2,2)))) + 
                      s1*(-(Power(-1 + t1 + t2,2)*
                       (1 + 2*t1 + 2*t2)) - 
                        Power(s2,2)*t2*
                       (-2 + 2*t1 + 3*t2 + 2*Power(t2,2)) + 
                        2*s*
                       (t1 + Power(t1,3) - 8*t1*t2 + 
                       (-1 + t2)*Power(t2,2) - 
                       2*Power(t1,2)*(1 + t2)) + 
                        Power(s,2)*
                       (Power(t1,2) - 2*Power(t1,3) + Power(t2,2) + 
                       2*t1*t2*(2 + t2)) + 
                        2*s2*
                        (1 + (-3 + s)*t2 + 5*s*Power(t2,2) + 
                        (2 + s)*Power(t2,3) + 
                        Power(t1,2)*(1 + t2 - s*t2) + 
                        t1*(-2 + (2 + 4*s)*t2 + 3*Power(t2,2)))) + 
                      Power(s1,2)*
                       (1 + (2 + 2*s + Power(s,2))*Power(t1,3) - 
                        3*t2 + 2*s*t2 + 3*Power(t2,2) - 
                        2*s*Power(t2,2) - Power(s,2)*Power(t2,2) - 
                        Power(t2,3) - Power(s,2)*Power(t2,3) + 
                        Power(s2,2)*
                        (-1 + t1 + 3*t2 + 5*Power(t2,2) + 
                        t1*Power(t2,2) - Power(t2,3)) - 
                        Power(t1,2)*
                        (-6*t2 - 4*s*(1 + t2) + Power(s,2)*(1 + t2)) \
+ t1*(Power(s,2)*(-4 + t2)*t2 + s*(4 + 4*t2 - 6*Power(t2,2)) + 
                        3*(-1 + Power(t2,2))) - 
                        2*s2*
                        (-1 - 2*t2 + 4*Power(t2,2) - Power(t2,3) + 
                        Power(t1,2)*(1 + t2) + 2*t1*t2*(2 + t2) + 
                        s*(Power(t1,2) + t1*Power(1 + t2,2) - 
                        t2*(-2 + Power(t2,2)))))))/
                  ((-1 + s1 + t1 - t2)*Power(-(s1*t1) + t2,2))))/
             ((s - s2 + t1)*(s - s1 + t2)) + 
            ((-1 + t2)*((2*(-1 + Power(s1,2)*(s2 - t1) + t1 - 
                      s*t1 + t2 - s*t2 - s2*t2 + 
                      s1*(-1 + s2 + t2 - s2*t2 + s*(t1 + t2)))*
                    (11 - 3*s2 + 6*t1 - 5*Power(t1,2) + 
                      Power(t1,3) + Power(s1,2)*(-6 + t1 - t2) + 
                      4*s2*t2 + 3*t1*t2 - s2*t1*t2 - 
                      Power(t1,2)*t2 + Power(t2,2) + 
                      s*(2 + 2*t1*(-1 + t2) - 5*t2 + s1*(2 + t2)) - 
                      s1*(-9 - 2*Power(t1,2) - 10*t2 + Power(t2,2) + 
                        s2*(3 + t2) + t1*(7 + 3*t2))))/(s1*t1 - t2) - 
                 4*(1 - 2*s2 + Power(s1,3)*(3 + s2) - 5*t1 - 
                    2*s2*t1 - Power(t1,2) + s2*Power(t1,2) - 
                    Power(t1,3) + 
                    Power(s,2)*(2 + s1 + t1*(-1 + t2) - 3*t2) + 
                    6*t2 + 6*s2*t2 - Power(s2,2)*t2 - 8*t1*t2 + 
                    4*Power(t1,2)*t2 - Power(t2,2) - t1*Power(t2,2) + 
                    Power(s1,2)*
                     (-6 + 5*t1 - Power(t1,2) + s2*(2 + 2*t1 - t2) - 
                       6*t2 + Power(t2,2)) + 
                    s1*(-14 + Power(s2,2) + 2*Power(t1,2) - 
                       Power(t1,3) - 3*t2 - Power(t2,2) + 
                       s2*(-6 + Power(t1,2) - 3*t2 + Power(t2,2)) + 
                       t1*(-4 - 2*t2 + Power(t2,2))) + 
                    s*(5 - 5*Power(s1,2) + 6*t1 - 4*Power(t1,2) + 
                       Power(t1,3) - 7*t2 + 5*t1*t2 - 
                       Power(t1,2)*t2 + s2*(-2 + t1 + 4*t2 - t1*t2) + 
                       s1*(6 - 2*s2 + Power(t1,2) + 10*t2 - 
                        Power(t2,2) - 2*t1*(3 + t2)))) - 
                 (2*(t1*(-1 + t2) - 2*t2 + s1*(1 + t2))*
                    (-1 + Power(s1,5)*Power(s2 - t1,2) + 
                      Power(-1 + s,2)*Power(t1,3) + 3*t2 - 2*s*t2 - 
                      2*s2*t2 + 2*s*s2*t2 - 3*Power(t2,2) + 
                      4*s*Power(t2,2) - Power(s,2)*Power(t2,2) + 
                      4*s2*Power(t2,2) - 6*s*s2*Power(t2,2) - 
                      Power(s2,2)*Power(t2,2) + Power(t2,3) - 
                      2*s*Power(t2,3) + Power(s,2)*Power(t2,3) - 
                      2*s2*Power(t2,3) - 4*s*s2*Power(t2,3) + 
                      Power(s2,2)*Power(t2,3) + 
                      Power(t1,2)*
                       (-3 + Power(s,2)*(-1 + t2) + (3 - 2*s2)*t2 + 
                        2*s*(2 + (-1 + s2)*t2)) + 
                      t1*(3 + (-6 + 4*s2)*t2 - 
                        3*Power(s,2)*Power(t2,2) + 
                        (3 - 4*s2 + Power(s2,2))*Power(t2,2) + 
                        2*s*
                        (-1 - 2*(-1 + s2)*t2 + (3 + s2)*Power(t2,2))) \
+ Power(s1,4)*(Power(s2,2)*(1 + t1 - 3*t2) + 
                        t1*(2 + Power(t1,2) - 2*s*(1 + t1) - 2*t2 - 
                       t1*(1 + t2)) + 
                        2*s2*
                        (-1 - Power(t1,2) + t2 + 2*t1*t2 + 
                        s*(t1 + t2))) + 
                      Power(s1,3)*
                       (1 - 2*Power(t1,3) - 2*t2 - 2*t1*t2 - 
                        2*Power(t1,2)*t2 + Power(t2,2) + 
                        2*t1*Power(t2,2) - 
                        2*s*t1*
                       (1 + t1 + Power(t1,2) - 4*t2 - t1*t2) + 
                        Power(s,2)*(Power(t1,2) + Power(t2,2)) + 
                        Power(s2,2)*
                       (-1 - 2*(2 + t1)*t2 + 3*Power(t2,2)) + 
                        2*s2*
                        (-1 - (-3 + s)*t2 - 2*(1 + s)*Power(t2,2) + 
                        Power(t1,2)*(1 + s + t2) + 
                        t1*(1 + 2*t2 - Power(t2,2)))) + 
                      s1*(-(Power(-1 + t1 + t2,2)*
                       (1 + 2*t1 + 2*t2)) - 
                        Power(s2,2)*t2*
                       (-2 + 2*t1 + 3*t2 + 2*Power(t2,2)) + 
                        2*s*
                       (t1 + Power(t1,3) - 8*t1*t2 + 
                       (-1 + t2)*Power(t2,2) - 
                       2*Power(t1,2)*(1 + t2)) + 
                        Power(s,2)*
                       (Power(t1,2) - 2*Power(t1,3) + Power(t2,2) + 
                       2*t1*t2*(2 + t2)) + 
                        2*s2*
                        (1 + (-3 + s)*t2 + 5*s*Power(t2,2) + 
                        (2 + s)*Power(t2,3) + 
                        Power(t1,2)*(1 + t2 - s*t2) + 
                        t1*(-2 + (2 + 4*s)*t2 + 3*Power(t2,2)))) + 
                      Power(s1,2)*
                       (1 + (2 + 2*s + Power(s,2))*Power(t1,3) - 
                        3*t2 + 2*s*t2 + 3*Power(t2,2) - 
                        2*s*Power(t2,2) - Power(s,2)*Power(t2,2) - 
                        Power(t2,3) - Power(s,2)*Power(t2,3) + 
                        Power(s2,2)*
                        (-1 + t1 + 3*t2 + 5*Power(t2,2) + 
                        t1*Power(t2,2) - Power(t2,3)) - 
                        Power(t1,2)*
                        (-6*t2 - 4*s*(1 + t2) + Power(s,2)*(1 + t2)) \
+ t1*(Power(s,2)*(-4 + t2)*t2 + s*(4 + 4*t2 - 6*Power(t2,2)) + 
                        3*(-1 + Power(t2,2))) - 
                        2*s2*
                        (-1 - 2*t2 + 4*Power(t2,2) - Power(t2,3) + 
                        Power(t1,2)*(1 + t2) + 2*t1*t2*(2 + t2) + 
                        s*(Power(t1,2) + t1*Power(1 + t2,2) - 
                        t2*(-2 + Power(t2,2)))))))/
                  ((-1 + s1 + t1 - t2)*Power(-(s1*t1) + t2,2))))/
             ((-1 + s2)*(-1 + t1)) + 
            ((-1 + t2)*((-2*(-1 + Power(s1,2)*(s2 - t1) + t1 - 
                      s*t1 + t2 - s*t2 - s2*t2 + 
                      s1*(-1 + s2 + t2 - s2*t2 + s*(t1 + t2)))*
                    (-9 + 3*s2 - 4*t1 - s2*t1 + 3*Power(t1,2) - 
                      3*t2 - 4*t1*t2 + 2*Power(t2,2) + 
                      Power(s1,2)*(5 + s2 + 3*t2) + 
                      s1*(-5 + 12*t1 - 17*t2 + Power(t2,2) + 
                        s2*(-5 + t1 + t2)) - 
                      s*(7 + Power(s1,2) + t1 - 4*t2 + 
                        s1*(-4 + t1 + 2*t2))))/(s1*t1 - t2) - 
                 4*(-2 - Power(s1,4) + 3*s2 - 11*t1 - s2*t1 + 
                    Power(t1,2) + Power(s,2)*(-3 + s1)*(-1 + t2) + 
                    3*t2 + 10*s2*t2 - 10*t1*t2 + 2*s2*t1*t2 + 
                    Power(t2,2) + Power(s1,3)*(2 - t1 + 2*t2) + 
                    s1*(-3 + t1 - Power(t1,2) + 
                       s2*(-12 + t1*(-2 + t2) - t2) - 7*t2 + 
                       3*t1*t2 - 2*Power(t2,2)) + 
                    Power(s1,2)*
                     (-5 + 2*t1 + 2*s2*(-1 + t2) - 5*t2 + Power(t2,2)) \
+ s*(-2 + Power(s1,3) + 8*t1 + Power(s1,2)*(-2 + t1 - 3*t2) + 
                       3*s2*(-1 + t2) - 11*t2 - t1*t2 + 
                       s1*(11 + s2 - 4*t1 + 11*t2 - s2*t2 - 
                        Power(t2,2)))) - 
                 (2*(Power(s1,2) + 2*(t1 - t2) + s1*(-3 + t1 + t2))*
                    (-1 + Power(s1,5)*Power(s2 - t1,2) + 
                      Power(-1 + s,2)*Power(t1,3) + 3*t2 - 2*s*t2 - 
                      2*s2*t2 + 2*s*s2*t2 - 3*Power(t2,2) + 
                      4*s*Power(t2,2) - Power(s,2)*Power(t2,2) + 
                      4*s2*Power(t2,2) - 6*s*s2*Power(t2,2) - 
                      Power(s2,2)*Power(t2,2) + Power(t2,3) - 
                      2*s*Power(t2,3) + Power(s,2)*Power(t2,3) - 
                      2*s2*Power(t2,3) - 4*s*s2*Power(t2,3) + 
                      Power(s2,2)*Power(t2,3) + 
                      Power(t1,2)*
                       (-3 + Power(s,2)*(-1 + t2) + (3 - 2*s2)*t2 + 
                        2*s*(2 + (-1 + s2)*t2)) + 
                      t1*(3 + (-6 + 4*s2)*t2 - 
                        3*Power(s,2)*Power(t2,2) + 
                        (3 - 4*s2 + Power(s2,2))*Power(t2,2) + 
                        2*s*(-1 - 2*(-1 + s2)*t2 + 
                        (3 + s2)*Power(t2,2))) + 
                      Power(s1,4)*
                       (Power(s2,2)*(1 + t1 - 3*t2) + 
                        t1*(2 + Power(t1,2) - 2*s*(1 + t1) - 2*t2 - 
                        t1*(1 + t2)) + 
                        2*s2*
                        (-1 - Power(t1,2) + t2 + 2*t1*t2 + 
                        s*(t1 + t2))) + 
                      Power(s1,3)*
                       (1 - 2*Power(t1,3) - 2*t2 - 2*t1*t2 - 
                        2*Power(t1,2)*t2 + Power(t2,2) + 
                        2*t1*Power(t2,2) - 
                        2*s*t1*
                       (1 + t1 + Power(t1,2) - 4*t2 - t1*t2) + 
                        Power(s,2)*(Power(t1,2) + Power(t2,2)) + 
                        Power(s2,2)*
                        (-1 - 2*(2 + t1)*t2 + 3*Power(t2,2)) + 
                        2*s2*
                        (-1 - (-3 + s)*t2 - 2*(1 + s)*Power(t2,2) + 
                        Power(t1,2)*(1 + s + t2) + 
                        t1*(1 + 2*t2 - Power(t2,2)))) + 
                      s1*(-(Power(-1 + t1 + t2,2)*
                       (1 + 2*t1 + 2*t2)) - 
                        Power(s2,2)*t2*
                        (-2 + 2*t1 + 3*t2 + 2*Power(t2,2)) + 
                        2*s*
                        (t1 + Power(t1,3) - 8*t1*t2 + 
                        (-1 + t2)*Power(t2,2) - 
                        2*Power(t1,2)*(1 + t2)) + 
                        Power(s,2)*
                        (Power(t1,2) - 2*Power(t1,3) + Power(t2,2) + 
                        2*t1*t2*(2 + t2)) + 
                        2*s2*
                        (1 + (-3 + s)*t2 + 5*s*Power(t2,2) + 
                        (2 + s)*Power(t2,3) + 
                        Power(t1,2)*(1 + t2 - s*t2) + 
                        t1*(-2 + (2 + 4*s)*t2 + 3*Power(t2,2)))) + 
                      Power(s1,2)*
                       (1 + (2 + 2*s + Power(s,2))*Power(t1,3) - 
                        3*t2 + 2*s*t2 + 3*Power(t2,2) - 
                        2*s*Power(t2,2) - Power(s,2)*Power(t2,2) - 
                        Power(t2,3) - Power(s,2)*Power(t2,3) + 
                        Power(s2,2)*
                        (-1 + t1 + 3*t2 + 5*Power(t2,2) + 
                        t1*Power(t2,2) - Power(t2,3)) - 
                        Power(t1,2)*
                        (-6*t2 - 4*s*(1 + t2) + Power(s,2)*(1 + t2)) + 
                        t1*(Power(s,2)*(-4 + t2)*t2 + 
                        s*(4 + 4*t2 - 6*Power(t2,2)) + 
                        3*(-1 + Power(t2,2))) - 
                        2*s2*
                        (-1 - 2*t2 + 4*Power(t2,2) - Power(t2,3) + 
                        Power(t1,2)*(1 + t2) + 2*t1*t2*(2 + t2) + 
                        s*(Power(t1,2) + t1*Power(1 + t2,2) - 
                        t2*(-2 + Power(t2,2)))))))/
                  ((-1 + s1 + t1 - t2)*Power(-(s1*t1) + t2,2))))/
             ((-1 + s1)*(-s + s1 - t2))))/Power(-1 + t2,2))*
     B1(s1,1 - s1 - t1 + t2,t2))/(16.*Power(Pi,2)) + 
  (((-4*(-1 + s2 - t1 + t2)*(1 + (-1 + s + s1)*t1 - t2 + s*t2 + 
            Power(s2,2)*(-s1 + t2) + 
            s2*(1 + s1*(-1 + t1) - (1 + s)*t1 - s*t2))*
          ((-3 - 4*s2 - 5*Power(s2,2) + 3*Power(s2,3) + 
               Power(s1,2)*(3 - 2*s2*(-2 + t1) - 3*t1) + 
               2*Power(s,3)*(-1 + t1) + 5*t1 + 16*s2*t1 + 
               4*Power(s2,2)*t1 - 2*Power(s2,3)*t1 + 2*Power(t1,2) - 
               2*s2*Power(t1,2) + 
               Power(s,2)*(3 - 4*s1*(-1 + t1) - 6*s2*(-1 + t1) + 
                  4*t1 - t2) - t2 + s2*t2 - Power(s2,2)*t2 + 
               s2*t1*t2 + s1*
                (2 + Power(s2,2)*(7 - 4*t1) - 2*Power(t1,2) + 
                  s2*(1 + t1 - t2) + t1*(11 + t2)) + 
               s*(2*Power(s1,2)*(-1 + t1) + Power(s2,2)*(-7 + 6*t1) + 
                  t1*(-15 + 2*t1 - t2) + 
                  s1*(-6 - t1 + 2*s2*(-5 + 4*t1) + t2) + 
                  s2*(2 - 9*t1 + 2*t2)))/((-1 + s2)*(-s + s2 - t1)) + 
            (5*s2 + 5*Power(s2,2) + 4*Power(s2,3) - 4*s2*t1 - 
               2*Power(s2,3)*t1 - 2*Power(t1,2) + 2*s2*Power(t1,2) + 
               Power(s1,2)*(-4 + s2 + t1 - 2*s2*t1) + 11*s2*t2 + 
               2*Power(s2,2)*t2 - s2*t1*t2 + 
               Power(s,2)*(3 - 2*s2*(-2 + t1) - t1 + 4*t2) + 
               s*(-5 + 4*Power(s2,2)*(-2 + t1) + 8*t1 - 
                  2*Power(t1,2) + s2*(-10 + t1 - 6*t2) + 
                  s1*(1 + s2*(-5 + 4*t1) - 4*t2) - 11*t2 + t1*t2) + 
               s1*(5 + Power(s2,2)*(5 - 4*t1) + 2*Power(t1,2) + 
                  11*t2 - t1*(8 + t2) + s2*(3 + t1 + 2*t2)))/
             ((-1 + s2)*(-1 + t1)) + 
            (-2 + s2 + 12*Power(s2,2) + 6*Power(s2,3) - 
               2*Power(s2,4) + t1 - 6*s2*t1 + Power(s2,2)*t1 - 
               Power(t1,2) - 2*t2 + 3*s2*t2 + 8*Power(s2,2)*t2 - 
               2*Power(s2,3)*t2 - 3*t1*t2 + 
               Power(s,2)*(2 + s2 - 2*Power(s2,2) + 2*t1 + t2 - 
                  2*s2*t2) + 
               Power(s1,2)*(1 + s2 - 2*Power(s2,2) + 2*t1 + 2*t2 - 
                  2*s2*t2) + 
               s1*(3 - 4*Power(s2,3) + Power(t1,2) + 
                  Power(s2,2)*(7 - 4*t2) - t2 - t1*(4 + t2) + 
                  s2*(7 + 2*t1 + 10*t2)) + 
               s*(-2 + 4*Power(s2,3) + 3*t1 - Power(t1,2) + t2 + 
                  t1*t2 + Power(s2,2)*(-7 + 4*t2) - 
                  s2*(8 + 2*t1 + 9*t2) + 
                  s1*(-3 + 4*Power(s2,2) - 4*t1 - 3*t2 + 
                     s2*(-2 + 4*t2))))/((-1 + t1)*(-1 + t2)) + 
            (-6 + 2*Power(s1,3) - 16*s2 + 5*Power(s2,2) + 
               4*Power(s2,3) - 2*Power(s2,4) + 5*t1 + 14*s2*t1 + 
               5*Power(s2,2)*t1 + Power(t1,2) - 2*t2 + 8*s2*t2 - 
               3*Power(s2,2)*t2 - 2*Power(s2,3)*t2 + 3*t1*t2 + 
               2*Power(s,3)*(-3 + s2 + t2) - 
               Power(s1,2)*(-1 + 2*Power(s2,2) + 3*t1 + 
                  2*s2*(-5 + t2) + 2*t2) + 
               s1*(-4*Power(s2,3) - Power(t1,2) + 
                  Power(s2,2)*(13 - 4*t2) + 2*s2*(3 + t1 - 2*t2) + 
                  7*(-1 + t2) + t1*(5 + t2)) - 
               2*Power(s,2)*
                (2 + 3*Power(s2,2) - 2*t1 + 3*s2*(-3 + t2) + t2 + 
                  2*s1*(-4 + s2 + t2)) + 
               s*(14 + 6*Power(s2,3) - 11*t1 - 7*t2 - t1*t2 + 
                  2*Power(s1,2)*(-6 + s2 + t2) + s2*(-10*t1 + 4*t2) + 
                  Power(s2,2)*(-17 + 6*t2) + 
                  s1*(2 + 8*Power(s2,2) + 4*t2 + s2*(-30 + 8*t2))))/
             ((s - s2 + t1)*(s - s1 + t2)) + 
            (-4 + 4*Power(s1,3) - 11*s2 + 10*Power(s2,2) + 
               10*Power(s2,3) - 2*Power(s2,4) + 3*t1 - s2*t1 - 
               2*Power(s2,2)*t1 + Power(t1,2) + 2*s2*Power(t1,2) - 
               5*t2 + 15*s2*t2 + 6*Power(s2,2)*t2 - 
               2*Power(s2,3)*t2 - 2*s2*t1*t2 + Power(t2,2) + 
               s2*Power(t2,2) + 2*Power(s,3)*(-3 + s2 + t2) - 
               2*Power(s1,2)*(4 + Power(s2,2) + s2*(-5 + t2) + t2) + 
               s1*(-3 - 4*Power(s2,3) + Power(t1,2) - 
                  s2*(4 + t1 - 5*t2) - 4*Power(s2,2)*(-4 + t2) + 
                  9*t2 + 2*Power(t2,2) - t1*(3 + 2*t2)) - 
               2*Power(s,2)*
                (-1 + 3*Power(s2,2) + 3*s2*(-3 + t2) - t2 + 
                  s1*(-7 + 2*s2 + 2*t2)) + 
               s*(6 + 6*Power(s2,3) + t1 - 2*Power(t1,2) + 
                  s2*(-7 + 2*t1 - 9*t2) - 9*t2 + 2*t1*t2 - 
                  2*Power(t2,2) + 2*Power(s1,2)*(-6 + s2 + t2) + 
                  Power(s2,2)*(-22 + 6*t2) + 
                  s1*(7 + 8*Power(s2,2) - t1 + s2*(-26 + 8*t2))))/
             ((-1 + s1)*(-s + s1 - t2)) + 
            (-2 + 2*Power(s1,3) - 2*s2 + 8*Power(s2,2) + 
               9*Power(s2,3) - 2*Power(s2,4) + t1 - s2*t1 + 
               2*Power(s2,2)*t1 - Power(t1,2) - 2*s2*Power(t1,2) - 
               5*t2 + 2*s2*t2 + 9*Power(s2,2)*t2 - 2*Power(s2,3)*t2 + 
               2*s2*t1*t2 - Power(t2,2) - s2*Power(t2,2) + 
               2*Power(s,3)*(-1 + s2 + t2) + 
               Power(s,2)*(2 - 6*Power(s2,2) + s2*(8 - 6*t2) + 
                  s1*(6 - 4*s2 - 4*t2) + 4*t2) + 
               Power(s1,2)*(-5 - 2*Power(s2,2) + s2*(7 - 2*t2) + 
                  4*t2) - s1*
                (-2 + 4*Power(s2,3) + t1 + Power(t1,2) + 
                  s2*(3 - 3*t1 - 10*t2) + 4*Power(s2,2)*(-3 + t2) - 
                  2*t1*t2 + 2*Power(t2,2)) + 
               s*(6*Power(s2,3) + 2*Power(s1,2)*(-3 + s2 + t2) + 
                  Power(s2,2)*(-13 + 6*t2) - s2*(7 + 2*t1 + 10*t2) + 
                  2*(Power(t1,2) - t1*t2 + Power(t2,2)) + 
                  s1*(3 + 8*Power(s2,2) - 8*t2 + s2*(-15 + 8*t2))))/
             ((-1 + s1)*(-1 + t2))))/
        (Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
          Power(s,2)*(-1 + s2)*(s2*t1 - Power(t1,2) - t2 + t1*t2) - 
          s*(-1 + s2)*(Power(s2,2)*(-1 + t1) + 2*t1*(1 - t1 + t2) + 
             s1*(1 + s2 - 2*t1 + s2*t1 - Power(t1,2) - t2 + t1*t2) - 
             s2*(1 + Power(t1,2) + 3*t2 - t1*(2 + t2)))) + 
       4*((-2 + 2*t2 - 2*s*t2 + Power(t1,2)*(5 + s - s1 + 4*t2) + 
             Power(s2,3)*(s1 + (3 - 4*t1)*t2) - 
             t1*(1 + (11 - 5*s + 2*s1)*t2 + 2*Power(t2,2)) + 
             Power(s2,2)*(-1 + 4*Power(t1,2) - 5*s*t2 + 
                2*Power(t2,2) + t1*(-3 + s + 3*t2 + 4*s*t2) + 
                s1*(3 - 4*(-1 + t1)*t2)) - 
             s2*(3 + 3*Power(t1,2) - 3*t2 + 6*t1*t2 - 12*Power(t2,2) + 
                4*t1*Power(t2,2) - 2*Power(t2,3) + 
                s1*(-2 + 6*t1 - 3*Power(t1,2) + 2*t2 - 
                   2*Power(t2,2)) + 
                s*(-7*t1 + 3*Power(t1,2) + t2 + t1*t2 + 6*Power(t2,2)))\
)/((-1 + s2)*(-1 + t1)*(-t1 + s2*t2)) + 
          (4*Power(s2,4)*t2 - 2*Power(t1,2)*(2 + 2*s - 2*s1 + t2) + 
             (1 + t2)*(1 + (-1 + s)*t2) - 
             t1*(3 + s + t2 + 3*s*t2 + 2*Power(t2,2) - 
                s1*(5 + 3*t2)) - 
             2*Power(s2,3)*(s1 - 2*s1*t2 + 2*(t1 + (2 + s - t2)*t2)) - 
             Power(s2,2)*(-2 - 2*t1*(4 + s - 2*t2) + 9*t2 - 4*s*t2 + 
                9*Power(t2,2) + 4*s*Power(t2,2) + 
                s1*(3 + 2*t1 + 3*t2 - 4*Power(t2,2))) + 
             s2*(3 + t2 + 3*s*t2 + 3*s*Power(t2,2) + 2*Power(t2,3) + 
                t1*(7 - 5*s + 13*t2 + 7*s*t2 + 2*Power(t2,2)) - 
                s1*(1 + 5*t2 + 2*Power(t2,2) + t1*(-5 + 7*t2))))/
           ((-1 + t1)*(-1 + t2)*(t1 - s2*t2)) + 
          (Power(t1,4)*(-17 + 6*Power(s,2) + s*(7 - 12*s1) - 7*s1 + 
                6*Power(s1,2) - 4*t2) + 
             2*(-1 + t2)*Power(1 + (-1 + s)*t2,2) + 
             Power(s2,5)*(2*Power(s1,2) + 3*(3 - 2*t1)*Power(t2,2) - 
                s1*(t2 + 2*t1*t2)) + 
             2*t1*(1 + 2*(-1 - s + Power(s,2))*t2 + 
                (1 + 2*s - Power(s,2))*Power(t2,2) + 
                s1*(-1 + t2)*(3 - 3*t2 + s*(-2 + 3*t2))) + 
             Power(t1,3)*(Power(s,2)*(4 - 8*t2) + 5*s1*(-5 + t2) - 
                6*Power(s1,2)*(-1 + t2) + 
                2*(-8 + 6*t2 + 3*Power(t2,2)) + 
                s*(7 - 8*t2 + 2*s1*(-5 + 7*t2))) + 
             Power(t1,2)*(2*Power(s,2)*(-3 + Power(t2,2)) + 
                (-1 + t2)*(11 + 4*Power(s1,2) + t2 - 2*Power(t2,2) + 
                   2*s1*(4 + t2)) - 
                s*(-16 + t2 - Power(t2,2) + 
                   2*s1*(1 + 3*t2 + Power(t2,2)))) + 
             Power(s2,3)*(2 - 4*t2 + 2*Power(t2,2) - 
                6*s*Power(t2,2) + 6*Power(s,2)*Power(t2,2) + 
                2*Power(t2,3) - 15*s*Power(t2,3) - 4*Power(t2,4) + 
                2*Power(t1,3)*(-3 + s - 6*t2 + s*t2) + 
                t1*(-1 + (13 + 37*s - 4*Power(s,2))*t2 + 
                   (1 + 3*s - 6*Power(s,2))*Power(t2,2) + 
                   (11 + 12*s)*Power(t2,3)) + 
                2*Power(s1,2)*
                 (-1 + Power(t2,2) + 2*Power(t1,2)*(2 + t2) - 
                   t1*(5 + 4*t2 + 3*Power(t2,2))) + 
                Power(t1,2)*
                 (7 - 8*t2 + 7*Power(t2,2) + 2*Power(s,2)*(1 + t2) - 
                   s*(3 + 22*t2 + 14*Power(t2,2))) - 
                s1*(4 + t2 + (-19 + 6*s)*Power(t2,2) - 
                   6*Power(t2,3) + 2*Power(t1,3)*(2 + t2) + 
                   Power(t1,2)*(1 + 2*s - 4*t2)*(5 + 3*t2) + 
                   2*t1*(-4 + (7 - 5*s)*t2 - 
                      2*(2 + 3*s)*Power(t2,2) + 5*Power(t2,3)))) + 
             Power(s2,4)*(Power(s1,2)*(2 + 4*t2 - 2*t1*(3 + t2)) + 
                s1*(-4 + 4*t2 + 5*Power(t2,2) + 
                   Power(t1,2)*(2 + 4*t2) + 
                   t1*(1 + 5*t2 - 12*Power(t2,2) + 2*s*(2 + t2))) + 
                t2*(1 - 9*t2 + 7*Power(t2,2) + 
                   6*Power(t1,2)*(2 + t2) - 
                   2*t1*(8 - t2 + 3*Power(t2,2)) + 
                   s*(-4 - 2*Power(t1,2) - 15*t2 + 3*t1*(1 + 4*t2)))) \
+ s2*((-7 - 9*s + 2*Power(s,2))*Power(t1,4) - 
                2*(1 - 2*t2 + Power(-1 + s,2)*Power(t2,2) + 
                   2*s*Power(t2,3)) + 
                t1*(-3 + (25 - 49*s + 12*Power(s,2))*t2 + 
                   (-17 + 18*s - 10*Power(s,2))*Power(t2,2) - 
                   (9 + s + 2*Power(s,2))*Power(t2,3) + 4*Power(t2,4)\
) + Power(t1,3)*(s - s*t2 - 2*Power(s,2)*(4 + 7*t2) + 
                   8*(1 + 6*t2 + Power(t2,2))) + 
                2*Power(s1,2)*t1*
                 (3 + Power(t1,3) - 2*t2 - Power(t2,2) - 
                   6*Power(t1,2)*(1 + t2) + 
                   t1*(-3 - 2*t2 + 5*Power(t2,2))) + 
                Power(t1,2)*
                 (2 + 30*t2 - 26*Power(t2,2) - 12*Power(t2,3) + 
                   2*Power(s,2)*(1 + 3*t2 + 7*Power(t2,2)) + 
                   s*(9 - 39*t2 + 11*Power(t2,2))) + 
                s1*((9 - 4*s)*Power(t1,4) + 
                   2*(2 - 3*t2 + s*Power(t2,2) - 
                      (-1 + s)*Power(t2,3)) + 
                   2*t1*(-4 + (17 + 3*s)*t2 + 
                      6*(-2 + s)*Power(t2,2) + (-1 + s)*Power(t2,3)) \
+ 2*Power(t1,3)*(8 + s*(10 + 13*t2)) + 
                   Power(t1,2)*
                    (8 + 59*t2 - 7*Power(t2,2) - 
                      2*s*(-7 + 2*t2 + 12*Power(t2,2))))) - 
             Power(s2,2)*(-2 + 2*(-3 + s)*Power(t1,4) + 3*t2 - 
                8*s*t2 + 8*Power(t2,2) - 17*s*Power(t2,2) + 
                2*Power(s,2)*Power(t2,2) - 5*Power(t2,3) + 
                9*s*Power(t2,3) - 6*Power(s,2)*Power(t2,3) - 
                6*Power(t2,4) + 2*Power(t2,5) + 
                t1*(2 + (-2 + 3*s + 12*Power(s,2))*t2 + 
                   (18 - 47*s + 6*Power(s,2))*Power(t2,2) + 
                   3*(-6 + s + 2*Power(s,2))*Power(t2,3) - 
                   6*Power(t2,4)) + 
                Power(t1,3)*(-6*(1 + t2) + 2*Power(s,2)*(2 + t2) - 
                   s*(10 + 23*t2)) - 
                2*Power(t1,2)*
                 (-2 - 8*t2 - 21*Power(t2,2) - 2*Power(t2,3) + 
                   Power(s,2)*(1 + 8*t2 + 4*Power(t2,2)) - 
                   s*(7 + 4*t2 + 9*Power(t2,2))) + 
                2*Power(s1,2)*
                 (1 - Power(t2,2) + Power(t1,3)*(3 + t2) - 
                   Power(t1,2)*(7 + 8*t2 + 3*Power(t2,2)) + 
                   t1*(-3 + t2 + 2*Power(t2,3))) + 
                s1*(-2*Power(t1,4) + 
                   2*(-2 + 2*t2 + (7 + 2*s)*Power(t2,2) + 
                      (-7 + 3*s)*Power(t2,3)) + 
                   Power(t1,3)*(5 + 21*t2 - 2*s*(5 + 2*t2)) + 
                   Power(t1,2)*
                    (3 + 20*t2 - 17*Power(t2,2) + 
                      2*s*(5 + 16*t2 + 7*Power(t2,2))) - 
                   t1*(13 - 39*t2 - 36*Power(t2,2) + 2*Power(t2,3) + 
                      2*s*(-4 + 4*t2 + 3*Power(t2,2) + 5*Power(t2,3))))\
))/((-1 + s2)*(-s + s2 - t1)*(-1 + s2 - t1 + t2)*Power(t1 - s2*t2,2)) + 
          (2*Power(t1,5) + Power(t1,4)*(2 - 4*t2) - 
             2*Power(s2,6)*t2*(s1 + 3*t2) + 
             2*(-1 + t2)*Power(1 + (-1 + s)*t2,2) + 
             Power(t1,3)*(3 + s1*(3 - 12*t2) - 4*t2 + 6*Power(t2,2) + 
                2*Power(s1,2)*(-5 + 3*t2) + Power(s,2)*(-4 + 6*t2) + 
                s*(-5 + 12*t2 - 2*s1*(-7 + 6*t2))) + 
             t1*((5 + 2*s1*(-4 + t2) - 2*t2)*Power(-1 + t2,2) + 
                2*Power(s,2)*t2*(3 - 3*t2 + Power(t2,2)) - 
                s*(-1 + t2)*(2 - 5*t2 + 2*s1*(2 - 4*t2 + Power(t2,2)))\
) + Power(t1,2)*(-4*Power(s,2)*(2 - 3*t2 + 2*Power(t2,2)) - 
                (-1 + t2)*(-6 + s1*(11 - 10*t2) - 4*t2 + 
                   4*Power(t2,2) + 2*Power(s1,2)*(-7 + 3*t2)) + 
                s*(13 + 16*t2 - 12*Power(t2,2) + 
                   2*s1*(5 - 17*t2 + 7*Power(t2,2)))) + 
             Power(s2,4)*(2*Power(s1,2)*(-1 + 2*t1 - 4*t2)*
                 (-1 + t2) + 2*Power(t1,2)*(-3 + s - 6*t2 + s*t2) + 
                t1*(-2 + (-58 - 13*s + 2*Power(s,2))*t2 + 
                   (9 - 16*s)*Power(t2,2) + 6*Power(t2,3)) - 
                t2*(5 + 12*t2 + 6*Power(s,2)*t2 - 39*Power(t2,2) + 
                   6*Power(t2,3) + s*(6 + 33*t2 - 24*Power(t2,2))) + 
                s1*(-4 + (5 - 2*s)*t2 + (39 + 14*s)*Power(t2,2) - 
                   22*Power(t2,3) - 2*Power(t1,2)*(2 + t2) + 
                   t1*(-5 + s*(2 - 6*t2) + 16*t2 + 16*Power(t2,2)))) + 
             Power(s2,5)*(-2*Power(s1,2)*(-1 + t2) + 
                s1*((5 + 2*s - 14*t2)*t2 + t1*(2 + 4*t2)) + 
                t2*(2 + 3*(9 + 4*s)*t2 - 12*Power(t2,2) - 
                   2*t1*(s - 3*(2 + t2)))) + 
             Power(s2,3)*(2 - 2*(-3 + s)*Power(t1,3) - 3*t2 + 
                2*s*t2 - 5*Power(t2,2) - 3*s*Power(t2,2) + 
                12*Power(s,2)*Power(t2,2) - 7*Power(t2,3) - 
                41*s*Power(t2,3) - 12*Power(s,2)*Power(t2,3) + 
                8*Power(t2,4) + 12*s*Power(t2,4) + 
                Power(t1,2)*
                 (31 + s + 20*t2 + 19*s*t2 - 2*Power(s,2)*t2 - 
                   20*Power(t2,2) + 2*s*Power(t2,2)) - 
                2*Power(s1,2)*
                 (1 + Power(t1,2)*(-1 + t2) - t2 - 7*Power(t2,2) + 
                   5*Power(t2,3) - t1*(-6 + t2 + 5*Power(t2,2))) + 
                t1*(5 + 36*t2 - 95*Power(t2,2) + 8*Power(t2,3) + 
                   2*Power(s,2)*t2*(2 + 5*t2) - 
                   2*s*(-1 - 36*t2 + 7*Power(t2,2) + 7*Power(t2,3))) \
+ s1*(-4 + 2*Power(t1,3) - (9 + 2*s)*t2 - 2*(5 + 12*s)*Power(t2,2) + 
                   (42 + 22*s)*Power(t2,3) - 10*Power(t2,4) - 
                   Power(t1,2)*
                    (2 + s*(2 - 4*t2) + 23*t2 + 2*Power(t2,2)) + 
                   t1*(7 - 86*t2 + 13*Power(t2,2) + 12*Power(t2,3) + 
                      s*(2 - 8*t2 - 20*Power(t2,2))))) + 
             Power(s2,2)*(2 + 3*t2 + 10*s*t2 - 20*Power(t2,2) + 
                17*s*Power(t2,2) + 13*Power(t2,3) - 2*s*Power(t2,3) + 
                10*Power(s,2)*Power(t2,3) + 6*Power(t2,4) - 
                8*s*Power(t2,4) - 6*Power(s,2)*Power(t2,4) - 
                4*Power(t2,5) + 
                Power(t1,3)*
                 (-17 + 22*t2 + 2*Power(t2,2) - s*(3 + 2*t2)) - 
                Power(t1,2)*
                 (24 - 64*t2 - 6*Power(t2,2) + 4*Power(t2,3) + 
                   2*Power(s,2)*t2*(4 + t2) + 
                   s*(31 + 40*t2 - 24*Power(t2,2))) - 
                2*Power(s1,2)*
                 (1 + t2 + Power(t2,2) - 5*Power(t2,3) + 
                   2*Power(t2,4) + 
                   Power(t1,2)*(-6 + 5*t2 + Power(t2,2)) - 
                   t1*(2 - 16*t2 + 5*Power(t2,2) + 3*Power(t2,3))) + 
                t1*(-3 + 14*t2 + 20*Power(t2,2) - 22*Power(t2,3) + 
                   6*Power(t2,4) + 
                   2*Power(s,2)*t2*(-12 + 5*t2 + 4*Power(t2,2)) - 
                   s*(2 + 20*t2 - 103*Power(t2,2) + 14*Power(t2,3))) \
+ s1*(4 + (-3 + 2*s)*t2 + (3 - 2*s)*Power(t2,2) - 
                   4*(3 + 5*s)*Power(t2,3) + 
                   2*(4 + 5*s)*Power(t2,4) + Power(t1,3)*(7 + 2*t2) + 
                   Power(t1,2)*
                    (35 + 41*t2 - 22*Power(t2,2) + 
                      2*s*(-3 + 9*t2 + 2*Power(t2,2))) + 
                   t1*(3*(7 + 14*t2 - 34*Power(t2,2) + 
                       4*Power(t2,3)) - 
                      2*s*(3 - 25*t2 + 10*Power(t2,2) + 7*Power(t2,3))\
))) + s2*(-2 + t2 - 2*s*t2 + 6*Power(t2,2) + 7*s*Power(t2,2) - 
                4*Power(s,2)*Power(t2,2) - 7*Power(t2,3) - 
                5*s*Power(t2,3) + 4*Power(s,2)*Power(t2,3) + 
                2*Power(t2,4) - 2*Power(s,2)*Power(t2,4) - 
                4*Power(t1,4)*(2 + t2) + 
                2*Power(t1,3)*
                 (-4 + s*(13 - 5*t2) - 5*t2 + Power(s,2)*t2 + 
                   4*Power(t2,2)) + 
                2*Power(s1,2)*t1*
                 (4 + Power(t1,2)*(-1 + t2) + 5*t2 - 14*Power(t2,2) + 
                   5*Power(t2,3) + t1*(6 + 6*t2 - 6*Power(t2,2))) - 
                Power(t1,2)*(2*Power(s,2)*
                    (-4 - 5*t2 + 7*Power(t2,2)) + 
                   s*(-23 + 57*t2 + 10*Power(t2,2)) + 
                   3*(1 + 6*t2 - 6*Power(t2,2) + 4*Power(t2,3))) + 
                t1*(-9 + 32*t2 - 17*Power(t2,2) - 14*Power(t2,3) + 
                   8*Power(t2,4) + 
                   2*Power(s,2)*t2*(6 - 13*t2 + 7*Power(t2,2)) + 
                   s*(-2 - 46*t2 - 6*Power(t2,2) + 20*Power(t2,3))) - 
                s1*(-2*(-1 + t2)*
                    (-2 + (3 - 2*s)*Power(t2,2) + 
                      (-1 + s)*Power(t2,3)) + 
                   2*Power(t1,3)*(14 - 5*t2 + s*(-1 + 2*t2)) + 
                   Power(t1,2)*
                    (44 - 61*t2 - 10*Power(t2,2) + 
                      s*(10 + 24*t2 - 26*Power(t2,2))) + 
                   t1*(9 + 2*t2 - 29*Power(t2,2) + 18*Power(t2,3) + 
                      s*(2 + 8*t2 - 54*Power(t2,2) + 24*Power(t2,3)))))\
)/((-1 + s1)*(-1 + t2)*(-1 + s2 - t1 + t2)*Power(t1 - s2*t2,2)) + 
          (2*(-41 - 2*Power(s,3) + 68*s2 + 26*Power(s2,2) - 
               14*Power(s2,3) + 13*t1 - 51*s2*t1 + 14*Power(s2,2)*t1 + 
               65*Power(t1,2) - 40*s2*Power(t1,2) + 95*t2 - 
               180*s2*t2 + 128*Power(s2,2)*t2 + 3*t1*t2 - 
               80*s2*t1*t2 - 86*Power(t2,2) + 106*s2*Power(t2,2) - 
               16*(-1 + s2)*(s1 - t2)*(-1 + s2 - t1 + t2) + 
               2*Power(s1,2)*(-3 + 7*s2 - 7*t1 + 5*t2) - 
               2*Power(s,2)*(-12 - 2*s1 + 2*s2 + t1 + 6*t2) - 
               2*s*(17 + Power(s1,2) + 3*Power(s2,2) + 21*t1 - 
                  20*Power(t1,2) + s1*(9 + 5*s2 - 8*t1 - t2) - 
                  24*t2 + 8*t1*t2 - 2*Power(t2,2) + 
                  s2*(-15 + 6*t1 + 16*t2)) + 
               ((-2 + 2*Power(s,2) + Power(s1,2) - 9*s2 + 
                    3*Power(s2,2) - t1 - 3*t2 + s2*t2 + 
                    s1*(-3 + 4*s2 + t2) - s*(-6 + 3*s1 + 5*s2 + t2))*
                  (-1 - (-1 + s + s1)*t1 + Power(s2,2)*(s1 - t2) + 
                    t2 - s*t2 + 
                    s2*(-1 + s1 + t1 + s*t1 - s1*t1 + s*t2)))/
                (t1 - s2*t2) - 
               s1*(48 + 120*Power(s2,2) + 40*Power(t1,2) - 70*t2 + 
                  14*Power(t2,2) - 5*t1*(-3 + 8*t2) + 
                  s2*(-120*t1 + 17*(-7 + 6*t2))) + 
               8*(2 - 2*s2 - 2*Power(s2,2) + Power(s2,3) + 2*t1 + 
                  2*s2*t1 - Power(s2,2)*t1 - 4*Power(t1,2) + 
                  2*s2*Power(t1,2) - 15*t2 + 27*s2*t2 - 
                  15*Power(s2,2)*t2 - 9*t1*t2 + 13*s2*t1*t2 + 
                  13*Power(t2,2) - 14*s2*Power(t2,2) - 
                  Power(s1,2)*(-1 + s2 - t1 + t2) + 
                  s*(1 + 3*t1 - 2*Power(t1,2) - t2 + t1*t2 - 
                     Power(t2,2) + s1*(-1 + s2 - t1 + t2) + 
                     s2*(-3 + t1 + t2)) + 
                  s1*(11 + 15*Power(s2,2) + 2*Power(t1,2) + 
                     t1*(9 - 2*t2) - 12*t2 + Power(t2,2) + 
                     s2*(-23 - 15*t1 + 14*t2))) + 
               2*(46 + 2*Power(s,3) - 62*s2 - 28*Power(s2,2) + 
                  18*Power(s2,3) + 10*t1 + 56*s2*t1 - 
                  19*Power(s2,2)*t1 - 80*Power(t1,2) + 
                  48*s2*Power(t1,2) + 
                  Power(s1,2)*(16 - 21*s2 + 19*t1 - 16*t2) - 
                  Power(s,2)*(14 + 4*s1 + s2 - 3*t1 - 8*t2) - 
                  140*t2 + 254*s2*t2 - 166*Power(s2,2)*t2 - 
                  24*t1*t2 + 112*s2*t1*t2 + 122*Power(t2,2) - 
                  141*s2*Power(t2,2) + 
                  s1*(76 + 157*Power(s2,2) + 48*Power(t1,2) - 
                     2*s2*(83 + 80*t1 - 69*t2) + t1*(34 - 48*t2) - 
                     102*t2 + 19*Power(t2,2)) + 
                  2*s*(13 + Power(s1,2) + 5*Power(s2,2) + 27*t1 - 
                     24*Power(t1,2) - 23*t2 + 11*t1*t2 - 
                     6*Power(t2,2) + s1*(-1 + 11*s2 - 11*t1 + 4*t2) + 
                     s2*(-28 + 8*t1 + 17*t2))) - 
               4*(19 - 20*s2 - 17*Power(s2,2) + 9*Power(s2,3) + 
                  11*t1 + 19*s2*t1 - 8*Power(s2,2)*t1 - 
                  32*Power(t1,2) + 18*s2*Power(t1,2) + 
                  Power(s1,2)*(7 - 8*s2 + 8*t1 - 7*t2) - 74*t2 + 
                  125*s2*t2 - 76*Power(s2,2)*t2 - 25*t1*t2 + 
                  58*s2*t1*t2 + 61*Power(t2,2) - 68*s2*Power(t2,2) + 
                  Power(s,2)*(-4 + s2 + t1 + 2*t2) + 
                  s*(7 + 23*t1 - 18*Power(t1,2) - 11*t2 + 9*t1*t2 - 
                     7*Power(t2,2) + s1*(-3 + 7*s2 - 9*t1 + 5*t2) + 
                     s2*(-19 + 7*t1 + 9*t2)) + 
                  s1*(45 + 77*Power(s2,2) + 18*Power(t1,2) - 54*t2 + 
                     8*Power(t2,2) - 9*t1*(-3 + 2*t2) + 
                     s2*(-96 - 76*t1 + 69*t2))) + 
               ((3 + s - s1 - s2)*
                  ((Power(s,2) + Power(-1 + s1,2) - 2*s*(1 + 2*s1))*
                     Power(t1,3) + Power(s2,5)*Power(s1 - t2,2) + 
                    (-1 + t2)*Power(1 + (-1 + s)*t2,2) + 
                    t1*(2*s*(-1 + s1)*Power(-1 + t2,2) - 
                       (-3 + 2*s1)*Power(-1 + t2,2) + 
                       Power(s,2)*Power(t2,2)) - 
                    Power(t1,2)*
                     (-((3 - 4*s1 + Power(s1,2))*(-1 + t2)) + 
                       Power(s,2)*(1 + 3*t2) - 
                       2*s*(2 + s1*(-3 + t2) + 3*t2)) + 
                    Power(s2,4)*
                     (Power(s1,2)*(1 - 3*t1 + t2) + 
                       t2*(2 - t2 + Power(t2,2) - 2*s*(1 + t2) - 
                       t1*(2 + t2)) + 
                       2*s1*
                        (-1 + s*t2 - Power(t2,2) + t1*(1 + s + 2*t2))\
) + Power(s2,3)*(1 - 2*s*t2 - 2*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) - 
                       2*s*Power(t2,3) + 
                       Power(t1,2)*(1 + Power(s,2) + 2*t2) + 
                       2*t1*
                       (-1 + (-1 + 4*s)*t2 + (-1 + s)*Power(t2,2)) + 
                       Power(s1,2)*
                       (-1 + 3*Power(t1,2) - 2*t1*(2 + t2)) + 
                       2*s1*
                        (-1 + t2 + (1 + s)*Power(t2,2) - 
                        Power(t1,2)*(2 + 2*s + t2) + 
                        t1*(3 - s + 2*t2 + Power(t2,2)))) + 
                    s2*(-1 + 2*(-1 + s)*Power(t1,3) + 2*s*t2 + 
                       3*Power(t2,2) - 4*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) + 
                       2*s*Power(t2,3) - 2*Power(s,2)*Power(t2,3) - 
                       Power(s1,2)*t1*
                       (-2 + 3*t1 + 2*Power(t1,2) + 2*t2) - 
                       2*t1*t2*
                       (-3 + 8*s - 2*Power(s,2) + 3*t2 + 2*s*t2) + 
                       Power(t1,2)*
                       (3 - 2*s + Power(s,2) - 6*t2 + 
                       2*Power(s,2)*t2) + 
                       2*s1*
                        ((2 + s)*Power(t1,3) + Power(-1 + t2,2) + 
                        Power(t1,2)*(5*s + 3*t2) + 
                        t1*(-3 + s + 2*t2 + 4*s*t2 + Power(t2,2) - 
                        s*Power(t2,2)))) - 
                    Power(s2,2)*
                     (-1 + (1 + Power(s,2))*Power(t1,3) + 3*t2 - 
                       4*s*t2 - 4*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) - 
                       2*s*Power(t2,3) - Power(s,2)*Power(t2,3) + 
                       Power(s1,2)*
                        (1 - 3*t1 + Power(t1,3) - t2 - 
                        Power(t1,2)*(5 + t2)) + 
                       Power(t1,2)*
                        (-(Power(s,2)*(-1 + t2)) - 3*(1 + t2) + 
                        s*(2 + 6*t2)) + 
                       t1*(3 - 6*Power(t2,2) + 
                        Power(s,2)*t2*(4 + t2) - 
                        2*s*(1 + 2*t2 + 2*Power(t2,2))) + 
                       2*s1*(-((1 + s)*Power(t1,3)) + 
                        (1 + t2)*(-1 + t2 + s*t2) + 
                        Power(t1,2)*(4 + (2 + s)*t2) + 
                        t1*(-2 + 4*t2 + Power(t2,2) + 2*s*(1 + t2)))))\
)/((-1 + s2 - t1 + t2)*Power(t1 - s2*t2,2))))/((-1 + s1)*(-s + s1 - t2)) \
+ (2*(6 - 2*Power(s,3) - 77*s2 + 14*Power(s2,2) + 14*Power(s2,3) + 
               2*Power(s1,2)*(3 + 2*s2 - 2*t1) + 36*t1 + 86*s2*t1 - 
               38*Power(s2,2)*t1 - 59*Power(t1,2) + 24*s2*Power(t1,2) - 
               27*t2 + 135*s2*t2 - 76*Power(s2,2)*t2 - 25*t1*t2 + 
               28*s2*t1*t2 + 31*Power(t2,2) - 46*s2*Power(t2,2) + 
               16*s*(-1 + s2)*(-1 + s2 - t1 + t2) + 
               2*Power(s,2)*(-2 + 2*s1 + 11*s2 - 14*t1 + 7*t2) + 
               s1*(-6 + 68*Power(s2,2) + 24*Power(t1,2) + 
                  t1*(45 - 22*t2) - 4*s2*(7 + 23*t1 - 6*t2) - 9*t2 + 
                  4*Power(t2,2)) - 
               2*s*(-54 + Power(s1,2) - 29*Power(s2,2) + 9*t1 + 
                  12*Power(t1,2) + s2*(76 - 4*t1 - 55*t2) + 60*t2 - 
                  5*t1*t2 - Power(t2,2) + s1*(1 + 13*s2 - 16*t1 + 7*t2)\
) + ((-3 + 2*Power(s,2) + Power(s1,2) - 8*s2 + 3*Power(s2,2) - 
                    2*t1 - 2*t2 + s2*t2 + s1*(-3 + 4*s2 + t2) - 
                    s*(-6 + 3*s1 + 5*s2 + t2))*
                  (-1 - (-1 + s + s1)*t1 + Power(s2,2)*(s1 - t2) + 
                    t2 - s*t2 + 
                    s2*(-1 + s1 + t1 + s*t1 - s1*t1 + s*t2)))/
                (t1 - s2*t2) + 
               2*(-26 + 2*Power(s,3) + 101*s2 + 8*Power(s2,2) - 
                  24*Power(s2,3) - 26*t1 - 107*s2*t1 + 
                  49*Power(s2,2)*t1 + 64*Power(t1,2) - 
                  26*s2*Power(t1,2) + 
                  Power(s,2)*(19 - 4*s1 - 32*s2 + 34*t1 - 23*t2) + 
                  53*t2 - 161*s2*t2 + 80*Power(s2,2)*t2 + 35*t1*t2 - 
                  30*s2*t1*t2 - 37*Power(t2,2) + 51*s2*Power(t2,2) - 
                  Power(s1,2)*(1 + 6*s2 - 4*t1 + t2) + 
                  s*(-150 + 2*Power(s1,2) - 78*Power(s2,2) + 10*t1 + 
                     26*Power(t1,2) + s2*(187 + 6*t1 - 139*t2) + 
                     154*t2 - 15*t1*t2 - Power(t2,2) + 
                     2*s1*(-9 + 19*s2 - 19*t1 + 12*t2)) + 
                  s1*(8 - 85*Power(s2,2) - 26*Power(t1,2) + 
                     s2*(57 + 108*t1 - 31*t2) + 6*t2 - 4*Power(t2,2) + 
                     t1*(-58 + 25*t2))) + 
               4*(12 - 36*s2 - 6*Power(s2,2) + 9*Power(s2,3) + 
                  Power(s1,2)*(s2 - t1) + 7*t1 + 44*s2*t1 - 
                  19*Power(s2,2)*t1 - 21*Power(t1,2) + 
                  9*s2*Power(t1,2) - 26*t2 + 63*s2*t2 - 
                  28*Power(s2,2)*t2 - 15*t1*t2 + 10*s2*t1*t2 + 
                  15*Power(t2,2) - 18*s2*Power(t2,2) + 
                  Power(s,2)*(-7 + 10*s2 - 12*t1 + 9*t2) + 
                  s1*(-1 + 29*Power(s2,2) + 9*Power(t1,2) + 
                     t1*(24 - 9*t2) - t2 + Power(t2,2) + 
                     s2*(-23 - 39*t1 + 10*t2)) + 
                  s*(71 + 48*Power(s2,2) + 8*t1 - 9*Power(t1,2) + 
                     s1*(7 - 11*s2 + 13*t1 - 9*t2) - 71*t2 + 7*t1*t2 + 
                     s2*(-103 - 19*t1 + 68*t2))) - 
               8*(2 + Power(s2,3) + 3*s1*t1 - 2*Power(t1,2) + 
                  s1*Power(t1,2) + 
                  Power(s2,2)*(-1 + 3*s1 - 2*t1 - 3*t2) - 4*t2 - 
                  2*t1*t2 - s1*t1*t2 + 2*Power(t2,2) + 
                  Power(s,2)*(-1 + s2 - t1 + t2) + 
                  s*(14 + 12*Power(s2,2) + 7*t1 - Power(t1,2) - 
                     14*t2 + t1*t2 - s1*(-1 + s2 - t1 + t2) + 
                     s2*(-24 - 9*t1 + 14*t2)) + 
                  s2*(Power(t1,2) + t1*(5 + t2) + 
                     s1*(-3 - 4*t1 + t2) - 2*(2 - 4*t2 + Power(t2,2)))) \
+ ((3 + s - s1 - s2)*((Power(s,2) + Power(-1 + s1,2) - 
                       2*s*(1 + 2*s1))*Power(t1,3) + 
                    Power(s2,5)*Power(s1 - t2,2) + 
                    (-1 + t2)*Power(1 + (-1 + s)*t2,2) + 
                    t1*(2*s*(-1 + s1)*Power(-1 + t2,2) - 
                       (-3 + 2*s1)*Power(-1 + t2,2) + 
                       Power(s,2)*Power(t2,2)) - 
                    Power(t1,2)*
                     (-((3 - 4*s1 + Power(s1,2))*(-1 + t2)) + 
                       Power(s,2)*(1 + 3*t2) - 
                       2*s*(2 + s1*(-3 + t2) + 3*t2)) + 
                    Power(s2,4)*
                     (Power(s1,2)*(1 - 3*t1 + t2) + 
                       t2*(2 - t2 + Power(t2,2) - 2*s*(1 + t2) - 
                        t1*(2 + t2)) + 
                       2*s1*(-1 + s*t2 - Power(t2,2) + 
                        t1*(1 + s + 2*t2))) + 
                    Power(s2,3)*
                     (1 - 2*s*t2 - 2*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) - 
                       2*s*Power(t2,3) + 
                       Power(t1,2)*(1 + Power(s,2) + 2*t2) + 
                       2*t1*
                        (-1 + (-1 + 4*s)*t2 + (-1 + s)*Power(t2,2)) + 
                       Power(s1,2)*
                        (-1 + 3*Power(t1,2) - 2*t1*(2 + t2)) + 
                       2*s1*(-1 + t2 + (1 + s)*Power(t2,2) - 
                        Power(t1,2)*(2 + 2*s + t2) + 
                        t1*(3 - s + 2*t2 + Power(t2,2)))) + 
                    s2*(-1 + 2*(-1 + s)*Power(t1,3) + 2*s*t2 + 
                       3*Power(t2,2) - 4*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) + 
                       2*s*Power(t2,3) - 2*Power(s,2)*Power(t2,3) - 
                       Power(s1,2)*t1*
                        (-2 + 3*t1 + 2*Power(t1,2) + 2*t2) - 
                       2*t1*t2*
                        (-3 + 8*s - 2*Power(s,2) + 3*t2 + 2*s*t2) + 
                       Power(t1,2)*
                        (3 - 2*s + Power(s,2) - 6*t2 + 
                        2*Power(s,2)*t2) + 
                       2*s1*((2 + s)*Power(t1,3) + Power(-1 + t2,2) + 
                        Power(t1,2)*(5*s + 3*t2) + 
                        t1*(-3 + s + 2*t2 + 4*s*t2 + Power(t2,2) - 
                        s*Power(t2,2)))) - 
                    Power(s2,2)*
                     (-1 + (1 + Power(s,2))*Power(t1,3) + 3*t2 - 
                       4*s*t2 - 4*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) - 
                       2*s*Power(t2,3) - Power(s,2)*Power(t2,3) + 
                       Power(s1,2)*
                        (1 - 3*t1 + Power(t1,3) - t2 - 
                        Power(t1,2)*(5 + t2)) + 
                       Power(t1,2)*
                        (-(Power(s,2)*(-1 + t2)) - 3*(1 + t2) + 
                        s*(2 + 6*t2)) + 
                       t1*(3 - 6*Power(t2,2) + 
                        Power(s,2)*t2*(4 + t2) - 
                        2*s*(1 + 2*t2 + 2*Power(t2,2))) + 
                       2*s1*(-((1 + s)*Power(t1,3)) + 
                        (1 + t2)*(-1 + t2 + s*t2) + 
                        Power(t1,2)*(4 + (2 + s)*t2) + 
                        t1*(-2 + 4*t2 + Power(t2,2) + 2*s*(1 + t2))))))/
                ((-1 + s2 - t1 + t2)*Power(t1 - s2*t2,2))))/
           ((s - s2 + t1)*(s - s1 + t2))) + 
       ((-4*(-2*Power(t1,3) + Power(s2,5)*(s1 + t2) + 
               (6 + t2)*(1 + (-1 + s)*t2) - 
               Power(s2,4)*(1 + t1 - s*t1 + s1*(1 + t1 - 4*t2) + 
                  4*t2 + s*t2 - 4*Power(t2,2)) - 
               Power(t1,2)*(4 - 3*s1*(-4 + t2) + 7*t2 + 
                  s*(2 + 3*t2)) - 
               Power(s2,3)*(-1 + 14*t2 + 3*s*t2 + 14*Power(t2,2) + 
                  4*s*Power(t2,2) - 3*Power(t2,3) + 
                  s1*(9 + 2*t1*t2 - 7*Power(t2,2)) + 
                  t1*(-5 + s + 7*t2 - 2*s*t2 + 2*Power(t2,2))) - 
               t1*(8 - 10*t2 + Power(t2,2) + 
                  s*(-2 + 9*t2 + Power(t2,2)) - 
                  s1*(8 + 5*t2 + 4*Power(t2,2))) - 
               Power(s2,2)*(-9 - (3 + s)*Power(t1,2) + 16*t2 - s*t2 + 
                  16*Power(t2,2) + 4*s*Power(t2,2) + 12*Power(t2,3) + 
                  3*s*Power(t2,3) - 
                  t1*(10 + 30*t2 - 2*Power(t2,2) + 
                     s*(-5 + 5*t2 + 3*Power(t2,2))) + 
                  s1*(13 + Power(t1,2) + 11*t2 + Power(t2,2) - 
                     4*Power(t2,3) + t1*(-16 + 9*t2 + 3*Power(t2,2)))) \
+ s2*(13 - 13*Power(t1,2) + 2*Power(t1,3) - 8*t2 + 27*t1*t2 + 
                  Power(t1,2)*t2 - 7*Power(t2,2) + 19*t1*Power(t2,2) + 
                  Power(t2,3) + 
                  s1*(-6 - 3*t2 - 4*Power(t2,2) - 4*Power(t2,3) + 
                     Power(t1,2)*(-5 + 3*t2) + 
                     t1*(27 + 10*t2 - 7*Power(t2,2))) + 
                  s*(-3*Power(t1,2)*(-1 + t2) + 
                     t2*(5 + 5*t2 + Power(t2,2)) + 
                     t1*(-5 + 2*t2 + 6*Power(t2,2))))))/
           ((-1 + t2)*(-t1 + s2*t2)) + 
          (4*(((-7 + Power(s2,2)*(-1 + t1) + 5*t1 - 3*t2 - 
                    Power(t2,2) + s2*(-2 + t1)*(1 + t2))*
                  (-1 - (-1 + s + s1)*t1 + Power(s2,2)*(s1 - t2) + 
                    t2 - s*t2 + 
                    s2*(-1 + s1 + t1 + s*t1 - s1*t1 + s*t2)))/
                (t1 - s2*t2) + 
               2*(5 + 15*s2 + 8*Power(s2,2) + 7*Power(s2,3) - 
                  Power(s2,4) - 11*t1 - 20*s2*t1 - 5*Power(s2,2)*t1 - 
                  Power(s2,3)*t1 + 8*Power(t1,2) + 4*s2*Power(t1,2) + 
                  Power(s2,2)*Power(t1,2) - 2*Power(t1,3) + 11*t2 + 
                  17*s2*t2 + 12*Power(s2,2)*t2 - Power(s2,3)*t2 - 
                  18*t1*t2 - 3*s2*t1*t2 - 2*Power(s2,2)*t1*t2 + 
                  4*Power(t1,2)*t2 + 10*Power(t2,2) + 
                  4*s2*Power(t2,2) - 2*t1*Power(t2,2) - 
                  s1*(2 + Power(s2,3) + 5*t2 - t1*(4 + t2) + 
                     Power(s2,2)*(-2 + t1 + t2) + 
                     s2*(5 - Power(t1,2) - 3*t2 + 2*t1*(1 + t2))) + 
                  s*(-5 + Power(s2,3) - 2*Power(t1,2) - 4*t2 - 
                     4*Power(t2,2) + Power(s2,2)*(-5 + t1 + t2) + 
                     t1*(7 + 3*t2) - 
                     s2*(5 + Power(t1,2) + 10*t2 - t1*(5 + 2*t2))))))/
           (-1 + s2) + (2*(-1 + t1)*
             ((2*(-1 - (-1 + s + s1)*t1 + Power(s2,2)*(s1 - t2) + 
                    t2 - s*t2 + 
                    s2*(-1 + s1 + t1 + s*t1 - s1*t1 + s*t2))*
                  (-4 - 11*s2 - 8*Power(s2,2) + 2*Power(s2,3) + 
                    5*t1 + Power(t1,2) - 2*t2 - 8*s2*t2 + 
                    5*Power(s2,2)*t2 - 2*t1*t2 - s2*t1*t2 - 
                    2*Power(t2,2) + 3*s2*Power(t2,2) + 
                    s*(6 - 2*Power(s2,2) + s2*(4 - 5*t2) + 
                       t1*(-3 + t2) + 2*t2 - 3*Power(t2,2)) + 
                    s1*(-2 + Power(s2,2) + t1 - 2*t2 - t1*t2 + 
                       2*Power(t2,2) + s2*(-4 + 3*t2))))/(t1 - s2*t2) \
- 4*(-3 - 4*s2 - 8*Power(s2,2) - 5*Power(s2,3) + Power(s2,4) + 6*t1 + 
                  2*s2*t1 - 3*Power(s2,2)*t1 + Power(t1,2) + 
                  3*s2*Power(t1,2) - t2 - 5*s2*t2 - 
                  8*Power(s2,2)*t2 + 3*Power(s2,3)*t2 + t1*t2 - 
                  7*s2*t1*t2 - Power(s2,2)*t1*t2 + 2*Power(t1,2)*t2 - 
                  2*s2*Power(t2,2) + 2*Power(s2,2)*Power(t2,2) - 
                  4*t1*Power(t2,2) + 2*Power(t2,3) - 
                  Power(s1,2)*(1 + s2 + 2*t2) + 
                  Power(s,2)*
                   (-1 + Power(s2,2) + 3*s2*t2 - t1*t2 + 
                     2*Power(t2,2)) + 
                  s1*(4 + Power(s2,3) + 4*t2 - 4*Power(t2,2) + 
                     t1*(-1 + 2*t2) + Power(s2,2)*(-4 + 3*t2) + 
                     s2*(-(t1*(-2 + t2)) + 2*(-4 + t2)*t2)) + 
                  s*(1 - 2*Power(s2,3) - 3*t1 - Power(t1,2) + 
                     Power(s2,2)*(4 - 6*t2) + t2 + 4*Power(t2,2) + 
                     s2*(5 + (7 + 2*t1)*t2 - 4*Power(t2,2)) + 
                     s1*(2 + s2 - Power(s2,2) - 3*s2*t2 + 
                        (2 + t1)*t2 - 2*Power(t2,2)))) - 
               (2*(-5 + Power(s2,2) + 2*t1 + 2*s2*(-2 + t2) - t2 + 
                    Power(t2,2))*
                  ((Power(s,2) + Power(-1 + s1,2) - 2*s*(1 + 2*s1))*
                     Power(t1,3) + Power(s2,5)*Power(s1 - t2,2) + 
                    (-1 + t2)*Power(1 + (-1 + s)*t2,2) + 
                    t1*(2*s*(-1 + s1)*Power(-1 + t2,2) - 
                       (-3 + 2*s1)*Power(-1 + t2,2) + 
                       Power(s,2)*Power(t2,2)) - 
                    Power(t1,2)*
                     (-((3 - 4*s1 + Power(s1,2))*(-1 + t2)) + 
                       Power(s,2)*(1 + 3*t2) - 
                       2*s*(2 + s1*(-3 + t2) + 3*t2)) + 
                    Power(s2,4)*
                     (Power(s1,2)*(1 - 3*t1 + t2) + 
                       t2*(2 - t2 + Power(t2,2) - 2*s*(1 + t2) - 
                       t1*(2 + t2)) + 
                       2*s1*
                        (-1 + s*t2 - Power(t2,2) + t1*(1 + s + 2*t2))\
) + Power(s2,3)*(1 - 2*s*t2 - 2*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) - 
                       2*s*Power(t2,3) + 
                       Power(t1,2)*(1 + Power(s,2) + 2*t2) + 
                       2*t1*
                       (-1 + (-1 + 4*s)*t2 + (-1 + s)*Power(t2,2)) + 
                       Power(s1,2)*
                       (-1 + 3*Power(t1,2) - 2*t1*(2 + t2)) + 
                       2*s1*
                        (-1 + t2 + (1 + s)*Power(t2,2) - 
                        Power(t1,2)*(2 + 2*s + t2) + 
                        t1*(3 - s + 2*t2 + Power(t2,2)))) + 
                    s2*(-1 + 2*(-1 + s)*Power(t1,3) + 2*s*t2 + 
                       3*Power(t2,2) - 4*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) + 
                       2*s*Power(t2,3) - 2*Power(s,2)*Power(t2,3) - 
                       Power(s1,2)*t1*
                       (-2 + 3*t1 + 2*Power(t1,2) + 2*t2) - 
                       2*t1*t2*
                       (-3 + 8*s - 2*Power(s,2) + 3*t2 + 2*s*t2) + 
                       Power(t1,2)*
                       (3 - 2*s + Power(s,2) - 6*t2 + 
                       2*Power(s,2)*t2) + 
                       2*s1*
                        ((2 + s)*Power(t1,3) + Power(-1 + t2,2) + 
                        Power(t1,2)*(5*s + 3*t2) + 
                        t1*(-3 + s + 2*t2 + 4*s*t2 + Power(t2,2) - 
                        s*Power(t2,2)))) - 
                    Power(s2,2)*
                     (-1 + (1 + Power(s,2))*Power(t1,3) + 3*t2 - 
                       4*s*t2 - 4*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) - 
                       2*s*Power(t2,3) - Power(s,2)*Power(t2,3) + 
                       Power(s1,2)*
                        (1 - 3*t1 + Power(t1,3) - t2 - 
                        Power(t1,2)*(5 + t2)) + 
                       Power(t1,2)*
                        (-(Power(s,2)*(-1 + t2)) - 3*(1 + t2) + 
                        s*(2 + 6*t2)) + 
                       t1*(3 - 6*Power(t2,2) + 
                        Power(s,2)*t2*(4 + t2) - 
                        2*s*(1 + 2*t2 + 2*Power(t2,2))) + 
                       2*s1*(-((1 + s)*Power(t1,3)) + 
                        (1 + t2)*(-1 + t2 + s*t2) + 
                        Power(t1,2)*(4 + (2 + s)*t2) + 
                        t1*(-2 + 4*t2 + Power(t2,2) + 2*s*(1 + t2)))))\
)/((-1 + s2 - t1 + t2)*Power(t1 - s2*t2,2))))/((-1 + s1)*(-1 + t2)) + 
          (2*(-1 + t1)*((2*(-1 - (-1 + s + s1)*t1 + 
                    Power(s2,2)*(s1 - t2) + t2 - s*t2 + 
                    s2*(-1 + s1 + t1 + s*t1 - s1*t1 + s*t2))*
                  (-11 - 3*s2 - 7*Power(s2,2) + Power(s2,3) - 2*t1 + 
                    2*s2*t1 + 2*Power(s2,2)*t1 + Power(t1,2) - 
                    s2*Power(t1,2) - 2*t2 - 7*s2*t2 + 
                    Power(s2,2)*t2 - 3*t1*t2 + 3*s2*t1*t2 + 
                    s1*(-1 + Power(s2,2) + t1 - Power(t1,2) - t2 + 
                       2*t1*t2 + s2*(-1 + t1 + t2)) - 
                    s*(-6 + Power(s2,2) + 2*t1 - Power(t1,2) - 
                       3*t2 + 3*t1*t2 + s2*(-3 + 2*t1 + t2))))/
                (t1 - s2*t2) - 
               4*(-3 - Power(s1,2)*(-1 + s2) - s2 + 3*Power(s2,2) - 
                  4*Power(s2,3) + Power(s2,4) - 3*t1 - 10*s2*t1 - 
                  5*Power(s2,2)*t1 + Power(s2,3)*t1 + 4*Power(t1,2) + 
                  6*s2*Power(t1,2) - Power(s2,2)*Power(t1,2) - 
                  2*Power(t1,3) + t2 - s2*t2 - 2*Power(s2,2)*t2 + 
                  Power(s2,3)*t2 - 8*t1*t2 - 11*s2*t1*t2 + 
                  2*Power(s2,2)*t1*t2 + 4*Power(t1,2)*t2 + 
                  2*s2*Power(t2,2) - 2*t1*Power(t2,2) + 
                  Power(s,2)*
                   (Power(s2,2) - Power(t1,2) - 3*t2 + 
                     2*t1*(1 + t2) + s2*(-3 + t1 + t2)) + 
                  s1*(-2 + Power(s2,3) - 2*Power(t1,2) - 3*t2 + 
                     Power(s2,2)*(-6 + t1 + t2) + t1*(2 + 3*t2) - 
                     s2*(1 + Power(t1,2) + 5*t2 - 2*t1*(3 + t2))) - 
                  s*(-3 + 2*Power(s2,3) - 2*t1 + Power(t1,2) - 2*t2 - 
                     4*t1*t2 + Power(t2,2) + 
                     Power(s2,2)*(-7 + 2*t1 + 2*t2) + 
                     s2*(1 + t1 - 2*Power(t1,2) - 6*t2 + 4*t1*t2) + 
                     s1*(1 + Power(s2,2) - Power(t1,2) - 3*t2 + 
                        2*t1*(1 + t2) + s2*(-4 + t1 + t2)))) - 
               (2*(-6 - t2 + t1*(1 + s2 + t2))*
                  ((Power(s,2) + Power(-1 + s1,2) - 2*s*(1 + 2*s1))*
                     Power(t1,3) + Power(s2,5)*Power(s1 - t2,2) + 
                    (-1 + t2)*Power(1 + (-1 + s)*t2,2) + 
                    t1*(2*s*(-1 + s1)*Power(-1 + t2,2) - 
                       (-3 + 2*s1)*Power(-1 + t2,2) + 
                       Power(s,2)*Power(t2,2)) - 
                    Power(t1,2)*
                     (-((3 - 4*s1 + Power(s1,2))*(-1 + t2)) + 
                       Power(s,2)*(1 + 3*t2) - 
                       2*s*(2 + s1*(-3 + t2) + 3*t2)) + 
                    Power(s2,4)*
                     (Power(s1,2)*(1 - 3*t1 + t2) + 
                       t2*(2 - t2 + Power(t2,2) - 2*s*(1 + t2) - 
                       t1*(2 + t2)) + 
                       2*s1*
                        (-1 + s*t2 - Power(t2,2) + t1*(1 + s + 2*t2))\
) + Power(s2,3)*(1 - 2*s*t2 - 2*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) - 
                       2*s*Power(t2,3) + 
                       Power(t1,2)*(1 + Power(s,2) + 2*t2) + 
                       2*t1*
                       (-1 + (-1 + 4*s)*t2 + (-1 + s)*Power(t2,2)) + 
                       Power(s1,2)*
                       (-1 + 3*Power(t1,2) - 2*t1*(2 + t2)) + 
                       2*s1*
                        (-1 + t2 + (1 + s)*Power(t2,2) - 
                        Power(t1,2)*(2 + 2*s + t2) + 
                        t1*(3 - s + 2*t2 + Power(t2,2)))) + 
                    s2*(-1 + 2*(-1 + s)*Power(t1,3) + 2*s*t2 + 
                       3*Power(t2,2) - 4*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) + 
                       2*s*Power(t2,3) - 2*Power(s,2)*Power(t2,3) - 
                       Power(s1,2)*t1*
                       (-2 + 3*t1 + 2*Power(t1,2) + 2*t2) - 
                       2*t1*t2*
                       (-3 + 8*s - 2*Power(s,2) + 3*t2 + 2*s*t2) + 
                       Power(t1,2)*
                       (3 - 2*s + Power(s,2) - 6*t2 + 
                       2*Power(s,2)*t2) + 
                       2*s1*
                        ((2 + s)*Power(t1,3) + Power(-1 + t2,2) + 
                        Power(t1,2)*(5*s + 3*t2) + 
                        t1*(-3 + s + 2*t2 + 4*s*t2 + Power(t2,2) - 
                        s*Power(t2,2)))) - 
                    Power(s2,2)*
                     (-1 + (1 + Power(s,2))*Power(t1,3) + 3*t2 - 
                       4*s*t2 - 4*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) - 
                       2*s*Power(t2,3) - Power(s,2)*Power(t2,3) + 
                       Power(s1,2)*
                        (1 - 3*t1 + Power(t1,3) - t2 - 
                        Power(t1,2)*(5 + t2)) + 
                       Power(t1,2)*
                        (-(Power(s,2)*(-1 + t2)) - 3*(1 + t2) + 
                        s*(2 + 6*t2)) + 
                       t1*(3 - 6*Power(t2,2) + 
                        Power(s,2)*t2*(4 + t2) - 
                        2*s*(1 + 2*t2 + 2*Power(t2,2))) + 
                       2*s1*(-((1 + s)*Power(t1,3)) + 
                        (1 + t2)*(-1 + t2 + s*t2) + 
                        Power(t1,2)*(4 + (2 + s)*t2) + 
                        t1*(-2 + 4*t2 + Power(t2,2) + 2*s*(1 + t2)))))\
)/((-1 + s2 - t1 + t2)*Power(t1 - s2*t2,2))))/((-1 + s2)*(-s + s2 - t1)) \
+ ((-1 + t1)*((-4*(-1 - (-1 + s + s1)*t1 + Power(s2,2)*(s1 - t2) + 
                    t2 - s*t2 + 
                    s2*(-1 + s1 + t1 + s*t1 - s1*t1 + s*t2))*
                  (2 + 16*s2 + 7*Power(s2,2) - 4*Power(s2,3) + 
                    Power(s1,2)*(-1 + t1) - 3*t1 + 
                    3*Power(s2,2)*t1 + Power(t1,2) + 7*t2 + 
                    10*s2*t2 - 7*Power(s2,2)*t2 - t1*t2 + s2*t1*t2 + 
                    3*Power(t2,2) - 3*s2*Power(t2,2) - 
                    Power(s,2)*(2 + s2 - 2*t1 + t2) + 
                    s1*(2 - 3*Power(s2,2) + s2*(1 + 4*t1 - 5*t2) + 
                       (2 + t1)*t2 - 2*Power(t2,2)) + 
                    s*(5*Power(s2,2) - (9 + t1 - 3*t2)*(1 + t2) + 
                       s1*(3 + s2 - 3*t1 + t2) + 
                       s2*(-4 - 5*t1 + 8*t2))))/(t1 - s2*t2) - 
               8*(3 - 2*s2 - 15*Power(s2,2) - 8*Power(s2,3) + 
                  2*Power(s2,4) + Power(s,3)*(-1 + t1) + 2*t1 + 
                  15*s2*t1 + 7*Power(s2,2)*t1 - 2*Power(s2,3)*t1 - 
                  5*Power(t1,2) - 3*s2*Power(t1,2) - 
                  Power(s,2)*
                   (-2*Power(s2,2) + 2*s1*(-1 + t1) + 
                     s2*(2 + 4*t1 - 4*t2) + 
                     (-1 + t1 - 2*t2)*(-3 + t2)) - 6*t2 - 19*s2*t2 - 
                  15*Power(s2,2)*t2 + 4*Power(s2,3)*t2 + 12*t1*t2 + 
                  9*s2*t1*t2 - Power(s2,2)*t1*t2 - 2*Power(t1,2)*t2 - 
                  11*Power(t2,2) - 9*s2*Power(t2,2) + 
                  2*Power(s2,2)*Power(t2,2) + 4*t1*Power(t2,2) - 
                  2*Power(t2,3) - 
                  Power(s1,2)*(2 - 3*t1 + s2*(4 + t1) + 4*t2) + 
                  s1*(7 + 2*Power(s2,3) + 7*t2 + 2*Power(t2,2) - 
                     t1*(1 + t2) + Power(s2,2)*(-3 - 3*t1 + 4*t2) - 
                     s2*(t1*(-2 + t2) + t2 - 2*Power(t2,2))) + 
                  s*(-3 - 4*Power(s2,3) + Power(s1,2)*(-1 + t1) - 
                     4*t1 + Power(t1,2) + 
                     Power(s2,2)*(7 + 5*t1 - 8*t2) + 3*t2 + 
                     2*Power(t2,2) + 
                     s1*(4 - 2*Power(s2,2) + s2*(4 + 5*t1 - 4*t2) + 
                        t1*(-5 + t2) + 7*t2 - 2*Power(t2,2)) + 
                     s2*(11 + 11*t2 - 4*Power(t2,2) + t1*(-3 + 2*t2)))\
) - (4*(-8 - 4*s2 + 3*Power(s2,2) - 2*t1 - s2*t1 + 
                    s*(1 - 2*s2 + t1 - 2*t2) - 3*t2 + 4*s2*t2 + 
                    Power(t2,2) + s1*(1 + s2 - t1 + t2))*
                  ((Power(s,2) + Power(-1 + s1,2) - 2*s*(1 + 2*s1))*
                     Power(t1,3) + Power(s2,5)*Power(s1 - t2,2) + 
                    (-1 + t2)*Power(1 + (-1 + s)*t2,2) + 
                    t1*(2*s*(-1 + s1)*Power(-1 + t2,2) - 
                       (-3 + 2*s1)*Power(-1 + t2,2) + 
                       Power(s,2)*Power(t2,2)) - 
                    Power(t1,2)*
                     (-((3 - 4*s1 + Power(s1,2))*(-1 + t2)) + 
                       Power(s,2)*(1 + 3*t2) - 
                       2*s*(2 + s1*(-3 + t2) + 3*t2)) + 
                    Power(s2,4)*
                     (Power(s1,2)*(1 - 3*t1 + t2) + 
                       t2*(2 - t2 + Power(t2,2) - 2*s*(1 + t2) - 
                       t1*(2 + t2)) + 
                       2*s1*
                       (-1 + s*t2 - Power(t2,2) + t1*(1 + s + 2*t2))\
) + Power(s2,3)*(1 - 2*s*t2 - 2*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) - 
                       2*s*Power(t2,3) + 
                       Power(t1,2)*(1 + Power(s,2) + 2*t2) + 
                       2*t1*
                       (-1 + (-1 + 4*s)*t2 + (-1 + s)*Power(t2,2)) \
+ Power(s1,2)*(-1 + 3*Power(t1,2) - 2*t1*(2 + t2)) + 
                       2*s1*
                       (-1 + t2 + (1 + s)*Power(t2,2) - 
                       Power(t1,2)*(2 + 2*s + t2) + 
                       t1*(3 - s + 2*t2 + Power(t2,2)))) + 
                    s2*(-1 + 2*(-1 + s)*Power(t1,3) + 2*s*t2 + 
                       3*Power(t2,2) - 4*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) + 
                       2*s*Power(t2,3) - 2*Power(s,2)*Power(t2,3) - 
                       Power(s1,2)*t1*
                       (-2 + 3*t1 + 2*Power(t1,2) + 2*t2) - 
                       2*t1*t2*
                       (-3 + 8*s - 2*Power(s,2) + 3*t2 + 2*s*t2) + 
                       Power(t1,2)*
                       (3 - 2*s + Power(s,2) - 6*t2 + 
                       2*Power(s,2)*t2) + 
                       2*s1*
                       ((2 + s)*Power(t1,3) + Power(-1 + t2,2) + 
                       Power(t1,2)*(5*s + 3*t2) + 
                       t1*(-3 + s + 2*t2 + 4*s*t2 + Power(t2,2) - 
                       s*Power(t2,2)))) - 
                    Power(s2,2)*
                     (-1 + (1 + Power(s,2))*Power(t1,3) + 3*t2 - 
                       4*s*t2 - 4*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) - 
                       2*s*Power(t2,3) - Power(s,2)*Power(t2,3) + 
                       Power(s1,2)*
                       (1 - 3*t1 + Power(t1,3) - t2 - 
                       Power(t1,2)*(5 + t2)) + 
                       Power(t1,2)*
                       (-(Power(s,2)*(-1 + t2)) - 3*(1 + t2) + 
                       s*(2 + 6*t2)) + 
                       t1*(3 - 6*Power(t2,2) + 
                       Power(s,2)*t2*(4 + t2) - 
                       2*s*(1 + 2*t2 + 2*Power(t2,2))) + 
                       2*s1*
                        (-((1 + s)*Power(t1,3)) + 
                        (1 + t2)*(-1 + t2 + s*t2) + 
                        Power(t1,2)*(4 + (2 + s)*t2) + 
                        t1*(-2 + 4*t2 + Power(t2,2) + 2*s*(1 + t2))))\
))/((-1 + s2 - t1 + t2)*Power(t1 - s2*t2,2)) - 
               ((-2 + s2 + t2)*
                  (-1 - (-1 + s + s1)*t1 + Power(s2,2)*(s1 - t2) + 
                    t2 - s*t2 + 
                    s2*(-1 + s1 + t1 + s*t1 - s1*t1 + s*t2))*
                  (4*(-1 + s2 - t1 + t2)*
                     Power(-1 - (-1 + s + s1)*t1 + 
                       Power(s2,2)*(s1 - t2) + t2 - s*t2 + 
                       s2*(-1 + s1 + t1 + s*t1 - s1*t1 + s*t2),2) + 
                    12*(t1 - s2*t2)*
                     (Power(-1 + s1*s2 + t1 - s1*t1 + t2 - s2*t2,2) + 
                       Power(s,2)*(-1 + s2)*
                        (s2*t1 - Power(t1,2) - t2 + t1*t2) - 
                       s*(-1 + s2)*
                        (Power(s2,2)*(-1 + t1) + 2*t1*(1 - t1 + t2) + 
                        s1*(1 + s2 - 2*t1 + s2*t1 - Power(t1,2) - 
                        t2 + t1*t2) - 
                        s2*(1 + Power(t1,2) + 3*t2 - t1*(2 + t2))))))/
                ((-1 + s2 - t1 + t2)*Power(-t1 + s2*t2,3))))/
           ((-1 + s1)*(-s + s1 - t2)) + 
          ((-1 + t1)*((-4*(-1 - (-1 + s + s1)*t1 + 
                    Power(s2,2)*(s1 - t2) + t2 - s*t2 + 
                    s2*(-1 + s1 + t1 + s*t1 - s1*t1 + s*t2))*
                  (4 + 14*s2 + 5*Power(s2,2) - 4*Power(s2,3) + 
                    Power(s1,2)*(-1 + t1) + 3*t1 + 2*s2*t1 + 
                    3*Power(s2,2)*t1 + Power(t1,2) + 
                    s1*(-3*Power(s2,2) + s2*(2 + 4*t1 - 5*t2) + 
                       (-1 + t1 - 2*t2)*(-2 + t2)) + 7*t2 + 5*s2*t2 - 
                    7*Power(s2,2)*t2 + t1*t2 + s2*t1*t2 - 
                    3*s2*Power(t2,2) - 
                    Power(s,2)*(2 + s2 - 2*t1 + t2) + 
                    s*(-9 + 5*Power(s2,2) - 3*t1 - 5*t2 - t1*t2 + 
                       3*Power(t2,2) + s1*(3 + s2 - 3*t1 + t2) + 
                       s2*(-3 - 5*t1 + 8*t2))))/(t1 - s2*t2) - 
               8*(4 - 7*Power(s2,2) - 2*Power(s2,3) + 2*Power(s2,4) + 
                  Power(s,3)*(-1 + t1) - 4*t1 + s2*t1 - 
                  3*Power(s2,2)*t1 - 2*Power(s2,3)*t1 + 
                  s2*Power(t1,2) + t2 - 11*s2*t2 - 3*Power(s2,2)*t2 + 
                  4*Power(s2,3)*t2 - t1*t2 - 3*s2*t1*t2 - 
                  Power(s2,2)*t1*t2 - 7*Power(t2,2) + s2*Power(t2,2) + 
                  2*Power(s2,2)*Power(t2,2) - 
                  Power(s1,2)*(1 + s2*t1 + 2*t2) - 
                  Power(s,2)*
                   (4 - 2*Power(s2,2) + 2*s1*(-1 + t1) - 4*t1 + 
                     s2*(2 + 4*t1 - 4*t2) + 5*t2 + t1*t2 - 
                     2*Power(t2,2)) + 
                  s1*(-1 + 2*Power(s2,3) - Power(t1,2) + 
                     2*Power(t2,2) + 2*t1*(2 + t2) + 
                     Power(s2,2)*(-3 - 3*t1 + 4*t2) - 
                     s2*(6 + t1*(-4 + t2) + 5*t2 - 2*Power(t2,2))) + 
                  s*(-3 - 4*Power(s2,3) + Power(s1,2)*(-1 + t1) - t1 + 
                     2*Power(t1,2) + Power(s2,2)*(5 + 5*t1 - 8*t2) - 
                     3*t2 + 4*t1*t2 - 2*Power(t2,2) + 
                     s2*(6 + 2*t1*(-1 + t2) + 7*t2 - 4*Power(t2,2)) + 
                     s1*(6 - 2*Power(s2,2) + s2*(4 + 5*t1 - 4*t2) + 
                        t1*(-5 + t2) + 9*t2 - 2*Power(t2,2)))) - 
               (4*(-10 - 2*s2 + 3*Power(s2,2) - 4*t1 - s2*t1 + 
                    s*(1 - 2*s2 + t1 - 2*t2) - t2 + 4*s2*t2 + 
                    Power(t2,2) + s1*(1 + s2 - t1 + t2))*
                  ((Power(s,2) + Power(-1 + s1,2) - 2*s*(1 + 2*s1))*
                     Power(t1,3) + Power(s2,5)*Power(s1 - t2,2) + 
                    (-1 + t2)*Power(1 + (-1 + s)*t2,2) + 
                    t1*(2*s*(-1 + s1)*Power(-1 + t2,2) - 
                       (-3 + 2*s1)*Power(-1 + t2,2) + 
                       Power(s,2)*Power(t2,2)) - 
                    Power(t1,2)*
                     (-((3 - 4*s1 + Power(s1,2))*(-1 + t2)) + 
                       Power(s,2)*(1 + 3*t2) - 
                       2*s*(2 + s1*(-3 + t2) + 3*t2)) + 
                    Power(s2,4)*
                     (Power(s1,2)*(1 - 3*t1 + t2) + 
                       t2*(2 - t2 + Power(t2,2) - 2*s*(1 + t2) - 
                       t1*(2 + t2)) + 
                       2*s1*
                        (-1 + s*t2 - Power(t2,2) + t1*(1 + s + 2*t2))\
) + Power(s2,3)*(1 - 2*s*t2 - 2*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) - 
                       2*s*Power(t2,3) + 
                       Power(t1,2)*(1 + Power(s,2) + 2*t2) + 
                       2*t1*
                       (-1 + (-1 + 4*s)*t2 + (-1 + s)*Power(t2,2)) + 
                       Power(s1,2)*
                       (-1 + 3*Power(t1,2) - 2*t1*(2 + t2)) + 
                       2*s1*
                        (-1 + t2 + (1 + s)*Power(t2,2) - 
                        Power(t1,2)*(2 + 2*s + t2) + 
                        t1*(3 - s + 2*t2 + Power(t2,2)))) + 
                    s2*(-1 + 2*(-1 + s)*Power(t1,3) + 2*s*t2 + 
                       3*Power(t2,2) - 4*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) + 
                       2*s*Power(t2,3) - 2*Power(s,2)*Power(t2,3) - 
                       Power(s1,2)*t1*
                       (-2 + 3*t1 + 2*Power(t1,2) + 2*t2) - 
                       2*t1*t2*
                       (-3 + 8*s - 2*Power(s,2) + 3*t2 + 2*s*t2) + 
                       Power(t1,2)*
                       (3 - 2*s + Power(s,2) - 6*t2 + 
                       2*Power(s,2)*t2) + 
                       2*s1*
                        ((2 + s)*Power(t1,3) + Power(-1 + t2,2) + 
                        Power(t1,2)*(5*s + 3*t2) + 
                        t1*(-3 + s + 2*t2 + 4*s*t2 + Power(t2,2) - 
                        s*Power(t2,2)))) - 
                    Power(s2,2)*
                     (-1 + (1 + Power(s,2))*Power(t1,3) + 3*t2 - 
                       4*s*t2 - 4*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) - 
                       2*s*Power(t2,3) - Power(s,2)*Power(t2,3) + 
                       Power(s1,2)*
                        (1 - 3*t1 + Power(t1,3) - t2 - 
                        Power(t1,2)*(5 + t2)) + 
                       Power(t1,2)*
                        (-(Power(s,2)*(-1 + t2)) - 3*(1 + t2) + 
                        s*(2 + 6*t2)) + 
                       t1*(3 - 6*Power(t2,2) + 
                        Power(s,2)*t2*(4 + t2) - 
                        2*s*(1 + 2*t2 + 2*Power(t2,2))) + 
                       2*s1*(-((1 + s)*Power(t1,3)) + 
                        (1 + t2)*(-1 + t2 + s*t2) + 
                        Power(t1,2)*(4 + (2 + s)*t2) + 
                        t1*(-2 + 4*t2 + Power(t2,2) + 2*s*(1 + t2)))))\
)/((-1 + s2 - t1 + t2)*Power(t1 - s2*t2,2)) - 
               ((-2 + s2 + t2)*
                  (-1 - (-1 + s + s1)*t1 + Power(s2,2)*(s1 - t2) + 
                    t2 - s*t2 + 
                    s2*(-1 + s1 + t1 + s*t1 - s1*t1 + s*t2))*
                  (4*(-1 + s2 - t1 + t2)*
                     Power(-1 - (-1 + s + s1)*t1 + 
                       Power(s2,2)*(s1 - t2) + t2 - s*t2 + 
                       s2*(-1 + s1 + t1 + s*t1 - s1*t1 + s*t2),2) + 
                    12*(t1 - s2*t2)*
                     (Power(-1 + s1*s2 + t1 - s1*t1 + t2 - s2*t2,2) + 
                       Power(s,2)*(-1 + s2)*
                        (s2*t1 - Power(t1,2) - t2 + t1*t2) - 
                       s*(-1 + s2)*
                        (Power(s2,2)*(-1 + t1) + 2*t1*(1 - t1 + t2) + 
                        s1*(1 + s2 - 2*t1 + s2*t1 - Power(t1,2) - 
                        t2 + t1*t2) - 
                        s2*(1 + Power(t1,2) + 3*t2 - t1*(2 + t2))))))/
                ((-1 + s2 - t1 + t2)*Power(-t1 + s2*t2,3))))/
           ((s - s2 + t1)*(s - s1 + t2)))/Power(-1 + t1,2))*
     B1(s2,1 - s2 + t1 - t2,t1))/(8.*Power(Pi,2)) + 
  (((-2*(1 - 11*s2 + 2*Power(s,3)*(-1 + s2)*s2 + 7*Power(s2,2) + 
            2*Power(s1,3)*Power(s2,2) + Power(s2,3) - t1 + 9*s2*t1 - 
            9*Power(s2,2)*t1 - Power(s2,3)*t1 + 2*Power(t1,2) + 
            2*s2*Power(t1,2) + 2*Power(s2,2)*Power(t1,2) - 
            2*Power(t1,3) + Power(s1,2)*s2*
             (4 + Power(s2,2) + t1 + s2*(-5 + 3*t1 - 4*t2)) - 5*t2 - 
            s2*t2 - 11*Power(s2,2)*t2 + 6*Power(s2,3)*t2 + 
            Power(s2,4)*t2 - s2*t1*t2 - 3*Power(s2,2)*t1*t2 - 
            2*Power(s2,3)*t1*t2 + 2*Power(t1,2)*t2 + 
            2*s2*Power(t1,2)*t2 - s2*Power(t2,2) + 
            2*Power(s2,2)*Power(t2,2) + Power(s2,3)*Power(t2,2) - 
            2*s2*t1*Power(t2,2) - 
            Power(s,2)*(8 + 4*Power(s2,3) + 
               Power(s2,2)*(6 + 2*s1 - 3*t1 - t2) + 
               s2*(-20 + 3*t1 + t2)) + 
            s1*(-6 - Power(s2,4) + 11*t1 - 2*Power(t1,2) + 
               s2*(4 + t1*(-4 + t2) - 3*t2) + 
               Power(s2,3)*(-7 + 3*t1 - 2*t2) + 
               Power(s2,2)*(10 + t1*(4 - 3*t2) + 3*t2 + 2*Power(t2,2))) \
- s*(-4 + 2*Power(s1,2)*(-1 + s2)*s2 - 2*Power(s2,4) + 11*t1 + 
               2*Power(t1,2) - 5*t2 - 2*t1*t2 + 
               Power(s2,3)*(-6 + 3*t1 + 2*t2) + 
               s2*(-3 + t1*(-10 + t2) - 10*t2 + Power(t2,2)) + 
               s1*(-3*Power(s2,3) + 2*(-3 + t1) + 
                  3*Power(s2,2)*(-5 + 2*t1 - t2) + s2*(20 - 2*t1 + t2)) \
+ Power(s2,2)*(19 + 7*t2 + Power(t2,2) - t1*(2 + 3*t2)))))/
        ((-1 + s2)*(-s + s2 - t1)*(1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2)) \
- (2*(3 + 2*Power(s,4)*(-1 + s2) - 5*s2 + 4*Power(s2,3) - 
            2*Power(s1,3)*s2*(-4 + t1) - 10*t1 + 2*s2*t1 - 
            2*Power(s2,2)*t1 - 4*Power(s2,3)*t1 + 9*Power(t1,2) + 
            3*s2*Power(t1,2) + 2*Power(s2,2)*Power(t1,2) - 
            2*Power(t1,3) + t2 + 7*s2*t2 - 8*Power(s2,2)*t2 + 
            4*Power(s2,4)*t2 - 2*t1*t2 - 16*s2*t1*t2 - 
            2*Power(s2,2)*t1*t2 - 2*Power(s2,3)*t1*t2 - 
            Power(t1,2)*t2 + 3*s2*Power(t1,2)*t2 + 5*s2*Power(t2,2) + 
            8*Power(s2,2)*Power(t2,2) - Power(s2,3)*Power(t2,2) + 
            t1*Power(t2,2) - Power(s2,2)*t1*Power(t2,2) - 
            2*s2*Power(t2,3) + 
            Power(s,3)*(2 + s1*(4 - 6*s2) + (-5 + s2)*t1 - t2 + 
               5*s2*t2) + Power(s,2)*
             (1 - 6*Power(s2,3) + Power(s1,2)*(-2 + 6*s2) - 6*t1 - 
               3*Power(t1,2) + Power(s2,2)*(6 + t1 - 2*t2) + 8*t2 - 
               4*t1*t2 + Power(t2,2) + 
               s1*(-14 - Power(s2,2) + 10*t1 + 
                  s2*(13 - 4*t1 - 9*t2) + t2) - 
               s2*(3 - 13*t1 + Power(t1,2) + 12*t2 - 4*t1*t2 - 
                  3*Power(t2,2))) + 
            Power(s1,2)*(-8 - 5*Power(s2,3) + 12*t1 - 3*Power(t1,2) + 
               Power(s2,2)*(5 + 2*t1) + 
               s2*(2 - Power(t1,2) - 18*t2 + 4*t1*(2 + t2))) + 
            s1*(-2 - 4*Power(s2,4) + 8*t2 + 3*Power(t1,2)*(2 + t2) + 
               2*Power(s2,3)*(t1 + 3*t2) - t1*(2 + 13*t2) - 
               Power(s2,2)*(-10 + (13 + t1)*t2) + 
               s2*(-8 - 7*t2 + Power(t1,2)*t2 + 12*Power(t2,2) - 
                  2*t1*(-7 + 4*t2 + Power(t2,2)))) - 
            s*(5 + 2*Power(s1,3)*s2 - 4*Power(s2,4) - 15*t1 + 
               12*Power(t1,2) + 7*t2 - 11*t1*t2 + 3*Power(t1,2)*t2 + 
               Power(t2,2) - t1*Power(t2,2) + 
               Power(s2,3)*(2 + 2*t1 + 7*t2) + 
               Power(s2,2)*(9 - 3*t1 - 11*t2 + 2*Power(t2,2)) - 
               Power(s1,2)*(Power(s2,2) - 5*(-2 + t1) + 
                  s2*(-19 + 5*t1 + 4*t2)) + 
               s2*(-8 + Power(t1,2)*(-2 + t2) + t2 + 15*Power(t2,2) + 
                  t1*(10 - 15*t2 - 3*Power(t2,2))) - 
               s1*(11 + 11*Power(s2,3) + 6*Power(t1,2) + 
                  4*t1*(-3 + t2) - 9*t2 + 
                  Power(s2,2)*(-10 - 3*t1 + t2) + 
                  s2*(2 + 2*Power(t1,2) + 34*t2 - 2*Power(t2,2) - 
                     t1*(17 + 8*t2))))))/
        ((s - s2 + t1)*(s - s1 + t2)*
          (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2)) - 
       (2*(-1 + s1 + t1 - t2)*
          (-2*Power(s,2)*(-1 + s2) - 
            (1 + s2)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2) + 
            s*(-4 + 2*s1*s2 + 2*Power(s2,2) + t1 + t2 - s2*(t1 + t2)))*
          ((-1 - Power(s,4) + 24*s2 + 24*Power(s2,2) + 
               3*Power(s2,3) - Power(s2,4) - 
               Power(s1,2)*(1 + s2)*(2 + s2 - t1) + 
               Power(s,3)*(-2 + 2*s1 + 4*s2 - t1) - 7*t1 - 29*s2*t1 - 
               11*Power(s2,2)*t1 + Power(s2,3)*t1 + 4*Power(t1,2) + 
               5*s2*Power(t1,2) - 
               Power(s,2)*(-21 + Power(s1,2) + 6*Power(s2,2) + 
                  s1*(3 + 6*s2 - 2*t1) + 9*t1 - s2*(10 + 3*t1) - 5*t2\
) - 4*t2 + 12*s2*t2 + 4*Power(s2,2)*t2 + t1*t2 - 5*s2*t1*t2 - 
               Power(t2,2) + s2*Power(t2,2) + 
               s*(-19 + 4*Power(s2,3) + 
                  Power(s1,2)*(4 + 2*s2 - t1) + 18*t1 - 
                  2*Power(t1,2) - Power(s2,2)*(11 + 3*t1) + 
                  s2*(-49 + 22*t1 - 9*t2) + 
                  s1*(-13 + 6*Power(s2,2) - 4*s2*(-1 + t1) + 6*t1 - 
                     5*t2) - 8*t2 + 2*t1*t2) + 
               s1*(10 - 2*Power(s2,3) + 2*Power(s2,2)*(-1 + t1) + 
                  Power(t1,2) + 9*t2 - 2*t1*(4 + t2) + 
                  s2*(11 - 4*t1 + t2)))/((-1 + s2)*(-1 + t1)) + 
            (-10 + Power(s,3)*(-2 + s2) + 10*s2 + 9*Power(s2,2) - 
               Power(s2,3) - Power(s2,4) + 12*t1 + 5*s2*t1 + 
               2*Power(s2,2)*t1 + Power(s2,3)*t1 - 4*Power(t1,2) - 
               5*s2*Power(t1,2) - 
               Power(s1,2)*(-5 + s2 + Power(s2,2) + 2*t1 - s2*t1) + 
               Power(s,2)*(8 - 2*s1*(-2 + s2) - 3*Power(s2,2) + 
                  s2*(2 + t1) - 2*t2) - t2 - 3*s2*t2 - 
               4*Power(s2,2)*t2 - t1*t2 + 5*s2*t1*t2 + Power(t2,2) - 
               s2*Power(t2,2) + 
               s1*(10 - 2*Power(s2,3) + t1 - Power(t1,2) + 
                  2*Power(s2,2)*(1 + t1) - 5*t2 + 2*t1*t2 - 
                  s2*(-16 + 8*t1 + t2)) + 
               s*(Power(s1,2)*(-2 + s2) + 3*Power(s2,3) + 
                  Power(s2,2)*(1 - 2*t1) - 7*t1 + 2*Power(t1,2) + 
                  4*t2 - 2*t1*t2 + 
                  s1*(-15 + 4*Power(s2,2) + 3*t1 - s2*(5 + 2*t1) + 
                     2*t2) + s2*(-17 - 2*t1 + 6*t2)))/
             ((-1 + s2)*(-s + s2 - t1)) + 
            (-6 - 2*Power(s1,3)*(-3 + s2) - 3*s2 + 20*Power(s2,2) + 
               8*Power(s2,3) - 2*Power(s2,4) - 5*t1 - s2*t1 - 
               2*Power(s2,2)*t1 + Power(t1,2) + 2*s2*Power(t1,2) + 
               6*t2 + 6*s2*t2 - 11*Power(s2,2)*t2 + 
               2*Power(s2,3)*t2 - 3*t1*t2 - 2*s2*t1*t2 - 
               2*Power(t2,2) + 4*s2*Power(t2,2) + 
               Power(s,3)*(-5 + 2*s2 + t2) - 
               Power(s1,2)*(1 + 6*Power(s2,2) + t1 + 8*t2 - 
                  2*s2*(9 + t2)) + 
               Power(s,2)*(6 + 15*s2 - 6*Power(s2,2) - 6*t2 + 
                  Power(t2,2) - 2*s1*(-7 + 3*s2 + t2)) + 
               s1*(-10 - 6*Power(s2,3) + Power(t1,2) + 7*t2 + 
                  2*Power(t2,2) + t1*(6 + t2) + 
                  Power(s2,2)*(19 + 4*t2) - s2*(-9 + t1 + 20*t2)) + 
               s*(3 + 6*Power(s2,3) + t1 - 10*t2 - 2*Power(t2,2) - 
                  3*Power(s2,2)*(6 + t2) + 
                  Power(s1,2)*(-16 + 6*s2 + t2) + 
                  s2*(-20 + 2*t1 + 15*t2 - Power(t2,2)) + 
                  s1*(1 + 12*Power(s2,2) + t1 + 15*t2 - Power(t2,2) - 
                     2*s2*(15 + t2))))/((-1 + s1)*(-s + s1 - t2)) + 
            (-11 + Power(s,4) - Power(s1,3)*(-3 + s2) - 2*s2 + 
               12*Power(s2,2) + Power(s2,3) - 2*Power(s2,4) + 10*t1 + 
               10*s2*t1 + 2*Power(s2,2)*t1 + Power(s2,3)*t1 - 
               Power(t1,2) + 6*t2 - 6*s2*t2 - 8*Power(s2,2)*t2 + 
               Power(s2,3)*t2 - 4*t1*t2 - s2*t1*t2 - Power(t2,2) + 
               3*s2*Power(t2,2) + 
               Power(s,3)*(-3 - 2*s1 - s2 + t1 + t2) + 
               Power(s1,2)*(5 - 4*Power(s2,2) - t1 - 5*t2 + 
                  s2*(9 + t1 + t2)) + 
               Power(s,2)*(-3 + Power(s1,2) - 3*Power(s2,2) + 3*t1 - 
                  8*t2 + t1*t2 - s2*(-8 + t1 + t2) - 
                  s1*(-15 + s2 + 2*t1 + t2)) + 
               s1*(-1 - 5*Power(s2,3) + Power(t1,2) - 4*t2 + 
                  2*Power(t2,2) + t1*(4 + t2) + 
                  2*Power(s2,2)*(4 + t1 + t2) - s2*(-20 + t1 + 15*t2)\
) + s*(13 + 5*Power(s2,3) - 12*t1 + Power(s1,2)*(-13 + 3*s2 + t1) - 
                  2*Power(t2,2) - Power(s2,2)*(6 + t1 + t2) - 
                  s2*(9 - 16*t2 + t1*(5 + t2)) + 
                  s1*(-9 + 8*Power(s2,2) + t1 + 14*t2 - t1*t2 - 
                     s2*(24 + t2))))/((s - s2 + t1)*(s - s1 + t2)) + 
            (-2*Power(s,4) - Power(s1,3)*(-1 + s2) + 17*s2 + 
               22*Power(s2,2) + 3*Power(s2,3) - 2*Power(s2,4) - 
               9*t1 - 13*s2*t1 - 2*Power(s2,2)*t1 + Power(s2,3)*t1 + 
               Power(t1,2) + 
               Power(s,3)*(-2 + 5*s1 + 8*s2 - t1 - 2*t2) - 3*t2 - 
               6*s2*t2 - 3*Power(s2,2)*t2 + Power(s2,3)*t2 + 
               4*t1*t2 + s2*t1*t2 + Power(t2,2) - 3*s2*Power(t2,2) + 
               Power(s1,2)*(3 - 4*Power(s2,2) + t2 + 
                  s2*(4 + t1 + t2)) - 
               s1*(-4 + 5*Power(s2,3) + Power(t1,2) + 3*t2 + t1*t2 + 
                  2*Power(t2,2) + s2*(-17 + 2*t2) - 
                  Power(s2,2)*(5 + 2*t1 + 2*t2)) + 
               Power(s,2)*(17 - 4*Power(s1,2) - 12*Power(s2,2) + 
                  t1 - 3*t2 - t1*t2 + s1*(4 - 15*s2 + 2*t1 + 3*t2) + 
                  s2*(10 + 3*t1 + 5*t2)) + 
               s*(-15 + Power(s1,3) + 8*Power(s2,3) + 9*t1 + 
                  Power(t1,2) + Power(s1,2)*(-4 + 8*s2 - t1 - t2) + 
                  7*t2 + 2*Power(t2,2) - 
                  Power(s2,2)*(11 + 3*t1 + 4*t2) + 
                  s2*(-41 + t1 + 8*t2 + t1*t2) + 
                  s1*(-12 + 15*Power(s2,2) + t1*(-2 + t2) + t2 - 
                     s2*(12 + 4*t1 + 5*t2))))/((-1 + t1)*(-1 + t2)) + 
            (-10 - 2*Power(s,4) - 2*Power(s1,3)*(-2 + s2) + 4*s2 + 
               18*Power(s2,2) + 7*Power(s2,3) - 2*Power(s2,4) - t1 + 
               s2*t1 + 2*Power(s2,2)*t1 - Power(t1,2) - 
               2*s2*Power(t1,2) + 6*t2 - 11*s2*t2 - 
               10*Power(s2,2)*t2 + 2*Power(s2,3)*t2 + 3*t1*t2 + 
               2*s2*t1*t2 + 2*Power(t2,2) - 4*s2*Power(t2,2) + 
               Power(s,3)*(6*s1 + 8*s2 - 3*(1 + t2)) + 
               Power(s1,2)*(2 - 6*Power(s2,2) + t1 - 2*t2 + 
                  s2*(13 + 2*t2)) - 
               Power(s,2)*(-12 + 6*Power(s1,2) + 12*Power(s2,2) + 
                  3*s1*(-3 + 6*s2 - 2*t2) + 4*t2 + Power(t2,2) - 
                  2*s2*(5 + 4*t2)) - 
               s1*(3 + 6*Power(s2,3) + Power(t1,2) - 
                  s2*(16 + t1 - 7*t2) + 4*t2 + t1*t2 + 
                  2*Power(t2,2) - Power(s2,2)*(17 + 4*t2)) + 
               s*(4 + 2*Power(s1,3) + 8*Power(s2,3) + 
                  Power(s1,2)*(-10 + 12*s2 - 3*t2) + 7*t2 + 
                  2*Power(t2,2) - 7*Power(s2,2)*(2 + t2) + 
                  s2*(-28 - 2*t1 + 12*t2 + Power(t2,2)) + 
                  s1*(-13 + 18*Power(s2,2) - 2*t1 + 5*t2 + 
                     Power(t2,2) - s2*(23 + 10*t2))))/
             ((-1 + s1)*(-1 + t2))))/
        (Power(s,3)*(-1 + s2)*(-1 + s1 + t1 - t2) - 
          Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
          Power(s,2)*(3 + Power(s1,2)*s2 - 3*t1 + Power(t1,2) + 
             2*Power(s2,2)*(-1 + t1 - t2) + 2*t2 - t1*t2 + 
             s2*t1*(1 - t1 + t2) + 
             s1*(-3 + 2*Power(s2,2) + t2 - s2*(1 + t2))) + 
          s*(Power(s1,2)*(1 + s2)*(s2 - t1) + 
             Power(s2,3)*(-1 + t1 - t2) - 
             Power(s2,2)*(-1 + t1)*(-1 + t1 - t2) - 
             s2*(-2 + t1 + Power(t1,2) - 4*t2 - t1*t2) + 
             2*(1 + Power(t1,2) - t2 - t1*(2 + t2)) + 
             s1*(-2 + Power(s2,3) + 5*t1 + Power(t1,2) + 2*t2 - 
                Power(s2,2)*t2 - t1*t2 - 
                s2*(3 + t1 + Power(t1,2) + t2 - t1*t2)))) + 
       (288 - 320*s - 96*Power(s,2) - 352*s1 + 480*s*s1 - 
          32*Power(s1,2) - 352*s2 + 480*s*s2 + 96*Power(s,2)*s2 + 
          352*s1*s2 - 544*s*s1*s2 + 32*Power(s1,2)*s2 - 
          96*s*Power(s2,2) + 64*s1*Power(s2,2) - 608*t1 + 320*s*t1 + 
          352*s1*t1 + 736*s2*t1 - 384*s*s2*t1 - 416*s1*s2*t1 + 
          320*Power(t1,2) - 384*s2*Power(t1,2) + 
          64*(-1 + s2)*(-1 + s + t1)*(-1 + s1 + t1 - t2) + 320*t2 - 
          448*s*t2 + 64*s1*t2 - 320*s2*t2 + 512*s*s2*t2 - 64*s1*s2*t2 - 
          64*Power(s2,2)*t2 - 320*t1*t2 + 384*s2*t1*t2 - 
          32*Power(t2,2) + 32*s2*Power(t2,2) - 
          16*(27 - 45*s2 + Power(s2,3) + 
             Power(s1,2)*(-7 + 8*s2 - t1) - 64*t1 + 97*s2*t1 + 
             Power(s2,2)*t1 + 36*Power(t1,2) - 53*s2*Power(t1,2) + 
             39*t2 - 39*s2*t2 - 18*Power(s2,2)*t2 - 37*t1*t2 + 
             53*s2*t1*t2 - 6*Power(t2,2) + 8*s2*Power(t2,2) + 
             Power(s,2)*(-23 + 26*s2 - 2*t1 + 2*t2) + 
             s1*(-45 + 19*Power(s2,2) + 44*t1 + 
                s2*(47 - 63*t1 - 16*t2) + 14*t2) - 
             s*(29 + 27*Power(s2,2) - 35*t1 + 
                s2*(-69 + 50*t1 - 84*t2) + 67*t2 + 
                s1*(-74 + 95*s2 - 3*t1 + t2))) - 
          (2*(-9 + Power(s,2) - s2 - Power(s2,2) - s1*(-6 + t1) + 
               5*t1 + s*(-3 - s1 + t1) - 3*t2)*
             (2*Power(s,2)*(-1 + s2) + 
               (1 + s2)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2) - 
               s*(-4 + 2*s1*s2 + 2*Power(s2,2) + t1 + t2 - 
                  s2*(t1 + t2))))/
           (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2) + 
          8*(30 + 2*Power(s,3) - 69*s2 + 2*Power(s2,2) + 
             2*Power(s2,3) + Power(s1,2)*(-19 + 21*s2 - 3*t1) - 
             94*t1 + 169*s2*t1 + 9*Power(s2,2)*t1 + 55*Power(t1,2) - 
             102*s2*Power(t1,2) + 71*t2 - 67*s2*t2 - 
             51*Power(s2,2)*t2 - 55*t1*t2 + 100*s2*t1*t2 - 
             14*Power(t2,2) + 21*s2*Power(t2,2) + 
             Power(s,2)*(-56 - 2*s1 + 62*s2 - 3*t1 + 5*t2) + 
             s1*(-78 + 52*Power(s2,2) - Power(t1,2) + 
                s2*(84 - 124*t1 - 43*t2) + 34*t2 + t1*(69 + t2)) + 
             s*(-38 - 66*Power(s2,2) + 65*t1 + Power(t1,2) + 
                s1*(155 - 209*s2 + 6*t1 - 2*t2) - 135*t2 - t1*t2 + 
                s2*(129 - 101*t1 + 184*t2))) - 
          4*(-1 + 2*Power(s,3) - 58*s2 - Power(s2,2) + Power(s2,3) + 
             Power(s1,2)*(-13 + 19*s2 - 3*t1) - 46*t1 + 117*s2*t1 + 
             13*Power(s2,2)*t1 + 29*Power(t1,2) - 72*s2*Power(t1,2) - 
             Power(s,2)*(44 + s1 - 49*s2 - 2*t2) + 48*t2 - 46*s2*t2 - 
             46*Power(s2,2)*t2 - 25*t1*t2 + 68*s2*t1*t2 - 
             12*Power(t2,2) + 18*s2*Power(t2,2) + 
             s1*(-58 + 46*Power(s2,2) - 2*Power(t1,2) + 
                s2*(59 - 86*t1 - 38*t2) + 25*t2 + 2*t1*(19 + t2)) + 
             s*(-Power(s1,2) - 52*Power(s2,2) + 
                s1*(-163*s2 + 3*(36 + t1)) + 
                s2*(92 - 79*t1 + 144*t2) + 
                2*(1 + 21*t1 + Power(t1,2) - 46*t2 - t1*t2))) + 
          ((-1 + s2)*(4*(-1 + s1 + t1 - t2)*
                Power(2*Power(s,2)*(-1 + s2) + 
                  (1 + s2)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2) - 
                  s*(-4 + 2*s1*s2 + 2*Power(s2,2) + t1 + t2 - 
                     s2*(t1 + t2)),2) - 
               2*(6*Power(s,4)*Power(-1 + s2,2)*(-1 + s1 + t1 - t2) + 
                  (-3 + s1*(1 + 4*s2 + Power(s2,2)) + 3*t1 + 
                     2*s2*(-1 + t1 - 2*t2) + 
                     Power(s2,2)*(-1 + t1 - t2) - t2)*
                   Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                  2*Power(s,3)*(-1 + s2)*
                   (12 + 6*Power(s1,2)*s2 - 15*t1 + 4*Power(t1,2) + 
                     6*Power(s2,2)*(-1 + t1 - t2) + 9*t2 - 2*t1*t2 - 
                     2*Power(t2,2) + 
                     3*s1*(-4 + 2*Power(s2,2) + t1 + 
                        s2*(-2 + t1 - 3*t2) + t2) + 
                     3*s2*(t1 - Power(t1,2) + t2 + Power(t2,2))) + 
                  Power(s,2)*
                   (-30 + 6*Power(s1,3)*Power(s2,2) + 52*t1 - 
                     25*Power(t1,2) + 3*Power(t1,3) + 
                     6*Power(s2,4)*(-1 + t1 - t2) - 8*t2 + 8*t1*t2 - 
                     Power(t1,2)*t2 + 11*Power(t2,2) - 
                     t1*Power(t2,2) - Power(t2,3) - 
                     6*Power(s2,3)*
                      (Power(t1,2) + t1*(-1 + t2) - 2*t2*(1 + t2)) + 
                     6*Power(s1,2)*
                      (3*Power(s2,3) + t1 + s2*(-5 + t1 + t2) - 
                        Power(s2,2)*(1 + t1 + 2*t2)) - 
                     2*s2*(2*Power(t1,3) + Power(t1,2)*(-10 + t2) + 
                        t1*(8 - 3*t2 - 2*Power(t2,2)) + 
                        t2*(15 + 7*t2 - Power(t2,2))) + 
                     Power(s2,2)*
                      (30 + Power(t1,3) + 26*t2 - 9*Power(t2,2) - 
                        Power(t2,3) + Power(t1,2)*(11 + 3*t2) - 
                        t1*(42 + 8*t2 + 3*Power(t2,2))) + 
                     s1*(30 + 6*Power(s2,4) + 3*Power(t1,2) + 
                        6*Power(s2,3)*(-3 + 2*t1 - 5*t2) + 
                        2*t1*(-16 + t2) - 18*t2 + Power(t2,2) + 
                        2*s2*
                        (15 + 4*Power(t1,2) + 23*t2 - 
                        4*Power(t2,2) - 6*t1*(1 + t2)) + 
                        Power(s2,2)*
                        (-30 - 11*Power(t1,2) + 14*t2 + 
                        7*Power(t2,2) + 10*t1*(2 + t2)))) - 
                  2*s*(-7 + 3*Power(s1,3)*s2*(1 + s2)*(s2 - t1) + 
                     17*t1 - 13*Power(t1,2) + 3*Power(t1,3) + 5*t2 - 
                     4*t1*t2 - Power(t1,2)*t2 + 3*Power(t2,2) - 
                     t1*Power(t2,2) - Power(t2,3) + 
                     3*Power(s2,4)*t2*(1 - t1 + t2) + 
                     Power(s2,2)*
                      (3 - 2*Power(t1,3) - 4*t2 - 13*Power(t2,2) + 
                        Power(t2,3) + Power(t1,2)*(7 + t2) + 
                        t1*(-8 + 3*t2)) + 
                     Power(s2,3)*
                      (3 + 4*t2 + Power(t2,2) - Power(t2,3) + 
                        Power(t1,2)*(3 + 2*t2) - 
                        t1*(6 + 6*t2 + Power(t2,2))) + 
                     s2*(-5 - Power(t1,3) - 14*t2 + 6*Power(t2,2) + 
                        Power(t2,3) - Power(t1,2)*(3 + 2*t2) + 
                        t1*(9 + 16*t2 + 2*Power(t2,2))) + 
                     Power(s1,2)*
                      (3*Power(s2,4) - t1*(-5 + 3*t1 + t2) - 
                        Power(s2,3)*(-1 + t1 + 7*t2) + 
                        s2*
                       (-9 - Power(t1,2) + 5*t2 + 2*t1*(9 + t2)) + 
                        Power(s2,2)*
                        (-13 - 2*Power(t1,2) - 4*t2 + t1*(2 + 5*t2))) \
+ s1*(6 + 3*Power(s2,4)*(-1 + t1 - 2*t2) - Power(t1,2)*(-12 + t2) - 
                        7*t2 + Power(t2,2) + 
                        t1*(-18 + 4*t2 + Power(t2,2)) + 
                        Power(s2,3)*
                        (-6 - 4*Power(t1,2) - 2*t2 + 5*Power(t2,2) + 
                        2*t1*(5 + t2)) + 
                        Power(s2,2)*
                        (4 + Power(t1,3) + Power(t1,2)*(-3 + t2) + 
                        25*t2 - 2*t1*(1 + Power(t2,2))) + 
                        s2*(17 + 7*Power(t1,2) - Power(t1,3) + 2*t2 - 
                        6*Power(t2,2) + t1*(-23 - 18*t2 + Power(t2,2)))\
)))))/((-1 + s1 + t1 - t2)*Power(-1 + s - s*s2 + s1*s2 + t1 - s2*t2,2)))/
        ((-1 + s2)*(-1 + t1)) + 
       (-64*(-1 + s2)*(-1 + s1 + t1 - t2)*(s - s1 + t2) + 
          (2*(6 + 2*Power(s,2) + Power(s1,2) + 4*s2 - 
               2*Power(s2,2) - 3*t2 - s2*t2 + Power(t2,2) - 
               s1*(-1 + s2 + 2*t2) + s*(-2 - 3*s1 + 3*t2))*
             (2*Power(s,2)*(-1 + s2) + 
               (1 + s2)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2) - 
               s*(-4 + 2*s1*s2 + 2*Power(s2,2) + t1 + t2 - 
                  s2*(t1 + t2))))/
           (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2) - 
          32*(3 + 2*Power(s,3) - Power(s1,3) + s2 - Power(s2,2) - 
             Power(s2,3) - 8*t1 + 2*s2*t1 + 2*Power(s2,2)*t1 + 
             4*Power(t1,2) - 2*s2*Power(t1,2) + 
             Power(s1,2)*(13*s2 - 2*(5 + t1 - t2)) - 7*t2 + 12*s2*t2 - 
             Power(s2,2)*t2 + 8*t1*t2 - 12*s2*t1*t2 - 12*Power(t2,2) + 
             13*s2*Power(t2,2) + 
             Power(s,2)*(-2 - 2*s1 - 2*s2 + 3*t1 + t2) + 
             s1*(5 + Power(s2,2) - 2*Power(t1,2) + 
                2*s2*(-6 + 6*t1 - 13*t2) + 2*t1*(-2 + t2) + 22*t2 - 
                Power(t2,2)) + 
             s*(-15 + Power(s1,2) + Power(s2,2) + 14*t1 + 
                2*Power(t1,2) - 16*t2 - t1*t2 + Power(t2,2) - 
                s1*(-18 + 15*s2 + t1 + 2*t2) + s2*(16 - 18*t1 + 15*t2))\
) + 16*(21 + 16*Power(s,3) - 7*Power(s1,3) + 6*s2 - 9*Power(s2,2) - 
             10*Power(s2,3) - 61*t1 + 22*s2*t1 + 19*Power(s2,2)*t1 + 
             30*Power(t1,2) - 18*s2*Power(t1,2) - 11*t2 + 49*s2*t2 - 
             9*Power(s2,2)*t2 + 22*t1*t2 - 51*s2*t1*t2 - 
             56*Power(t2,2) + 62*s2*Power(t2,2) + 
             Power(s,2)*(-23 - 14*s1 - 17*s2 + 24*t1 + 6*t2) + 
             2*Power(s1,2)*(-19 + 30*s2 - 8*t1 + 7*t2) + 
             s1*(-3 + 7*Power(s2,2) - 16*Power(t1,2) + 
                s2*(-49 + 51*t1 - 122*t2) + 94*t2 - 7*Power(t2,2) + 
                2*t1*(5 + 8*t2)) + 
             s*(-70 + 5*Power(s1,2) + 11*Power(s2,2) + 63*t1 + 
                16*Power(t1,2) - 91*t2 - 8*t1*t2 + 6*Power(t2,2) - 
                s1*(-105 + 79*s2 + 8*t1 + 11*t2) + 
                s2*(88 - 104*t1 + 81*t2))) + 
          4*(48 + 29*Power(s,3) - 12*Power(s1,3) + 5*s2 - 
             20*Power(s2,2) - 22*Power(s2,3) - 123*t1 + 67*s2*t1 + 
             44*Power(s2,2)*t1 + 61*Power(t1,2) - 40*s2*Power(t1,2) + 
             32*t2 + 61*s2*t2 - 24*Power(s2,2)*t2 + 2*t1*t2 - 
             72*s2*t1*t2 - 81*Power(t2,2) + 99*s2*Power(t2,2) - 
             Power(s,2)*(77 + 10*s1 + 28*s2 - 58*t1 + 3*t2) + 
             Power(s1,2)*(-38 + 88*s2 - 36*t1 + 25*t2) + 
             s1*(-58 + 18*Power(s2,2) - 36*Power(t1,2) + 
                s2*(-53 + 64*t1 - 185*t2) + 121*t2 - 13*Power(t2,2) + 
                6*t1*(11 + 6*t2)) + 
             s*(-82 - 7*Power(s1,2) + 21*Power(s2,2) + 66*t1 + 
                36*Power(t1,2) - 2*s1*(-93 + 75*s2 + 11*t1 - t2) - 
                163*t2 - 14*t1*t2 + 4*Power(t2,2) + 
                s2*(172 - 192*t1 + 151*t2))) - 
          8*(54 + 39*Power(s,3) - 16*Power(s1,3) - 29*Power(s2,2) - 
             25*Power(s2,3) - 153*t1 + 71*s2*t1 + 52*Power(s2,2)*t1 + 
             74*Power(t1,2) - 48*s2*Power(t1,2) + 14*t2 + 94*s2*t2 - 
             27*Power(s2,2)*t2 + 21*t1*t2 - 98*s2*t1*t2 - 
             113*Power(t2,2) + 129*s2*Power(t2,2) + 
             Power(s,2)*(-80 - 26*s1 - 39*s2 + 65*t1 + 8*t2) + 
             Power(s1,2)*(-65 + 122*s2 - 42*t1 + 33*t2) + 
             s1*(-50 + 23*Power(s2,2) - 42*Power(t1,2) + 
                s2*(-90 + 94*t1 - 249*t2) + 180*t2 - 17*Power(t2,2) + 
                t1*(61 + 42*t2)) + 
             s*(-122 + 3*Power(s1,2) + 25*Power(s2,2) + 111*t1 + 
                42*Power(t1,2) - 212*t2 - 19*t1*t2 + 11*Power(t2,2) - 
                s1*(-244 + 184*s2 + 23*t1 + 15*t2) + 
                s2*(214 - 240*t1 + 185*t2))) - 
          ((s - s1 + t2)*(4*(-1 + s1 + t1 - t2)*
                Power(2*Power(s,2)*(-1 + s2) + 
                  (1 + s2)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2) - 
                  s*(-4 + 2*s1*s2 + 2*Power(s2,2) + t1 + t2 - 
                     s2*(t1 + t2)),2) - 
               2*(6*Power(s,4)*Power(-1 + s2,2)*(-1 + s1 + t1 - t2) + 
                  (-3 + s1*(1 + 4*s2 + Power(s2,2)) + 3*t1 + 
                     2*s2*(-1 + t1 - 2*t2) + 
                     Power(s2,2)*(-1 + t1 - t2) - t2)*
                   Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                  2*Power(s,3)*(-1 + s2)*
                   (12 + 6*Power(s1,2)*s2 - 15*t1 + 4*Power(t1,2) + 
                     6*Power(s2,2)*(-1 + t1 - t2) + 9*t2 - 2*t1*t2 - 
                     2*Power(t2,2) + 
                     3*s1*(-4 + 2*Power(s2,2) + t1 + 
                        s2*(-2 + t1 - 3*t2) + t2) + 
                     3*s2*(t1 - Power(t1,2) + t2 + Power(t2,2))) + 
                  Power(s,2)*
                   (-30 + 6*Power(s1,3)*Power(s2,2) + 52*t1 - 
                     25*Power(t1,2) + 3*Power(t1,3) + 
                     6*Power(s2,4)*(-1 + t1 - t2) - 8*t2 + 8*t1*t2 - 
                     Power(t1,2)*t2 + 11*Power(t2,2) - 
                     t1*Power(t2,2) - Power(t2,3) - 
                     6*Power(s2,3)*
                      (Power(t1,2) + t1*(-1 + t2) - 2*t2*(1 + t2)) + 
                     6*Power(s1,2)*
                      (3*Power(s2,3) + t1 + s2*(-5 + t1 + t2) - 
                        Power(s2,2)*(1 + t1 + 2*t2)) - 
                     2*s2*(2*Power(t1,3) + Power(t1,2)*(-10 + t2) + 
                        t1*(8 - 3*t2 - 2*Power(t2,2)) + 
                        t2*(15 + 7*t2 - Power(t2,2))) + 
                     Power(s2,2)*
                      (30 + Power(t1,3) + 26*t2 - 9*Power(t2,2) - 
                        Power(t2,3) + Power(t1,2)*(11 + 3*t2) - 
                        t1*(42 + 8*t2 + 3*Power(t2,2))) + 
                     s1*(30 + 6*Power(s2,4) + 3*Power(t1,2) + 
                        6*Power(s2,3)*(-3 + 2*t1 - 5*t2) + 
                        2*t1*(-16 + t2) - 18*t2 + Power(t2,2) + 
                        2*s2*
                        (15 + 4*Power(t1,2) + 23*t2 - 
                        4*Power(t2,2) - 6*t1*(1 + t2)) + 
                        Power(s2,2)*
                        (-30 - 11*Power(t1,2) + 14*t2 + 
                        7*Power(t2,2) + 10*t1*(2 + t2)))) - 
                  2*s*(-7 + 3*Power(s1,3)*s2*(1 + s2)*(s2 - t1) + 
                     17*t1 - 13*Power(t1,2) + 3*Power(t1,3) + 5*t2 - 
                     4*t1*t2 - Power(t1,2)*t2 + 3*Power(t2,2) - 
                     t1*Power(t2,2) - Power(t2,3) + 
                     3*Power(s2,4)*t2*(1 - t1 + t2) + 
                     Power(s2,2)*
                      (3 - 2*Power(t1,3) - 4*t2 - 13*Power(t2,2) + 
                        Power(t2,3) + Power(t1,2)*(7 + t2) + 
                        t1*(-8 + 3*t2)) + 
                     Power(s2,3)*
                      (3 + 4*t2 + Power(t2,2) - Power(t2,3) + 
                        Power(t1,2)*(3 + 2*t2) - 
                        t1*(6 + 6*t2 + Power(t2,2))) + 
                     s2*(-5 - Power(t1,3) - 14*t2 + 6*Power(t2,2) + 
                        Power(t2,3) - Power(t1,2)*(3 + 2*t2) + 
                        t1*(9 + 16*t2 + 2*Power(t2,2))) + 
                     Power(s1,2)*
                      (3*Power(s2,4) - t1*(-5 + 3*t1 + t2) - 
                        Power(s2,3)*(-1 + t1 + 7*t2) + 
                        s2*
                       (-9 - Power(t1,2) + 5*t2 + 2*t1*(9 + t2)) + 
                        Power(s2,2)*
                        (-13 - 2*Power(t1,2) - 4*t2 + t1*(2 + 5*t2))) \
+ s1*(6 + 3*Power(s2,4)*(-1 + t1 - 2*t2) - Power(t1,2)*(-12 + t2) - 
                        7*t2 + Power(t2,2) + 
                        t1*(-18 + 4*t2 + Power(t2,2)) + 
                        Power(s2,3)*
                        (-6 - 4*Power(t1,2) - 2*t2 + 5*Power(t2,2) + 
                        2*t1*(5 + t2)) + 
                        Power(s2,2)*
                        (4 + Power(t1,3) + Power(t1,2)*(-3 + t2) + 
                        25*t2 - 2*t1*(1 + Power(t2,2))) + 
                        s2*(17 + 7*Power(t1,2) - Power(t1,3) + 2*t2 - 
                        6*Power(t2,2) + t1*(-23 - 18*t2 + Power(t2,2)))\
)))))/((-1 + s1 + t1 - t2)*Power(-1 + s - s*s2 + s1*s2 + t1 - s2*t2,2)))/
        ((-1 + s1)*(-s + s1 - t2)) + 
       (2*(-32*s*(-1 + s2)*(-1 + s1 + t1 - t2) - 
            ((-12 + 6*Power(s,2) + 3*Power(s1,2) - 13*s2 + 
                 6*Power(s2,2) + 2*t1 - s2*t1 + 
                 s1*(-4 + 9*s2 - t1 - 3*t2) + 2*t2 - 5*s2*t2 + 
                 t1*t2 + s*(10 - 9*s1 - 12*s2 + t1 + 6*t2))*
               (2*Power(s,2)*(-1 + s2) + 
                 (1 + s2)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2) - 
                 s*(-4 + 2*s1*s2 + 2*Power(s2,2) + t1 + t2 - 
                    s2*(t1 + t2))))/
             (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2) - 
            16*(2 + 2*Power(s,2)*(-1 + s2) - Power(s2,2) + 
               Power(s1,2)*(s2 - t1) - 4*t1 + s2*t1 + Power(s2,2)*t1 + 
               2*Power(t1,2) - s2*Power(t1,2) - 3*t2 + 3*s2*t2 - 
               Power(s2,2)*t2 + 2*t1*t2 - s2*t1*t2 - Power(t2,2) + 
               s2*Power(t2,2) + 
               s1*(1 + Power(s2,2) + t1 - Power(t1,2) + t2 + t1*t2 - 
                  2*s2*(1 + t2)) + 
               s*(-14 - 2*Power(s2,2) + 12*t1 + Power(t1,2) + 
                  s1*(13 - 16*s2 + t1) - 12*t2 - t1*t2 + 
                  s2*(17 - 14*t1 + 15*t2))) - 
            8*(-17 + Power(s,3) + 9*Power(s2,2) + 36*t1 - 10*s2*t1 - 
               10*Power(s2,2)*t1 - 18*Power(t1,2) + 
               10*s2*Power(t1,2) - 
               Power(s,2)*(-17 + 3*s1 + 18*s2 + t1) + 
               Power(s1,2)*(-10*s2 + 8*t1) + 20*t2 - 23*s2*t2 + 
               11*Power(s2,2)*t2 - 17*t1*t2 + 9*s2*t1*t2 + 
               8*Power(t2,2) - 10*s2*Power(t2,2) + 
               s*(69 + 2*Power(s1,2) + 17*Power(s2,2) - 58*t1 - 
                  8*Power(t1,2) + s2*(-92 + 72*t1 - 80*t2) + 
                  s1*(-63 + 91*s2 - 7*t1 - 2*t2) + 57*t2 + 5*t1*t2 + 
                  Power(t2,2)) + 
               s1*(-7 - 11*Power(s2,2) + 8*Power(t1,2) - 9*t2 - 
                  t1*(6 + 7*t2) + s2*(15 + t1 + 20*t2))) + 
            4*(-47 + 5*Power(s,3) + 6*s2 + 23*Power(s2,2) - 
               2*Power(s2,3) + 110*t1 - 31*s2*t1 - 31*Power(s2,2)*t1 - 
               53*Power(t1,2) + 31*s2*Power(t1,2) + 36*t2 - 56*s2*t2 + 
               39*Power(s2,2)*t2 - 43*t1*t2 + 23*s2*t1*t2 + 
               21*Power(t2,2) - 31*s2*Power(t2,2) - 
               Power(s1,2)*(35*s2 - 21*t1 + t2) + 
               Power(s,2)*(46 - 17*s1 - 58*s2 - 7*t1 + 4*t2) + 
               s*(135 + 12*Power(s1,2) + 55*Power(s2,2) - 129*t1 - 
                  21*Power(t1,2) + s2*(-200 + 161*t1 - 189*t2) + 
                  s1*(-134 + 228*s2 - 14*t1 - 11*t2) + 126*t2 + 
                  10*t1*t2 + 3*Power(t2,2)) + 
               s1*(-8 - 42*Power(s2,2) + 21*Power(t1,2) - 25*t2 + 
                  Power(t2,2) - 6*t1*(2 + 3*t2) + 
                  s2*(33 + 5*t1 + 65*t2))) + 
            2*(49 - Power(s1,3) + 21*s2 - 6*Power(s2,2) - 
               2*Power(s2,3) - 113*t1 + 27*s2*t1 + 32*Power(s2,2)*t1 + 
               48*Power(t1,2) - 30*s2*Power(t1,2) + 
               2*Power(s,2)*(-16 + 5*s1 + 21*s2 + 6*t1 - t2) - 22*t2 + 
               37*s2*t2 - 38*Power(s2,2)*t2 + 34*t1*t2 - 19*s2*t1*t2 - 
               19*Power(t2,2) + 30*s2*Power(t2,2) + 
               Power(s1,2)*(3 + 31*s2 - 17*t1 + 3*t2) + 
               s1*(7 + 36*Power(s2,2) - 18*Power(t1,2) + 20*t2 - 
                  2*Power(t2,2) + 3*t1*(3 + 5*t2) - 
                  s2*(11 + 3*t1 + 61*t2)) + 
               s*(-119 - 9*Power(s1,2) - 40*Power(s2,2) + 106*t1 + 
                  18*Power(t1,2) - 98*t2 - 6*t1*t2 - 2*Power(t2,2) + 
                  s1*(96 - 180*s2 + 5*t1 + 8*t2) + 
                  s2*(129 - 134*t1 + 152*t2))) + 
            ((2 + s - s1 - s2 + t2)*
               (4*(-1 + s1 + t1 - t2)*
                  Power(2*Power(s,2)*(-1 + s2) + 
                    (1 + s2)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2) - 
                    s*(-4 + 2*s1*s2 + 2*Power(s2,2) + t1 + t2 - 
                       s2*(t1 + t2)),2) - 
                 2*(6*Power(s,4)*Power(-1 + s2,2)*
                     (-1 + s1 + t1 - t2) + 
                    (-3 + s1*(1 + 4*s2 + Power(s2,2)) + 3*t1 + 
                       2*s2*(-1 + t1 - 2*t2) + 
                       Power(s2,2)*(-1 + t1 - t2) - t2)*
                     Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                    2*Power(s,3)*(-1 + s2)*
                     (12 + 6*Power(s1,2)*s2 - 15*t1 + 4*Power(t1,2) + 
                       6*Power(s2,2)*(-1 + t1 - t2) + 9*t2 - 
                       2*t1*t2 - 2*Power(t2,2) + 
                       3*s1*
                        (-4 + 2*Power(s2,2) + t1 + 
                        s2*(-2 + t1 - 3*t2) + t2) + 
                       3*s2*(t1 - Power(t1,2) + t2 + Power(t2,2))) + 
                    Power(s,2)*
                     (-30 + 6*Power(s1,3)*Power(s2,2) + 52*t1 - 
                       25*Power(t1,2) + 3*Power(t1,3) + 
                       6*Power(s2,4)*(-1 + t1 - t2) - 8*t2 + 
                       8*t1*t2 - Power(t1,2)*t2 + 11*Power(t2,2) - 
                       t1*Power(t2,2) - Power(t2,3) - 
                       6*Power(s2,3)*
                        (Power(t1,2) + t1*(-1 + t2) - 2*t2*(1 + t2)) \
+ 6*Power(s1,2)*(3*Power(s2,3) + t1 + s2*(-5 + t1 + t2) - 
                        Power(s2,2)*(1 + t1 + 2*t2)) - 
                       2*s2*
                        (2*Power(t1,3) + Power(t1,2)*(-10 + t2) + 
                        t1*(8 - 3*t2 - 2*Power(t2,2)) + 
                        t2*(15 + 7*t2 - Power(t2,2))) + 
                       Power(s2,2)*
                        (30 + Power(t1,3) + 26*t2 - 9*Power(t2,2) - 
                        Power(t2,3) + Power(t1,2)*(11 + 3*t2) - 
                        t1*(42 + 8*t2 + 3*Power(t2,2))) + 
                       s1*(30 + 6*Power(s2,4) + 3*Power(t1,2) + 
                        6*Power(s2,3)*(-3 + 2*t1 - 5*t2) + 
                        2*t1*(-16 + t2) - 18*t2 + Power(t2,2) + 
                        2*s2*
                        (15 + 4*Power(t1,2) + 23*t2 - 
                        4*Power(t2,2) - 6*t1*(1 + t2)) + 
                        Power(s2,2)*
                        (-30 - 11*Power(t1,2) + 14*t2 + 
                        7*Power(t2,2) + 10*t1*(2 + t2)))) - 
                    2*s*(-7 + 3*Power(s1,3)*s2*(1 + s2)*(s2 - t1) + 
                       17*t1 - 13*Power(t1,2) + 3*Power(t1,3) + 5*t2 - 
                       4*t1*t2 - Power(t1,2)*t2 + 3*Power(t2,2) - 
                       t1*Power(t2,2) - Power(t2,3) + 
                       3*Power(s2,4)*t2*(1 - t1 + t2) + 
                       Power(s2,2)*
                        (3 - 2*Power(t1,3) - 4*t2 - 13*Power(t2,2) + 
                        Power(t2,3) + Power(t1,2)*(7 + t2) + 
                        t1*(-8 + 3*t2)) + 
                       Power(s2,3)*
                        (3 + 4*t2 + Power(t2,2) - Power(t2,3) + 
                        Power(t1,2)*(3 + 2*t2) - 
                        t1*(6 + 6*t2 + Power(t2,2))) + 
                       s2*(-5 - Power(t1,3) - 14*t2 + 6*Power(t2,2) + 
                        Power(t2,3) - Power(t1,2)*(3 + 2*t2) + 
                        t1*(9 + 16*t2 + 2*Power(t2,2))) + 
                       Power(s1,2)*
                        (3*Power(s2,4) - t1*(-5 + 3*t1 + t2) - 
                        Power(s2,3)*(-1 + t1 + 7*t2) + 
                        s2*
                       (-9 - Power(t1,2) + 5*t2 + 2*t1*(9 + t2)) + 
                        Power(s2,2)*
                        (-13 - 2*Power(t1,2) - 4*t2 + t1*(2 + 5*t2))) \
+ s1*(6 + 3*Power(s2,4)*(-1 + t1 - 2*t2) - Power(t1,2)*(-12 + t2) - 
                        7*t2 + Power(t2,2) + 
                        t1*(-18 + 4*t2 + Power(t2,2)) + 
                        Power(s2,3)*
                        (-6 - 4*Power(t1,2) - 2*t2 + 5*Power(t2,2) + 
                        2*t1*(5 + t2)) + 
                        Power(s2,2)*
                        (4 + Power(t1,3) + Power(t1,2)*(-3 + t2) + 
                        25*t2 - 2*t1*(1 + Power(t2,2))) + 
                        s2*(17 + 7*Power(t1,2) - Power(t1,3) + 2*t2 - 
                        6*Power(t2,2) + t1*(-23 - 18*t2 + Power(t2,2)))\
)))))/((-1 + s1 + t1 - t2)*Power(-1 + s - s*s2 + s1*s2 + t1 - s2*t2,2)))\
)/((-1 + t1)*(-1 + t2)) + (2*
          (32*(-1 + s2)*(-1 + s1 + t1 - t2)*(-1 + s + t2) - 
            ((-11 + 6*Power(s,2) + 4*Power(s1,2) - 14*s2 + 
                 6*Power(s2,2) + t1 + s1*(-6 + 10*s2 - 5*t2) + 
                 4*t2 - 6*s2*t2 + Power(t2,2) + 
                 s*(11 - 10*s1 - 12*s2 + 7*t2))*
               (2*Power(s,2)*(-1 + s2) + 
                 (1 + s2)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2) - 
                 s*(-4 + 2*s1*s2 + 2*Power(s2,2) + t1 + t2 - 
                    s2*(t1 + t2))))/
             (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2) - 
            2*(33 + 38*Power(s,3) - 6*Power(s1,3) - 67*s2 - 
               49*Power(s2,2) - 4*Power(s2,3) - 117*t1 + 79*s2*t1 + 
               64*Power(s2,2)*t1 + 61*Power(t1,2) - 
               40*s2*Power(t1,2) - 
               2*Power(s,2)*(50 + 15*s1 + 18*s2 - 31*t1 - t2) + 
               58*t2 + 89*s2*t2 - 64*Power(s2,2)*t2 + 2*t1*t2 - 
               72*s2*t1*t2 - 81*Power(t2,2) + 99*s2*Power(t2,2) + 
               Power(s1,2)*(-30 + 50*s2 - 26*t1 + 19*t2) + 
               s*(-45 - 2*Power(s1,2) + 2*Power(s2,2) + 92*t1 + 
                  36*Power(t1,2) - 6*s2*(-40 + 36*t1 - 29*t2) - 
                  203*t2 - 14*t1*t2 + 4*Power(t2,2) + 
                  s1*(213 - 154*s2 - 36*t1 + t2)) + 
               s1*(-48 + 50*Power(s2,2) - 36*Power(t1,2) + 
                  s2*(-65 + 36*t1 - 145*t2) + 95*t2 - 13*Power(t2,2) + 
                  4*t1*(13 + 9*t2))) + 
            16*(9 + 2*Power(s,3) - 8*s2 - 2*Power(s2,2) - 14*t1 + 
               10*s2*t1 + 3*Power(s2,2)*t1 + 4*Power(t1,2) - 
               2*s2*Power(t1,2) + t2 + 5*s2*t2 - 3*Power(s2,2)*t2 + 
               8*t1*t2 - 12*s2*t1*t2 - 12*Power(t2,2) + 
               13*s2*Power(t2,2) + Power(s1,2)*(-2 + 2*s2 - t1 + t2) + 
               Power(s,2)*(-2 - 2*s1 - 2*s2 + 3*t1 + t2) + 
               s1*(-8 + 2*Power(s2,2) - 2*Power(t1,2) + 
                  s2*(5 + 2*t1 - 15*t2) + 13*t2 - Power(t2,2) + 
                  2*t1*(1 + t2)) + 
               s*(-15 + 16*t1 + 2*Power(t1,2) - 18*t2 - t1*t2 + 
                  Power(t2,2) - s1*(-18 + 15*s2 + 2*t1 + t2) + 
                  s2*(18 - 19*t1 + 16*t2))) - 
            8*(35 + 17*Power(s,3) - Power(s1,3) - 27*s2 - 
               20*Power(s2,2) - 72*t1 + 44*s2*t1 + 28*Power(s2,2)*t1 + 
               30*Power(t1,2) - 18*s2*Power(t1,2) + 12*t2 + 40*s2*t2 - 
               28*Power(s2,2)*t2 + 22*t1*t2 - 51*s2*t1*t2 - 
               56*Power(t2,2) + 62*s2*Power(t2,2) + 
               Power(s,2)*(-22 - 17*s1 - 18*s2 + 25*t1 + 6*t2) + 
               Power(s1,2)*(-15 + 19*s2 - 9*t1 + 8*t2) + 
               s1*(-27 + 21*Power(s2,2) - 16*Power(t1,2) + 
                  s2*(-7 + 16*t1 - 80*t2) + 63*t2 - 7*Power(t2,2) + 
                  2*t1*(9 + 8*t2)) + 
               s*(-71 + Power(s1,2) + Power(s2,2) + 79*t1 + 
                  16*Power(t1,2) - 106*t2 - 8*t1*t2 + 6*Power(t2,2) - 
                  s1*(-104 + 81*s2 + 16*t1 + 6*t2) + 
                  s2*(108 - 114*t1 + 91*t2))) + 
            4*(57 + 48*Power(s,3) - 6*Power(s1,3) - 50*s2 - 
               48*Power(s2,2) - 5*Power(s2,3) - 154*t1 + 96*s2*t1 + 
               76*Power(s2,2)*t1 + 74*Power(t1,2) - 
               48*s2*Power(t1,2) + 46*t2 + 98*s2*t2 - 
               73*Power(s2,2)*t2 + 21*t1*t2 - 98*s2*t1*t2 - 
               113*Power(t2,2) + 129*s2*Power(t2,2) + 
               Power(s,2)*(-81 - 48*s1 - 55*s2 + 69*t1 + 13*t2) + 
               Power(s1,2)*(-34 + 50*s2 - 27*t1 + 23*t2) + 
               s1*(-48 + 52*Power(s2,2) - 42*Power(t1,2) + 
                  2*s2*(-24 + 21*t1 - 88*t2) + 127*t2 - 
                  17*Power(t2,2) + 6*t1*(9 + 7*t2)) + 
               s*(-124 + 6*Power(s1,2) + 12*Power(s2,2) + 148*t1 + 
                  42*Power(t1,2) - 245*t2 - 19*t1*t2 + 
                  11*Power(t2,2) - 
                  s1*(-243 + 173*s2 + 42*t1 + 13*t2) + 
                  s2*(260 - 268*t1 + 206*t2))) + 
            ((2 + s - s1 - s2 + t2)*
               (4*(-1 + s1 + t1 - t2)*
                  Power(2*Power(s,2)*(-1 + s2) + 
                    (1 + s2)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2) - 
                    s*(-4 + 2*s1*s2 + 2*Power(s2,2) + t1 + t2 - 
                       s2*(t1 + t2)),2) - 
                 2*(6*Power(s,4)*Power(-1 + s2,2)*
                     (-1 + s1 + t1 - t2) + 
                    (-3 + s1*(1 + 4*s2 + Power(s2,2)) + 3*t1 + 
                       2*s2*(-1 + t1 - 2*t2) + 
                       Power(s2,2)*(-1 + t1 - t2) - t2)*
                     Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                    2*Power(s,3)*(-1 + s2)*
                     (12 + 6*Power(s1,2)*s2 - 15*t1 + 4*Power(t1,2) + 
                       6*Power(s2,2)*(-1 + t1 - t2) + 9*t2 - 
                       2*t1*t2 - 2*Power(t2,2) + 
                       3*s1*
                        (-4 + 2*Power(s2,2) + t1 + 
                        s2*(-2 + t1 - 3*t2) + t2) + 
                       3*s2*(t1 - Power(t1,2) + t2 + Power(t2,2))) + 
                    Power(s,2)*
                     (-30 + 6*Power(s1,3)*Power(s2,2) + 52*t1 - 
                       25*Power(t1,2) + 3*Power(t1,3) + 
                       6*Power(s2,4)*(-1 + t1 - t2) - 8*t2 + 
                       8*t1*t2 - Power(t1,2)*t2 + 11*Power(t2,2) - 
                       t1*Power(t2,2) - Power(t2,3) - 
                       6*Power(s2,3)*
                        (Power(t1,2) + t1*(-1 + t2) - 2*t2*(1 + t2)) \
+ 6*Power(s1,2)*(3*Power(s2,3) + t1 + s2*(-5 + t1 + t2) - 
                        Power(s2,2)*(1 + t1 + 2*t2)) - 
                       2*s2*
                        (2*Power(t1,3) + Power(t1,2)*(-10 + t2) + 
                        t1*(8 - 3*t2 - 2*Power(t2,2)) + 
                        t2*(15 + 7*t2 - Power(t2,2))) + 
                       Power(s2,2)*
                        (30 + Power(t1,3) + 26*t2 - 9*Power(t2,2) - 
                        Power(t2,3) + Power(t1,2)*(11 + 3*t2) - 
                        t1*(42 + 8*t2 + 3*Power(t2,2))) + 
                       s1*(30 + 6*Power(s2,4) + 3*Power(t1,2) + 
                        6*Power(s2,3)*(-3 + 2*t1 - 5*t2) + 
                        2*t1*(-16 + t2) - 18*t2 + Power(t2,2) + 
                        2*s2*
                        (15 + 4*Power(t1,2) + 23*t2 - 
                        4*Power(t2,2) - 6*t1*(1 + t2)) + 
                        Power(s2,2)*
                        (-30 - 11*Power(t1,2) + 14*t2 + 
                        7*Power(t2,2) + 10*t1*(2 + t2)))) - 
                    2*s*(-7 + 3*Power(s1,3)*s2*(1 + s2)*(s2 - t1) + 
                       17*t1 - 13*Power(t1,2) + 3*Power(t1,3) + 5*t2 - 
                       4*t1*t2 - Power(t1,2)*t2 + 3*Power(t2,2) - 
                       t1*Power(t2,2) - Power(t2,3) + 
                       3*Power(s2,4)*t2*(1 - t1 + t2) + 
                       Power(s2,2)*
                        (3 - 2*Power(t1,3) - 4*t2 - 13*Power(t2,2) + 
                        Power(t2,3) + Power(t1,2)*(7 + t2) + 
                        t1*(-8 + 3*t2)) + 
                       Power(s2,3)*
                        (3 + 4*t2 + Power(t2,2) - Power(t2,3) + 
                        Power(t1,2)*(3 + 2*t2) - 
                        t1*(6 + 6*t2 + Power(t2,2))) + 
                       s2*(-5 - Power(t1,3) - 14*t2 + 6*Power(t2,2) + 
                        Power(t2,3) - Power(t1,2)*(3 + 2*t2) + 
                        t1*(9 + 16*t2 + 2*Power(t2,2))) + 
                       Power(s1,2)*
                        (3*Power(s2,4) - t1*(-5 + 3*t1 + t2) - 
                        Power(s2,3)*(-1 + t1 + 7*t2) + 
                        s2*
                       (-9 - Power(t1,2) + 5*t2 + 2*t1*(9 + t2)) + 
                        Power(s2,2)*
                        (-13 - 2*Power(t1,2) - 4*t2 + t1*(2 + 5*t2))) \
+ s1*(6 + 3*Power(s2,4)*(-1 + t1 - 2*t2) - Power(t1,2)*(-12 + t2) - 
                        7*t2 + Power(t2,2) + 
                        t1*(-18 + 4*t2 + Power(t2,2)) + 
                        Power(s2,3)*
                        (-6 - 4*Power(t1,2) - 2*t2 + 5*Power(t2,2) + 
                        2*t1*(5 + t2)) + 
                        Power(s2,2)*
                        (4 + Power(t1,3) + Power(t1,2)*(-3 + t2) + 
                        25*t2 - 2*t1*(1 + Power(t2,2))) + 
                        s2*(17 + 7*Power(t1,2) - Power(t1,3) + 2*t2 - 
                        6*Power(t2,2) + t1*(-23 - 18*t2 + Power(t2,2)))\
)))))/((-1 + s1 + t1 - t2)*Power(-1 + s - s*s2 + s1*s2 + t1 - s2*t2,2)))\
)/((-1 + s1)*(-1 + t2)) + ((2*
             (((-3 - Power(s1,2) + Power(s,2)*s2 + 6*Power(s2,2) + 
                    Power(s2,3) - 4*t1 - 7*s2*t1 - Power(s2,2)*t1 - 
                    s1*(4 + s*(-1 + s2) + s2 - Power(s2,2) + t1 + 
                       s2*t1 - 2*t2) + 4*t2 + s2*t2 - 
                    Power(s2,2)*t2 + t1*t2 + s2*t1*t2 - 
                    Power(t2,2) + 
                    s*(-2*Power(s2,2) + t1 - t2 + s2*(-6 + t1 + t2))\
)*(2*Power(s,2)*(-1 + s2) + 
                    (1 + s2)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2) - 
                    s*(-4 + 2*s1*s2 + 2*Power(s2,2) + t1 + t2 - 
                       s2*(t1 + t2))))/
                (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2) - 
               2*(1 + Power(s,3)*(-2 + s2) + 13*s2 + 9*Power(s2,2) - 
                  2*Power(s2,3) - Power(s2,4) - 
                  Power(s1,2)*(-2 + s2)*(2 + s2 - t1) + 14*t1 + 
                  7*s2*t1 + 4*Power(s2,2)*t1 + Power(s2,3)*t1 - 
                  7*Power(t1,2) - 3*s2*Power(t1,2) - 13*t2 - 
                  17*s2*t2 - 6*Power(s2,2)*t2 + Power(s2,3)*t2 + 
                  3*t1*t2 + 5*s2*t1*t2 - Power(s2,2)*t1*t2 + 
                  2*Power(t1,2)*t2 + 4*Power(t2,2) + 
                  2*s2*Power(t2,2) - 2*t1*Power(t2,2) + 
                  Power(s,2)*
                   (7 - 2*s1*(-1 + s2) - 3*Power(s2,2) - 4*t2 + 
                     s2*(2 + t1 + t2)) - 
                  s1*(-13 + 2*Power(s2,3) + 4*Power(t1,2) + 
                     t1*(2 - 4*t2) + 8*t2 - 
                     Power(s2,2)*(2 + 2*t1 + t2) + 
                     s2*(-14 + t1 + 2*t2 + t1*t2)) + 
                  s*(-5 + Power(s1,2)*(-2 + s2) + 3*Power(s2,3) - 
                     8*t1 + 2*Power(t1,2) + 15*t2 - 2*t1*t2 - 
                     2*Power(t2,2) - 2*Power(s2,2)*(-1 + t1 + t2) + 
                     s2*(-13 + t1*(-3 + t2) + 9*t2) + 
                     s1*(4*Power(s2,2) - 2*(7 + t1 - 2*t2) - 
                        s2*(3 + 2*t1 + t2))))))/
           ((-1 + s2)*(-s + s2 - t1)) + 
          (2*(15 + 4*Power(s,4)*(-1 + s2) - 17*s2 + 4*Power(s1,4)*s2 - 
               24*Power(s2,2) - 3*Power(s2,3) + Power(s2,4) - 53*t1 + 
               7*s2*t1 + 26*Power(s2,2)*t1 + Power(s2,3)*t1 - 
               Power(s2,4)*t1 + 45*Power(t1,2) + 10*s2*Power(t1,2) - 
               Power(s2,2)*Power(t1,2) + 2*Power(s2,3)*Power(t1,2) - 
               3*Power(t1,3) - Power(s2,2)*Power(t1,3) - 
               4*Power(t1,4) + 5*t2 + 58*s2*t2 - Power(s2,2)*t2 - 
               17*Power(s2,3)*t2 - 2*Power(s2,4)*t2 + Power(s2,5)*t2 + 
               2*t1*t2 - 95*s2*t1*t2 - 27*Power(s2,2)*t1*t2 - 
               2*Power(s2,4)*t1*t2 - 16*Power(t1,2)*t2 + 
               7*s2*Power(t1,2)*t2 + Power(s2,3)*Power(t1,2)*t2 + 
               Power(t1,3)*t2 + 5*s2*Power(t1,3)*t2 - 5*Power(t2,2) - 
               2*s2*Power(t2,2) + 39*Power(s2,2)*Power(t2,2) + 
               8*Power(s2,3)*Power(t2,2) + 12*t1*Power(t2,2) + 
               34*s2*t1*Power(t2,2) + 3*Power(s2,2)*t1*Power(t2,2) + 
               Power(s2,3)*t1*Power(t2,2) + 
               2*Power(t1,2)*Power(t2,2) - 
               3*s2*Power(t1,2)*Power(t2,2) - 
               Power(s2,2)*Power(t1,2)*Power(t2,2) + Power(t2,3) - 
               15*s2*Power(t2,3) - 13*Power(s2,2)*Power(t2,3) - 
               Power(s2,3)*Power(t2,3) - 3*t1*Power(t2,3) - 
               2*s2*t1*Power(t2,3) + Power(s2,2)*t1*Power(t2,3) + 
               4*s2*Power(t2,4) - 
               Power(s1,3)*(4 + Power(s2,3) - 7*t1 - 
                  2*Power(s2,2)*t1 + Power(t1,2) + 
                  s2*(-13 - 9*t1 + Power(t1,2) + 16*t2)) + 
               Power(s,3)*(14 - 10*t1 - Power(t1,2) + 
                  Power(s2,2)*(-18 + t1 - t2) - 6*t2 + 2*t1*t2 - 
                  Power(t2,2) - 
                  s2*(-8 - 11*t1 + Power(t1,2) - 5*t2 - 2*t1*t2 + 
                     Power(t2,2)) + 
                  s1*(6 + 2*Power(s2,2) - t1 + t2 + 
                     s2*(-12 - t1 + t2))) + 
               Power(s1,2)*(-13 - 2*Power(s2,4) - Power(t1,3) + 
                  t1*(14 - 17*t2) + 9*t2 + 
                  Power(s2,3)*(-2 + 3*t1 + t2) + 
                  Power(t1,2)*(9 + 2*t2) + 
                  s2*(10 - Power(t1,3) + t1*(11 - 20*t2) - 41*t2 + 
                     24*Power(t2,2) + Power(t1,2)*(9 + 2*t2)) - 
                  Power(s2,2)*(13*(-3 + t2) + t1*(4 + 3*t2))) + 
               Power(s,2)*(8 + t1 - 13*Power(t1,2) - Power(t1,3) + 
                  21*t2 - 6*t1*t2 + Power(t1,2)*t2 - 5*Power(t2,2) + 
                  t1*Power(t2,2) - Power(t2,3) + 
                  Power(s2,3)*(-2*t1 + 3*(8 + t2)) - 
                  Power(s1,2)*
                   (6 + 3*Power(s2,2) - 2*t1 + t2 + 
                     s2*(-15 - 2*t1 + t2)) + 
                  Power(s2,2)*
                   (25 + 3*Power(t1,2) - 18*t2 + 2*Power(t2,2) - 
                     t1*(29 + 5*t2)) + 
                  s2*(-79 - Power(t1,3) + 12*t2 + 17*Power(t2,2) - 
                     Power(t2,3) + Power(t1,2)*(22 + t2) + 
                     t1*(30 - 15*t2 + Power(t2,2))) + 
                  s1*(-29 - 5*Power(s2,3) + 5*t1 + Power(t1,2) + 
                     11*t2 - 3*t1*t2 + 2*Power(t2,2) + 
                     Power(s2,2)*(25 + 3*t1 + t2) + 
                     s2*(25 - 8*t1 + Power(t1,2) - 32*t2 - 3*t1*t2 + 
                        2*Power(t2,2)))) + 
               s1*(-Power(s2,5) + Power(t1,3)*(2 + t2) + 
                  Power(s2,4)*(3 + t1 + 2*t2) - 
                  Power(t1,2)*(-13 + 11*t2 + Power(t2,2)) - 
                  2*(5 - 9*t2 + 3*Power(t2,2)) + 
                  t1*(3 - 26*t2 + 13*Power(t2,2)) + 
                  Power(s2,2)*
                   (19 - Power(t1,3) + t1*(-3 + t2) - 78*t2 + 
                     26*Power(t2,2) + Power(t1,2)*(13 + t2)) + 
                  Power(s2,3)*
                   (26 + Power(t1,2) - 6*t2 + Power(t2,2) - 
                     t1*(11 + 4*t2)) + 
                  s2*(-53 + Power(t1,3)*(-7 + t2) - 8*t2 + 
                     43*Power(t2,2) - 16*Power(t2,3) - 
                     Power(t1,2)*(-9 + 6*t2 + Power(t2,2)) + 
                     t1*(76 - 45*t2 + 13*Power(t2,2)))) + 
               s*(-34 + 65*t1 - 21*Power(t1,2) - 14*Power(t1,3) + 
                  Power(s1,3)*(4 + Power(s2,2) - t1 - s2*(11 + t1)) + 
                  Power(s2,4)*(-10 + t1 - 3*t2) - 15*t2 + 14*t1*t2 + 
                  10*Power(t1,2)*t2 - Power(t1,3)*t2 + 9*Power(t2,2) - 
                  12*t1*Power(t2,2) + 2*Power(t1,2)*Power(t2,2) - 
                  t1*Power(t2,3) - 
                  Power(s2,3)*
                   (24 + 2*Power(t1,2) - 13*t2 + Power(t2,2) - 
                     t1*(21 + 5*t2)) + 
                  Power(s1,2)*
                   (17 + 5*Power(s2,3) + Power(t1,2) - 
                     Power(s2,2)*(7 + 6*t1) + t1*(-6 + t2) - 8*t2 + 
                     s2*(-39 + Power(t1,2) + t1*(-12 + t2) + 40*t2)) + 
                  Power(s2,2)*
                   (26 + Power(t1,3) + 43*t2 - 26*Power(t2,2) + 
                     2*Power(t2,3) - Power(t1,2)*(15 + t2) + 
                     t1*(-7 + t2 - 2*Power(t2,2))) - 
                  s2*(-70 + Power(t1,3)*(-7 + t2) + 88*t2 + 
                     20*Power(t2,2) - 18*Power(t2,3) + 
                     Power(t1,2)*(2 - 13*t2 - 2*Power(t2,2)) + 
                     t1*(78 - 82*t2 + 22*Power(t2,2) + Power(t2,3))) + 
                  s1*(28 + 4*Power(s2,4) + 2*Power(t1,3) - 26*t2 + 
                     4*Power(t2,2) - Power(t1,2)*(8 + 3*t2) - 
                     Power(s2,3)*(14 + 3*t1 + 4*t2) + 
                     t1*(-22 + 18*t2 + Power(t2,2)) + 
                     Power(s2,2)*
                      (-73 - 3*Power(t1,2) + 33*t2 - 3*Power(t2,2) + 
                        8*t1*(4 + t2)) + 
                     s2*(57 + 2*Power(t1,3) + 59*t2 - 47*Power(t2,2) - 
                        Power(t1,2)*(31 + 3*t2) + 
                        t1*(-27 + 34*t2 + Power(t2,2)))))))/
           ((s - s2 + t1)*(s - s1 + t2)*
             (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2)) + 
          ((-2*(2*Power(s,2)*(-1 + s2) + 
                  (1 + s2)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2) - 
                  s*(-4 + 2*s1*s2 + 2*Power(s2,2) + t1 + t2 - 
                     s2*(t1 + t2)))*
                (15 + 2*Power(s,3) + 17*s2 + 4*Power(s2,2) - 
                  3*Power(s2,3) - 9*t1 - 10*s2*t1 + 
                  3*Power(s2,2)*t1 - Power(t1,2) + 
                  Power(s1,2)*(-4 - s2 + t1) - 3*t2 + 5*s2*t2 + 
                  2*Power(s2,2)*t2 + 8*t1*t2 - 3*s2*t1*t2 - 
                  4*Power(t2,2) + s2*Power(t2,2) - 
                  s1*(-2 + 4*Power(s2,2) + s2*(2 - 4*t1) - 8*t2 + 
                     t1*(9 + t2)) + 
                  Power(s,2)*(-3*s1 - 7*s2 + 2*(1 + t1 + t2)) + 
                  s*(Power(s1,2) + 8*Power(s2,2) + 
                     s1*(1 + 7*s2 - 3*t1 - t2) + 
                     (-1 + t1)*(7 + 2*t2) - s2*(8 + 5*t1 + 4*t2))))/
              (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2) + 
             4*(-4 + Power(s,4) - 25*s2 - 17*Power(s2,2) - 
                2*Power(s2,3) + Power(s2,4) + 
                Power(s1,2)*s2*(4 + s2 - t1) + 7*t1 + 15*s2*t1 + 
                6*Power(s2,2)*t1 - Power(s2,3)*t1 - 3*Power(t1,2) - 
                s2*Power(t1,2) + t2 + 10*s2*t2 - 2*Power(s2,2)*t2 - 
                Power(s2,3)*t2 - 10*t1*t2 - 5*s2*t1*t2 + 
                Power(s2,2)*t1*t2 + 2*Power(t1,2)*t2 + 
                7*Power(t2,2) + 3*s2*Power(t2,2) - 2*t1*Power(t2,2) + 
                Power(s,3)*(2 - 2*s1 - 4*s2 + t1 + t2) + 
                Power(s,2)*(-9 + Power(s1,2) + 6*Power(s2,2) + 
                   6*t1 + s1*(-1 + 6*s2 - 2*t1 - t2) - t2 + t1*t2 - 
                   s2*(8 + 3*t1 + 3*t2)) + 
                s1*(-5 + 2*Power(s2,3) - 2*Power(t1,2) - 6*t2 - 
                   Power(s2,2)*(2*t1 + t2) + t1*(5 + 3*t2) + 
                   s2*(-9*(1 + t2) + t1*(7 + t2))) - 
                s*(-20 + 4*Power(s2,3) + 
                   Power(s1,2)*(3 + 2*s2 - t1) + 10*t1 + 
                   Power(t1,2) + 9*t2 - 9*t1*t2 + 5*Power(t2,2) - 
                   Power(s2,2)*(8 + 3*t1 + 3*t2) + 
                   s1*(6*Power(s2,2) - 9*(1 + t2) + t1*(11 + t2) - 
                      2*s2*(2*t1 + t2)) + 
                   s2*(-28 - 3*t2 + t1*(13 + 2*t2)))) + 
             ((1 + Power(s,2) + s1 - 3*s2 + 2*s1*s2 + 
                  2*Power(s2,2) + t1 - s1*t1 - s2*t1 - t2 - 
                  2*s2*t2 + t1*t2 + s*(-s1 - 3*s2 + t1 + t2))*
                (4*(-1 + s1 + t1 - t2)*
                   Power(2*Power(s,2)*(-1 + s2) + 
                     (1 + s2)*
                      (-1 + s1*(s2 - t1) + t1 + t2 - s2*t2) - 
                     s*(-4 + 2*s1*s2 + 2*Power(s2,2) + t1 + t2 - 
                       s2*(t1 + t2)),2) - 
                  2*(6*Power(s,4)*Power(-1 + s2,2)*
                      (-1 + s1 + t1 - t2) + 
                     (-3 + s1*(1 + 4*s2 + Power(s2,2)) + 3*t1 + 
                       2*s2*(-1 + t1 - 2*t2) + 
                       Power(s2,2)*(-1 + t1 - t2) - t2)*
                      Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                     2*Power(s,3)*(-1 + s2)*
                      (12 + 6*Power(s1,2)*s2 - 15*t1 + 
                        4*Power(t1,2) + 
                        6*Power(s2,2)*(-1 + t1 - t2) + 9*t2 - 
                        2*t1*t2 - 2*Power(t2,2) + 
                        3*s1*
                       (-4 + 2*Power(s2,2) + t1 + 
                       s2*(-2 + t1 - 3*t2) + t2) + 
                        3*s2*(t1 - Power(t1,2) + t2 + Power(t2,2))) + 
                     Power(s,2)*
                      (-30 + 6*Power(s1,3)*Power(s2,2) + 52*t1 - 
                        25*Power(t1,2) + 3*Power(t1,3) + 
                        6*Power(s2,4)*(-1 + t1 - t2) - 8*t2 + 
                        8*t1*t2 - Power(t1,2)*t2 + 11*Power(t2,2) - 
                        t1*Power(t2,2) - Power(t2,3) - 
                        6*Power(s2,3)*
                       (Power(t1,2) + t1*(-1 + t2) - 2*t2*(1 + t2)) \
+ 6*Power(s1,2)*(3*Power(s2,3) + t1 + s2*(-5 + t1 + t2) - 
                       Power(s2,2)*(1 + t1 + 2*t2)) - 
                        2*s2*
                       (2*Power(t1,3) + Power(t1,2)*(-10 + t2) + 
                       t1*(8 - 3*t2 - 2*Power(t2,2)) + 
                       t2*(15 + 7*t2 - Power(t2,2))) + 
                        Power(s2,2)*
                       (30 + Power(t1,3) + 26*t2 - 9*Power(t2,2) - 
                       Power(t2,3) + Power(t1,2)*(11 + 3*t2) - 
                       t1*(42 + 8*t2 + 3*Power(t2,2))) + 
                        s1*(30 + 6*Power(s2,4) + 3*Power(t1,2) + 
                        6*Power(s2,3)*(-3 + 2*t1 - 5*t2) + 
                        2*t1*(-16 + t2) - 18*t2 + Power(t2,2) + 
                        2*s2*
                       (15 + 4*Power(t1,2) + 23*t2 - 
                       4*Power(t2,2) - 6*t1*(1 + t2)) + 
                        Power(s2,2)*
                        (-30 - 11*Power(t1,2) + 14*t2 + 
                        7*Power(t2,2) + 10*t1*(2 + t2)))) - 
                     2*s*(-7 + 3*Power(s1,3)*s2*(1 + s2)*(s2 - t1) + 
                        17*t1 - 13*Power(t1,2) + 3*Power(t1,3) + 
                        5*t2 - 4*t1*t2 - Power(t1,2)*t2 + 
                        3*Power(t2,2) - t1*Power(t2,2) - 
                        Power(t2,3) + 
                        3*Power(s2,4)*t2*(1 - t1 + t2) + 
                        Power(s2,2)*
                        (3 - 2*Power(t1,3) - 4*t2 - 13*Power(t2,2) + 
                        Power(t2,3) + Power(t1,2)*(7 + t2) + 
                        t1*(-8 + 3*t2)) + 
                        Power(s2,3)*
                        (3 + 4*t2 + Power(t2,2) - Power(t2,3) + 
                        Power(t1,2)*(3 + 2*t2) - 
                        t1*(6 + 6*t2 + Power(t2,2))) + 
                        s2*(-5 - Power(t1,3) - 14*t2 + 
                        6*Power(t2,2) + Power(t2,3) - 
                        Power(t1,2)*(3 + 2*t2) + 
                        t1*(9 + 16*t2 + 2*Power(t2,2))) + 
                        Power(s1,2)*
                        (3*Power(s2,4) - t1*(-5 + 3*t1 + t2) - 
                        Power(s2,3)*(-1 + t1 + 7*t2) + 
                        s2*
                       (-9 - Power(t1,2) + 5*t2 + 2*t1*(9 + t2)) + 
                        Power(s2,2)*
                        (-13 - 2*Power(t1,2) - 4*t2 + t1*(2 + 5*t2))) \
+ s1*(6 + 3*Power(s2,4)*(-1 + t1 - 2*t2) - Power(t1,2)*(-12 + t2) - 
                        7*t2 + Power(t2,2) + 
                        t1*(-18 + 4*t2 + Power(t2,2)) + 
                        Power(s2,3)*
                        (-6 - 4*Power(t1,2) - 2*t2 + 5*Power(t2,2) + 
                        2*t1*(5 + t2)) + 
                        Power(s2,2)*
                        (4 + Power(t1,3) + Power(t1,2)*(-3 + t2) + 
                        25*t2 - 2*t1*(1 + Power(t2,2))) + 
                        s2*(17 + 7*Power(t1,2) - Power(t1,3) + 2*t2 - 
                        6*Power(t2,2) + t1*(-23 - 18*t2 + Power(t2,2))\
))))))/((-1 + s1 + t1 - t2)*Power(-1 + s - s*s2 + s1*s2 + t1 - s2*t2,
                 2)))/((-1 + s2)*(-1 + t1)) + 
          ((2*(2*Power(s,2)*(-1 + s2) + 
                  (1 + s2)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2) - 
                  s*(-4 + 2*s1*s2 + 2*Power(s2,2) + t1 + t2 - 
                     s2*(t1 + t2)))*
                (-12 + 2*Power(s,3) + Power(s1,3) - 17*s2 + 
                  4*Power(s2,2) - 2*Power(s2,3) - 3*t1 - 3*s2*t1 + 
                  5*Power(s2,2)*t1 - Power(t1,2) - s2*Power(t1,2) + 
                  Power(s1,2)*(-1 + 3*s2 + t1 - 2*t2) + 12*t2 + 
                  5*s2*t2 - 2*Power(s2,2)*t2 - 4*s2*t1*t2 + 
                  Power(t1,2)*t2 - 4*Power(t2,2) + 
                  4*s2*Power(t2,2) - t1*Power(t2,2) - 
                  Power(s,2)*(-4 + 6*s2 - 5*t1 + t2) + 
                  s1*(-10 - 3*t1 - Power(t1,2) + 
                     s2*(2 + 6*t1 - 7*t2) + 5*t2 + Power(t2,2)) + 
                  s*(11 - 3*Power(s1,2) + 6*Power(s2,2) + 3*t1 + 
                     Power(t1,2) - 4*t2 + 4*t1*t2 - 3*Power(t2,2) + 
                     s2*(-8 - 10*t1 + 3*t2) + s1*(-2 - 6*t1 + 6*t2)))\
)/(1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2) - 
             4*(5 + Power(s,4) - Power(s1,3)*(-1 + s2) + 17*s2 + 
                12*Power(s2,2) - 2*Power(s2,3) + Power(s2,4) - 4*t1 + 
                2*s2*t1 + Power(s2,2)*t1 - 3*Power(s2,3)*t1 + 
                3*Power(t1,2) + 5*s2*Power(t1,2) + 
                Power(s2,2)*Power(t1,2) - 2*Power(t1,3) - 8*t2 - 
                16*s2*t2 - 8*Power(s2,2)*t2 + 2*Power(s2,3)*t2 - 
                t1*t2 + s2*t1*t2 + 2*Power(s2,2)*t1*t2 - 
                Power(t1,2)*t2 - s2*Power(t1,2)*t2 + Power(t2,2) + 
                4*s2*Power(t2,2) - 3*Power(s2,2)*Power(t2,2) - 
                t1*Power(t2,2) + s2*t1*Power(t2,2) + 2*Power(t2,3) - 
                Power(s,3)*(-1 + 4*s2 - 3*t1 + t2) - 
                Power(s,2)*(-7 + 2*Power(s1,2) - 6*Power(s2,2) + 
                   t1 - Power(t1,2) + s1*(1 + 4*t1 - 4*t2) + 
                   s2*(4 + 9*t1 - 4*t2) + 5*t2 - 2*t1*t2 + 
                   2*Power(t2,2)) + 
                Power(s1,2)*
                 (4 - 2*Power(s2,2) + t1 - t2 + s2*(2 - t1 + 3*t2)) + 
                s1*(6 - 5*t2 - 2*Power(t2,2) + 2*t1*(2 + t2) + 
                   Power(s2,2)*(-2 - 4*t1 + 6*t2) + 
                   s2*(15 + 6*t1 + Power(t1,2) - 9*t2 - 
                      2*Power(t2,2))) + 
                s*(-16 + Power(s1,3) - 4*Power(s2,3) - t1 - 
                   5*Power(t1,2) + Power(s2,2)*(5 + 9*t1 - 5*t2) + 
                   Power(s1,2)*(-2 + 4*s2 + t1 - 2*t2) + 14*t2 - 
                   t1*t2 + Power(t1,2)*t2 - 3*Power(t2,2) - 
                   t1*Power(t2,2) - 
                   s2*(19 + t1 + 2*Power(t1,2) - 13*t2 + 4*t1*t2 - 
                      5*Power(t2,2)) + 
                   s1*(-11 - 4*t1 - Power(t1,2) + 
                      2*s2*(1 + 4*t1 - 5*t2) + 8*t2 + Power(t2,2)))) - 
             ((5 + Power(s,2) - Power(s1,2) - 2*s2 + Power(s2,2) + 
                  2*t1 - 2*s2*t1 + s*(3 - 2*s2 + 2*t1) + t2 + 
                  2*t1*t2 - Power(t2,2) + s1*(-1 - 2*t1 + 2*t2))*
                (4*(-1 + s1 + t1 - t2)*
                   Power(2*Power(s,2)*(-1 + s2) + 
                     (1 + s2)*
                      (-1 + s1*(s2 - t1) + t1 + t2 - s2*t2) - 
                     s*(-4 + 2*s1*s2 + 2*Power(s2,2) + t1 + t2 - 
                       s2*(t1 + t2)),2) - 
                  2*(6*Power(s,4)*Power(-1 + s2,2)*
                      (-1 + s1 + t1 - t2) + 
                     (-3 + s1*(1 + 4*s2 + Power(s2,2)) + 3*t1 + 
                       2*s2*(-1 + t1 - 2*t2) + 
                       Power(s2,2)*(-1 + t1 - t2) - t2)*
                      Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                     2*Power(s,3)*(-1 + s2)*
                      (12 + 6*Power(s1,2)*s2 - 15*t1 + 
                        4*Power(t1,2) + 
                        6*Power(s2,2)*(-1 + t1 - t2) + 9*t2 - 
                        2*t1*t2 - 2*Power(t2,2) + 
                        3*s1*
                       (-4 + 2*Power(s2,2) + t1 + 
                       s2*(-2 + t1 - 3*t2) + t2) + 
                        3*s2*(t1 - Power(t1,2) + t2 + Power(t2,2))) + 
                     Power(s,2)*
                      (-30 + 6*Power(s1,3)*Power(s2,2) + 52*t1 - 
                        25*Power(t1,2) + 3*Power(t1,3) + 
                        6*Power(s2,4)*(-1 + t1 - t2) - 8*t2 + 
                        8*t1*t2 - Power(t1,2)*t2 + 11*Power(t2,2) - 
                        t1*Power(t2,2) - Power(t2,3) - 
                        6*Power(s2,3)*
                       (Power(t1,2) + t1*(-1 + t2) - 2*t2*(1 + t2)) \
+ 6*Power(s1,2)*(3*Power(s2,3) + t1 + s2*(-5 + t1 + t2) - 
                       Power(s2,2)*(1 + t1 + 2*t2)) - 
                        2*s2*
                       (2*Power(t1,3) + Power(t1,2)*(-10 + t2) + 
                       t1*(8 - 3*t2 - 2*Power(t2,2)) + 
                       t2*(15 + 7*t2 - Power(t2,2))) + 
                        Power(s2,2)*
                       (30 + Power(t1,3) + 26*t2 - 9*Power(t2,2) - 
                       Power(t2,3) + Power(t1,2)*(11 + 3*t2) - 
                       t1*(42 + 8*t2 + 3*Power(t2,2))) + 
                        s1*(30 + 6*Power(s2,4) + 3*Power(t1,2) + 
                        6*Power(s2,3)*(-3 + 2*t1 - 5*t2) + 
                        2*t1*(-16 + t2) - 18*t2 + Power(t2,2) + 
                        2*s2*
                       (15 + 4*Power(t1,2) + 23*t2 - 
                       4*Power(t2,2) - 6*t1*(1 + t2)) + 
                        Power(s2,2)*
                        (-30 - 11*Power(t1,2) + 14*t2 + 
                        7*Power(t2,2) + 10*t1*(2 + t2)))) - 
                     2*s*(-7 + 3*Power(s1,3)*s2*(1 + s2)*(s2 - t1) + 
                        17*t1 - 13*Power(t1,2) + 3*Power(t1,3) + 
                        5*t2 - 4*t1*t2 - Power(t1,2)*t2 + 
                        3*Power(t2,2) - t1*Power(t2,2) - 
                        Power(t2,3) + 
                        3*Power(s2,4)*t2*(1 - t1 + t2) + 
                        Power(s2,2)*
                        (3 - 2*Power(t1,3) - 4*t2 - 13*Power(t2,2) + 
                        Power(t2,3) + Power(t1,2)*(7 + t2) + 
                        t1*(-8 + 3*t2)) + 
                        Power(s2,3)*
                        (3 + 4*t2 + Power(t2,2) - Power(t2,3) + 
                        Power(t1,2)*(3 + 2*t2) - 
                        t1*(6 + 6*t2 + Power(t2,2))) + 
                        s2*(-5 - Power(t1,3) - 14*t2 + 
                        6*Power(t2,2) + Power(t2,3) - 
                        Power(t1,2)*(3 + 2*t2) + 
                        t1*(9 + 16*t2 + 2*Power(t2,2))) + 
                        Power(s1,2)*
                        (3*Power(s2,4) - t1*(-5 + 3*t1 + t2) - 
                        Power(s2,3)*(-1 + t1 + 7*t2) + 
                        s2*
                       (-9 - Power(t1,2) + 5*t2 + 2*t1*(9 + t2)) + 
                        Power(s2,2)*
                        (-13 - 2*Power(t1,2) - 4*t2 + t1*(2 + 5*t2))) \
+ s1*(6 + 3*Power(s2,4)*(-1 + t1 - 2*t2) - Power(t1,2)*(-12 + t2) - 
                        7*t2 + Power(t2,2) + 
                        t1*(-18 + 4*t2 + Power(t2,2)) + 
                        Power(s2,3)*
                        (-6 - 4*Power(t1,2) - 2*t2 + 5*Power(t2,2) + 
                        2*t1*(5 + t2)) + 
                        Power(s2,2)*
                        (4 + Power(t1,3) + Power(t1,2)*(-3 + t2) + 
                        25*t2 - 2*t1*(1 + Power(t2,2))) + 
                        s2*(17 + 7*Power(t1,2) - Power(t1,3) + 2*t2 - 
                        6*Power(t2,2) + t1*(-23 - 18*t2 + Power(t2,2))\
))))))/((-1 + s1 + t1 - t2)*Power(-1 + s - s*s2 + s1*s2 + t1 - s2*t2,
                 2)))/((-1 + t1)*(-1 + t2)) + 
          ((2*(2*Power(s,2)*(-1 + s2) + 
                  (1 + s2)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2) - 
                  s*(-4 + 2*s1*s2 + 2*Power(s2,2) + t1 + t2 - 
                     s2*(t1 + t2)))*
                (-12 + 2*Power(s,3) + 2*Power(s1,3) - 16*s2 + 
                  5*Power(s2,2) - 2*Power(s2,3) - 4*t1 - 5*s2*t1 + 
                  4*Power(s2,2)*t1 + 
                  Power(s1,2)*(1 + 3*s2 + 3*t1 - 5*t2) + 15*t2 + 
                  6*s2*t2 - Power(s2,2)*t2 - 3*t1*t2 - 5*s2*t1*t2 - 
                  2*Power(t2,2) + 4*s2*Power(t2,2) + 
                  t1*Power(t2,2) - Power(t2,3) - 
                  Power(s,2)*(-4 + 6*s2 - 5*t1 + t2) + 
                  s1*(-13 + s2 - Power(s2,2) + 7*s2*t1 + t2 - 
                     7*s2*t2 - 4*t1*t2 + 4*Power(t2,2)) + 
                  s*(13 - 4*Power(s1,2) + 6*Power(s2,2) + 2*t1 - 
                     3*t2 + 6*t1*t2 - 4*Power(t2,2) + 
                     s2*(-9 - 9*t1 + 2*t2) + 
                     s1*(-3 + s2 - 8*t1 + 8*t2))))/
              (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2) - 
             4*(-6 + Power(s,4) - Power(s1,3)*(-3 + s2) + 17*s2 + 
                11*Power(s2,2) - 3*Power(s2,3) + Power(s2,4) + 4*t1 + 
                7*s2*t1 + 3*Power(s2,2)*t1 - 2*Power(s2,3)*t1 - 
                3*t2 - 24*s2*t2 - 8*Power(s2,2)*t2 + Power(s2,3)*t2 - 
                5*t1*t2 + 4*s2*t1*t2 + 3*Power(s2,2)*t1*t2 + 
                4*Power(t2,2) + 5*s2*Power(t2,2) - 
                3*Power(s2,2)*Power(t2,2) - t1*Power(t2,2) - 
                s2*t1*Power(t2,2) + Power(t2,3) + s2*Power(t2,3) - 
                Power(s,3)*(-1 + 4*s2 - 3*t1 + t2) + 
                Power(s1,2)*
                 (8 + s2 - Power(s2,2) + 6*t1 - 2*s2*t1 - 6*t2 + 
                   4*s2*t2) + 
                Power(s,2)*(9 - 3*Power(s1,2) + 6*Power(s2,2) + t1 - 
                   5*t2 + 4*t1*t2 - 3*Power(t2,2) + 
                   s2*(-5 - 8*t1 + 3*t2) + s1*(-2 + s2 - 6*t1 + 6*t2)\
) + s*(-8 + 2*Power(s1,3) - 4*Power(s2,3) - 7*t1 + 
                   Power(s1,2)*(4*s2 + 3*t1 - 5*t2) + 
                   Power(s2,2)*(7 + 7*t1 - 3*t2) + 18*t2 - t1*t2 - 
                   4*Power(t2,2) + t1*Power(t2,2) - Power(t2,3) + 
                   s1*(-21 - 2*Power(s2,2) - 3*t1 + 
                      s2*(6 + 10*t1 - 11*t2) + 7*t2 - 4*t1*t2 + 
                      4*Power(t2,2)) + 
                   s2*(-20 - 3*t1 + 13*t2 - 7*t1*t2 + 6*Power(t2,2))) \
+ s1*(3 + Power(s2,3) + 2*Power(t1,2) - 3*t1*(-1 + t2) - 12*t2 + 
                   2*Power(t2,2) + Power(s2,2)*(-3 - 4*t1 + 5*t2) + 
                   s2*(25 - 9*t2 - 4*Power(t2,2) + t1*(5 + 3*t2)))) - 
             ((5 + Power(s,2) - Power(s1,2) - 2*s2 + Power(s2,2) + 
                  2*t1 - 2*s2*t1 + s*(3 - 2*s2 + 2*t1) + t2 + 
                  2*t1*t2 - Power(t2,2) + s1*(-1 - 2*t1 + 2*t2))*
                (4*(-1 + s1 + t1 - t2)*
                   Power(2*Power(s,2)*(-1 + s2) + 
                     (1 + s2)*
                      (-1 + s1*(s2 - t1) + t1 + t2 - s2*t2) - 
                     s*(-4 + 2*s1*s2 + 2*Power(s2,2) + t1 + t2 - 
                       s2*(t1 + t2)),2) - 
                  2*(6*Power(s,4)*Power(-1 + s2,2)*
                      (-1 + s1 + t1 - t2) + 
                     (-3 + s1*(1 + 4*s2 + Power(s2,2)) + 3*t1 + 
                       2*s2*(-1 + t1 - 2*t2) + 
                       Power(s2,2)*(-1 + t1 - t2) - t2)*
                      Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                     2*Power(s,3)*(-1 + s2)*
                      (12 + 6*Power(s1,2)*s2 - 15*t1 + 
                        4*Power(t1,2) + 
                        6*Power(s2,2)*(-1 + t1 - t2) + 9*t2 - 
                        2*t1*t2 - 2*Power(t2,2) + 
                        3*s1*
                       (-4 + 2*Power(s2,2) + t1 + 
                       s2*(-2 + t1 - 3*t2) + t2) + 
                        3*s2*(t1 - Power(t1,2) + t2 + Power(t2,2))) + 
                     Power(s,2)*
                      (-30 + 6*Power(s1,3)*Power(s2,2) + 52*t1 - 
                        25*Power(t1,2) + 3*Power(t1,3) + 
                        6*Power(s2,4)*(-1 + t1 - t2) - 8*t2 + 
                        8*t1*t2 - Power(t1,2)*t2 + 11*Power(t2,2) - 
                        t1*Power(t2,2) - Power(t2,3) - 
                        6*Power(s2,3)*
                       (Power(t1,2) + t1*(-1 + t2) - 2*t2*(1 + t2)) \
+ 6*Power(s1,2)*(3*Power(s2,3) + t1 + s2*(-5 + t1 + t2) - 
                       Power(s2,2)*(1 + t1 + 2*t2)) - 
                        2*s2*
                       (2*Power(t1,3) + Power(t1,2)*(-10 + t2) + 
                       t1*(8 - 3*t2 - 2*Power(t2,2)) + 
                       t2*(15 + 7*t2 - Power(t2,2))) + 
                        Power(s2,2)*
                       (30 + Power(t1,3) + 26*t2 - 9*Power(t2,2) - 
                       Power(t2,3) + Power(t1,2)*(11 + 3*t2) - 
                       t1*(42 + 8*t2 + 3*Power(t2,2))) + 
                        s1*(30 + 6*Power(s2,4) + 3*Power(t1,2) + 
                        6*Power(s2,3)*(-3 + 2*t1 - 5*t2) + 
                        2*t1*(-16 + t2) - 18*t2 + Power(t2,2) + 
                        2*s2*
                       (15 + 4*Power(t1,2) + 23*t2 - 
                       4*Power(t2,2) - 6*t1*(1 + t2)) + 
                        Power(s2,2)*
                        (-30 - 11*Power(t1,2) + 14*t2 + 
                        7*Power(t2,2) + 10*t1*(2 + t2)))) - 
                     2*s*(-7 + 3*Power(s1,3)*s2*(1 + s2)*(s2 - t1) + 
                        17*t1 - 13*Power(t1,2) + 3*Power(t1,3) + 
                        5*t2 - 4*t1*t2 - Power(t1,2)*t2 + 
                        3*Power(t2,2) - t1*Power(t2,2) - 
                        Power(t2,3) + 
                        3*Power(s2,4)*t2*(1 - t1 + t2) + 
                        Power(s2,2)*
                        (3 - 2*Power(t1,3) - 4*t2 - 13*Power(t2,2) + 
                        Power(t2,3) + Power(t1,2)*(7 + t2) + 
                        t1*(-8 + 3*t2)) + 
                        Power(s2,3)*
                        (3 + 4*t2 + Power(t2,2) - Power(t2,3) + 
                        Power(t1,2)*(3 + 2*t2) - 
                        t1*(6 + 6*t2 + Power(t2,2))) + 
                        s2*(-5 - Power(t1,3) - 14*t2 + 
                        6*Power(t2,2) + Power(t2,3) - 
                        Power(t1,2)*(3 + 2*t2) + 
                        t1*(9 + 16*t2 + 2*Power(t2,2))) + 
                        Power(s1,2)*
                        (3*Power(s2,4) - t1*(-5 + 3*t1 + t2) - 
                        Power(s2,3)*(-1 + t1 + 7*t2) + 
                        s2*
                       (-9 - Power(t1,2) + 5*t2 + 2*t1*(9 + t2)) + 
                        Power(s2,2)*
                        (-13 - 2*Power(t1,2) - 4*t2 + t1*(2 + 5*t2))) \
+ s1*(6 + 3*Power(s2,4)*(-1 + t1 - 2*t2) - Power(t1,2)*(-12 + t2) - 
                        7*t2 + Power(t2,2) + 
                        t1*(-18 + 4*t2 + Power(t2,2)) + 
                        Power(s2,3)*
                        (-6 - 4*Power(t1,2) - 2*t2 + 5*Power(t2,2) + 
                        2*t1*(5 + t2)) + 
                        Power(s2,2)*
                        (4 + Power(t1,3) + Power(t1,2)*(-3 + t2) + 
                        25*t2 - 2*t1*(1 + Power(t2,2))) + 
                        s2*(17 + 7*Power(t1,2) - Power(t1,3) + 2*t2 - 
                        6*Power(t2,2) + t1*(-23 - 18*t2 + Power(t2,2))\
))))))/((-1 + s1 + t1 - t2)*Power(-1 + s - s*s2 + s1*s2 + t1 - s2*t2,
                 2)))/((-1 + s1)*(-1 + t2)) + 
          ((-2*(2*Power(s,2)*(-1 + s2) + 
                  (1 + s2)*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2) - 
                  s*(-4 + 2*s1*s2 + 2*Power(s2,2) + t1 + t2 - 
                     s2*(t1 + t2)))*
                (11 + 4*Power(s,3) + 19*s2 - Power(s2,3) - 2*t1 - 
                  8*s2*t1 + Power(s2,2)*t1 - 4*t2 + 5*s2*t2 + 
                  3*Power(s2,2)*t2 - 2*t1*t2 - 2*s2*t1*t2 + 
                  2*Power(t2,2) - s2*Power(t2,2) + t1*Power(t2,2) - 
                  Power(t2,3) - Power(s1,2)*(2*s2 - 2*t1 + t2) + 
                  Power(s,2)*(4 - 7*s1 - 9*s2 + 3*t1 + 5*t2) + 
                  s1*(4 - 3*Power(s2,2) - (2 + 3*t1)*t2 + 
                     2*Power(t2,2) + s2*(-8 + 3*t1 + 3*t2)) + 
                  s*(-12 + 3*Power(s1,2) + 6*Power(s2,2) + 3*t1 + 
                     s1*(-1 + 10*s2 - 5*t1 - 3*t2) + t2 + 4*t1*t2 - 
                     2*s2*(1 + 2*t1 + 4*t2))))/
              (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2) + 
             4*(9 + 2*Power(s,4) - 4*Power(s1,3) - 15*s2 - 
                17*Power(s2,2) - Power(s2,3) - 7*t1 + 4*s2*t1 + 
                5*Power(s2,2)*t1 - 2*Power(t1,2) - 
                Power(s1,2)*(9 + 2*s2 + 3*t1 - 10*t2) + t2 + 
                12*s2*t2 - 3*Power(s2,2)*t2 - Power(s2,3)*t2 + 
                6*s2*t1*t2 + Power(s2,2)*t1*t2 + Power(t2,2) - 
                5*s2*Power(t2,2) - t1*Power(t2,2) - 
                s2*t1*Power(t2,2) + Power(t2,3) + s2*Power(t2,3) + 
                2*Power(s,3)*(1 - 2*s1 - 3*s2 + t1 + t2) + 
                Power(s,2)*(-8 + 2*Power(s1,2) + 6*Power(s2,2) + t1 + 
                   s1*(-3 + 8*s2 - 4*t1 - t2) + t2 + 3*t1*t2 - 
                   Power(t2,2) - s2*(3 + 4*t1 + 5*t2)) + 
                s1*(3 - Power(s2,2)*(-5 + t2) + 7*t2 - 
                   7*Power(t2,2) + t1*(-11 + 6*t2) + 
                   s2*(-17 + t1*(-4 + t2) + 7*t2 - Power(t2,2))) + 
                s*(3 - 2*Power(s2,3) + 6*t1 + 2*Power(t1,2) - 12*t2 - 
                   7*t1*t2 + 5*Power(t2,2) + t1*Power(t2,2) - 
                   Power(t2,3) - Power(s1,2)*(-6 + 2*s2 - 2*t1 + t2) + 
                   2*Power(s2,2)*(1 + t1 + 2*t2) + 
                   s2*(25 + 2*t2 + Power(t2,2) - 4*t1*(1 + t2)) + 
                   s1*(-4*Power(s2,2) + t1*(5 - 3*t2) + 
                      2*s2*(-2 + 2*t1 + t2) + 
                      2*(6 - 5*t2 + Power(t2,2))))) + 
             ((-3 + 2*Power(s,2) + Power(s1,2) + s2 + Power(s2,2) + 
                  2*t1 - s2*t1 + s1*(2*s2 - t1 - 2*t2) - 2*s2*t2 + 
                  t1*t2 + Power(t2,2) + 
                  s*(2 - 3*s1 - 3*s2 + t1 + 3*t2))*
                (4*(-1 + s1 + t1 - t2)*
                   Power(2*Power(s,2)*(-1 + s2) + 
                     (1 + s2)*
                      (-1 + s1*(s2 - t1) + t1 + t2 - s2*t2) - 
                     s*(-4 + 2*s1*s2 + 2*Power(s2,2) + t1 + t2 - 
                        s2*(t1 + t2)),2) - 
                  2*(6*Power(s,4)*Power(-1 + s2,2)*
                      (-1 + s1 + t1 - t2) + 
                     (-3 + s1*(1 + 4*s2 + Power(s2,2)) + 3*t1 + 
                        2*s2*(-1 + t1 - 2*t2) + 
                        Power(s2,2)*(-1 + t1 - t2) - t2)*
                      Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                     2*Power(s,3)*(-1 + s2)*
                      (12 + 6*Power(s1,2)*s2 - 15*t1 + 
                        4*Power(t1,2) + 
                        6*Power(s2,2)*(-1 + t1 - t2) + 9*t2 - 
                        2*t1*t2 - 2*Power(t2,2) + 
                        3*s1*
                        (-4 + 2*Power(s2,2) + t1 + 
                        s2*(-2 + t1 - 3*t2) + t2) + 
                        3*s2*(t1 - Power(t1,2) + t2 + Power(t2,2))) + 
                     Power(s,2)*
                      (-30 + 6*Power(s1,3)*Power(s2,2) + 52*t1 - 
                        25*Power(t1,2) + 3*Power(t1,3) + 
                        6*Power(s2,4)*(-1 + t1 - t2) - 8*t2 + 
                        8*t1*t2 - Power(t1,2)*t2 + 11*Power(t2,2) - 
                        t1*Power(t2,2) - Power(t2,3) - 
                        6*Power(s2,3)*
                        (Power(t1,2) + t1*(-1 + t2) - 2*t2*(1 + t2)) \
+ 6*Power(s1,2)*(3*Power(s2,3) + t1 + s2*(-5 + t1 + t2) - 
                        Power(s2,2)*(1 + t1 + 2*t2)) - 
                        2*s2*
                        (2*Power(t1,3) + Power(t1,2)*(-10 + t2) + 
                        t1*(8 - 3*t2 - 2*Power(t2,2)) + 
                        t2*(15 + 7*t2 - Power(t2,2))) + 
                        Power(s2,2)*
                        (30 + Power(t1,3) + 26*t2 - 9*Power(t2,2) - 
                        Power(t2,3) + Power(t1,2)*(11 + 3*t2) - 
                        t1*(42 + 8*t2 + 3*Power(t2,2))) + 
                        s1*(30 + 6*Power(s2,4) + 3*Power(t1,2) + 
                        6*Power(s2,3)*(-3 + 2*t1 - 5*t2) + 
                        2*t1*(-16 + t2) - 18*t2 + Power(t2,2) + 
                        2*s2*
                        (15 + 4*Power(t1,2) + 23*t2 - 
                        4*Power(t2,2) - 6*t1*(1 + t2)) + 
                        Power(s2,2)*
                        (-30 - 11*Power(t1,2) + 14*t2 + 
                        7*Power(t2,2) + 10*t1*(2 + t2)))) - 
                     2*s*(-7 + 3*Power(s1,3)*s2*(1 + s2)*(s2 - t1) + 
                        17*t1 - 13*Power(t1,2) + 3*Power(t1,3) + 
                        5*t2 - 4*t1*t2 - Power(t1,2)*t2 + 
                        3*Power(t2,2) - t1*Power(t2,2) - Power(t2,3) + 
                        3*Power(s2,4)*t2*(1 - t1 + t2) + 
                        Power(s2,2)*
                        (3 - 2*Power(t1,3) - 4*t2 - 13*Power(t2,2) + 
                        Power(t2,3) + Power(t1,2)*(7 + t2) + 
                        t1*(-8 + 3*t2)) + 
                        Power(s2,3)*
                        (3 + 4*t2 + Power(t2,2) - Power(t2,3) + 
                        Power(t1,2)*(3 + 2*t2) - 
                        t1*(6 + 6*t2 + Power(t2,2))) + 
                        s2*(-5 - Power(t1,3) - 14*t2 + 
                        6*Power(t2,2) + Power(t2,3) - 
                        Power(t1,2)*(3 + 2*t2) + 
                        t1*(9 + 16*t2 + 2*Power(t2,2))) + 
                        Power(s1,2)*
                        (3*Power(s2,4) - t1*(-5 + 3*t1 + t2) - 
                        Power(s2,3)*(-1 + t1 + 7*t2) + 
                        s2*
                       (-9 - Power(t1,2) + 5*t2 + 2*t1*(9 + t2)) + 
                        Power(s2,2)*
                        (-13 - 2*Power(t1,2) - 4*t2 + t1*(2 + 5*t2))) \
+ s1*(6 + 3*Power(s2,4)*(-1 + t1 - 2*t2) - Power(t1,2)*(-12 + t2) - 
                        7*t2 + Power(t2,2) + 
                        t1*(-18 + 4*t2 + Power(t2,2)) + 
                        Power(s2,3)*
                        (-6 - 4*Power(t1,2) - 2*t2 + 5*Power(t2,2) + 
                        2*t1*(5 + t2)) + 
                        Power(s2,2)*
                        (4 + Power(t1,3) + Power(t1,2)*(-3 + t2) + 
                        25*t2 - 2*t1*(1 + Power(t2,2))) + 
                        s2*(17 + 7*Power(t1,2) - Power(t1,3) + 2*t2 - 
                        6*Power(t2,2) + t1*(-23 - 18*t2 + Power(t2,2)))\
)))))/((-1 + s1 + t1 - t2)*Power(-1 + s - s*s2 + s1*s2 + t1 - s2*t2,2)))/
           ((-1 + s1)*(-s + s1 - t2)))/(s - s2 + t1))*
     B1(s2,1 - s1 - t1 + t2,1 - s + s2 - t1))/(4.*Power(Pi,2)) + 
  (((-64*s*(2 - t1 - Power(t1,2) - s1*Power(t1,2) + 
            Power(s,2)*(t1 - t2) - s2*(-1 + t1)*(-1 + 2*t1 - t2) + 
            Power(s2,2)*(-2 + s1 + 2*t1 - t2) - 2*t2 + 3*t1*t2 - 
            s*(1 - Power(t1,2) + s2*(-2 + s1 + 3*t1 - 2*t2) - 3*t2 + 
               t1*(1 + s1 + t2)))*
          ((-1 - Power(s,4) + 24*s2 + 24*Power(s2,2) + 
               3*Power(s2,3) - Power(s2,4) - 
               Power(s1,2)*(1 + s2)*(2 + s2 - t1) + 
               Power(s,3)*(-2 + 2*s1 + 4*s2 - t1) - 7*t1 - 29*s2*t1 - 
               11*Power(s2,2)*t1 + Power(s2,3)*t1 + 4*Power(t1,2) + 
               5*s2*Power(t1,2) - 
               Power(s,2)*(-21 + Power(s1,2) + 6*Power(s2,2) + 
                  s1*(3 + 6*s2 - 2*t1) + 9*t1 - s2*(10 + 3*t1) - 5*t2\
) - 4*t2 + 12*s2*t2 + 4*Power(s2,2)*t2 + t1*t2 - 5*s2*t1*t2 - 
               Power(t2,2) + s2*Power(t2,2) + 
               s*(-19 + 4*Power(s2,3) + 
                  Power(s1,2)*(4 + 2*s2 - t1) + 18*t1 - 
                  2*Power(t1,2) - Power(s2,2)*(11 + 3*t1) + 
                  s2*(-49 + 22*t1 - 9*t2) + 
                  s1*(-13 + 6*Power(s2,2) - 4*s2*(-1 + t1) + 6*t1 - 
                     5*t2) - 8*t2 + 2*t1*t2) + 
               s1*(10 - 2*Power(s2,3) + 2*Power(s2,2)*(-1 + t1) + 
                  Power(t1,2) + 9*t2 - 2*t1*(4 + t2) + 
                  s2*(11 - 4*t1 + t2)))/((-1 + s2)*(-1 + t1)) + 
            (-10 + Power(s,3)*(-2 + s2) + 10*s2 + 9*Power(s2,2) - 
               Power(s2,3) - Power(s2,4) + 12*t1 + 5*s2*t1 + 
               2*Power(s2,2)*t1 + Power(s2,3)*t1 - 4*Power(t1,2) - 
               5*s2*Power(t1,2) - 
               Power(s1,2)*(-5 + s2 + Power(s2,2) + 2*t1 - s2*t1) + 
               Power(s,2)*(8 - 2*s1*(-2 + s2) - 3*Power(s2,2) + 
                  s2*(2 + t1) - 2*t2) - t2 - 3*s2*t2 - 
               4*Power(s2,2)*t2 - t1*t2 + 5*s2*t1*t2 + Power(t2,2) - 
               s2*Power(t2,2) + 
               s1*(10 - 2*Power(s2,3) + t1 - Power(t1,2) + 
                  2*Power(s2,2)*(1 + t1) - 5*t2 + 2*t1*t2 - 
                  s2*(-16 + 8*t1 + t2)) + 
               s*(Power(s1,2)*(-2 + s2) + 3*Power(s2,3) + 
                  Power(s2,2)*(1 - 2*t1) - 7*t1 + 2*Power(t1,2) + 
                  4*t2 - 2*t1*t2 + 
                  s1*(-15 + 4*Power(s2,2) + 3*t1 - s2*(5 + 2*t1) + 
                     2*t2) + s2*(-17 - 2*t1 + 6*t2)))/
             ((-1 + s2)*(-s + s2 - t1)) + 
            (-6 - 2*Power(s1,3)*(-3 + s2) - 3*s2 + 20*Power(s2,2) + 
               8*Power(s2,3) - 2*Power(s2,4) - 5*t1 - s2*t1 - 
               2*Power(s2,2)*t1 + Power(t1,2) + 2*s2*Power(t1,2) + 
               6*t2 + 6*s2*t2 - 11*Power(s2,2)*t2 + 
               2*Power(s2,3)*t2 - 3*t1*t2 - 2*s2*t1*t2 - 
               2*Power(t2,2) + 4*s2*Power(t2,2) + 
               Power(s,3)*(-5 + 2*s2 + t2) - 
               Power(s1,2)*(1 + 6*Power(s2,2) + t1 + 8*t2 - 
                  2*s2*(9 + t2)) + 
               Power(s,2)*(6 + 15*s2 - 6*Power(s2,2) - 6*t2 + 
                  Power(t2,2) - 2*s1*(-7 + 3*s2 + t2)) + 
               s1*(-10 - 6*Power(s2,3) + Power(t1,2) + 7*t2 + 
                  2*Power(t2,2) + t1*(6 + t2) + 
                  Power(s2,2)*(19 + 4*t2) - s2*(-9 + t1 + 20*t2)) + 
               s*(3 + 6*Power(s2,3) + t1 - 10*t2 - 2*Power(t2,2) - 
                  3*Power(s2,2)*(6 + t2) + 
                  Power(s1,2)*(-16 + 6*s2 + t2) + 
                  s2*(-20 + 2*t1 + 15*t2 - Power(t2,2)) + 
                  s1*(1 + 12*Power(s2,2) + t1 + 15*t2 - Power(t2,2) - 
                     2*s2*(15 + t2))))/((-1 + s1)*(-s + s1 - t2)) + 
            (-11 + Power(s,4) - Power(s1,3)*(-3 + s2) - 2*s2 + 
               12*Power(s2,2) + Power(s2,3) - 2*Power(s2,4) + 10*t1 + 
               10*s2*t1 + 2*Power(s2,2)*t1 + Power(s2,3)*t1 - 
               Power(t1,2) + 6*t2 - 6*s2*t2 - 8*Power(s2,2)*t2 + 
               Power(s2,3)*t2 - 4*t1*t2 - s2*t1*t2 - Power(t2,2) + 
               3*s2*Power(t2,2) + 
               Power(s,3)*(-3 - 2*s1 - s2 + t1 + t2) + 
               Power(s1,2)*(5 - 4*Power(s2,2) - t1 - 5*t2 + 
                  s2*(9 + t1 + t2)) + 
               Power(s,2)*(-3 + Power(s1,2) - 3*Power(s2,2) + 3*t1 - 
                  8*t2 + t1*t2 - s2*(-8 + t1 + t2) - 
                  s1*(-15 + s2 + 2*t1 + t2)) + 
               s1*(-1 - 5*Power(s2,3) + Power(t1,2) - 4*t2 + 
                  2*Power(t2,2) + t1*(4 + t2) + 
                  2*Power(s2,2)*(4 + t1 + t2) - s2*(-20 + t1 + 15*t2)\
) + s*(13 + 5*Power(s2,3) - 12*t1 + Power(s1,2)*(-13 + 3*s2 + t1) - 
                  2*Power(t2,2) - Power(s2,2)*(6 + t1 + t2) - 
                  s2*(9 - 16*t2 + t1*(5 + t2)) + 
                  s1*(-9 + 8*Power(s2,2) + t1 + 14*t2 - t1*t2 - 
                     s2*(24 + t2))))/((s - s2 + t1)*(s - s1 + t2)) + 
            (-2*Power(s,4) - Power(s1,3)*(-1 + s2) + 17*s2 + 
               22*Power(s2,2) + 3*Power(s2,3) - 2*Power(s2,4) - 
               9*t1 - 13*s2*t1 - 2*Power(s2,2)*t1 + Power(s2,3)*t1 + 
               Power(t1,2) + 
               Power(s,3)*(-2 + 5*s1 + 8*s2 - t1 - 2*t2) - 3*t2 - 
               6*s2*t2 - 3*Power(s2,2)*t2 + Power(s2,3)*t2 + 
               4*t1*t2 + s2*t1*t2 + Power(t2,2) - 3*s2*Power(t2,2) + 
               Power(s1,2)*(3 - 4*Power(s2,2) + t2 + 
                  s2*(4 + t1 + t2)) - 
               s1*(-4 + 5*Power(s2,3) + Power(t1,2) + 3*t2 + t1*t2 + 
                  2*Power(t2,2) + s2*(-17 + 2*t2) - 
                  Power(s2,2)*(5 + 2*t1 + 2*t2)) + 
               Power(s,2)*(17 - 4*Power(s1,2) - 12*Power(s2,2) + 
                  t1 - 3*t2 - t1*t2 + s1*(4 - 15*s2 + 2*t1 + 3*t2) + 
                  s2*(10 + 3*t1 + 5*t2)) + 
               s*(-15 + Power(s1,3) + 8*Power(s2,3) + 9*t1 + 
                  Power(t1,2) + Power(s1,2)*(-4 + 8*s2 - t1 - t2) + 
                  7*t2 + 2*Power(t2,2) - 
                  Power(s2,2)*(11 + 3*t1 + 4*t2) + 
                  s2*(-41 + t1 + 8*t2 + t1*t2) + 
                  s1*(-12 + 15*Power(s2,2) + t1*(-2 + t2) + t2 - 
                     s2*(12 + 4*t1 + 5*t2))))/((-1 + t1)*(-1 + t2)) + 
            (-10 - 2*Power(s,4) - 2*Power(s1,3)*(-2 + s2) + 4*s2 + 
               18*Power(s2,2) + 7*Power(s2,3) - 2*Power(s2,4) - t1 + 
               s2*t1 + 2*Power(s2,2)*t1 - Power(t1,2) - 
               2*s2*Power(t1,2) + 6*t2 - 11*s2*t2 - 
               10*Power(s2,2)*t2 + 2*Power(s2,3)*t2 + 3*t1*t2 + 
               2*s2*t1*t2 + 2*Power(t2,2) - 4*s2*Power(t2,2) + 
               Power(s,3)*(6*s1 + 8*s2 - 3*(1 + t2)) + 
               Power(s1,2)*(2 - 6*Power(s2,2) + t1 - 2*t2 + 
                  s2*(13 + 2*t2)) - 
               Power(s,2)*(-12 + 6*Power(s1,2) + 12*Power(s2,2) + 
                  3*s1*(-3 + 6*s2 - 2*t2) + 4*t2 + Power(t2,2) - 
                  2*s2*(5 + 4*t2)) - 
               s1*(3 + 6*Power(s2,3) + Power(t1,2) - 
                  s2*(16 + t1 - 7*t2) + 4*t2 + t1*t2 + 
                  2*Power(t2,2) - Power(s2,2)*(17 + 4*t2)) + 
               s*(4 + 2*Power(s1,3) + 8*Power(s2,3) + 
                  Power(s1,2)*(-10 + 12*s2 - 3*t2) + 7*t2 + 
                  2*Power(t2,2) - 7*Power(s2,2)*(2 + t2) + 
                  s2*(-28 - 2*t1 + 12*t2 + Power(t2,2)) + 
                  s1*(-13 + 18*Power(s2,2) - 2*t1 + 5*t2 + 
                     Power(t2,2) - s2*(23 + 10*t2))))/
             ((-1 + s1)*(-1 + t2))))/
        (Power(s,3)*(-1 + s2)*(-1 + s1 + t1 - t2) - 
          Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
          Power(s,2)*(3 + Power(s1,2)*s2 - 3*t1 + Power(t1,2) + 
             2*Power(s2,2)*(-1 + t1 - t2) + 2*t2 - t1*t2 + 
             s2*t1*(1 - t1 + t2) + 
             s1*(-3 + 2*Power(s2,2) + t2 - s2*(1 + t2))) + 
          s*(Power(s1,2)*(1 + s2)*(s2 - t1) + 
             Power(s2,3)*(-1 + t1 - t2) - 
             Power(s2,2)*(-1 + t1)*(-1 + t1 - t2) - 
             s2*(-2 + t1 + Power(t1,2) - 4*t2 - t1*t2) + 
             2*(1 + Power(t1,2) - t2 - t1*(2 + t2)) + 
             s1*(-2 + Power(s2,3) + 5*t1 + Power(t1,2) + 2*t2 - 
                Power(s2,2)*t2 - t1*t2 - 
                s2*(3 + t1 + Power(t1,2) + t2 - t1*t2)))) + 
       64*((-8 - 2*Power(s2,4)*(-1 + t1) + 22*t1 - 2*Power(s,4)*t1 + 
             4*s1*t1 - 35*Power(t1,2) + 27*Power(t1,3) - 
             s1*Power(t1,3) - 6*Power(t1,4) + 
             Power(s,3)*(t1 + 2*s1*t1 - 4*Power(t1,2) + 
                s2*(-2 + 8*t1) - t2) + 8*t2 - 12*t1*t2 - 
             Power(t1,2)*t2 + 2*Power(t1,3)*t2 + 
             Power(s2,3)*(2 + s1 - 6*t1 - 2*s1*t1 + 4*Power(t1,2) + 
                t2) + Power(s2,2)*
              (-15 + (-6 + 4*s1)*Power(t1,2) - 2*Power(t1,3) + t2 + 
                t1*(23 - 7*s1 + 2*t2)) + 
             Power(s,2)*(-1 + Power(s2,2)*(6 - 12*t1) + 29*t1 - 
                16*Power(t1,2) - 2*Power(t1,3) + 
                s1*(s2 - 6*s2*t1 + t1*(-9 + 4*t1)) + 7*t2 + 4*t1*t2 + 
                s2*(-8*t1 + 12*Power(t1,2) + 3*t2)) + 
             s2*(-14 + 50*t1 + 16*Power(t1,3) + 
                s1*(-4 + 7*Power(t1,2) - 2*Power(t1,3)) + 4*t2 - 
                Power(t1,2)*(52 + 5*t2)) + 
             s*(6 + (-23 + 2*s1)*Power(t1,3) + 
                Power(s2,3)*(-6 + 8*t1) + 
                Power(s2,2)*(-2 + 13*t1 - 12*Power(t1,2) + 
                   s1*(-2 + 6*t1) - 3*t2) - 14*t2 + 
                2*t1*(-23 + 3*s1 + t2) + 
                Power(t1,2)*(66 - 10*s1 + 7*t2) - 
                2*s2*(-12 + s1 + 29*t1 - 8*s1*t1 - 11*Power(t1,2) + 
                   4*s1*Power(t1,2) - 2*Power(t1,3) + 3*t2 + 3*t1*t2)))/
           ((-1 + s2)*(-1 + t1)*(s2*(-1 + t1) - t1*(-1 + s + t1))) + 
          (-8 - 2*Power(s2,4)*(-1 + t1) + 10*t1 + 10*s1*t1 + 
             13*Power(t1,2) - 12*s1*Power(t1,2) - 21*Power(t1,3) + 
             5*s1*Power(t1,3) + 6*Power(t1,4) + 
             2*Power(s,3)*((2 + s2)*t1 - t2) + 8*t2 - 20*t1*t2 + 
             11*Power(t1,2)*t2 - 2*Power(t1,3)*t2 + 
             Power(s2,3)*(6 + s1 - 10*t1 - 2*s1*t1 + 4*Power(t1,2) + 
                t2) + Power(s2,2)*
              (3 + 12*Power(t1,2) - 2*Power(t1,3) + 
                s1*(-4 - 3*t1 + 4*Power(t1,2)) + 5*t2 - t1*(13 + 2*t2)\
) - s2*(2 + 10*Power(t1,3) + 
                s1*(10 - 16*t1 + 3*Power(t1,2) + 2*Power(t1,3)) - 
                12*t2 + 16*t1*(1 + t2) - Power(t1,2)*(28 + 3*t2)) + 
             Power(s,2)*(-2 + Power(s2,2)*(2 - 6*t1) + 
                13*Power(t1,2) + t1*(-6 + 4*s1 - 7*t2) + 10*t2 + 
                s2*(6 - 17*t1 + 4*Power(t1,2) - 2*s1*(1 + t1) + 5*t2)) \
+ s*(8 + 15*Power(t1,3) + Power(s2,3)*(-4 + 6*t1) + 
                Power(t1,2)*(-27 + 9*s1 - 7*t2) - 16*t2 + 
                t1*(-5 - 12*s1 + 23*t2) + 
                Power(s2,2)*(s1 + 23*t1 + 4*s1*t1 - 8*Power(t1,2) - 
                   4*(3 + t2)) + 
                s2*(-7 - 26*Power(t1,2) + 2*Power(t1,3) - 
                   2*s1*(-5 + 2*t1 + 2*Power(t1,2)) - 17*t2 + 
                   9*t1*(3 + t2))))/
           ((-1 + s2)*(-s + s2 - t1)*(s2*(-1 + t1) - t1*(-1 + s + t1))) \
- (Power(s,6)*t1*(t1 - t2) - 
             4*(s2 - t1)*(-1 + t1)*
              Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
             Power(s,5)*(3*Power(t1,3) + 
                Power(t1,2)*(-8 + s2 - 3*t2) - t2*(s2 + 2*t2) + 
                t1*(-1 + (6 + s1)*t2 + s2*(3 - s1 + 6*t2))) + 
             Power(s,4)*(3*Power(t1,4) - 3*Power(t1,3)*(4 + t2) + 
                4*t2*(-1 + 3*t2) + 
                Power(s2,2)*
                 (2 - 13*Power(t1,2) + s1*(-1 + 5*t1) + 5*t2 - 
                   14*t1*t2) + 
                t1*(1 - 5*s1*(-1 + t2) - 14*t2 - 6*Power(t2,2)) + 
                Power(t1,2)*
                 (11 - Power(s1,2) + 12*t2 + 3*s1*(6 + t2)) + 
                s2*(-1 + 3*Power(t1,3) - 
                   4*Power(t1,2)*(-6 + 3*s1 - 4*t2) + 9*t2 - 
                   3*s1*t2 + 8*Power(t2,2) + 
                   t1*(-10 - s1 + Power(s1,2) - 25*t2 - 5*s1*t2))) - 
             s*(Power(t1,5)*(5 - 2*s1 + Power(s1,2) - 4*t2) + 
                2*Power(s2,5)*(-1 + t1)*(2 + s1 - 2*t1 - t2) + 
                8*Power(-1 + t2,2) - 
                2*t1*(-1 + t2)*(-13 + s1 + 16*t2) + 
                Power(t1,4)*
                 (-11 + 4*Power(s1,2) - 5*s1*(-3 + t2) + 21*t2) + 
                Power(t1,2)*
                 (27 + 12*Power(s1,2) + 5*s1*(-5 + t2) - 23*t2 + 
                   28*Power(t2,2)) - 
                Power(t1,3)*
                 (3 + 19*Power(s1,2) + 36*t2 + 6*Power(t2,2) - 
                   2*s1*(5 + 3*t2)) + 
                Power(s2,4)*
                 (4 + 2*Power(s1,2)*t1 + 10*Power(t1,3) + 
                   t1*(2 - 10*t2) + 3*t2 + 2*Power(t2,2) - 
                   s1*(9 + 13*Power(t1,2) + 2*t1*(-11 + t2) + 
                      2*t2) + Power(t1,2)*(-16 + 7*t2)) - 
                Power(s2,3)*(-1 + t1)*
                 (15 + 8*Power(t1,3) + Power(s1,2)*(1 + 9*t1) + 
                   s1*(15 - 23*Power(t1,2) + t1*(8 - 7*t2) - 9*t2) - 
                   6*t2 + 6*Power(t2,2) - 2*t1*(10 + t2) + 
                   Power(t1,2)*(-3 + 8*t2)) + 
                s2*(3*Power(t1,5) + 
                   Power(s1,2)*t1*(-24 + 39*t1 - 7*Power(t1,3)) + 
                   Power(t1,4)*(1 + 2*t2) + 
                   6*(1 - 3*t2 + 2*Power(t2,2)) - 
                   Power(t1,3)*(5 + 32*t2 + 2*Power(t2,2)) - 
                   2*t1*(3 + 17*t2 + 12*Power(t2,2)) + 
                   Power(t1,2)*(1 + 82*t2 + 22*Power(t2,2)) + 
                   s1*(3*Power(t1,5) + Power(t1,2)*(17 - 43*t2) + 
                      2*(-1 + t2) + 8*t1*(3 + 2*t2) + 
                      Power(t1,4)*(11 + 3*t2) + 
                      Power(t1,3)*(-53 + 6*t2))) + 
                Power(s2,2)*(-5 - 21*Power(t1,3) + 2*Power(t1,5) + 
                   Power(s1,2)*
                    (12 - 21*t1 - 12*Power(t1,2) + 13*Power(t1,3)) + 
                   3*Power(t1,4)*(-2 + t2) + 25*t2 + 12*Power(t2,2) + 
                   Power(t1,2)*(47 + 4*t2 + 6*Power(t2,2)) - 
                   t1*(17 + 32*t2 + 26*Power(t2,2)) - 
                   s1*(-1 + 15*Power(t1,4) + t1*(42 - 46*t2) + 
                      Power(t1,2)*(-54 + t2) + 21*t2 + 
                      Power(t1,3)*(-2 + 8*t2)))) - 
             Power(s,3)*(2 - Power(t1,5) - 20*t2 + Power(t1,4)*t2 + 
                26*Power(t2,2) + 
                Power(s2,3)*
                 (-23*Power(t1,2) + s1*(-4 + 9*t1) + 
                   t1*(17 - 16*t2) + 9*t2) + 
                Power(t1,3)*
                 (1 + 3*Power(s1,2) - 10*t2 - s1*(38 + 3*t2)) + 
                Power(t1,2)*
                 (3 + 4*Power(s1,2) + 34*t2 + 6*Power(t2,2) + 
                   s1*(37 + 5*t2)) - 
                t1*(11 + 2*t2 + 32*Power(t2,2) + s1*(-20 + 9*t2)) + 
                Power(s2,2)*
                 (4 + 25*Power(t1,3) + Power(s1,2)*(1 + 4*t1) + 
                   21*t2 + 12*Power(t2,2) - 2*t1*(13 + 21*t2) + 
                   Power(t1,2)*(7 + 30*t2) - 
                   s1*(2 + 37*Power(t1,2) + 8*t2 + t1*(-20 + 9*t2))) \
- s2*(4 + 3*Power(t1,4) + Power(s1,2)*t1*(-1 + 9*t1) - 17*t2 - 
                   30*Power(t2,2) + Power(t1,3)*(19 + 14*t2) - 
                   2*Power(t1,2)*(23 + 16*t2) + 
                   2*t1*(-3 + 22*t2 + 9*Power(t2,2)) + 
                   s1*(1 - 24*Power(t1,3) + 12*t2 + t1*(31 + 5*t2) - 
                      Power(t1,2)*(31 + 13*t2)))) + 
             Power(s,2)*(4*Power(t1,5) + 
                Power(s2,4)*(-6 - 16*Power(t1,2) + s1*(-5 + 7*t1) + 
                   t1*(22 - 9*t2) + 7*t2) + 
                t1*(-29 - 7*s1*(-3 + t2) + 49*t2 - 54*Power(t2,2)) - 
                Power(t1,3)*(-19 + 8*Power(s1,2) + s1*(57 - 5*t2) + 
                   41*t2 + 2*Power(t2,2)) + 
                8*(1 - 4*t2 + 3*Power(t2,2)) + 
                Power(t1,4)*(-18 - 3*Power(s1,2) + 8*t2 + 
                   s1*(22 + t2)) + 
                Power(t1,2)*(16 + 15*Power(s1,2) + 32*t2 + 
                   26*Power(t2,2) + s1*(-2 + 3*t2)) + 
                Power(s2,3)*(9 + 29*Power(t1,3) + 
                   Power(s1,2)*(1 + 5*t1) + 15*t2 + 8*Power(t2,2) + 
                   Power(t1,2)*(-25 + 24*t2) - t1*(13 + 33*t2) - 
                   s1*(11 + 38*Power(t1,2) + 7*t2 + t1*(-43 + 7*t2))) \
+ s2*(-1 + Power(t1,5) + Power(s1,2)*t1*(-20 - t1 + 15*Power(t1,2)) + 
                   t2 + 34*Power(t2,2) - 10*Power(t1,3)*(2 + t2) + 
                   Power(t1,4)*(-5 + 4*t2) + 
                   t1*(19 - 71*t2 - 54*Power(t2,2)) + 
                   6*Power(t1,2)*(1 + 9*t2 + 2*Power(t2,2)) - 
                   s1*(3 + 16*Power(t1,4) + t1*(26 - 34*t2) + 11*t2 - 
                      2*Power(t1,2)*(54 + t2) + 
                      Power(t1,3)*(41 + 11*t2))) + 
                Power(s2,2)*(-9 - 14*Power(t1,4) + 
                   Power(s1,2)*(5 + 8*t1 - 17*Power(t1,2)) + 
                   Power(t1,3)*(4 - 19*t2) + 17*t2 + 24*Power(t2,2) + 
                   5*Power(t1,2)*(11 + 6*t2) - 
                   2*t1*(18 + 13*t2 + 9*Power(t2,2)) + 
                   s1*(12 + 47*Power(t1,3) + 2*t1*(-21 + t2) - 21*t2 + 
                      Power(t1,2)*(-19 + 17*t2)))))/
           (s*(s - s2 + t1)*Power(s2 - s2*t1 + t1*(-1 + s + t1),2)*
             (s - s1 + t2)) + 
          (2*Power(s,6)*t1*(3*t1 + t2) + 
             4*(s2 - t1)*(-1 + t1)*
              Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
             Power(s,5)*(15*Power(t1,3) - 
                Power(t1,2)*(1 + 10*s1 + 28*s2 - 7*t2) + 
                2*t2*(s2 + t2) + 
                t1*(2 + 2*s2*(5 + s1 - 5*t2) - 5*t2 - 2*s1*t2)) + 
             Power(s,4)*(12*Power(t1,4) + 4*(1 - 3*t2)*t2 + 
                Power(t1,3)*(9 - 24*s1 + 8*t2) + 
                Power(t1,2)*
                 (-51 + 10*s1 + 4*Power(s1,2) - 5*t2 - 7*s1*t2) + 
                Power(s2,2)*
                 (4 + s1*(2 - 8*t1) - 38*t1 + 52*Power(t1,2) - 
                   8*t2 + 20*t1*t2) + 
                t1*(1 - 6*s1 - 2*t2 + 9*s1*t2 + 6*Power(t2,2)) + 
                s2*(2 - 55*Power(t1,3) + 
                   Power(t1,2)*(17 + 43*s1 - 28*t2) - 7*t2 + 
                   2*s1*t2 - 8*Power(t2,2) + 
                   t1*(1 - 15*s1 - 2*Power(s1,2) + 24*t2 + 8*s1*t2))) \
+ Power(s,3)*(2 + 3*Power(t1,5) + 
                Power(t1,3)*
                 (9*Power(s1,2) + s1*(4 - 8*t2) + 12*(-10 + t2)) - 
                20*t2 + 26*Power(t2,2) + 
                Power(t1,4)*(14 - 18*s1 + 3*t2) + 
                2*Power(s2,3)*
                 (-6 + 27*t1 - 24*Power(t1,2) + s1*(-3 + 6*t1) + 
                   6*t2 - 10*t1*t2) + 
                t1*(-25 + s1*(23 - 14*t2) + 25*t2 - 32*Power(t2,2)) + 
                Power(t1,2)*
                 (113 - 7*Power(s1,2) - 30*t2 + 6*Power(t2,2) + 
                   8*s1*(2 + 3*t2)) + 
                Power(s2,2)*
                 (6*Power(s1,2)*t1 + 75*Power(t1,3) - 
                   21*t1*(1 + 2*t2) + Power(t1,2)*(-43 + 42*t2) + 
                   4*(1 + 4*t2 + 3*Power(t2,2)) - 
                   s1*(7 + 69*Power(t1,2) + 6*t2 + 3*t1*(-17 + 4*t2))\
) - s2*(1 + 32*Power(t1,4) + Power(s1,2)*t1*(-9 + 17*t1) - 3*t2 - 
                   30*Power(t2,2) + 24*Power(t1,3)*(1 + t2) - 
                   Power(t1,2)*(157 + 22*t2) + 
                   3*t1*(24 + t2 + 6*Power(t2,2)) + 
                   s1*(2 - 71*Power(t1,3) + 
                      Power(t1,2)*(44 - 21*t2) + 9*t2 + t1*(7 + 16*t2)\
))) + s*(2*Power(s2,5)*(-1 + t1)*(2 + s1 - 2*t1 - t2) + 
                8*Power(-1 + t2,2) + 
                Power(t1,5)*(-19 - 2*s1 + Power(s1,2) + 4*t2) - 
                2*t1*(-1 + t2)*(-14 + s1 + 16*t2) + 
                Power(t1,4)*(58 - 15*t2 + 3*s1*(6 + t2)) - 
                Power(t1,3)*
                 (71 + 11*Power(s1,2) - 13*t2 + 6*Power(t2,2) + 
                   5*s1*(-1 + 2*t2)) - 
                Power(s2,3)*(-1 + t1)*
                 (26 + 8*Power(t1,3) + Power(s1,2)*(1 + 9*t1) + 
                   s1*(12 - 23*Power(t1,2) + t1*(9 - 7*t2) - 9*t2) - 
                   2*t2 + 6*Power(t2,2) - t1*(47 + 4*t2) + 
                   Power(t1,2)*(13 + 8*t2)) + 
                Power(t1,2)*
                 (52 + 8*Power(s1,2) - 46*t2 + 28*Power(t2,2) + 
                   s1*(-23 + 13*t2)) + 
                Power(s2,4)*
                 (2*Power(s1,2)*t1 + 10*Power(t1,3) - 
                   3*t1*(2 + 3*t2) + Power(t1,2)*(-12 + 7*t2) - 
                   s1*(8 - 21*t1 + 13*Power(t1,2) + 2*t2 + 
                      2*t1*t2) + 2*(4 + t2 + Power(t2,2))) + 
                s2*(-5*Power(t1,5) + 
                   Power(s1,2)*t1*
                    (-16 + 23*t1 + 8*Power(t1,2) - 7*Power(t1,3)) + 
                   Power(t1,4)*(76 - 10*t2) + 
                   t1*(-50 + 6*t2 - 24*Power(t2,2)) + 
                   Power(t1,3)*(-171 + 30*t2 - 2*Power(t2,2)) + 
                   4*(2 - 5*t2 + 3*Power(t2,2)) + 
                   2*Power(t1,2)*(71 - 3*t2 + 11*Power(t2,2)) + 
                   s1*(20*t1 + 3*Power(t1,5) + 
                      Power(t1,2)*(24 - 11*t2) + 2*(-1 + t2) + 
                      Power(t1,4)*(8 + 3*t2) - 
                      Power(t1,3)*(53 + 10*t2))) + 
                Power(s2,2)*(2*Power(t1,5) + 
                   Power(s1,2)*
                    (8 - 13*t1 - 16*Power(t1,2) + 13*Power(t1,3)) + 
                   Power(t1,3)*(-107 + t2) + 
                   Power(t1,4)*(14 + 3*t2) + 
                   t1*(-101 + 3*t2 - 26*Power(t2,2)) + 
                   Power(t1,2)*(178 - 15*t2 + 6*Power(t2,2)) + 
                   2*(7 + 4*t2 + 6*Power(t2,2)) + 
                   s1*(3 - 15*Power(t1,4) + Power(t1,3)*(7 - 8*t2) - 
                      13*t2 + Power(t1,2)*(46 + 7*t2) + 
                      t1*(-41 + 30*t2)))) - 
             Power(s,2)*(4*(-1 + s1)*Power(t1,5) + 
                t1*(-42 + s1*(23 - 9*t2) + 64*t2 - 54*Power(t2,2)) + 
                8*(1 - 4*t2 + 3*Power(t2,2)) + 
                Power(t1,3)*(-158 + 7*Power(s1,2) + 45*t2 - 
                   2*Power(t2,2) - 6*s1*(8 + 3*t2)) + 
                Power(t1,4)*(82 - 6*Power(s1,2) - 16*t2 + 
                   s1*(8 + 3*t2)) - 
                2*Power(s2,4)*
                 (6 + s1*(3 - 4*t1) + 11*Power(t1,2) - 4*t2 + 
                   t1*(-17 + 5*t2)) + 
                Power(t1,2)*(114 + 3*Power(s1,2) - 44*t2 + 
                   26*Power(t2,2) + s1*(-4 + 26*t2)) + 
                Power(s2,2)*(23 - 28*Power(t1,4) + 
                   Power(s1,2)*(1 + 17*t1 - 22*Power(t1,2)) + 3*t2 + 
                   24*Power(t2,2) - 4*Power(t1,3)*(5 + 6*t2) + 
                   Power(t1,2)*(173 + 29*t2) - 
                   t1*(148 + 7*t2 + 18*Power(t2,2)) + 
                   s1*(12 + 70*Power(t1,3) - 18*t2 - t1*(17 + 5*t2) + 
                      3*Power(t1,2)*(-22 + 7*t2))) + 
                Power(s2,3)*(18 + 6*Power(s1,2)*t1 + 45*Power(t1,3) + 
                   11*t2 + 8*Power(t2,2) - 8*t1*(3 + 4*t2) + 
                   Power(t1,2)*(-39 + 28*t2) - 
                   s1*(15 + 49*Power(t1,2) + 6*t2 + t1*(-57 + 8*t2))) \
+ s2*(5*Power(t1,5) + 2*Power(s1,2)*t1*(-2 - 12*t1 + 11*Power(t1,2)) + 
                   13*Power(t1,3)*(-17 + t2) + 
                   Power(t1,4)*(29 + 6*t2) + 
                   t1*(-107 + 15*t2 - 54*Power(t2,2)) + 
                   Power(t1,2)*(284 - 45*t2 + 12*Power(t2,2)) + 
                   2*(5 - 6*t2 + 17*Power(t2,2)) + 
                   s1*(-5 - 33*Power(t1,4) + Power(t1,3)*(7 - 16*t2) + 
                      8*t1*(-3 + t2) - 9*t2 + Power(t1,2)*(78 + 31*t2))\
)))/(s*(-1 + t1)*Power(s2 - s2*t1 + t1*(-1 + s + t1),2)*(-1 + t2)) + 
          (2*(5 + 12*Power(s,3) + 2*Power(s1,3) + 4*s2 - 
               4*Power(s2,2) - 40*Power(s2,3) - 11*t1 + 82*s2*t1 + 
               80*Power(s2,2)*t1 - 10*Power(t1,2) - 
               80*s2*Power(t1,2) + 
               Power(s,2)*(-92 + 60*s1 - 18*s2 + 166*t1 - 119*t2) + 
               16*s*(s - s2 + t1)*(-1 + s1 + t1 - t2) - 38*t2 - 
               65*s2*t2 - 18*Power(s2,2)*t2 + 75*t1*t2 + 40*s2*t1*t2 - 
               2*Power(s1,2)*(-4 + 2*s2 - 2*t1 + t2) + 
               s1*(9 - 10*Power(s2,2) - 18*t1 - 40*Power(t1,2) - 
                  4*t2 + s2*(65 + 8*t1 + 2*t2)) + 
               s*(1 + 16*Power(s1,2) + 46*Power(s2,2) - 147*t1 + 
                  166*Power(t1,2) + 113*t2 - 130*t1*t2 - 
                  s1*(61 + 84*s2 - 46*t1 + 15*t2) + 
                  s2*(132 - 234*t1 + 137*t2)) - 
               ((-5 + 2*Power(s,2) + 2*Power(s1,2) - 9*s2 + 
                    3*Power(s2,2) + s1*(-5 + 5*s2 - t2) + 3*t2 - 
                    s2*t2 + s*(5 - 4*s1 - 5*s2 + t2))*
                  (-2 + t1 + Power(t1,2) + s1*Power(t1,2) + 
                    s2*(-1 + t1)*(-1 + 2*t1 - t2) + 2*t2 - 3*t1*t2 + 
                    Power(s2,2)*(2 - s1 - 2*t1 + t2) + 
                    Power(s,2)*(-t1 + t2) + 
                    s*(1 - Power(t1,2) + 
                       s2*(-2 + s1 + 3*t1 - 2*t2) - 3*t2 + 
                       t1*(1 + s1 + t2))))/
                (s2*(-1 + t1) - t1*(-1 + s + t1)) - 
               8*(1 + Power(s,3) - 2*Power(s2,3) - 2*t1 + 
                  2*Power(t1,2) - 2*s1*Power(t1,2) + 
                  Power(s,2)*(-12 + 11*s1 - 2*s2 + 17*t1 - 14*t2) + 
                  Power(s2,2)*(1 + 4*t1 - t2) - t2 + 2*t1*t2 + 
                  2*s2*(-1 + s1 + t1 - 2*Power(t1,2) - t2 + t1*t2) + 
                  s*(-1 + Power(s1,2) + 3*Power(s2,2) - 13*t1 + 
                     17*Power(t1,2) - 3*s2*(-5 + 7*t1 - 5*t2) + 
                     5*t2 - 15*t1*t2 - s1*(3 + 13*s2 - 11*t1 + t2))) + 
               4*(5 + 8*Power(s,3) - 14*s2 + 5*Power(s2,2) - 
                  19*Power(s2,3) - 11*t1 + 26*s2*t1 + 
                  36*Power(s2,2)*t1 + 8*Power(t1,2) - 
                  36*s2*Power(t1,2) + Power(s1,2)*(1 - 2*s2 + t1) + 
                  Power(s,2)*(-55 + 45*s1 - 16*s2 + 95*t1 - 72*t2) - 
                  11*t2 - 23*s2*t2 - 8*Power(s2,2)*t2 + 24*t1*t2 + 
                  18*s2*t1*t2 + 
                  s1*(1 - 4*Power(s2,2) - 4*t1 - 18*Power(t1,2) - 
                     t2 + s2*(23 + 2*t1 + t2)) + 
                  s*(-5 + 8*Power(s1,2) + 27*Power(s2,2) - 70*t1 + 
                     95*Power(t1,2) + 45*t2 - 79*t1*t2 - 
                     s1*(26 + 58*s2 - 42*t1 + 7*t2) + 
                     s2*(77 - 129*t1 + 80*t2))) - 
               2*(13 + 17*Power(s,3) + 2*Power(s1,3) - 32*s2 - 
                  6*Power(s2,2) - 47*Power(s2,3) - 20*t1 + 86*s2*t1 + 
                  96*Power(s2,2)*t1 + 2*Power(t1,2) - 
                  96*s2*Power(t1,2) + 
                  s1*(4 - 6*Power(s2,2) - 16*t1 - 48*Power(t1,2) + 
                     8*s2*(7 + t1)) + 
                  Power(s,2)*
                   (-121 + 87*s1 - 28*s2 + 213*t1 - 156*t2) - 39*t2 - 
                  62*s2*t2 - 23*Power(s2,2)*t2 + 78*t1*t2 + 
                  48*s2*t1*t2 - Power(s1,2)*(-1 + s2 - 4*t1 + 2*t2) + 
                  s*(-1 + 17*Power(s1,2) + 58*Power(s2,2) - 171*t1 + 
                     213*Power(t1,2) + 123*t2 - 171*t1*t2 - 
                     s1*(61 + 124*s2 - 71*t1 + 16*t2) + 
                     s2*(182 - 299*t1 + 179*t2))) + 
               ((3 + s - s1 - s2)*
                  (Power(s,5)*Power(t1 - t2,2) + 
                    2*(s2 - t1)*(-1 + t1)*
                     Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                    2*Power(s,4)*
                     (-2*Power(t1,2) + Power(t1,3) + t2 + 5*t1*t2 - 
                       2*Power(t1,2)*t2 - 3*Power(t2,2) + 
                       t1*Power(t2,2) - 
                       s1*(t1 + Power(t1,2) - s2*t2 - t1*t2) + 
                       s2*(t1 - 2*Power(t1,2) + 4*t1*t2 - 
                        2*t2*(1 + t2))) + 
                    Power(s,3)*
                     (1 + Power(t1,4) - 10*t2 + 13*Power(t2,2) - 
                       2*Power(t1,3)*(4 + 2*s1 + t2) + 
                       Power(t1,2)*
                       (7 + Power(s1,2) + 20*t2 + 4*s1*t2 + 
                       Power(t2,2)) + 
                       Power(s2,2)*
                       (2 + Power(s1,2) + 7*Power(t1,2) + 
                       2*s1*(-1 + t1 - 3*t2) + 10*t2 + 
                       6*Power(t2,2) - 4*t1*(2 + 3*t2)) - 
                       2*t1*
                       (1 + 6*t2 + 6*Power(t2,2) + s1*(-5 + 4*t2)) - 
                       2*s2*
                        (1 + 3*Power(t1,3) + 3*(-2 + s1)*t2 - 
                        7*Power(t2,2) - 
                        3*Power(t1,2)*(2 + s1 + 2*t2) + 
                        t1*(1 + 17*t2 + 3*Power(t2,2)))) + 
                    s*(4 - 10*t1 + 13*Power(t1,2) - 12*Power(t1,3) + 
                       5*Power(t1,4) + 
                       Power(s1,2)*Power(s2 - t1,2)*
                       (2 + 2*s2 + Power(s2,2) - 4*t1 + Power(t1,2)) \
- 8*t2 + 24*t1*t2 - 6*Power(t1,2)*t2 - 10*Power(t1,3)*t2 + 
                       4*Power(t2,2) - 14*t1*Power(t2,2) + 
                       9*Power(t1,2)*Power(t2,2) + 
                       Power(s2,4)*
                       (2 + 2*Power(t1,2) + 2*t2 + Power(t2,2) - 
                       2*t1*(2 + t2)) - 
                       2*Power(s2,3)*(-1 + t1)*
                       (2*Power(t1,2) + Power(1 + t2,2) - 
                       t1*(3 + 2*t2)) + 
                       2*s1*(s2 - t1)*
                       (-2*Power(t1,3) + t1*(11 - 7*t2) + 
                       Power(s2,3)*(-1 + t1 - t2) + 2*(-1 + t2) + 
                       Power(t1,2)*(-7 + 4*t2) - 
                       s2*(1 - 2*Power(t1,2) + Power(t1,3) + 
                       3*t2 - 6*t1*t2) + 
                       Power(s2,2)*(-1 + t1 - 2*t2 + t1*t2)) + 
                       2*s2*
                       (Power(t1,4) + 2*(-1 + t2)*t2 - 
                       3*Power(t1,3)*(1 + 2*t2) - 
                       t1*(1 + 12*t2 + 3*Power(t2,2)) + 
                       Power(t1,2)*(3 + 20*t2 + 3*Power(t2,2))) + 
                       Power(s2,2)*
                        (-3 + 2*Power(t1,4) + 14*t2 + 
                        5*Power(t2,2) - 2*Power(t1,3)*(4 + t2) + 
                        Power(t1,2)*(7 + 18*t2 + Power(t2,2)) - 
                        2*t1*(-1 + 15*t2 + 5*Power(t2,2)))) - 
                    2*Power(s,2)*
                     (2 + (2 + s1)*Power(t1,4) - 8*t2 + 
                       6*Power(t2,2) + 
                       t1*(-4 + s1*(7 - 5*t2) + 6*t2 - 
                       11*Power(t2,2)) - 
                       Power(t1,3)*
                        (6 + Power(s1,2) + 5*t2 + s1*(3 + t2)) + 
                       Power(s2,3)*
                        (Power(s1,2) + 3*Power(t1,2) + 
                        s1*(-2 + 2*t1 - 3*t2) + 2*Power(1 + t2,2) - 
                        t1*(5 + 4*t2)) + 
                       Power(t1,2)*
                        (Power(s1,2) + 2*s1*(-5 + 4*t2) + 
                        3*(2 + 4*t2 + Power(t2,2))) + 
                       Power(s2,2)*
                        (-(Power(s1,2)*(-1 + t1)) - 4*Power(t1,3) + 
                        t2*(9 + 5*t2) + Power(t1,2)*(9 + 6*t2) + 
                        s1*
                       (-1 + t1 + Power(t1,2) - 5*t2 + 3*t1*t2) - 
                        t1*(5 + 16*t2 + 3*Power(t2,2))) + 
                       s2*(-2 + Power(s1,2)*(-2 + t1)*t1 + 
                        Power(t1,4) + 4*t2 + 7*Power(t2,2) - 
                        2*Power(t1,3)*(3 + t2) + 
                        t1*(2 - 27*t2 - 10*Power(t2,2)) + 
                        Power(t1,2)*(5 + 21*t2 + Power(t2,2)) + 
                        s1*(-1 - 4*Power(t1,3) - t2 + 
                        Power(t1,2)*(4 + t2) + t1*(5 + 3*t2))))))/
                (s*Power(s2 - s2*t1 + t1*(-1 + s + t1),2))))/
           ((-1 + s1)*(-s + s1 - t2)) + 
          (2*(-40 - 12*Power(s,3) + 2*Power(s1,3) + 16*s2 + 
               61*Power(s2,2) + 57*t1 - 78*s2*t1 - 80*Power(s2,2)*t1 + 
               10*Power(t1,2) + 80*s2*Power(t1,2) - 
               Power(s,2)*(-62 + 10*s1 - 2*s2 + 166*t1 - 121*t2) - 
               16*s*(s - s2 + t1)*(t1 - t2) + 49*t2 + 18*s2*t2 + 
               42*Power(s2,2)*t2 - 75*t1*t2 - 40*s2*t1*t2 - 
               2*Power(s1,2)*(-5 + 2*s2 - 2*t1 + t2) - 
               2*s1*(3 + 15*Power(s2,2) + 18*t1 - 20*Power(t1,2) + 
                  s2*(-30 + 6*t1 - t2) + 4*t2) - 
               s*(-71 + 12*Power(s1,2) - 10*Power(s2,2) - 82*t1 + 
                  166*Power(t1,2) + 90*t2 - 130*t1*t2 - 
                  s1*(-58 + 62*s2 + 4*t1 + 21*t2) + 
                  s2*(175 - 246*t1 + 163*t2)) - 
               ((-5 + 3*Power(s,2) + 2*Power(s1,2) - 9*s2 + 
                    3*Power(s2,2) + s1*(-5 + 5*s2 - t2) + 3*t2 - 
                    s2*t2 + s*(7 - 5*s1 - 6*s2 + t2))*
                  (-2 + t1 + Power(t1,2) + s1*Power(t1,2) + 
                    s2*(-1 + t1)*(-1 + 2*t1 - t2) + 2*t2 - 3*t1*t2 + 
                    Power(s2,2)*(2 - s1 - 2*t1 + t2) + 
                    Power(s,2)*(-t1 + t2) + 
                    s*(1 - Power(t1,2) + s2*(-2 + s1 + 3*t1 - 2*t2) - 
                       3*t2 + t1*(1 + s1 + t2))))/
                (s2*(-1 + t1) - t1*(-1 + s + t1)) + 
               8*(2 + Power(s,3) - 5*t1 + 3*s1*t1 + 2*Power(t1,2) - 
                  2*s1*Power(t1,2) + 
                  Power(s2,2)*(-3 + s1 + 4*t1 - 2*t2) - 2*t2 + 
                  2*t1*t2 - Power(s,2)*(2 + s2 - 17*t1 + 14*t2) + 
                  s2*(3 + s1*(-3 + t1) + t1 - 4*Power(t1,2) + 
                     2*t1*t2) - 
                  s*(7 + t1 - 17*Power(t1,2) + 
                     s2*(-7 + 21*t1 - 16*t2) - 4*t2 + 15*t1*t2 + 
                     s1*(-4 + 2*s2 + t1 + t2))) - 
               4*(19 + 8*Power(s,3) + 17*s2 - 28*Power(s2,2) + 
                  Power(s2,3) + Power(s1,2)*(-1 + 2*s2 - t1) - 38*t1 + 
                  21*s2*t1 + 36*Power(s2,2)*t1 + 8*Power(t1,2) - 
                  36*s2*Power(t1,2) - 20*t2 - 3*s2*t2 - 
                  19*Power(s2,2)*t2 + 24*t1*t2 + 18*s2*t1*t2 + 
                  s1*(1 + 13*Power(s2,2) + 25*t1 - 18*Power(t1,2) + 
                     s2*(-28 + 7*t1 - t2) + t2) - 
                  Power(s,2)*(19 + 7*s2 - 95*t1 + 72*t2) + 
                  s*(-52 + 2*Power(s1,2) - 2*Power(s2,2) - 19*t1 + 
                     95*Power(t1,2) + 36*t2 - 79*t1*t2 - 
                     s1*(-31 + 21*s2 + 7*t1 + 9*t2) + 
                     s2*(65 - 131*t1 + 91*t2))) + 
               2*(49 + 21*Power(s,3) - 2*Power(s1,3) + 27*s2 - 
                  65*Power(s2,2) - Power(s2,3) - 84*t1 + 78*s2*t1 + 
                  96*Power(s2,2)*t1 + 2*Power(t1,2) - 
                  96*s2*Power(t1,2) - 61*t2 - 22*s2*t2 - 
                  49*Power(s2,2)*t2 + 78*t1*t2 + 48*s2*t1*t2 + 
                  Power(s1,2)*(-3 + s2 - 4*t1 + 2*t2) + 
                  s1*(9 + 30*Power(s2,2) + 56*t1 - 48*Power(t1,2) + 
                     s2*(-61 + 16*t1) + 4*t2) - 
                  Power(s,2)*(56 + 2*s1 + 20*s2 - 213*t1 + 154*t2) + 
                  s*(-121 + 15*Power(s1,2) - 76*t1 + 213*Power(t1,2) + 
                     107*t2 - 171*t1*t2 - 
                     s1*(-66 + 51*s2 + 12*t1 + 28*t2) + 
                     s2*(176 - 309*t1 + 203*t2))) + 
               ((3 + s - s1 - s2)*
                  (Power(s,5)*Power(t1 - t2,2) + 
                    2*(s2 - t1)*(-1 + t1)*
                     Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                    2*Power(s,4)*
                     (-2*Power(t1,2) + Power(t1,3) + t2 + 5*t1*t2 - 
                       2*Power(t1,2)*t2 - 3*Power(t2,2) + 
                       t1*Power(t2,2) - 
                       s1*(t1 + Power(t1,2) - s2*t2 - t1*t2) + 
                       s2*(t1 - 2*Power(t1,2) + 4*t1*t2 - 
                        2*t2*(1 + t2))) + 
                    Power(s,3)*
                     (1 + Power(t1,4) - 10*t2 + 13*Power(t2,2) - 
                       2*Power(t1,3)*(4 + 2*s1 + t2) + 
                       Power(t1,2)*
                        (7 + Power(s1,2) + 20*t2 + 4*s1*t2 + 
                        Power(t2,2)) + 
                       Power(s2,2)*
                        (2 + Power(s1,2) + 7*Power(t1,2) + 
                        2*s1*(-1 + t1 - 3*t2) + 10*t2 + 
                        6*Power(t2,2) - 4*t1*(2 + 3*t2)) - 
                       2*t1*
                        (1 + 6*t2 + 6*Power(t2,2) + s1*(-5 + 4*t2)) - 
                       2*s2*(1 + 3*Power(t1,3) + 3*(-2 + s1)*t2 - 
                        7*Power(t2,2) - 
                        3*Power(t1,2)*(2 + s1 + 2*t2) + 
                        t1*(1 + 17*t2 + 3*Power(t2,2)))) + 
                    s*(4 - 10*t1 + 13*Power(t1,2) - 12*Power(t1,3) + 
                       5*Power(t1,4) + 
                       Power(s1,2)*Power(s2 - t1,2)*
                        (2 + 2*s2 + Power(s2,2) - 4*t1 + Power(t1,2)) \
- 8*t2 + 24*t1*t2 - 6*Power(t1,2)*t2 - 10*Power(t1,3)*t2 + 
                       4*Power(t2,2) - 14*t1*Power(t2,2) + 
                       9*Power(t1,2)*Power(t2,2) + 
                       Power(s2,4)*
                        (2 + 2*Power(t1,2) + 2*t2 + Power(t2,2) - 
                        2*t1*(2 + t2)) - 
                       2*Power(s2,3)*(-1 + t1)*
                        (2*Power(t1,2) + Power(1 + t2,2) - 
                        t1*(3 + 2*t2)) + 
                       2*s1*(s2 - t1)*
                        (-2*Power(t1,3) + t1*(11 - 7*t2) + 
                        Power(s2,3)*(-1 + t1 - t2) + 2*(-1 + t2) + 
                        Power(t1,2)*(-7 + 4*t2) - 
                        s2*(1 - 2*Power(t1,2) + Power(t1,3) + 
                       3*t2 - 6*t1*t2) + 
                        Power(s2,2)*(-1 + t1 - 2*t2 + t1*t2)) + 
                       2*s2*
                        (Power(t1,4) + 2*(-1 + t2)*t2 - 
                        3*Power(t1,3)*(1 + 2*t2) - 
                        t1*(1 + 12*t2 + 3*Power(t2,2)) + 
                        Power(t1,2)*(3 + 20*t2 + 3*Power(t2,2))) + 
                       Power(s2,2)*
                        (-3 + 2*Power(t1,4) + 14*t2 + 5*Power(t2,2) - 
                        2*Power(t1,3)*(4 + t2) + 
                        Power(t1,2)*(7 + 18*t2 + Power(t2,2)) - 
                        2*t1*(-1 + 15*t2 + 5*Power(t2,2)))) - 
                    2*Power(s,2)*
                     (2 + (2 + s1)*Power(t1,4) - 8*t2 + 
                       6*Power(t2,2) + 
                       t1*(-4 + s1*(7 - 5*t2) + 6*t2 - 
                       11*Power(t2,2)) - 
                       Power(t1,3)*
                        (6 + Power(s1,2) + 5*t2 + s1*(3 + t2)) + 
                       Power(s2,3)*
                        (Power(s1,2) + 3*Power(t1,2) + 
                        s1*(-2 + 2*t1 - 3*t2) + 2*Power(1 + t2,2) - 
                        t1*(5 + 4*t2)) + 
                       Power(t1,2)*
                        (Power(s1,2) + 2*s1*(-5 + 4*t2) + 
                        3*(2 + 4*t2 + Power(t2,2))) + 
                       Power(s2,2)*
                        (-(Power(s1,2)*(-1 + t1)) - 4*Power(t1,3) + 
                        t2*(9 + 5*t2) + Power(t1,2)*(9 + 6*t2) + 
                        s1*(-1 + t1 + Power(t1,2) - 5*t2 + 3*t1*t2) - 
                        t1*(5 + 16*t2 + 3*Power(t2,2))) + 
                       s2*(-2 + Power(s1,2)*(-2 + t1)*t1 + 
                        Power(t1,4) + 4*t2 + 7*Power(t2,2) - 
                        2*Power(t1,3)*(3 + t2) + 
                        t1*(2 - 27*t2 - 10*Power(t2,2)) + 
                        Power(t1,2)*(5 + 21*t2 + Power(t2,2)) + 
                        s1*(-1 - 4*Power(t1,3) - t2 + 
                        Power(t1,2)*(4 + t2) + t1*(5 + 3*t2))))))/
                (s*Power(s2 - s2*t1 + t1*(-1 + s + t1),2))))/
           ((-1 + s1)*(-1 + t2))) + 
       (((-64*(-13 + Power(s,2)*s2 + Power(s2,3) + s2*(3 - 4*t1) - 
                  Power(s2,2)*(-3 + t1) + t1 + 2*Power(t1,2) + 
                  s*(4 - 2*Power(s2,2) + s2*(-4 + t1) + t1))*
                (-2 + t1 + Power(t1,2) + s1*Power(t1,2) + 
                  s2*(-1 + t1)*(-1 + 2*t1 - t2) + 2*t2 - 3*t1*t2 + 
                  Power(s2,2)*(2 - s1 - 2*t1 + t2) + 
                  Power(s,2)*(-t1 + t2) + 
                  s*(1 - Power(t1,2) + s2*(-2 + s1 + 3*t1 - 2*t2) - 
                     3*t2 + t1*(1 + s1 + t2))))/
              (s2*(-1 + t1) - t1*(-1 + s + t1)) + 
             128*(19 + Power(s,3)*(-4 + s2) + 15*s2 + 11*Power(s2,2) + 
                Power(s2,3) - Power(s2,4) + 2*t1 - 7*s2*t1 - 
                2*Power(s2,2)*t1 + Power(s2,3)*t1 - 3*Power(t1,2) - 
                s1*(-7 + Power(s2,3) - Power(s2,2)*(-1 + t1) + 2*t1 + 
                   s2*t1 - Power(t1,2)) + 
                Power(s,2)*(17 + s1 - s1*s2 - 3*Power(s2,2) - 6*t1 + 
                   s2*(8 + t1) - t2) - 9*t2 - s2*t2 - 
                2*Power(s2,2)*t2 + 3*t1*t2 + s2*t1*t2 + 
                s*(-25 + 3*Power(s2,3) + 3*t1 - 2*Power(t1,2) - 
                   Power(s2,2)*(5 + 2*t1) + 
                   s1*(-7 + 2*Power(s2,2) - s2*(-2 + t1) + 2*t1) + 
                   7*t2 - t1*t2 + s2*(-25 + 11*t1 + 2*t2))))/
           (-s + s2 - t1) + (64*
             (((17 + Power(s,3) - Power(s2,3) + s2*(3 - 2*t1) + 
                    Power(s2,2)*(-1 + t1) - 11*t1 + 2*Power(t1,2) + 
                    Power(s,2)*(-3*s2 + t1) + 
                    s*(-11 + s2 + 3*Power(s2,2) + 3*t1 - 2*s2*t1))*
                  (-2 + t1 + Power(t1,2) + s1*Power(t1,2) + 
                    s2*(-1 + t1)*(-1 + 2*t1 - t2) + 2*t2 - 3*t1*t2 + 
                    Power(s2,2)*(2 - s1 - 2*t1 + t2) + 
                    Power(s,2)*(-t1 + t2) + 
                    s*(1 - Power(t1,2) + 
                       s2*(-2 + s1 + 3*t1 - 2*t2) - 3*t2 + 
                       t1*(1 + s1 + t2))))/
                (s2*(-1 + t1) - t1*(-1 + s + t1)) - 
               2*(-17 + Power(s,4) - 25*s2 - 13*Power(s2,2) - 
                  3*Power(s2,3) + Power(s2,4) + 18*t1 + 13*s2*t1 + 
                  6*Power(s2,2)*t1 - Power(s2,3)*t1 - 3*Power(t1,2) + 
                  Power(s,3)*(-s1 - 4*s2 + t1) + 
                  s1*(5 + Power(s2,3) - Power(s2,2)*(-3 + t1) - 
                     s2*(-2 + t1) - 2*t1 + Power(t1,2)) + 
                  Power(s,2)*
                   (-12 + 6*Power(s2,2) + s1*(4 + 3*s2 - t1) + 3*t1 - 
                     3*s2*(2 + t1) - 4*t2) - 15*t2 - 3*s2*t2 - 
                  2*Power(s2,2)*t2 + 3*t1*t2 + s2*t1*t2 + 
                  s*(31 - 4*Power(s2,3) - 16*t1 - 2*Power(t1,2) + 
                     3*Power(s2,2)*(3 + t1) + 
                     s1*(-8 - 5*s2 - 3*Power(s2,2) + 2*t1 + 2*s2*t1) + 
                     15*t2 - t1*t2 + s2*(29 - 8*t1 + 4*t2)))))/(-1 + t1) \
+ (32*(-1 + s2)*((-2*(15 + Power(s,3) + 13*s2 - 5*Power(s2,2) - 
                    2*Power(s2,3) - 
                    Power(s,2)*(-3 + s1 + 4*s2 - 2*t1) + 3*t1 + 
                    2*s2*t1 + 2*Power(s2,2)*t1 + 
                    s1*(7 - Power(s2,2) + s2*(-5 + t1) + 4*t1) - 
                    9*t2 + s2*t2 - 2*t1*t2 + 
                    s*(-13 + s1 + 3*s1*s2 + 5*Power(s2,2) - 2*t1 - 
                       s1*t1 - 3*s2*t1 + Power(t1,2) + t2))*
                  (2 - t1 - Power(t1,2) - s1*Power(t1,2) + 
                    Power(s,2)*(t1 - t2) - 
                    s2*(-1 + t1)*(-1 + 2*t1 - t2) + 
                    Power(s2,2)*(-2 + s1 + 2*t1 - t2) - 2*t2 + 
                    3*t1*t2 - 
                    s*(1 - Power(t1,2) + 
                       s2*(-2 + s1 + 3*t1 - 2*t2) - 3*t2 + 
                       t1*(1 + s1 + t2))))/
                (s2*(-1 + t1) - t1*(-1 + s + t1)) + 
               4*(1 + Power(s,4) - Power(s1,2)*(-5 + s2) + 14*s2 + 
                  11*Power(s2,2) - Power(s2,3) - Power(s2,4) + 
                  10*t1 + 6*s2*t1 + Power(s2,3)*t1 - Power(t1,2) + 
                  s1*(7 + 15*s2 - Power(s2,3) + 
                     Power(s2,2)*(-3 + t1) + Power(t1,2) + 
                     t1*(-4 + t2) - 9*t2) - 4*t2 - 19*s2*t2 - 
                  2*Power(s2,2)*t2 + s2*t1*t2 + 3*Power(t2,2) + 
                  s2*Power(t2,2) + Power(s,3)*(-8 - s1 + t1 + t2) - 
                  Power(s,2)*
                   (-12 + 4*Power(s2,2) + 3*t1 + 
                     s1*(-13 + s2 + t1) + 12*t2 - t1*t2 + 
                     s2*(-13 - t1 + t2)) + 
                  s*(-2 + Power(s1,2)*(-3 + s2) + 4*Power(s2,3) - 
                     6*t1 - Power(t1,2) - Power(s2,2)*(4 + 3*t1) + 
                     18*t2 - 4*t1*t2 - 2*Power(t2,2) + 
                     s2*(-27 + t1 + Power(t1,2) + 14*t2 - t1*t2) + 
                     s1*(-23 + 3*Power(s2,2) + 5*t1 + 6*t2 - 
                        s2*(8 + t2)))) + 
               (2*(6 + s - 4*s2 + s*s2 - Power(s2,2) + 2*t1 + s2*t1)*
                  (Power(s,5)*Power(t1 - t2,2) + 
                    2*(s2 - t1)*(-1 + t1)*
                     Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                    2*Power(s,4)*
                     (-2*Power(t1,2) + Power(t1,3) + t2 + 5*t1*t2 - 
                       2*Power(t1,2)*t2 - 3*Power(t2,2) + 
                       t1*Power(t2,2) - 
                       s1*(t1 + Power(t1,2) - s2*t2 - t1*t2) + 
                       s2*(t1 - 2*Power(t1,2) + 4*t1*t2 - 
                        2*t2*(1 + t2))) + 
                    Power(s,3)*
                     (1 + Power(t1,4) - 10*t2 + 13*Power(t2,2) - 
                       2*Power(t1,3)*(4 + 2*s1 + t2) + 
                       Power(t1,2)*
                       (7 + Power(s1,2) + 20*t2 + 4*s1*t2 + 
                       Power(t2,2)) + 
                       Power(s2,2)*
                       (2 + Power(s1,2) + 7*Power(t1,2) + 
                       2*s1*(-1 + t1 - 3*t2) + 10*t2 + 
                       6*Power(t2,2) - 4*t1*(2 + 3*t2)) - 
                       2*t1*
                       (1 + 6*t2 + 6*Power(t2,2) + s1*(-5 + 4*t2)) - 
                       2*s2*
                        (1 + 3*Power(t1,3) + 3*(-2 + s1)*t2 - 
                        7*Power(t2,2) - 
                        3*Power(t1,2)*(2 + s1 + 2*t2) + 
                        t1*(1 + 17*t2 + 3*Power(t2,2)))) + 
                    s*(4 - 10*t1 + 13*Power(t1,2) - 12*Power(t1,3) + 
                       5*Power(t1,4) + 
                       Power(s1,2)*Power(s2 - t1,2)*
                       (2 + 2*s2 + Power(s2,2) - 4*t1 + Power(t1,2)) \
- 8*t2 + 24*t1*t2 - 6*Power(t1,2)*t2 - 10*Power(t1,3)*t2 + 
                       4*Power(t2,2) - 14*t1*Power(t2,2) + 
                       9*Power(t1,2)*Power(t2,2) + 
                       Power(s2,4)*
                       (2 + 2*Power(t1,2) + 2*t2 + Power(t2,2) - 
                       2*t1*(2 + t2)) - 
                       2*Power(s2,3)*(-1 + t1)*
                       (2*Power(t1,2) + Power(1 + t2,2) - 
                       t1*(3 + 2*t2)) + 
                       2*s1*(s2 - t1)*
                       (-2*Power(t1,3) + t1*(11 - 7*t2) + 
                       Power(s2,3)*(-1 + t1 - t2) + 2*(-1 + t2) + 
                       Power(t1,2)*(-7 + 4*t2) - 
                       s2*(1 - 2*Power(t1,2) + Power(t1,3) + 
                       3*t2 - 6*t1*t2) + 
                       Power(s2,2)*(-1 + t1 - 2*t2 + t1*t2)) + 
                       2*s2*
                       (Power(t1,4) + 2*(-1 + t2)*t2 - 
                       3*Power(t1,3)*(1 + 2*t2) - 
                       t1*(1 + 12*t2 + 3*Power(t2,2)) + 
                       Power(t1,2)*(3 + 20*t2 + 3*Power(t2,2))) + 
                       Power(s2,2)*
                        (-3 + 2*Power(t1,4) + 14*t2 + 
                        5*Power(t2,2) - 2*Power(t1,3)*(4 + t2) + 
                        Power(t1,2)*(7 + 18*t2 + Power(t2,2)) - 
                        2*t1*(-1 + 15*t2 + 5*Power(t2,2)))) - 
                    2*Power(s,2)*
                     (2 + (2 + s1)*Power(t1,4) - 8*t2 + 
                       6*Power(t2,2) + 
                       t1*(-4 + s1*(7 - 5*t2) + 6*t2 - 
                       11*Power(t2,2)) - 
                       Power(t1,3)*
                        (6 + Power(s1,2) + 5*t2 + s1*(3 + t2)) + 
                       Power(s2,3)*
                        (Power(s1,2) + 3*Power(t1,2) + 
                        s1*(-2 + 2*t1 - 3*t2) + 2*Power(1 + t2,2) - 
                        t1*(5 + 4*t2)) + 
                       Power(t1,2)*
                        (Power(s1,2) + 2*s1*(-5 + 4*t2) + 
                        3*(2 + 4*t2 + Power(t2,2))) + 
                       Power(s2,2)*
                        (-(Power(s1,2)*(-1 + t1)) - 4*Power(t1,3) + 
                        t2*(9 + 5*t2) + Power(t1,2)*(9 + 6*t2) + 
                        s1*
                       (-1 + t1 + Power(t1,2) - 5*t2 + 3*t1*t2) - 
                        t1*(5 + 16*t2 + 3*Power(t2,2))) + 
                       s2*(-2 + Power(s1,2)*(-2 + t1)*t1 + 
                        Power(t1,4) + 4*t2 + 7*Power(t2,2) - 
                        2*Power(t1,3)*(3 + t2) + 
                        t1*(2 - 27*t2 - 10*Power(t2,2)) + 
                        Power(t1,2)*(5 + 21*t2 + Power(t2,2)) + 
                        s1*(-1 - 4*Power(t1,3) - t2 + 
                        Power(t1,2)*(4 + t2) + t1*(5 + 3*t2))))))/
                (s*Power(s2 - s2*t1 + t1*(-1 + s + t1),2))))/
           ((s - s2 + t1)*(s - s1 + t2)) + 
          (32*(-1 + s2)*((-2*
                  (16 + 2*Power(s,3) + 19*s2 - 2*Power(s2,2) - 
                    2*Power(s2,3) - 5*t1 - 4*s2*t1 + 
                    2*Power(s2,2)*t1 + 
                    Power(s,2)*(3 - 2*s1 - 6*s2 + t1) - 
                    s1*(-7 + 3*s2 + Power(s2,2) + 2*t1 - s2*t1) - 
                    5*t2 + s2*t2 + 2*t1*t2 + 
                    s*(-22 + 6*Power(s2,2) + 3*s2*(-1 + s1 - t1) + 
                       2*t1 - s1*t1 + 2*t2))*
                  (2 - t1 - Power(t1,2) - s1*Power(t1,2) + 
                    Power(s,2)*(t1 - t2) - 
                    s2*(-1 + t1)*(-1 + 2*t1 - t2) + 
                    Power(s2,2)*(-2 + s1 + 2*t1 - t2) - 2*t2 + 
                    3*t1*t2 - 
                    s*(1 - Power(t1,2) + 
                       s2*(-2 + s1 + 3*t1 - 2*t2) - 3*t2 + 
                       t1*(1 + s1 + t2))))/
                (s2*(-1 + t1) - t1*(-1 + s + t1)) - 
               4*(-12 + Power(s,4) + Power(s1,2)*(-1 + s2) - 17*s2 - 
                  17*Power(s2,2) + Power(s2,4) + 4*t1 + 4*s2*t1 + 
                  2*Power(s2,2)*t1 - Power(s2,3)*t1 + 9*t2 + 
                  2*Power(s2,2)*t2 + t1*t2 + 3*Power(t2,2) + 
                  s2*Power(t2,2) + Power(s,3)*(-2*s1 - 4*s2 + t2) + 
                  Power(s,2)*
                   (-11 + Power(s1,2) + 6*Power(s2,2) - t1 + 
                     s1*(2 + 5*s2 - t2) - t2 - s2*(3 + t1 + 2*t2)) + 
                  s1*(-4 + Power(s2,3) - Power(s2,2)*(-2 + t1) + 
                     t1*(-3 + t2) - 3*t2 + s2*(t1 - 2*(2 + t2))) + 
                  s*(21 - 4*Power(s2,3) - Power(s1,2)*(1 + s2) + t1 - 
                     Power(t1,2) - 5*t2 - 2*t1*t2 - 2*Power(t2,2) - 
                     s2*(-30 + t1 + t2) + 
                     Power(s2,2)*(3 + 2*t1 + t2) + 
                     s1*(4 - 4*Power(s2,2) + t1 + 4*t2 + 
                        s2*(-3 + t1 + t2)))) - 
               (2*(Power(s,2) + (-2 + s2)*(4 + s2 - t1) + 
                    s*(2 - 2*s2 + t1))*
                  (Power(s,5)*Power(t1 - t2,2) + 
                    2*(s2 - t1)*(-1 + t1)*
                     Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                    2*Power(s,4)*
                     (-2*Power(t1,2) + Power(t1,3) + t2 + 5*t1*t2 - 
                       2*Power(t1,2)*t2 - 3*Power(t2,2) + 
                       t1*Power(t2,2) - 
                       s1*(t1 + Power(t1,2) - s2*t2 - t1*t2) + 
                       s2*(t1 - 2*Power(t1,2) + 4*t1*t2 - 
                        2*t2*(1 + t2))) + 
                    Power(s,3)*
                     (1 + Power(t1,4) - 10*t2 + 13*Power(t2,2) - 
                       2*Power(t1,3)*(4 + 2*s1 + t2) + 
                       Power(t1,2)*
                       (7 + Power(s1,2) + 20*t2 + 4*s1*t2 + 
                       Power(t2,2)) + 
                       Power(s2,2)*
                       (2 + Power(s1,2) + 7*Power(t1,2) + 
                       2*s1*(-1 + t1 - 3*t2) + 10*t2 + 
                       6*Power(t2,2) - 4*t1*(2 + 3*t2)) - 
                       2*t1*
                       (1 + 6*t2 + 6*Power(t2,2) + s1*(-5 + 4*t2)) - 
                       2*s2*
                        (1 + 3*Power(t1,3) + 3*(-2 + s1)*t2 - 
                        7*Power(t2,2) - 
                        3*Power(t1,2)*(2 + s1 + 2*t2) + 
                        t1*(1 + 17*t2 + 3*Power(t2,2)))) + 
                    s*(4 - 10*t1 + 13*Power(t1,2) - 12*Power(t1,3) + 
                       5*Power(t1,4) + 
                       Power(s1,2)*Power(s2 - t1,2)*
                       (2 + 2*s2 + Power(s2,2) - 4*t1 + Power(t1,2)) \
- 8*t2 + 24*t1*t2 - 6*Power(t1,2)*t2 - 10*Power(t1,3)*t2 + 
                       4*Power(t2,2) - 14*t1*Power(t2,2) + 
                       9*Power(t1,2)*Power(t2,2) + 
                       Power(s2,4)*
                       (2 + 2*Power(t1,2) + 2*t2 + Power(t2,2) - 
                       2*t1*(2 + t2)) - 
                       2*Power(s2,3)*(-1 + t1)*
                       (2*Power(t1,2) + Power(1 + t2,2) - 
                       t1*(3 + 2*t2)) + 
                       2*s1*(s2 - t1)*
                       (-2*Power(t1,3) + t1*(11 - 7*t2) + 
                       Power(s2,3)*(-1 + t1 - t2) + 2*(-1 + t2) + 
                       Power(t1,2)*(-7 + 4*t2) - 
                       s2*(1 - 2*Power(t1,2) + Power(t1,3) + 
                       3*t2 - 6*t1*t2) + 
                       Power(s2,2)*(-1 + t1 - 2*t2 + t1*t2)) + 
                       2*s2*
                       (Power(t1,4) + 2*(-1 + t2)*t2 - 
                       3*Power(t1,3)*(1 + 2*t2) - 
                       t1*(1 + 12*t2 + 3*Power(t2,2)) + 
                       Power(t1,2)*(3 + 20*t2 + 3*Power(t2,2))) + 
                       Power(s2,2)*
                        (-3 + 2*Power(t1,4) + 14*t2 + 
                        5*Power(t2,2) - 2*Power(t1,3)*(4 + t2) + 
                        Power(t1,2)*(7 + 18*t2 + Power(t2,2)) - 
                        2*t1*(-1 + 15*t2 + 5*Power(t2,2)))) - 
                    2*Power(s,2)*
                     (2 + (2 + s1)*Power(t1,4) - 8*t2 + 
                       6*Power(t2,2) + 
                       t1*(-4 + s1*(7 - 5*t2) + 6*t2 - 
                       11*Power(t2,2)) - 
                       Power(t1,3)*
                        (6 + Power(s1,2) + 5*t2 + s1*(3 + t2)) + 
                       Power(s2,3)*
                        (Power(s1,2) + 3*Power(t1,2) + 
                        s1*(-2 + 2*t1 - 3*t2) + 2*Power(1 + t2,2) - 
                        t1*(5 + 4*t2)) + 
                       Power(t1,2)*
                        (Power(s1,2) + 2*s1*(-5 + 4*t2) + 
                        3*(2 + 4*t2 + Power(t2,2))) + 
                       Power(s2,2)*
                        (-(Power(s1,2)*(-1 + t1)) - 4*Power(t1,3) + 
                        t2*(9 + 5*t2) + Power(t1,2)*(9 + 6*t2) + 
                        s1*
                       (-1 + t1 + Power(t1,2) - 5*t2 + 3*t1*t2) - 
                        t1*(5 + 16*t2 + 3*Power(t2,2))) + 
                       s2*(-2 + Power(s1,2)*(-2 + t1)*t1 + 
                        Power(t1,4) + 4*t2 + 7*Power(t2,2) - 
                        2*Power(t1,3)*(3 + t2) + 
                        t1*(2 - 27*t2 - 10*Power(t2,2)) + 
                        Power(t1,2)*(5 + 21*t2 + Power(t2,2)) + 
                        s1*(-1 - 4*Power(t1,3) - t2 + 
                        Power(t1,2)*(4 + t2) + t1*(5 + 3*t2))))))/
                (s*Power(s2 - s2*t1 + t1*(-1 + s + t1),2))))/
           ((-1 + t1)*(-1 + t2)) + 
          (16*(-1 + s2)*((4*
                  (2 - t1 - Power(t1,2) - s1*Power(t1,2) + 
                    Power(s,2)*(t1 - t2) - 
                    s2*(-1 + t1)*(-1 + 2*t1 - t2) + 
                    Power(s2,2)*(-2 + s1 + 2*t1 - t2) - 2*t2 + 
                    3*t1*t2 - 
                    s*(1 - Power(t1,2) + 
                       s2*(-2 + s1 + 3*t1 - 2*t2) - 3*t2 + 
                       t1*(1 + s1 + t2)))*
                  (-10 - 4*Power(s,3) - Power(s1,2)*(-3 + s2) - 
                    21*s2 - 2*Power(s2,2) + Power(s2,3) + 
                    Power(s,2)*(3 + 5*s1 + 9*s2 + t1 - 2*t2) + 
                    2*s1*(-7 + 2*s2 - t2) + 10*t2 + s2*t2 - 
                    Power(s2,2)*t2 + 
                    s*(18 - Power(s1,2) - 6*Power(s2,2) - 2*t1 + 
                       Power(t1,2) - s2*(2 + t1 - 3*t2) + 3*t2 + 
                       t1*t2 + s1*(-6 - 5*s2 + t2))))/
                (s2*(-1 + t1) - t1*(-1 + s + t1)) - 
               8*(-5 + 2*Power(s,4) - Power(s1,3)*(-1 + s2) - 11*s2 - 
                  11*Power(s2,2) - Power(s2,3) - 4*t1 - 3*s2*t1 + 
                  3*Power(t1,2) - 
                  Power(s,3)*(2 + 3*s1 + 6*s2 - 3*t2) + 4*t2 + 
                  14*s2*t2 + Power(s2,2)*t2 - 2*t1*t2 - s2*t1*t2 + 
                  3*Power(t2,2) + s2*Power(t2,2) + 
                  Power(s,2)*
                   (-12 + Power(s1,2) + 6*Power(s2,2) - 2*t1 + 
                     s2*(3 + t1 - 5*t2) + 
                     s1*(1 + 6*s2 + t1 - 2*t2) - 3*t2 + Power(t2,2)) \
+ Power(s1,2)*(-8 - 2*Power(s2,2) + 2*t1 - t2 + s2*(4 + t2)) + 
                  s1*(-2 - Power(s2,3) + t1 + Power(t1,2) + 4*t2 - 
                     t1*t2 + Power(s2,2)*(3 + t2) - 2*s2*(9 + 2*t2)) \
+ s*(14 - 2*Power(s2,3) + Power(s1,2)*(s2 - t1) + 7*t1 - 
                     Power(t1,2) - Power(s2,2)*(t1 - 2*t2) - 11*t2 - 
                     2*t1*t2 - 3*Power(t2,2) + 
                     s1*(16 - 2*Power(s2,2) + s2*(-5 + t2) + 
                        (3 + t1)*t2) + 
                     s2*(22 + Power(t1,2) + t1*(-1 + t2) + t2 - 
                        Power(t2,2)))) - 
               (4*(-10 + 3*Power(s,2) + s1 - 2*s2 + s1*s2 + 
                    2*Power(s2,2) + t2 - s2*t2 + 
                    s*(-2 - 2*s1 - 5*s2 + t2))*
                  (Power(s,5)*Power(t1 - t2,2) + 
                    2*(s2 - t1)*(-1 + t1)*
                     Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                    2*Power(s,4)*
                     (-2*Power(t1,2) + Power(t1,3) + t2 + 5*t1*t2 - 
                       2*Power(t1,2)*t2 - 3*Power(t2,2) + 
                       t1*Power(t2,2) - 
                       s1*(t1 + Power(t1,2) - s2*t2 - t1*t2) + 
                       s2*(t1 - 2*Power(t1,2) + 4*t1*t2 - 
                       2*t2*(1 + t2))) + 
                    Power(s,3)*
                     (1 + Power(t1,4) - 10*t2 + 13*Power(t2,2) - 
                       2*Power(t1,3)*(4 + 2*s1 + t2) + 
                       Power(t1,2)*
                       (7 + Power(s1,2) + 20*t2 + 4*s1*t2 + 
                       Power(t2,2)) + 
                       Power(s2,2)*
                       (2 + Power(s1,2) + 7*Power(t1,2) + 
                       2*s1*(-1 + t1 - 3*t2) + 10*t2 + 
                       6*Power(t2,2) - 4*t1*(2 + 3*t2)) - 
                       2*t1*
                       (1 + 6*t2 + 6*Power(t2,2) + s1*(-5 + 4*t2)) \
- 2*s2*(1 + 3*Power(t1,3) + 3*(-2 + s1)*t2 - 7*Power(t2,2) - 
                       3*Power(t1,2)*(2 + s1 + 2*t2) + 
                       t1*(1 + 17*t2 + 3*Power(t2,2)))) + 
                    s*(4 - 10*t1 + 13*Power(t1,2) - 
                       12*Power(t1,3) + 5*Power(t1,4) + 
                       Power(s1,2)*Power(s2 - t1,2)*
                       (2 + 2*s2 + Power(s2,2) - 4*t1 + 
                       Power(t1,2)) - 8*t2 + 24*t1*t2 - 
                       6*Power(t1,2)*t2 - 10*Power(t1,3)*t2 + 
                       4*Power(t2,2) - 14*t1*Power(t2,2) + 
                       9*Power(t1,2)*Power(t2,2) + 
                       Power(s2,4)*
                       (2 + 2*Power(t1,2) + 2*t2 + Power(t2,2) - 
                       2*t1*(2 + t2)) - 
                       2*Power(s2,3)*(-1 + t1)*
                       (2*Power(t1,2) + Power(1 + t2,2) - 
                       t1*(3 + 2*t2)) + 
                       2*s1*(s2 - t1)*
                       (-2*Power(t1,3) + t1*(11 - 7*t2) + 
                       Power(s2,3)*(-1 + t1 - t2) + 2*(-1 + t2) + 
                       Power(t1,2)*(-7 + 4*t2) - 
                       s2*(1 - 2*Power(t1,2) + Power(t1,3) + 
                       3*t2 - 6*t1*t2) + 
                       Power(s2,2)*(-1 + t1 - 2*t2 + t1*t2)) + 
                       2*s2*
                       (Power(t1,4) + 2*(-1 + t2)*t2 - 
                       3*Power(t1,3)*(1 + 2*t2) - 
                       t1*(1 + 12*t2 + 3*Power(t2,2)) + 
                       Power(t1,2)*(3 + 20*t2 + 3*Power(t2,2))) + 
                       Power(s2,2)*
                       (-3 + 2*Power(t1,4) + 14*t2 + 
                       5*Power(t2,2) - 2*Power(t1,3)*(4 + t2) + 
                       Power(t1,2)*(7 + 18*t2 + Power(t2,2)) - 
                       2*t1*(-1 + 15*t2 + 5*Power(t2,2)))) - 
                    2*Power(s,2)*
                     (2 + (2 + s1)*Power(t1,4) - 8*t2 + 
                       6*Power(t2,2) + 
                       t1*(-4 + s1*(7 - 5*t2) + 6*t2 - 
                       11*Power(t2,2)) - 
                       Power(t1,3)*
                       (6 + Power(s1,2) + 5*t2 + s1*(3 + t2)) + 
                       Power(s2,3)*
                       (Power(s1,2) + 3*Power(t1,2) + 
                       s1*(-2 + 2*t1 - 3*t2) + 2*Power(1 + t2,2) - 
                       t1*(5 + 4*t2)) + 
                       Power(t1,2)*
                       (Power(s1,2) + 2*s1*(-5 + 4*t2) + 
                       3*(2 + 4*t2 + Power(t2,2))) + 
                       Power(s2,2)*
                       (-(Power(s1,2)*(-1 + t1)) - 4*Power(t1,3) + 
                       t2*(9 + 5*t2) + Power(t1,2)*(9 + 6*t2) + 
                       s1*(-1 + t1 + Power(t1,2) - 5*t2 + 
                      3*t1*t2) - t1*(5 + 16*t2 + 3*Power(t2,2))) + 
                       s2*(-2 + Power(s1,2)*(-2 + t1)*t1 + 
                        Power(t1,4) + 4*t2 + 7*Power(t2,2) - 
                        2*Power(t1,3)*(3 + t2) + 
                        t1*(2 - 27*t2 - 10*Power(t2,2)) + 
                        Power(t1,2)*(5 + 21*t2 + Power(t2,2)) + 
                        s1*(-1 - 4*Power(t1,3) - t2 + 
                        Power(t1,2)*(4 + t2) + t1*(5 + 3*t2))))))/
                (s*Power(s2 - s2*t1 + t1*(-1 + s + t1),2)) + 
               ((1 + s - s2)*
                  (2 - t1 - Power(t1,2) - s1*Power(t1,2) + 
                    Power(s,2)*(t1 - t2) - 
                    s2*(-1 + t1)*(-1 + 2*t1 - t2) + 
                    Power(s2,2)*(-2 + s1 + 2*t1 - t2) - 2*t2 + 
                    3*t1*t2 - 
                    s*(1 - Power(t1,2) + 
                       s2*(-2 + s1 + 3*t1 - 2*t2) - 3*t2 + 
                       t1*(1 + s1 + t2)))*
                  (4*s*Power(-2 + t1 + Power(t1,2) + 
                       s1*Power(t1,2) + 
                       s2*(-1 + t1)*(-1 + 2*t1 - t2) + 2*t2 - 
                       3*t1*t2 + Power(s2,2)*(2 - s1 - 2*t1 + t2) + 
                       Power(s,2)*(-t1 + t2) + 
                       s*(1 - Power(t1,2) + 
                       s2*(-2 + s1 + 3*t1 - 2*t2) - 3*t2 + 
                       t1*(1 + s1 + t2)),2) - 
                    12*(s2*(-1 + t1) - t1*(-1 + s + t1))*
                     (Power(s,3)*(-1 + s2)*(-1 + s1 + t1 - t2) - 
                       Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                       Power(s,2)*
                        (3 + Power(s1,2)*s2 - 3*t1 + Power(t1,2) + 
                        2*Power(s2,2)*(-1 + t1 - t2) + 2*t2 - 
                        t1*t2 + s2*t1*(1 - t1 + t2) + 
                        s1*(-3 + 2*Power(s2,2) + t2 - s2*(1 + t2))) + 
                       s*(Power(s1,2)*(1 + s2)*(s2 - t1) + 
                        Power(s2,3)*(-1 + t1 - t2) - 
                        Power(s2,2)*(-1 + t1)*(-1 + t1 - t2) - 
                        s2*(-2 + t1 + Power(t1,2) - 4*t2 - t1*t2) + 
                        2*(1 + Power(t1,2) - t2 - t1*(2 + t2)) + 
                        s1*(-2 + Power(s2,3) + 5*t1 + Power(t1,2) + 
                        2*t2 - Power(s2,2)*t2 - t1*t2 - 
                        s2*(3 + t1 + Power(t1,2) + t2 - t1*t2))))))/
                (s*Power(s2 - s2*t1 + t1*(-1 + s + t1),3))))/
           ((-1 + s1)*(-1 + t2)) + 
          (16*(-1 + s2)*((-4*
                  (2 - t1 - Power(t1,2) - s1*Power(t1,2) + 
                    Power(s,2)*(t1 - t2) - 
                    s2*(-1 + t1)*(-1 + 2*t1 - t2) + 
                    Power(s2,2)*(-2 + s1 + 2*t1 - t2) - 2*t2 + 
                    3*t1*t2 - 
                    s*(1 - Power(t1,2) + 
                       s2*(-2 + s1 + 3*t1 - 2*t2) - 3*t2 + 
                       t1*(1 + s1 + t2)))*
                  (9 + 4*Power(s,3) + Power(s1,2)*(-3 + s2) + 22*s2 + 
                    2*Power(s2,2) - Power(s2,3) + 
                    Power(s,2)*(1 - 5*s1 - 11*s2 + 5*t1) - 8*t2 - 
                    3*s2*t2 + Power(s2,2)*t2 + 
                    s1*(13 - 3*s2 + 2*t2) + 
                    s*(-17 + Power(s1,2) + 8*Power(s2,2) + 2*t1 + 
                       Power(t1,2) + s1*(3 + 5*s2 - 4*t1 - t2) + 
                       3*t2 + t1*t2 - s2*(-2 + 7*t1 + t2))))/
                (s2*(-1 + t1) - t1*(-1 + s + t1)) + 
               8*(2 + Power(s1,3)*(-1 + s2) + 13*s2 + 12*Power(s2,2) + 
                  Power(s2,3) - 5*t1 - 2*s2*t1 + 3*Power(t1,2) + 
                  7*t2 - 8*s2*t2 - 2*t1*t2 - s2*t1*t2 + 
                  3*Power(t2,2) + s2*Power(t2,2) + 
                  Power(s,3)*(-8 + s1 + 4*s2 + t2) + 
                  Power(s,2)*
                   (20 - Power(s1,2) - 8*Power(s2,2) - 6*t1 + 
                     s1*(9 - 4*s2 + t1) + s2*(7 + 5*t1 - t2) - 7*t2 + 
                     Power(t2,2)) + 
                  Power(s1,2)*
                   (11 + 2*Power(s2,2) - 2*t1 + t2 - s2*(3 + t2)) + 
                  s1*(-5 + Power(s2,3) + Power(t1,2) - t1*(-6 + t2) - 
                     12*t2 - Power(s2,2)*(3 + t2) + 
                     s2*(13 - t1 + 2*t2)) + 
                  s*(-17 + 4*Power(s2,3) + 4*t1 - 5*Power(s2,2)*t1 - 
                     Power(t1,2) - Power(s1,2)*(4 + s2 + t1) + 6*t2 - 
                     2*t1*t2 - 3*Power(t2,2) + 
                     s2*(-27 + Power(t1,2) + 4*t2 - Power(t2,2) + 
                        t1*(4 + t2)) + 
                     s1*(-13 + 2*Power(s2,2) + 5*t2 + t1*(6 + t2) + 
                        s2*(-4*t1 + t2)))) - 
               (4*(-10 + 3*Power(s,2) + s1 - 2*s2 + s1*s2 + 
                    2*Power(s2,2) + t2 - s2*t2 + 
                    s*(-2 - 2*s1 - 5*s2 + 2*t1 + t2))*
                  (Power(s,5)*Power(t1 - t2,2) + 
                    2*(s2 - t1)*(-1 + t1)*
                     Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                    2*Power(s,4)*
                     (-2*Power(t1,2) + Power(t1,3) + t2 + 5*t1*t2 - 
                       2*Power(t1,2)*t2 - 3*Power(t2,2) + 
                       t1*Power(t2,2) - 
                       s1*(t1 + Power(t1,2) - s2*t2 - t1*t2) + 
                       s2*(t1 - 2*Power(t1,2) + 4*t1*t2 - 
                        2*t2*(1 + t2))) + 
                    Power(s,3)*
                     (1 + Power(t1,4) - 10*t2 + 13*Power(t2,2) - 
                       2*Power(t1,3)*(4 + 2*s1 + t2) + 
                       Power(t1,2)*
                       (7 + Power(s1,2) + 20*t2 + 4*s1*t2 + 
                       Power(t2,2)) + 
                       Power(s2,2)*
                       (2 + Power(s1,2) + 7*Power(t1,2) + 
                       2*s1*(-1 + t1 - 3*t2) + 10*t2 + 
                       6*Power(t2,2) - 4*t1*(2 + 3*t2)) - 
                       2*t1*
                       (1 + 6*t2 + 6*Power(t2,2) + s1*(-5 + 4*t2)) - 
                       2*s2*
                        (1 + 3*Power(t1,3) + 3*(-2 + s1)*t2 - 
                        7*Power(t2,2) - 
                        3*Power(t1,2)*(2 + s1 + 2*t2) + 
                        t1*(1 + 17*t2 + 3*Power(t2,2)))) + 
                    s*(4 - 10*t1 + 13*Power(t1,2) - 12*Power(t1,3) + 
                       5*Power(t1,4) + 
                       Power(s1,2)*Power(s2 - t1,2)*
                       (2 + 2*s2 + Power(s2,2) - 4*t1 + Power(t1,2)) \
- 8*t2 + 24*t1*t2 - 6*Power(t1,2)*t2 - 10*Power(t1,3)*t2 + 
                       4*Power(t2,2) - 14*t1*Power(t2,2) + 
                       9*Power(t1,2)*Power(t2,2) + 
                       Power(s2,4)*
                       (2 + 2*Power(t1,2) + 2*t2 + Power(t2,2) - 
                       2*t1*(2 + t2)) - 
                       2*Power(s2,3)*(-1 + t1)*
                       (2*Power(t1,2) + Power(1 + t2,2) - 
                       t1*(3 + 2*t2)) + 
                       2*s1*(s2 - t1)*
                       (-2*Power(t1,3) + t1*(11 - 7*t2) + 
                       Power(s2,3)*(-1 + t1 - t2) + 2*(-1 + t2) + 
                       Power(t1,2)*(-7 + 4*t2) - 
                       s2*(1 - 2*Power(t1,2) + Power(t1,3) + 
                       3*t2 - 6*t1*t2) + 
                       Power(s2,2)*(-1 + t1 - 2*t2 + t1*t2)) + 
                       2*s2*
                       (Power(t1,4) + 2*(-1 + t2)*t2 - 
                       3*Power(t1,3)*(1 + 2*t2) - 
                       t1*(1 + 12*t2 + 3*Power(t2,2)) + 
                       Power(t1,2)*(3 + 20*t2 + 3*Power(t2,2))) + 
                       Power(s2,2)*
                        (-3 + 2*Power(t1,4) + 14*t2 + 
                        5*Power(t2,2) - 2*Power(t1,3)*(4 + t2) + 
                        Power(t1,2)*(7 + 18*t2 + Power(t2,2)) - 
                        2*t1*(-1 + 15*t2 + 5*Power(t2,2)))) - 
                    2*Power(s,2)*
                     (2 + (2 + s1)*Power(t1,4) - 8*t2 + 
                       6*Power(t2,2) + 
                       t1*(-4 + s1*(7 - 5*t2) + 6*t2 - 
                       11*Power(t2,2)) - 
                       Power(t1,3)*
                        (6 + Power(s1,2) + 5*t2 + s1*(3 + t2)) + 
                       Power(s2,3)*
                        (Power(s1,2) + 3*Power(t1,2) + 
                        s1*(-2 + 2*t1 - 3*t2) + 2*Power(1 + t2,2) - 
                        t1*(5 + 4*t2)) + 
                       Power(t1,2)*
                        (Power(s1,2) + 2*s1*(-5 + 4*t2) + 
                        3*(2 + 4*t2 + Power(t2,2))) + 
                       Power(s2,2)*
                        (-(Power(s1,2)*(-1 + t1)) - 4*Power(t1,3) + 
                        t2*(9 + 5*t2) + Power(t1,2)*(9 + 6*t2) + 
                        s1*
                       (-1 + t1 + Power(t1,2) - 5*t2 + 3*t1*t2) - 
                        t1*(5 + 16*t2 + 3*Power(t2,2))) + 
                       s2*(-2 + Power(s1,2)*(-2 + t1)*t1 + 
                        Power(t1,4) + 4*t2 + 7*Power(t2,2) - 
                        2*Power(t1,3)*(3 + t2) + 
                        t1*(2 - 27*t2 - 10*Power(t2,2)) + 
                        Power(t1,2)*(5 + 21*t2 + Power(t2,2)) + 
                        s1*(-1 - 4*Power(t1,3) - t2 + 
                        Power(t1,2)*(4 + t2) + t1*(5 + 3*t2))))))/
                (s*Power(s2 - s2*t1 + t1*(-1 + s + t1),2)) + 
               ((1 + s - s2)*
                  (2 - t1 - Power(t1,2) - s1*Power(t1,2) + 
                    Power(s,2)*(t1 - t2) - 
                    s2*(-1 + t1)*(-1 + 2*t1 - t2) + 
                    Power(s2,2)*(-2 + s1 + 2*t1 - t2) - 2*t2 + 
                    3*t1*t2 - 
                    s*(1 - Power(t1,2) + s2*(-2 + s1 + 3*t1 - 2*t2) - 
                       3*t2 + t1*(1 + s1 + t2)))*
                  (4*s*Power(-2 + t1 + Power(t1,2) + 
                       s1*Power(t1,2) + 
                       s2*(-1 + t1)*(-1 + 2*t1 - t2) + 2*t2 - 
                       3*t1*t2 + Power(s2,2)*(2 - s1 - 2*t1 + t2) + 
                       Power(s,2)*(-t1 + t2) + 
                       s*(1 - Power(t1,2) + 
                        s2*(-2 + s1 + 3*t1 - 2*t2) - 3*t2 + 
                        t1*(1 + s1 + t2)),2) - 
                    12*(s2*(-1 + t1) - t1*(-1 + s + t1))*
                     (Power(s,3)*(-1 + s2)*(-1 + s1 + t1 - t2) - 
                       Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                       Power(s,2)*
                        (3 + Power(s1,2)*s2 - 3*t1 + Power(t1,2) + 
                        2*Power(s2,2)*(-1 + t1 - t2) + 2*t2 - t1*t2 + 
                        s2*t1*(1 - t1 + t2) + 
                        s1*(-3 + 2*Power(s2,2) + t2 - s2*(1 + t2))) + 
                       s*(Power(s1,2)*(1 + s2)*(s2 - t1) + 
                        Power(s2,3)*(-1 + t1 - t2) - 
                        Power(s2,2)*(-1 + t1)*(-1 + t1 - t2) - 
                        s2*(-2 + t1 + Power(t1,2) - 4*t2 - t1*t2) + 
                        2*(1 + Power(t1,2) - t2 - t1*(2 + t2)) + 
                        s1*(-2 + Power(s2,3) + 5*t1 + Power(t1,2) + 
                        2*t2 - Power(s2,2)*t2 - t1*t2 - 
                        s2*(3 + t1 + Power(t1,2) + t2 - t1*t2))))))/
                (s*Power(s2 - s2*t1 + t1*(-1 + s + t1),3))))/
           ((-1 + s1)*(-s + s1 - t2)))/Power(-1 + s2,2))*
     B1(1 - s + s2 - t1,s,s2))/(128.*Power(Pi,2)) + 
  (((-64*(-1 + s2 - t1 + t2)*
          (2*(s - s2 + t1)*Power(-1 + s1 + t1 - t2,2) + 
            (-4 + s2 + t1 + 
               (-1 + s - s2 + t1)*(-2 + 2*s - s2 + 3*t1 - 2*t2))*
             (1 - s1 - t1 + t2) + 
            (-2 + s - s2 + t1)*
             (-1 + s2 + s*t1 - s2*t1 + Power(t1,2) + t2 - s*t2 - t1*t2)\
)*((-Power(s,4) + 22*s2 - 5*Power(s2,2) - 4*Power(s2,3) + 
               Power(s2,4) - 2*Power(s1,3)*(s2 - t1) - t1 + 6*s2*t1 + 
               9*Power(s2,2)*t1 - 2*Power(s2,3)*t1 - 4*Power(t1,2) - 
               6*s2*Power(t1,2) + Power(s2,2)*Power(t1,2) + 
               Power(t1,3) + Power(s,3)*(-1 + 4*s1 + 2*s2 - 2*t2) - 
               17*s2*t2 - 6*Power(s2,2)*t2 + Power(s2,3)*t2 + 
               2*t1*t2 + 7*s2*t1*t2 - s2*Power(t1,2)*t2 + 
               2*Power(t2,2) + s2*Power(t2,2) - 
               Power(s2,2)*Power(t2,2) + t1*Power(t2,2) + 
               s2*t1*Power(t2,2) - 2*Power(t2,3) + 
               Power(s1,2)*(13 - 2*Power(s2,2) + Power(t1,2) - 
                  3*t2 - t1*(2 + 3*t2) + s2*(2 + t1 + 3*t2)) + 
               Power(s,2)*(16 - 5*Power(s1,2) + 4*t1 + Power(t1,2) - 
                  8*t2 - t1*t2 - Power(t2,2) + 
                  s2*(-2 - 2*t1 + 5*t2) + s1*(3 - 7*s2 + 2*t1 + 5*t2)\
) + s1*(21 + Power(s2,3) - 18*t2 + 5*Power(t2,2) - 
                  Power(t1,2)*(4 + t2) + 
                  Power(s2,2)*(-1 - 3*t1 + 4*t2) + 
                  s2*(14 + 2*Power(t1,2) + t1*(5 - 3*t2) - 6*t2 - 
                     Power(t2,2)) + t1*(-5 + 5*t2 + Power(t2,2))) + 
               s*(-22 + 2*Power(s1,3) - 2*Power(s2,3) - 2*t1 + 
                  6*Power(t1,2) + Power(s2,2)*(7 + 4*t1 - 4*t2) + 
                  Power(s1,2)*(-2 + 7*s2 - 4*t1 - 3*t2) + 20*t2 - 
                  5*t1*t2 + Power(t1,2)*t2 - 4*Power(t2,2) - 
                  t1*Power(t2,2) + 
                  s1*(-32 + 2*Power(s2,2) - t1 - 2*Power(t1,2) + 
                     s2*(-1 + t1 - 9*t2) + 11*t2 + 4*t1*t2 + 
                     Power(t2,2)) + 
                  s2*(-10 - 2*Power(t1,2) + t1*(-12 + t2) + 13*t2 + 
                     2*Power(t2,2))))/((-1 + s2)*(-s + s2 - t1)) - 
            (Power(s,4) + 2*Power(s1,4) - 23*s2 - 11*Power(s2,2) + 
               4*Power(s2,3) - Power(s2,4) + 11*t1 + 11*s2*t1 + 
               Power(s2,3)*t1 - Power(t1,2) - 3*s2*Power(t1,2) + 
               Power(s1,3)*(4*s2 + t1 - 5*t2) - 
               Power(s,3)*(-2 + 5*s1 + 2*s2 + t1 - 3*t2) + 6*t2 + 
               8*s2*t2 + 2*Power(s2,2)*t2 + t1*t2 + 2*s2*t1*t2 - 
               2*Power(s2,2)*t1*t2 - 6*Power(t2,2) - 
               3*s2*Power(t2,2) + 2*Power(s2,2)*Power(t2,2) - 
               2*t1*Power(t2,2) + s2*t1*Power(t2,2) + 2*Power(t2,3) - 
               s2*Power(t2,3) + 
               Power(s,2)*(-27 + 9*Power(s1,2) + 
                  s1*(-3 + 8*s2 + 3*t1 - 11*t2) + 
                  t1*(3 + 3*s2 - 2*t2) + 4*t2 - 6*s2*t2 + 
                  3*Power(t2,2)) + 
               Power(s1,2)*(-19 + Power(s2,2) + 
                  3*s2*(1 + t1 - 3*t2) + 5*t2 + 4*Power(t2,2) - 
                  2*t1*(1 + t2)) + 
               s1*(-18 - 2*Power(s2,3) - 4*Power(t1,2) + 
                  Power(s2,2)*(6 + 3*t1 - 4*t2) + 17*t2 - 
                  7*Power(t2,2) - Power(t2,3) + 
                  t1*(12 + 5*t2 + Power(t2,2)) + 
                  s2*(-30 - 4*(-1 + t1)*t2 + 6*Power(t2,2))) - 
               s*(-29 + 7*Power(s1,3) - 2*Power(s2,3) + 15*t1 - 
                  4*Power(t1,2) + 
                  Power(s1,2)*(-1 + 10*s2 + 3*t1 - 13*t2) + 
                  3*Power(s2,2)*(2 + t1 - t2) + 16*t2 + 4*t1*t2 - 
                  5*Power(t2,2) + t1*Power(t2,2) - Power(t2,3) + 
                  s2*(-36 + 3*t1 + 6*t2 - 4*t1*t2 + 5*Power(t2,2)) + 
                  s1*(-44 + Power(s2,2) + t1 + 
                     3*s2*(1 + 2*t1 - 5*t2) + 8*t2 - 4*t1*t2 + 
                     7*Power(t2,2))))/((-1 + s2)*(-1 + t1)) + 
            (3*Power(s1,3) + 36*s2 + 15*Power(s2,2) - 4*Power(s2,3) + 
               Power(s,3)*(-3 + t1) - 24*t1 - 19*s2*t1 + 
               5*Power(s2,2)*t1 + 4*Power(t1,2) - s2*Power(t1,2) + 
               Power(s1,2)*(-2 + s2 + 4*t1 - s2*t1 + Power(t1,2) - 
                  8*t2) + 12*t2 - 3*s2*t2 - 7*Power(s2,2)*t2 + 
               9*t1*t2 + 8*s2*t1*t2 + Power(s2,2)*t1*t2 - 
               3*Power(t1,2)*t2 - s2*Power(t1,2)*t2 - 
               12*Power(t2,2) - 2*s2*Power(t2,2) + 5*t1*Power(t2,2) - 
               2*Power(t2,3) + 
               Power(s,2)*(18 - 2*s1*(-3 + t1) - 2*s2*(-1 + t1) - 
                  t1 + Power(t1,2) - 7*t2 + t1*t2) + 
               s1*(-(Power(s2,2)*(1 + t1)) - Power(t1,2)*(-2 + t2) - 
                  5*t1*(4 + t2) + t2*(18 + 7*t2) + 
                  s2*(18 + t1 + Power(t1,2) - 3*t2 + t1*t2)) + 
               s*(-24 + Power(s1,2)*(-6 + t1) + 28*t1 - 
                  2*Power(t1,2) + Power(s2,2)*(5 + t1) - 9*t2 - 
                  t1*t2 + Power(t1,2)*t2 - 2*Power(t2,2) - 
                  s2*(33 + Power(t1,2) - 12*t2 + 2*t1*(1 + t2)) + 
                  s1*(3*s2*(-1 + t1) - 2*Power(t1,2) + 12*(-1 + t2) - 
                     t1*(3 + t2))))/((-1 + s1)*(-s + s1 - t2)) + 
            (Power(s1,3) + 20*s2 - 4*Power(s2,3) + 
               Power(s,3)*(-2 + t1) + 4*s2*t1 + 8*Power(s2,2)*t1 - 
               4*Power(t1,2) - 4*s2*Power(t1,2) + 
               Power(s1,2)*(2 - s2*(-7 + t1) - 5*t1 + Power(t1,2) - 
                  t2) - 10*s2*t2 - 6*Power(s2,2)*t2 + 3*t1*t2 + 
               2*s2*t1*t2 + Power(s2,2)*t1*t2 + 3*Power(t1,2)*t2 - 
               s2*Power(t1,2)*t2 - Power(t2,2) + 5*s2*Power(t2,2) - 
               4*t1*Power(t2,2) + Power(t2,3) + 
               Power(s,2)*(9 - 2*s1*(-4 + t1) + Power(t1,2) - 
                  10*t2 + t1*(5 - 2*s2 + t2)) + 
               s*(-20 + Power(s1,2)*(-7 + t1) - t1 + 7*Power(t1,2) + 
                  Power(s2,2)*(6 + t1) + 8*t2 - 7*t1*t2 + 
                  Power(t1,2)*t2 - 3*Power(t2,2) - 
                  s1*(15 + s2*(7 - 3*t1) + 2*Power(t1,2) - 14*t2 + 
                     t1*t2) - 
                  s2*(10 + Power(t1,2) - 17*t2 + 2*t1*(7 + t2))) - 
               s1*(-20 + Power(s2,2)*(2 + t1) + t1*(7 - 13*t2) + 
                  5*t2 + Power(t2,2) + Power(t1,2)*(7 + t2) - 
                  s2*(10 + Power(t1,2) - 16*t2 + t1*(10 + t2))))/
             ((s - s2 + t1)*(s - s1 + t2)) + 
            (-2*Power(s,4) + 31*s2 + 20*Power(s2,2) + Power(s2,3) - 
               Power(s2,4) - 18*t1 - 30*s2*t1 - 4*Power(s2,2)*t1 + 
               2*Power(s2,3)*t1 + 9*Power(t1,2) + 4*s2*Power(t1,2) - 
               Power(s2,2)*Power(t1,2) - Power(t1,3) + 
               Power(s1,3)*(4 - s2 + t1) + 
               Power(s,3)*(1 + 5*s1 + 7*s2 - 4*t1 - t2) + 
               Power(s1,2)*(1 - 4*Power(s2,2) + t1 - 2*Power(t1,2) + 
                  s2*(7 + 6*t1) - t2) + 10*t2 - 6*s2*t2 - 
               Power(s2,2)*t2 + 2*Power(s2,3)*t2 + 3*t1*t2 + 
               2*s2*t1*t2 - 4*Power(s2,2)*t1*t2 + 
               2*s2*Power(t1,2)*t2 - 5*Power(t2,2) + 
               Power(s2,2)*Power(t2,2) - 2*t1*Power(t2,2) - 
               s2*t1*Power(t2,2) + 3*Power(t2,3) + 
               Power(s,2)*(16 - 4*Power(s1,2) - 9*Power(s2,2) - 
                  3*t1 - 2*Power(t1,2) - t2 - 3*t1*t2 + 
                  Power(t2,2) + s1*(1 - 14*s2 + 9*t1 + t2) + 
                  2*s2*(-1 + 5*t1 + 2*t2)) + 
               s1*(3 - 4*Power(s2,3) + 3*t2 - 6*Power(t2,2) + 
                  2*Power(t1,2)*(1 + t2) + 
                  Power(s2,2)*(3 + 7*t1 + 2*t2) - 
                  s2*(-18 + t1 + 3*Power(t1,2) + 6*t2 + 4*t1*t2 - 
                     Power(t2,2)) - t1*(17 - 5*t2 + Power(t2,2))) + 
               s*(-21 + Power(s1,3) + 5*Power(s2,3) + 29*t1 - 
                  4*Power(t1,2) + Power(s1,2)*(8*s2 - 6*(1 + t1)) + 
                  2*t2 - 5*t1*t2 - 2*Power(t1,2)*t2 + 4*Power(t2,2) + 
                  t1*Power(t2,2) - Power(s2,2)*(8*t1 + 5*t2) + 
                  s2*(-35 + 6*t1 + 3*Power(t1,2) + 3*t2 + 7*t1*t2 - 
                     2*Power(t2,2)) + 
                  s1*(-14 + 13*Power(s2,2) + t1 + 4*Power(t1,2) + 
                     2*t2 + 3*t1*t2 - Power(t2,2) - 
                     s2*(4 + 16*t1 + 3*t2))))/((-1 + s1)*(-1 + t2)) + 
            (-2*Power(s,4) - Power(s1,4) + 18*s2 + 15*Power(s2,2) + 
               Power(s2,3) - Power(s2,4) - 6*t1 - 12*s2*t1 - 
               3*Power(s2,2)*t1 + Power(s2,3)*t1 + s2*Power(t1,2) + 
               Power(s,3)*(-1 + 7*s1 + 7*s2 - 2*t1 - 3*t2) - 5*t2 - 
               21*s2*t2 - 3*Power(s2,2)*t2 + 3*Power(s2,3)*t2 + 
               9*t1*t2 + 5*s2*t1*t2 - 3*Power(s2,2)*t1*t2 - 
               Power(s2,2)*Power(t2,2) - t1*Power(t2,2) + 
               2*s2*t1*Power(t2,2) + Power(t2,3) - s2*Power(t2,3) + 
               Power(s1,3)*(2 - 5*s2 + 2*t1 + t2) + 
               Power(s1,2)*(13 - 8*Power(s2,2) - 2*t2 + Power(t2,2) - 
                  t1*(3 + 4*t2) + s2*(6 + 5*t1 + 6*t2)) + 
               Power(s,2)*(21 - 9*Power(s1,2) - 9*Power(s2,2) - 
                  2*t1 - 4*t2 - 4*t1*t2 + 
                  s1*(5 - 19*s2 + 6*t1 + 7*t2) + s2*(4 + 5*t1 + 9*t2)) \
+ s*(-23 + 5*Power(s1,3) + 5*Power(s2,3) + 16*t1 + 
                  Power(s1,2)*(-6 + 17*s2 - 6*t1 - 5*t2) + 20*t2 - 
                  5*t1*t2 - 2*t1*Power(t2,2) + Power(t2,3) - 
                  Power(s2,2)*(4 + 4*t1 + 9*t2) + 
                  s2*(-37 + 6*t1 + 6*t2 + 7*t1*t2 + Power(t2,2)) + 
                  s1*(-36 + 17*Power(s2,2) + 5*t1 + 7*t2 + 8*t1*t2 - 
                     Power(t2,2) - 5*s2*(2 + 2*t1 + 3*t2))) - 
               s1*(-17 + 5*Power(s2,3) + 13*t2 + Power(t2,2) + 
                  Power(t2,3) - Power(s2,2)*(5 + 4*t1 + 8*t2) - 
                  2*t1*(-9 + 3*t2 + Power(t2,2)) + 
                  s2*(-33 + 5*t2 + t1*(5 + 7*t2))))/
             ((-1 + t1)*(-1 + t2))))/
        ((s - s2 + t1)*Power(-1 + s1 + t1 - t2,3)*(-1 + s2 - t1 + t2) - 
          Power(-1 + s2 + s*t1 - s2*t1 + Power(t1,2) + t2 - s*t2 - 
            t1*t2,2) - Power(-1 + s1 + t1 - t2,2)*
           (3 + Power(s2,2) + 2*t1 - s2*(3 + t1) + (-3 + t1)*t2 + 
             2*Power(-1 + s - s2 + t1,2)*(-1 + s2 - t1 + t2) + 
             (-1 + s - s2 + t1)*(s2 - t2)*(-1 + s2 - t1 + t2)) + 
          (1 - s1 - t1 + t2)*(2*(-1 + t1)*(-1 + t2) + 
             Power(s2,2)*(2 + t2) - 
             Power(-1 + s - s2 + t1,3)*(-1 + s2 - t1 + t2) + 
             Power(-1 + s - s2 + t1,2)*(1 - s2 + t2)*
              (-1 + s2 - t1 + t2) + 
             (1 - s + s2 - t1)*
              (2 - s2 + 4*t1 - Power(s2,2)*(1 + t2) + 
                t2*(-3 - t1 + t2) - s2*(1 + t2)*(-t1 + t2)) - 
             s2*(2*(2 + t1) + t2*(-5 + t1 + t2)))) + 
       8*(((-128*(-7 + Power(s,2) + s1 - s2 + 2*s1*s2 + 
                  Power(s2,2) + t1 - 2*s1*t1 - 2*s2*t1 + 
                  Power(t1,2) + s*(3 - 2*s1 - 2*s2 + 2*t1) + 2*t2)*
                (-s2 + 2*Power(s1,2)*(s2 - t1) + t1 + s2*t1 - 
                  Power(t1,2) + Power(s,2)*(-2 + 2*s1 + t1 - t2) + 
                  2*t2 - s2*t2 - Power(s2,2)*t2 + t1*t2 + s2*t1*t2 - 
                  2*Power(t2,2) + 
                  s1*(-2 + Power(s2,2) - Power(t1,2) + 2*t2 - 
                     2*s2*t2 + 2*t1*t2) + 
                  s*(3 - 2*Power(s1,2) - 3*t1 + Power(t1,2) - t2 - 
                     t1*t2 + s2*(2 - t1 + 2*t2) + 
                     s1*(-3*s2 + t1 + 2*t2))))/
              (1 + s*(-1 + s1) + s1*(-s2 + t1) - t2) + 
             256*(15 + 2*Power(s,3) + 8*s2 - Power(s2,2) - t1 + 
                3*s2*t1 - Power(s2,2)*t1 - 2*Power(t1,2) + 
                2*s2*Power(t1,2) - Power(t1,3) + 
                Power(s1,2)*(-1 - 4*s2 + 4*t1) - 11*t2 - 3*s2*t2 + 
                2*Power(s2,2)*t2 + t1*t2 - 3*s2*t1*t2 + 
                Power(t1,2)*t2 + 4*Power(t2,2) + 
                2*Power(s,2)*(2 - 3*s1 - 2*s2 + t1 + t2) + 
                s1*(15 - 3*Power(s2,2) - 4*t2 - 3*t1*(1 + t2) + 
                   s2*(4 + 3*t1 + 3*t2)) + 
                s*(-18 + 4*Power(s1,2) + 2*Power(s2,2) - t1 - 
                   Power(t1,2) + s1*(-4 + 9*s2 - 6*t1 - 3*t2) + 8*t2 + 
                   3*t1*t2 - s2*(2 + t1 + 4*t2))))/
           (16.*(-1 + s2)*(-s + s2 - t1)) + 
          (8*(30 + 4*s2 - 7*Power(s2,2) - Power(s2,3) - 
               2*Power(s1,3)*(s2 - t1) - 14*t1 + 12*s2*t1 + 
               6*Power(s2,2)*t1 + Power(s2,3)*t1 - 5*Power(t1,2) - 
               9*s2*Power(t1,2) - 3*Power(s2,2)*Power(t1,2) + 
               4*Power(t1,3) + 3*s2*Power(t1,3) - Power(t1,4) + 
               Power(s1,2)*(2 + Power(s2,3) - 6*t1 - 6*Power(t1,2) + 
                  Power(t1,3) - Power(s2,2)*(4 + t1) + 
                  s2*(6 + 10*t1 - Power(t1,2)) - 2*t2) + 
               Power(s,4)*(t1 - t2) - 30*t2 - 10*s2*t2 + 
               7*Power(s2,2)*t2 - 2*Power(s2,3)*t2 - Power(s2,4)*t2 + 
               20*t1*t2 - 10*s2*t1*t2 + 3*Power(s2,2)*t1*t2 + 
               3*Power(s2,3)*t1*t2 + 3*Power(t1,2)*t2 - 
               3*Power(s2,2)*Power(t1,2)*t2 - Power(t1,3)*t2 + 
               s2*Power(t1,3)*t2 - 2*Power(t2,2) + 10*s2*Power(t2,2) - 
               10*t1*Power(t2,2) - 2*s2*t1*Power(t2,2) + 
               2*Power(t1,2)*Power(t2,2) + 2*Power(t2,3) + 
               Power(s,3)*(3 - 4*t1 + 3*Power(t1,2) - 
                  s1*(4 + s2 + 2*t1 - t2) + 4*t2 - 3*t1*t2 + 
                  s2*(2 - 3*t1 + 4*t2)) + 
               Power(s,2)*(-1 + 10*t1 - 11*Power(t1,2) + 
                  3*Power(t1,3) + Power(s1,2)*(-4 + s2 + t1) + 
                  Power(s2,2)*(-4 + 3*t1 - 6*t2) - 5*t2 + 11*t1*t2 - 
                  3*Power(t1,2)*t2 + 
                  s1*(5 + 3*Power(s2,2) - 6*t1 - 3*Power(t1,2) + 
                     s2*(11 - 3*t2) + t2) + 
                  s2*(-11 + 15*t1 - 6*Power(t1,2) - 10*t2 + 9*t1*t2)) \
+ s*(-34 + 2*Power(s1,3) + 10*t1 + 13*Power(t1,2) - 8*Power(t1,3) + 
                  Power(t1,4) - 
                  2*Power(s1,2)*
                   (6 - 4*s2 + Power(s2,2) + 5*t1 - Power(t1,2)) + 
                  2*t2 - 8*t1*t2 + 6*Power(t1,2)*t2 - 
                  Power(t1,3)*t2 - 4*Power(t2,2) + 2*t1*Power(t2,2) + 
                  Power(s2,3)*(2 - t1 + 4*t2) + 
                  Power(s2,2)*
                   (11 + 3*Power(t1,2) + 8*t2 - 3*t1*(4 + 3*t2)) + 
                  s1*(36 - 3*Power(s2,3) + 
                     s2*(4 + 20*t1 - 3*Power(t1,2) - 8*t2) + 12*t2 - 
                     2*Power(t2,2) - 3*Power(t1,2)*(2 + t2) + 
                     Power(s2,2)*(-14 + 6*t1 + 3*t2) + 
                     2*t1*(-7 + 5*t2)) + 
                  s2*(2 - 3*Power(t1,3) + 4*t2 + 
                     6*Power(t1,2)*(3 + t2) - 2*t1*(12 + 7*t2))) + 
               s1*(Power(s2,4) + Power(t1,4) + 6*(-1 + t2) - 
                  2*Power(t1,3)*(2 + t2) - 
                  Power(s2,3)*(-7 + 4*t1 + t2) + 
                  Power(s2,2)*(-5 - 18*t1 + 6*Power(t1,2) + 7*t2) + 
                  Power(t1,2)*(-13 + 9*t2) + 
                  t1*(28 + 6*t2 - 2*Power(t2,2)) + 
                  s2*(-4*Power(t1,3) + 3*Power(t1,2)*(5 + t2) - 
                     2*t1*(-9 + 8*t2) + 2*(-13 - 4*t2 + Power(t2,2)))))\
)/((s - s2 + t1)*(1 + s*(-1 + s1) + s1*(-s2 + t1) - t2)*(s - s1 + t2)) + 
          (8*(18 + 34*s2 - 5*Power(s2,2) - Power(s2,3) - 
               6*Power(s1,3)*(s2 - t1) - 40*t1 - 2*s2*t1 + 
               4*Power(s2,2)*t1 + Power(s2,3)*t1 + 7*Power(t1,2) - 
               5*s2*Power(t1,2) - 3*Power(s2,2)*Power(t1,2) + 
               2*Power(t1,3) + 3*s2*Power(t1,3) - Power(t1,4) + 
               Power(s,4)*(t1 - t2) + 2*t2 - 38*s2*t2 + 
               7*Power(s2,2)*t2 - Power(s2,4)*t2 + 34*t1*t2 - 
               2*s2*t1*t2 - Power(s2,2)*t1*t2 + 3*Power(s2,3)*t1*t2 - 
               5*Power(t1,2)*t2 + 2*s2*Power(t1,2)*t2 - 
               3*Power(s2,2)*Power(t1,2)*t2 - Power(t1,3)*t2 + 
               s2*Power(t1,3)*t2 - 16*Power(t2,2) + 8*s2*Power(t2,2) + 
               2*Power(s2,2)*Power(t2,2) + 2*t1*Power(t2,2) - 
               4*s2*t1*Power(t2,2) + 2*Power(t1,2)*Power(t2,2) - 
               4*Power(t2,3) + 
               Power(s,3)*(1 - 4*t1 + 3*Power(t1,2) - 
                  s1*(2 + s2 + 2*t1 - t2) + 4*t2 - 3*t1*t2 + 
                  s2*(2 - 3*t1 + 4*t2)) + 
               Power(s,2)*(13 + 6*t1 - 9*Power(t1,2) + 
                  3*Power(t1,3) + Power(s1,2)*(2 + s2 + t1) + 
                  Power(s2,2)*(-4 + 3*t1 - 6*t2) - 9*t2 + 7*t1*t2 - 
                  3*Power(t1,2)*t2 + 2*Power(t2,2) + 
                  s1*(-13 + 3*Power(s2,2) + 2*t1 - 3*Power(t1,2) + 
                     s2*(5 - 3*t2) + t2) + 
                  s2*(-7 + 13*t1 - 6*Power(t1,2) - 8*t2 + 9*t1*t2)) + 
               s1*(Power(s2,4) + Power(t1,4) - 
                  2*Power(t1,3)*(-3 + t2) - 
                  Power(s2,3)*(-5 + 4*t1 + t2) + 
                  Power(s2,2)*(-35 - 4*t1 + 6*Power(t1,2) + t2) - 
                  Power(t1,2)*(43 + 3*t2) + 
                  2*t1*(15 + 6*t2 + 2*Power(t2,2)) + 
                  2*(-7 + 2*t2 + 5*Power(t2,2)) - 
                  s2*(22 + 4*Power(t1,3) + Power(t1,2)*(7 - 3*t2) + 
                     20*t2 + 4*Power(t2,2) - 2*t1*(39 + t2))) + 
               Power(s1,2)*(6 + Power(s2,3) - Power(s2,2)*(-4 + t1) + 
                  8*Power(t1,2) + Power(t1,3) - 6*t2 - 
                  2*t1*(7 + 5*t2) - 
                  s2*(12*t1 + Power(t1,2) - 2*(7 + 5*t2))) + 
               s*(-32 + 6*Power(s1,3) + 40*t1 - 3*Power(t1,2) - 
                  6*Power(t1,3) + Power(t1,4) - 2*t2 + 10*t1*t2 + 
                  2*Power(t1,2)*t2 - Power(t1,3)*t2 - 16*Power(t2,2) + 
                  4*t1*Power(t2,2) + Power(s2,3)*(2 - t1 + 4*t2) - 
                  2*Power(s1,2)*
                   (12 + 3*s2 + Power(s2,2) - 5*t1 - Power(t1,2) + 
                     5*t2) + 
                  s1*(34 - 3*Power(s2,3) + Power(t1,2)*(10 - 3*t2) + 
                     40*t2 + 4*Power(t2,2) - 2*t1*(37 + t2) - 
                     s2*(-58 + 2*t1 + 3*Power(t1,2) + 2*t2) + 
                     Power(s2,2)*(-8 + 6*t1 + 3*t2)) + 
                  Power(s2,2)*
                   (7 + 3*Power(t1,2) + 4*t2 - t1*(10 + 9*t2)) - 
                  s2*(3*Power(t1,3) - 2*Power(t1,2)*(7 + 3*t2) + 
                     t1*(4 + 6*t2) + 2*(16 + t2 + 2*Power(t2,2))))))/
           ((-1 + s1)*(-s + s1 - t2)*
             (1 + s*(-1 + s1) + s1*(-s2 + t1) - t2)) + 
          (8*(16 + 47*s2 + 8*Power(s2,2) - 4*Power(s2,3) - 39*t1 - 
               13*s2*t1 + 6*Power(s2,2)*t1 - 2*Power(s2,3)*t1 + 
               5*Power(t1,2) + 6*Power(s2,2)*Power(t1,2) - 
               2*Power(t1,3) - 6*s2*Power(t1,3) + 2*Power(t1,4) + 
               2*Power(s1,3)*
                (Power(s2,2) - 2*s2*(2 + t1) + t1*(4 + t1)) + 
               Power(s1,2)*(8 + 7*Power(s2,3) - 5*Power(t1,3) + 
                  Power(s2,2)*(-19*t1 + 2*(-5 + t2)) + 
                  s2*(8 + 17*Power(t1,2) - 4*t1*(-3 + t2) - 4*t2) + 
                  4*t1*(-2 + t2) + 2*Power(t1,2)*(-1 + t2) - 8*t2) - 
               10*t2 - 52*s2*t2 - 11*Power(s2,2)*t2 + 
               8*Power(s2,3)*t2 + 2*Power(s2,4)*t2 + 52*t1*t2 + 
               20*s2*t1*t2 - 18*Power(s2,2)*t1*t2 - 
               6*Power(s2,3)*t1*t2 - 9*Power(t1,2)*t2 + 
               12*s2*Power(t1,2)*t2 + 6*Power(s2,2)*Power(t1,2)*t2 - 
               2*Power(t1,3)*t2 - 2*s2*Power(t1,3)*t2 - 
               12*Power(t2,2) + s2*Power(t2,2) + 
               Power(s2,2)*Power(t2,2) - 9*t1*Power(t2,2) - 
               3*s2*t1*Power(t2,2) + 2*Power(t1,2)*Power(t2,2) + 
               6*Power(t2,3) + 2*Power(s,4)*(-2 + 2*s1 - t1 + t2) - 
               Power(s,3)*(-4 + 6*Power(s1,2) + t1 + 6*Power(t1,2) - 
                  2*s2*(5 + 3*t1 - 4*t2) + 9*t2 - 6*t1*t2 + 
                  s1*(-8 + 16*s2 - 15*t1 + 3*t2)) + 
               s1*(4*Power(s2,4) + 2*Power(t1,4) + 
                  Power(t1,3)*(5 - 2*t2) + 
                  Power(s2,3)*(-8 - 14*t1 + 3*t2) + 
                  2*Power(t1,2)*(-23 + 7*t2) + 
                  Power(s2,2)*
                   (-54 + 18*Power(t1,2) + t1*(21 - 8*t2) + 10*t2) + 
                  t1*(18 + 10*t2 - 8*Power(t2,2)) - 
                  2*(5 - 8*t2 + 3*Power(t2,2)) - 
                  s2*(6 + 10*Power(t1,3) + Power(t1,2)*(18 - 7*t2) + 
                     22*t2 - 8*Power(t2,2) + 4*t1*(-25 + 6*t2))) + 
               Power(s,2)*(17 + 2*Power(s1,3) - 8*t1 + 
                  11*Power(t1,2) - 6*Power(t1,3) - 
                  2*Power(s2,2)*(4 + 3*t1 - 6*t2) + 2*t2 - 18*t1*t2 + 
                  6*Power(t1,2)*t2 + Power(t2,2) + 
                  Power(s1,2)*(-8 + 19*s2 - 17*t1 + 2*t2) + 
                  s2*(-4 + 12*Power(t1,2) + 26*t2 - 
                     3*t1*(1 + 6*t2)) + 
                  s1*(-21 + 24*Power(s2,2) + 22*t1 + 20*Power(t1,2) + 
                     5*t2 - 8*t1*t2 + s2*(-21 - 44*t1 + 9*t2))) + 
               s*(-33 + 46*t1 - 8*Power(t1,2) + 10*Power(t1,3) - 
                  2*Power(t1,4) + Power(s1,3)*(8 - 4*s2 + 4*t1) + 
                  2*Power(s2,3)*(1 + t1 - 4*t2) + 16*t2 - 25*t1*t2 - 
                  11*Power(t1,2)*t2 + 2*Power(t1,3)*t2 + 
                  9*Power(t2,2) + 3*t1*Power(t2,2) - 
                  2*Power(s1,2)*
                   (10 + 10*Power(s2,2) + 5*t1 + 8*Power(t1,2) - 
                     2*t2 - 2*t1*t2 + s2*(-9 - 18*t1 + 2*t2)) + 
                  Power(s2,2)*
                   (2 - 6*Power(t1,2) - 25*t2 + 6*t1*(1 + 3*t2)) + 
                  s2*(-49 + 6*Power(t1,3) + 11*t2 - 2*Power(t2,2) - 
                     6*Power(t1,2)*(3 + 2*t2) + 6*t1*(1 + 6*t2)) + 
                  s1*(22 - 16*Power(s2,3) + 11*Power(t1,3) + 
                     Power(s2,2)*(21 + 43*t1 - 9*t2) + 
                     Power(t1,2)*(19 - 7*t2) + 10*t2 - 8*Power(t2,2) + 
                     t1*(-73 + 19*t2) + 
                     s2*(69 - 38*Power(t1,2) - 15*t2 + 
                        8*t1*(-5 + 2*t2))))))/
           ((-1 + s1)*(1 + s*(-1 + s1) + s1*(-s2 + t1) - t2)*(-1 + t2)) \
+ (4*(64*(s - s2 + t1)*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) - 
               (2*(-19 + 3*Power(s,2) + 6*Power(s1,2) - 3*s2 + 
                    Power(s2,2) - t1 - Power(t1,2) + 
                    s1*(-8 + 6*s2 - t1 - 5*t2) + 10*t2 - 3*s2*t2 + 
                    2*t1*t2 + s*(11 - 9*s1 - 4*s2 + t1 + 3*t2))*
                  (-s2 + 2*Power(s1,2)*(s2 - t1) + t1 + s2*t1 - 
                    Power(t1,2) + Power(s,2)*(-2 + 2*s1 + t1 - t2) + 
                    2*t2 - s2*t2 - Power(s2,2)*t2 + t1*t2 + 
                    s2*t1*t2 - 2*Power(t2,2) + 
                    s1*(-2 + Power(s2,2) - Power(t1,2) + 2*t2 - 
                       2*s2*t2 + 2*t1*t2) + 
                    s*(3 - 2*Power(s1,2) - 3*t1 + Power(t1,2) - t2 - 
                       t1*t2 + s2*(2 - t1 + 2*t2) + 
                       s1*(-3*s2 + t1 + 2*t2))))/
                (1 + s*(-1 + s1) + s1*(-s2 + t1) - t2) + 
               32*(2 - 2*Power(s1,3) + 6*s2 - 9*Power(s2,2) - 9*t1 + 
                  13*s2*t1 + 11*Power(s2,2)*t1 - 3*Power(t1,2) - 
                  23*s2*Power(t1,2) + 12*Power(t1,3) - 3*t2 + 
                  4*s2*t2 - 15*Power(s2,2)*t2 + 39*s2*t1*t2 - 
                  25*Power(t1,2)*t2 + Power(t2,2) - 
                  14*s2*Power(t2,2) + 13*t1*Power(t2,2) + 
                  Power(s,2)*(-2 + 2*s2 - t1 + t2) + 
                  Power(s1,2)*(2 + s2 - 3*t1 + 4*t2) + 
                  s1*(15*Power(s2,2) + 13*Power(t1,2) + 
                     t1*(15 - 11*t2) - 2*t2*(1 + t2) + 
                     s2*(-15 - 28*t1 + 13*t2)) + 
                  s*(-11 + 2*Power(s1,2) - 2*Power(s2,2) - 5*t1 + 
                     12*Power(t1,2) + 
                     s1*(13 - 17*s2 + 17*t1 - 17*t2) + t2 - 
                     26*t1*t2 + 14*Power(t2,2) + 
                     s2*(15 - 10*t1 + 14*t2))) + 
               16*(-9 + 17*Power(s1,3) - 14*s2 + 29*Power(s2,2) - 
                  Power(s2,3) + 29*t1 - 50*s2*t1 - 
                  43*Power(s2,2)*t1 + 15*Power(t1,2) + 
                  96*s2*Power(t1,2) - 52*Power(t1,3) + 
                  Power(s,2)*(16 + s1 - 18*s2 + 13*t1 - 13*t2) + 
                  17*t2 - 26*s2*t2 + 76*Power(s2,2)*t2 + 2*t1*t2 - 
                  186*s2*t1*t2 + 115*Power(t1,2)*t2 - 2*Power(t2,2) + 
                  72*s2*Power(t2,2) - 63*t1*Power(t2,2) - 
                  Power(s1,2)*(9 + 7*s2 - 25*t1 + 34*t2) + 
                  s1*(-2 - 79*Power(s2,2) - 63*Power(t1,2) + 
                     s2*(72 + 143*t1 - 67*t2) + 5*t2 + 
                     17*Power(t2,2) + 23*t1*(-3 + 2*t2)) + 
                  s*(40 - 18*Power(s1,2) + 19*Power(s2,2) + 28*t1 - 
                     48*Power(t1,2) + 3*s2*(-22 + 9*t1 - 21*t2) + 
                     3*t2 + 120*t1*t2 - 72*Power(t2,2) + 
                     s1*(-63 + 92*s2 - 97*t1 + 98*t2))) - 
               4*(-2 + 3*Power(s,3) - 38*Power(s1,3) + 10*s2 - 
                  23*Power(s2,2) + 3*Power(s2,3) - 36*t1 + 35*s2*t1 + 
                  52*Power(s2,2)*t1 - 11*Power(t1,2) - 
                  119*s2*Power(t1,2) + 64*Power(t1,3) + 33*t2 + 
                  33*s2*t2 - 121*Power(s2,2)*t2 + 6*t1*t2 + 
                  281*s2*t1*t2 - 159*Power(t1,2)*t2 - 
                  36*Power(t2,2) - 117*s2*Power(t2,2) + 
                  95*t1*Power(t2,2) + 
                  Power(s,2)*(-45 - 4*s1 + 35*s2 - 27*t1 + 34*t2) + 
                  Power(s1,2)*(-20 + 17*s2 - 63*t1 + 78*t2) + 
                  s*(-8 + 39*Power(s1,2) - 41*Power(s2,2) - 29*t1 + 
                     57*Power(t1,2) - 57*t2 - 170*t1*t2 + 
                     117*Power(t2,2) - 
                     2*s1*(-62 + 75*s2 - 84*t1 + 85*t2) + 
                     s2*(80 - 9*t1 + 87*t2)) + 
                  s1*(-38 + 132*Power(s2,2) + 85*Power(t1,2) + 
                     t1*(75 - 44*t2) + 55*t2 - 40*Power(t2,2) + 
                     s2*(-110 - 225*t1 + 108*t2))) + 
               8*(19 + 3*Power(s,3) - 48*Power(s1,3) + 33*s2 - 
                  43*Power(s2,2) + 4*Power(s2,3) - 56*t1 + 72*s2*t1 + 
                  75*Power(s2,2)*t1 - 21*Power(t1,2) - 
                  175*s2*Power(t1,2) + 96*Power(t1,3) + 
                  Power(s,2)*(-10*s1 + 43*s2 - 6*(7 + 6*t1 - 7*t2)) - 
                  17*t2 + 50*s2*t2 - 159*Power(s2,2)*t2 - 3*t1*t2 + 
                  381*s2*t1*t2 - 226*Power(t1,2)*t2 - 
                  14*Power(t2,2) - 156*s2*Power(t2,2) + 
                  130*t1*Power(t2,2) + 
                  2*Power(s1,2)*(1 + 7*s2 - 35*t1 + 48*t2) + 
                  s*(-72 + 55*Power(s1,2) - 50*Power(s2,2) - 49*t1 + 
                     84*Power(t1,2) + 
                     s1*(138 - 192*s2 + 219*t1 - 229*t2) - 31*t2 - 
                     237*t1*t2 + 156*Power(t2,2) + 
                     s2*(118 - 25*t1 + 117*t2)) + 
                  s1*(4 + 171*Power(s2,2) + 125*Power(t1,2) + 
                     t1*(125 - 78*t2) + 18*t2 - 48*Power(t2,2) + 
                     s2*(-147 - 304*t1 + 150*t2))) + 
               (2*(4 + s - 2*s1 - s2 + t1)*
                  (Power(s2,2) + Power(s2,3) - 2*s2*t1 - 
                    5*Power(s2,2)*t1 - 2*Power(s2,3)*t1 + 
                    Power(t1,2) + 7*s2*Power(t1,2) + 
                    7*Power(s2,2)*Power(t1,2) + 
                    Power(s2,3)*Power(t1,2) - 3*Power(t1,3) - 
                    8*s2*Power(t1,3) - 3*Power(s2,2)*Power(t1,3) + 
                    3*Power(t1,4) + 3*s2*Power(t1,4) - Power(t1,5) + 
                    4*s2*t2 - 5*Power(s2,2)*t2 + 2*Power(s2,3)*t2 + 
                    2*Power(s2,4)*t2 - 4*t1*t2 + 4*s2*t1*t2 - 
                    8*Power(s2,3)*t1*t2 - 2*Power(s2,4)*t1*t2 + 
                    Power(t1,2)*t2 - 4*s2*Power(t1,2)*t2 + 
                    11*Power(s2,2)*Power(t1,2)*t2 + 
                    6*Power(s2,3)*Power(t1,2)*t2 + 2*Power(t1,3)*t2 - 
                    6*s2*Power(t1,3)*t2 - 
                    6*Power(s2,2)*Power(t1,3)*t2 + Power(t1,4)*t2 + 
                    2*s2*Power(t1,4)*t2 - 2*Power(t2,2) - 
                    4*s2*Power(t2,2) + 7*Power(s2,2)*Power(t2,2) - 
                    3*Power(s2,3)*Power(t2,2) + 
                    Power(s2,4)*Power(t2,2) + 
                    Power(s2,5)*Power(t2,2) + 4*t1*Power(t2,2) - 
                    Power(s2,2)*t1*Power(t2,2) - 
                    4*Power(s2,3)*t1*Power(t2,2) - 
                    3*Power(s2,4)*t1*Power(t2,2) - 
                    7*Power(t1,2)*Power(t2,2) + 
                    7*s2*Power(t1,2)*Power(t2,2) + 
                    5*Power(s2,2)*Power(t1,2)*Power(t2,2) + 
                    3*Power(s2,3)*Power(t1,2)*Power(t2,2) - 
                    3*Power(t1,3)*Power(t2,2) - 
                    2*s2*Power(t1,3)*Power(t2,2) - 
                    Power(s2,2)*Power(t1,3)*Power(t2,2) + 
                    6*Power(t2,3) - 4*s2*Power(t2,3) - 
                    5*Power(s2,2)*Power(t2,3) + 
                    4*Power(s2,3)*Power(t2,3) + 
                    Power(s2,4)*Power(t2,3) + 4*t1*Power(t2,3) - 
                    6*Power(s2,2)*t1*Power(t2,3) - 
                    2*Power(s2,3)*t1*Power(t2,3) + 
                    5*Power(t1,2)*Power(t2,3) + 
                    2*s2*Power(t1,2)*Power(t2,3) + 
                    Power(s2,2)*Power(t1,2)*Power(t2,3) - 
                    6*Power(t2,4) + 4*s2*Power(t2,4) + 
                    2*Power(s2,2)*Power(t2,4) - 4*t1*Power(t2,4) - 
                    2*s2*t1*Power(t2,4) + 2*Power(t2,5) + 
                    2*Power(s1,4)*Power(s2 - t1,2)*
                     (-1 + s2 - t1 + t2) + 
                    Power(s,4)*(-1 + s2 - t1 + t2)*
                     (2 + 2*Power(s1,2) + Power(t1,2) + 
                       2*s1*(-2 + t1 - t2) + 2*t2 + Power(t2,2) - 
                       2*t1*(1 + t2)) + 
                    2*Power(s1,3)*(s2 - t1)*
                     (Power(s2,3) + Power(t1,3) + 
                       2*Power(-1 + t2,2) + 2*t1*Power(-1 + t2,2) - 
                       3*Power(t1,2)*t2 - 
                       Power(s2,2)*(2 + t1 + t2) - 
                       s2*(Power(t1,2) + 2*Power(-1 + t2,2) - 
                        2*t1*(1 + 2*t2))) + 
                    Power(s1,2)*
                     (Power(s2,5) - Power(t1,5) + 
                       Power(s2,4)*(1 - 3*t1 - 3*t2) + 
                       2*Power(-1 + t2,3) + 
                       Power(t1,4)*(-1 + 3*t2) + 
                       2*t1*Power(-1 + t2,2)*(-1 + 4*t2) + 
                       Power(t1,3)*(-6 + 6*t2 - 4*Power(t2,2)) + 
                       2*Power(t1,2)*
                       (2 + 4*t2 - 7*Power(t2,2) + Power(t2,3)) + 
                       2*Power(s2,3)*
                       (-2 + 2*Power(t1,2) + 5*t2 - Power(t2,2) + 
                       t1*(-2 + 4*t2)) - 
                       2*Power(s2,2)*
                       (-4 + 2*Power(t1,3) + 
                       2*Power(t1,2)*(-1 + t2) + 5*Power(t2,2) - 
                       Power(t2,3) + t1*(-1 + 7*t2)) + 
                       s2*(3*Power(t1,4) - 4*Power(t1,3)*t2 - 
                        2*Power(-1 + t2,2)*(-1 + 4*t2) + 
                        Power(t1,2)*(8 - 2*t2 + 6*Power(t2,2)) - 
                        4*t1*(3 + 2*t2 - 6*Power(t2,2) + Power(t2,3))\
)) - 2*s1*(-2*Power(t1,4)*(-1 + t2) + Power(s2,5)*t2 + 
                       2*Power(-1 + t2,3)*t2 + 
                       t1*Power(-1 + t2,2)*
                       (-2 - 3*t2 + 2*Power(t2,2)) + 
                       Power(t1,3)*(-2 - 3*t2 + 4*Power(t2,2)) - 
                       2*Power(t1,2)*
                       (-1 + t2 - 2*Power(t2,2) + 2*Power(t2,3)) + 
                       Power(s2,4)*(1 + t2 - t1*(1 + 3*t2)) + 
                       Power(s2,3)*
                       (1 - 4*t2 + 5*Power(t2,2) - Power(t2,3) + 
                       3*Power(t1,2)*(1 + t2) + 
                       t1*(-4 - 3*t2 + Power(t2,2))) - 
                       Power(s2,2)*
                       (1 - 7*t2 + 5*Power(t2,2) + Power(t2,3) + 
                       Power(t1,3)*(3 + t2) + 
                       Power(t1,2)*(-7 - t2 + 2*Power(t2,2)) + 
                       t1*(3 - 3*t2 + 5*Power(t2,2) - 2*Power(t2,3))\
) + s2*(Power(t1,4) + Power(t1,3)*(-6 + 3*t2 + Power(t2,2)) - 
                        Power(-1 + t2,2)*
                       (-2 - 3*t2 + 2*Power(t2,2)) - 
                        Power(t1,2)*
                       (-4 - 4*t2 + 4*Power(t2,2) + Power(t2,3)) + 
                        t1*(-1 - 5*t2 + Power(t2,2) + 5*Power(t2,3)))\
) - 2*Power(s,3)*(-3 + 2*t1 + 2*Power(t1,2) - 3*Power(t1,3) + 
                       Power(t1,4) + 3*t2 - 2*t1*t2 + 
                       7*Power(t1,2)*t2 - 3*Power(t1,3)*t2 - 
                       5*t1*Power(t2,2) + 
                       3*Power(t1,2)*Power(t2,2) + Power(t2,3) - 
                       t1*Power(t2,3) + 
                       2*Power(s1,3)*(-1 + s2 - t1 + t2) + 
                       Power(s1,2)*(-2 + 3*s2 - 3*t2)*
                       (-1 + s2 - t1 + t2) + 
                       Power(s2,2)*
                       (Power(t1,2) + 2*Power(1 + t2,2) - 
                       t1*(2 + 3*t2)) + 
                       s2*(1 - 2*Power(t1,3) - 2*t2 + 
                       3*Power(t2,2) + 2*Power(t2,3) + 
                       Power(t1,2)*(5 + 6*t2) - 
                       t1*(5 + 8*t2 + 6*Power(t2,2))) + 
                       s1*(3 + Power(t1,3) + 
                        Power(s2,2)*(-5 + 2*t1 - 4*t2) - 6*t2 + 
                        Power(t2,2) + Power(t2,3) - 
                        Power(t1,2)*(5 + t2) - 
                        t1*(2 - 4*t2 + Power(t2,2)) + 
                        s2*(2 - 3*Power(t1,2) + 2*t2 - 
                        3*Power(t2,2) + t1*(8 + 6*t2)))) + 
                    Power(s,2)*
                     (-5 + 7*t1 + 3*Power(t1,2) - 11*Power(t1,3) + 
                       7*Power(t1,4) - Power(t1,5) + 13*t2 - 
                       6*t1*t2 + 15*Power(t1,2)*t2 - 
                       18*Power(t1,3)*t2 + 3*Power(t1,4)*t2 - 
                       15*Power(t2,2) - 9*t1*Power(t2,2) + 
                       17*Power(t1,2)*Power(t2,2) - 
                       3*Power(t1,3)*Power(t2,2) + 5*Power(t2,3) - 
                       8*t1*Power(t2,3) + Power(t1,2)*Power(t2,3) + 
                       2*Power(t2,4) + 
                       2*Power(s1,4)*(-1 + s2 - t1 + t2) + 
                       Power(s2,3)*
                       (2 + Power(t1,2) + 10*t2 + 6*Power(t2,2) - 
                       2*t1*(1 + 3*t2)) + 
                       2*Power(s1,3)*
                       (5*Power(s2,2) + 3*Power(t1,2) - 
                       t1*(-3 + t2) - 2*(-1 + t2)*t2 + 
                       s2*(-5 - 8*t1 + 3*t2)) + 
                       Power(s2,2)*
                       (6 - 3*Power(t1,3) + 6*t2 + 10*Power(t2,2) + 
                       6*Power(t2,3) + Power(t1,2)*(11 + 15*t2) - 
                       2*t1*(7 + 17*t2 + 9*Power(t2,2))) + 
                       Power(s1,2)*
                       (7*Power(s2,3) + 5*Power(t1,3) - 
                       3*Power(t1,2)*(3 + 5*t2) - 
                       Power(s2,2)*(15 + 7*t1 + 9*t2) + 
                       s2*(-2 + 22*t1 - 5*Power(t1,2) + 18*t2 + 
                       26*t1*t2 - 14*Power(t2,2)) + 
                       t1*(-8 - 2*t2 + 8*Power(t2,2)) + 
                       2*(5 - 8*t2 + 2*Power(t2,2) + Power(t2,3))) + 
                       s2*(-3 + 3*Power(t1,4) - 10*t2 + 
                       21*Power(t2,2) + 8*Power(t2,3) - 
                       4*Power(t1,3)*(4 + 3*t2) + 
                       Power(t1,2)*(23 + 42*t2 + 15*Power(t2,2)) - 
                       2*t1*
                       (5 + 7*t2 + 17*Power(t2,2) + 3*Power(t2,3))) \
+ 2*s1*(Power(t1,4) + Power(s2,3)*(-4 + t1 - 6*t2) - 
                        4*Power(t1,3)*(-1 + t2) - 
                        2*Power(-1 + t2,2)*(1 + 2*t2) + 
                        Power(t1,2)*(-4 - t2 + 5*Power(t2,2)) + 
                        t1*(-2 + 10*t2 + Power(t2,2) - 
                       2*Power(t2,3)) - 
                        Power(s2,2)*
                       (2 + Power(t1,2) - 5*t2 + 3*Power(t2,2) - 
                       t1*(10 + 9*t2)) - 
                        s2*(-8 + Power(t1,3) - 
                        Power(t1,2)*(-10 + t2) + 11*t2 + 
                        7*Power(t2,2) - 3*Power(t2,3) + 
                        3*t1*(-1 + Power(t2,2))))) - 
                    2*s*(t1 - 4*Power(t1,3) + 4*Power(t1,4) - 
                       Power(t1,5) + 4*t2 - 2*t1*t2 + 
                       Power(t1,2)*t2 - 6*Power(t1,3)*t2 + 
                       3*Power(t1,4)*t2 - 9*Power(t2,2) + 
                       4*t1*Power(t2,2) + 5*Power(t1,2)*Power(t2,2) - 
                       4*Power(t1,3)*Power(t2,2) + 6*Power(t2,3) - 
                       2*t1*Power(t2,3) + 3*Power(t1,2)*Power(t2,3) - 
                       Power(t2,4) - t1*Power(t2,4) + 
                       2*Power(s1,4)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                       Power(s2,4)*t2*(2 - t1 + 2*t2) + 
                       Power(s2,3)*
                        (1 + 6*t2 + 3*Power(t2,2) + 2*Power(t2,3) + 
                        Power(t1,2)*(1 + 3*t2) - 
                        t1*(2 + 11*t2 + 6*Power(t2,2))) + 
                       Power(s2,2)*
                        (1 - 5*t2 + 8*Power(t2,2) + 5*Power(t2,3) - 
                        3*Power(t1,3)*(1 + t2) + 
                        Power(t1,2)*(8 + 19*t2 + 6*Power(t2,2)) - 
                        t1*(6 + 15*t2 + 14*Power(t2,2) + 
                        3*Power(t2,3))) + 
                       s2*(-1 - t2 + 3*Power(t2,2) - 3*Power(t2,3) + 
                        2*Power(t2,4) + Power(t1,4)*(3 + t2) - 
                        Power(t1,3)*(10 + 13*t2 + 2*Power(t2,2)) + 
                        Power(t1,2)*
                       (9 + 15*t2 + 15*Power(t2,2) + Power(t2,3)) - 
                        t1*(1 - 5*t2 + 15*Power(t2,2) + 
                       7*Power(t2,3))) + 
                       Power(s1,3)*
                        (4*Power(s2,3) - Power(s2,2)*(5 + 8*t1) + 
                        2*Power(-1 + t2,2) - 
                        Power(t1,2)*(1 + 4*t2) + 
                        t1*(2 - 6*t2 + 4*Power(t2,2)) + 
                        s2*(-2 + 4*Power(t1,2) + 6*t2 - 
                        4*Power(t2,2) + t1*(6 + 4*t2))) + 
                       Power(s1,2)*
                        (2*Power(s2,4) - Power(t1,4) - 
                        4*Power(-1 + t2,2)*t2 + 
                        2*Power(t1,3)*(2 + t2) - 
                        Power(s2,3)*(3 + 3*t1 + 5*t2) + 
                        Power(t1,2)*(5 - 3*t2 + Power(t2,2)) + 
                        t1*(-7 + 4*t2 + 5*Power(t2,2) - 
                       2*Power(t2,3)) - 
                        Power(s2,2)*
                       (5 + Power(t1,2) - 13*t2 + 5*Power(t2,2) - 
                       2*t1*(4 + 7*t2)) + 
                        s2*(9 + 3*Power(t1,3) - 8*t2 - 
                        3*Power(t2,2) + 2*Power(t2,3) + 
                        2*t1*t2*(-5 + 2*t2) - Power(t1,2)*(9 + 11*t2)\
)) - s1*(Power(t1,5) - 3*Power(t1,4)*t2 + Power(s2,4)*(1 + 4*t2) - 
                        Power(-1 + t2,2)*(-4 + t2 + 2*Power(t2,2)) + 
                        Power(t1,3)*(-6 + 4*t2 + 3*Power(t2,2)) + 
                        t1*(-2 - 2*t2 + 4*Power(t2,2)) - 
                        Power(t1,2)*
                        (-3 - 7*t2 + 2*Power(t2,2) + Power(t2,3)) - 
                        Power(s2,3)*
                        (-3 + Power(t1,2) + t2 - Power(t2,2) + 
                        t1*(4 + 8*t2)) + 
                        Power(s2,2)*
                        (-4 + 3*Power(t1,3) + t2 + 13*Power(t2,2) - 
                        3*Power(t2,3) + Power(t1,2)*(5 + t2) + 
                        t1*(-8 + 3*Power(t2,2))) - 
                        s2*(1 + 3*Power(t1,4) + 
                        Power(t1,3)*(2 - 6*t2) - 11*t2 + 
                        13*Power(t2,2) - 3*Power(t2,3) + 
                        Power(t1,2)*(-11 + 3*t2 + 7*Power(t2,2)) - 
                        2*t1*
                        (1 - 5*t2 - 5*Power(t2,2) + 2*Power(t2,3))))))\
)/(Power(-1 + s - s*s1 + s1*(s2 - t1) + t2,2)*(-1 + s2 - t1 + t2))))/
           ((-1 + s2)*(-1 + t1)) + 
          (8*(-32*(s - s2 + t1)*(s1 - t2)*(-1 + s2 - t1 + t2) - 
               ((-14 + 6*Power(s,2) + 3*Power(s1,2) - 8*s2 + 
                    5*Power(s2,2) + 3*t1 - 5*s2*t1 + 
                    s1*(-8 + 9*s2 - 6*t1 - t2) + 6*t2 - 3*s2*t2 + 
                    3*t1*t2 + s*(9 - 9*s1 - 11*s2 + 6*t1 + 3*t2))*
                  (-s2 + 2*Power(s1,2)*(s2 - t1) + t1 + s2*t1 - 
                    Power(t1,2) + Power(s,2)*(-2 + 2*s1 + t1 - t2) + 
                    2*t2 - s2*t2 - Power(s2,2)*t2 + t1*t2 + 
                    s2*t1*t2 - 2*Power(t2,2) + 
                    s1*(-2 + Power(s2,2) - Power(t1,2) + 2*t2 - 
                       2*s2*t2 + 2*t1*t2) + 
                    s*(3 - 2*Power(s1,2) - 3*t1 + Power(t1,2) - t2 - 
                       t1*t2 + s2*(2 - t1 + 2*t2) + 
                       s1*(-3*s2 + t1 + 2*t2))))/
                (1 + s*(-1 + s1) + s1*(-s2 + t1) - t2) + 
               8*(4 + Power(s1,3) + 20*s2 - 10*Power(s2,2) + 
                  Power(s2,3) + Power(s1,2)*(19*s2 - 17*t1) - 24*t1 + 
                  10*s2*t1 + 5*Power(s2,2)*t1 + 2*Power(t1,2) - 
                  15*s2*Power(t1,2) + 9*Power(t1,3) + 
                  Power(s,2)*(-19 + 17*s1 + 9*s2 - 6*t2) + 9*t2 + 
                  40*s2*t2 - 68*Power(s2,2)*t2 - 35*t1*t2 + 
                  136*s2*t1*t2 - 71*Power(t1,2)*t2 - 14*Power(t2,2) - 
                  54*s2*Power(t2,2) + 55*t1*Power(t2,2) + 
                  s1*(-15 + 71*Power(s2,2) + 55*Power(t1,2) + 
                     t1*(70 - 40*t2) + 16*t2 - Power(t2,2) + 
                     s2*(-71 - 125*t1 + 37*t2)) - 
                  s*(22 + 18*Power(s1,2) + 10*Power(s2,2) + 21*t1 - 
                     13*Power(t1,2) + s2*(-43 + 4*t1 - 74*t2) + 
                     40*t2 + 80*t1*t2 - 54*Power(t2,2) + 
                     s1*(-75 + 95*s2 - 76*t1 + 45*t2))) + 
               2*(19 + 2*Power(s,3) + 3*Power(s1,3) + 35*s2 - 
                  9*Power(s2,2) + 6*Power(s2,3) - 36*t1 + 18*s2*t1 + 
                  6*Power(s2,2)*t1 - 2*Power(t1,2) - 
                  32*s2*Power(t1,2) + 20*Power(t1,3) + 
                  Power(s,2)*(-35 + 39*s1 + 22*s2 - 8*t2) - 7*t2 + 
                  42*s2*t2 - 106*Power(s2,2)*t2 - 24*t1*t2 + 
                  208*s2*t1*t2 - 110*Power(t1,2)*t2 - 24*Power(t2,2) - 
                  78*s2*Power(t2,2) + 77*t1*Power(t2,2) + 
                  Power(s1,2)*(7 + 54*s2 - 44*t1 + 6*t2) + 
                  s1*(-1 + 123*Power(s2,2) + 76*Power(t1,2) + 
                     t1*(91 - 39*t2) + 24*t2 - 9*Power(t2,2) + 
                     s2*(-97 - 192*t1 + 35*t2)) + 
                  s*(-52 - 44*Power(s1,2) - 30*Power(s2,2) - 42*t1 + 
                     22*Power(t1,2) + 
                     s1*(87 - 174*s2 + 122*t1 - 56*t2) - 39*t2 - 
                     123*t1*t2 + 78*Power(t2,2) + 
                     s2*(71 + 2*t1 + 114*t2))) + 
               4*(-9 + 2*Power(s,3) - 5*Power(s1,3) - 27*s2 + 
                  23*Power(s2,2) - 7*Power(s2,3) + 42*t1 - 31*s2*t1 - 
                  7*Power(s2,2)*t1 + Power(t1,2) + 38*s2*Power(t1,2) - 
                  24*Power(t1,3) + 
                  Power(s1,2)*(-63*s2 + 51*t1 - 2*t2) - 17*t2 - 
                  74*s2*t2 + 144*Power(s2,2)*t2 + 54*t1*t2 - 
                  283*s2*t1*t2 + 149*Power(t1,2)*t2 + 33*Power(t2,2) + 
                  107*s2*Power(t2,2) - 108*t1*Power(t2,2) + 
                  Power(s,2)*(-53*s1 + 3*(17 - 11*s2 + t1 + 5*t2)) + 
                  s1*(37 - 162*Power(s2,2) - 108*Power(t1,2) + 
                     s2*(150 + 263*t1 - 55*t2) - 40*t2 + 
                     7*Power(t2,2) + t1*(-141 + 64*t2)) + 
                  s*(36 + 56*Power(s1,2) + 38*Power(s2,2) + 58*t1 - 
                     29*Power(t1,2) + 73*t2 + 169*t1*t2 - 
                     107*Power(t2,2) + 
                     s1*(-149 + 231*s2 - 167*t1 + 78*t2) - 
                     s2*(108 + 2*t1 + 159*t2))) - 
               16*(4*s2 - Power(s2,2) + 2*Power(s1,2)*(s2 - t1) - 
                  4*t1 + Power(s2,2)*t1 + Power(t1,2) - 
                  2*s2*Power(t1,2) + Power(t1,3) + 
                  Power(s,2)*(-2 + 2*s1 + s2 - t2) + 2*t2 + 10*s2*t2 - 
                  14*Power(s2,2)*t2 - 10*t1*t2 + 28*s2*t1*t2 - 
                  14*Power(t1,2)*t2 - 2*Power(t2,2) - 
                  12*s2*Power(t2,2) + 12*t1*Power(t2,2) + 
                  2*s1*(-1 + 7*Power(s2,2) + 6*Power(t1,2) + 
                     t1*(7 - 5*t2) + t2 + s2*(-7 - 13*t1 + 5*t2)) - 
                  s*(2*Power(s1,2) + Power(s2,2) + 
                     s2*(-5 + t1 - 15*t2) + 
                     s1*(-15 + 17*s2 - 15*t1 + 11*t2) + 
                     2*(2 + t1 - Power(t1,2) + 5*t2 + 8*t1*t2 - 
                        6*Power(t2,2)))) + 
               ((3 + 2*s - s1 - 2*s2 + 2*t1)*
                  (Power(s2,2) + Power(s2,3) - 2*s2*t1 - 
                    5*Power(s2,2)*t1 - 2*Power(s2,3)*t1 + 
                    Power(t1,2) + 7*s2*Power(t1,2) + 
                    7*Power(s2,2)*Power(t1,2) + 
                    Power(s2,3)*Power(t1,2) - 3*Power(t1,3) - 
                    8*s2*Power(t1,3) - 3*Power(s2,2)*Power(t1,3) + 
                    3*Power(t1,4) + 3*s2*Power(t1,4) - Power(t1,5) + 
                    4*s2*t2 - 5*Power(s2,2)*t2 + 2*Power(s2,3)*t2 + 
                    2*Power(s2,4)*t2 - 4*t1*t2 + 4*s2*t1*t2 - 
                    8*Power(s2,3)*t1*t2 - 2*Power(s2,4)*t1*t2 + 
                    Power(t1,2)*t2 - 4*s2*Power(t1,2)*t2 + 
                    11*Power(s2,2)*Power(t1,2)*t2 + 
                    6*Power(s2,3)*Power(t1,2)*t2 + 2*Power(t1,3)*t2 - 
                    6*s2*Power(t1,3)*t2 - 
                    6*Power(s2,2)*Power(t1,3)*t2 + Power(t1,4)*t2 + 
                    2*s2*Power(t1,4)*t2 - 2*Power(t2,2) - 
                    4*s2*Power(t2,2) + 7*Power(s2,2)*Power(t2,2) - 
                    3*Power(s2,3)*Power(t2,2) + 
                    Power(s2,4)*Power(t2,2) + 
                    Power(s2,5)*Power(t2,2) + 4*t1*Power(t2,2) - 
                    Power(s2,2)*t1*Power(t2,2) - 
                    4*Power(s2,3)*t1*Power(t2,2) - 
                    3*Power(s2,4)*t1*Power(t2,2) - 
                    7*Power(t1,2)*Power(t2,2) + 
                    7*s2*Power(t1,2)*Power(t2,2) + 
                    5*Power(s2,2)*Power(t1,2)*Power(t2,2) + 
                    3*Power(s2,3)*Power(t1,2)*Power(t2,2) - 
                    3*Power(t1,3)*Power(t2,2) - 
                    2*s2*Power(t1,3)*Power(t2,2) - 
                    Power(s2,2)*Power(t1,3)*Power(t2,2) + 
                    6*Power(t2,3) - 4*s2*Power(t2,3) - 
                    5*Power(s2,2)*Power(t2,3) + 
                    4*Power(s2,3)*Power(t2,3) + 
                    Power(s2,4)*Power(t2,3) + 4*t1*Power(t2,3) - 
                    6*Power(s2,2)*t1*Power(t2,3) - 
                    2*Power(s2,3)*t1*Power(t2,3) + 
                    5*Power(t1,2)*Power(t2,3) + 
                    2*s2*Power(t1,2)*Power(t2,3) + 
                    Power(s2,2)*Power(t1,2)*Power(t2,3) - 
                    6*Power(t2,4) + 4*s2*Power(t2,4) + 
                    2*Power(s2,2)*Power(t2,4) - 4*t1*Power(t2,4) - 
                    2*s2*t1*Power(t2,4) + 2*Power(t2,5) + 
                    2*Power(s1,4)*Power(s2 - t1,2)*
                     (-1 + s2 - t1 + t2) + 
                    Power(s,4)*(-1 + s2 - t1 + t2)*
                     (2 + 2*Power(s1,2) + Power(t1,2) + 
                       2*s1*(-2 + t1 - t2) + 2*t2 + Power(t2,2) - 
                       2*t1*(1 + t2)) + 
                    2*Power(s1,3)*(s2 - t1)*
                     (Power(s2,3) + Power(t1,3) + 
                       2*Power(-1 + t2,2) + 2*t1*Power(-1 + t2,2) - 
                       3*Power(t1,2)*t2 - Power(s2,2)*(2 + t1 + t2) - 
                       s2*(Power(t1,2) + 2*Power(-1 + t2,2) - 
                        2*t1*(1 + 2*t2))) + 
                    Power(s1,2)*
                     (Power(s2,5) - Power(t1,5) + 
                       Power(s2,4)*(1 - 3*t1 - 3*t2) + 
                       2*Power(-1 + t2,3) + Power(t1,4)*(-1 + 3*t2) + 
                       2*t1*Power(-1 + t2,2)*(-1 + 4*t2) + 
                       Power(t1,3)*(-6 + 6*t2 - 4*Power(t2,2)) + 
                       2*Power(t1,2)*
                        (2 + 4*t2 - 7*Power(t2,2) + Power(t2,3)) + 
                       2*Power(s2,3)*
                        (-2 + 2*Power(t1,2) + 5*t2 - Power(t2,2) + 
                        t1*(-2 + 4*t2)) - 
                       2*Power(s2,2)*
                        (-4 + 2*Power(t1,3) + 
                        2*Power(t1,2)*(-1 + t2) + 5*Power(t2,2) - 
                        Power(t2,3) + t1*(-1 + 7*t2)) + 
                       s2*(3*Power(t1,4) - 4*Power(t1,3)*t2 - 
                        2*Power(-1 + t2,2)*(-1 + 4*t2) + 
                        Power(t1,2)*(8 - 2*t2 + 6*Power(t2,2)) - 
                        4*t1*(3 + 2*t2 - 6*Power(t2,2) + Power(t2,3)))\
) - 2*s1*(-2*Power(t1,4)*(-1 + t2) + Power(s2,5)*t2 + 
                       2*Power(-1 + t2,3)*t2 + 
                       t1*Power(-1 + t2,2)*
                        (-2 - 3*t2 + 2*Power(t2,2)) + 
                       Power(t1,3)*(-2 - 3*t2 + 4*Power(t2,2)) - 
                       2*Power(t1,2)*
                        (-1 + t2 - 2*Power(t2,2) + 2*Power(t2,3)) + 
                       Power(s2,4)*(1 + t2 - t1*(1 + 3*t2)) + 
                       Power(s2,3)*
                        (1 - 4*t2 + 5*Power(t2,2) - Power(t2,3) + 
                        3*Power(t1,2)*(1 + t2) + 
                        t1*(-4 - 3*t2 + Power(t2,2))) - 
                       Power(s2,2)*
                        (1 - 7*t2 + 5*Power(t2,2) + Power(t2,3) + 
                        Power(t1,3)*(3 + t2) + 
                        Power(t1,2)*(-7 - t2 + 2*Power(t2,2)) + 
                        t1*(3 - 3*t2 + 5*Power(t2,2) - 2*Power(t2,3))\
) + s2*(Power(t1,4) + Power(t1,3)*(-6 + 3*t2 + Power(t2,2)) - 
                        Power(-1 + t2,2)*
                       (-2 - 3*t2 + 2*Power(t2,2)) - 
                        Power(t1,2)*
                        (-4 - 4*t2 + 4*Power(t2,2) + Power(t2,3)) + 
                        t1*(-1 - 5*t2 + Power(t2,2) + 5*Power(t2,3)))) \
- 2*Power(s,3)*(-3 + 2*t1 + 2*Power(t1,2) - 3*Power(t1,3) + 
                       Power(t1,4) + 3*t2 - 2*t1*t2 + 
                       7*Power(t1,2)*t2 - 3*Power(t1,3)*t2 - 
                       5*t1*Power(t2,2) + 3*Power(t1,2)*Power(t2,2) + 
                       Power(t2,3) - t1*Power(t2,3) + 
                       2*Power(s1,3)*(-1 + s2 - t1 + t2) + 
                       Power(s1,2)*(-2 + 3*s2 - 3*t2)*
                        (-1 + s2 - t1 + t2) + 
                       Power(s2,2)*
                        (Power(t1,2) + 2*Power(1 + t2,2) - 
                        t1*(2 + 3*t2)) + 
                       s2*(1 - 2*Power(t1,3) - 2*t2 + 
                        3*Power(t2,2) + 2*Power(t2,3) + 
                        Power(t1,2)*(5 + 6*t2) - 
                        t1*(5 + 8*t2 + 6*Power(t2,2))) + 
                       s1*(3 + Power(t1,3) + 
                        Power(s2,2)*(-5 + 2*t1 - 4*t2) - 6*t2 + 
                        Power(t2,2) + Power(t2,3) - 
                        Power(t1,2)*(5 + t2) - 
                        t1*(2 - 4*t2 + Power(t2,2)) + 
                        s2*(2 - 3*Power(t1,2) + 2*t2 - 
                        3*Power(t2,2) + t1*(8 + 6*t2)))) + 
                    Power(s,2)*
                     (-5 + 7*t1 + 3*Power(t1,2) - 11*Power(t1,3) + 
                       7*Power(t1,4) - Power(t1,5) + 13*t2 - 
                       6*t1*t2 + 15*Power(t1,2)*t2 - 
                       18*Power(t1,3)*t2 + 3*Power(t1,4)*t2 - 
                       15*Power(t2,2) - 9*t1*Power(t2,2) + 
                       17*Power(t1,2)*Power(t2,2) - 
                       3*Power(t1,3)*Power(t2,2) + 5*Power(t2,3) - 
                       8*t1*Power(t2,3) + Power(t1,2)*Power(t2,3) + 
                       2*Power(t2,4) + 
                       2*Power(s1,4)*(-1 + s2 - t1 + t2) + 
                       Power(s2,3)*
                        (2 + Power(t1,2) + 10*t2 + 6*Power(t2,2) - 
                        2*t1*(1 + 3*t2)) + 
                       2*Power(s1,3)*
                        (5*Power(s2,2) + 3*Power(t1,2) - 
                        t1*(-3 + t2) - 2*(-1 + t2)*t2 + 
                        s2*(-5 - 8*t1 + 3*t2)) + 
                       Power(s2,2)*
                        (6 - 3*Power(t1,3) + 6*t2 + 10*Power(t2,2) + 
                        6*Power(t2,3) + Power(t1,2)*(11 + 15*t2) - 
                        2*t1*(7 + 17*t2 + 9*Power(t2,2))) + 
                       Power(s1,2)*
                        (7*Power(s2,3) + 5*Power(t1,3) - 
                        3*Power(t1,2)*(3 + 5*t2) - 
                        Power(s2,2)*(15 + 7*t1 + 9*t2) + 
                        s2*(-2 + 22*t1 - 5*Power(t1,2) + 18*t2 + 
                       26*t1*t2 - 14*Power(t2,2)) + 
                        t1*(-8 - 2*t2 + 8*Power(t2,2)) + 
                        2*(5 - 8*t2 + 2*Power(t2,2) + Power(t2,3))) + 
                       s2*(-3 + 3*Power(t1,4) - 10*t2 + 
                        21*Power(t2,2) + 8*Power(t2,3) - 
                        4*Power(t1,3)*(4 + 3*t2) + 
                        Power(t1,2)*(23 + 42*t2 + 15*Power(t2,2)) - 
                        2*t1*
                        (5 + 7*t2 + 17*Power(t2,2) + 3*Power(t2,3))) \
+ 2*s1*(Power(t1,4) + Power(s2,3)*(-4 + t1 - 6*t2) - 
                        4*Power(t1,3)*(-1 + t2) - 
                        2*Power(-1 + t2,2)*(1 + 2*t2) + 
                        Power(t1,2)*(-4 - t2 + 5*Power(t2,2)) + 
                        t1*(-2 + 10*t2 + Power(t2,2) - 
                       2*Power(t2,3)) - 
                        Power(s2,2)*
                        (2 + Power(t1,2) - 5*t2 + 3*Power(t2,2) - 
                        t1*(10 + 9*t2)) - 
                        s2*(-8 + Power(t1,3) - 
                        Power(t1,2)*(-10 + t2) + 11*t2 + 
                        7*Power(t2,2) - 3*Power(t2,3) + 
                        3*t1*(-1 + Power(t2,2))))) - 
                    2*s*(t1 - 4*Power(t1,3) + 4*Power(t1,4) - 
                       Power(t1,5) + 4*t2 - 2*t1*t2 + Power(t1,2)*t2 - 
                       6*Power(t1,3)*t2 + 3*Power(t1,4)*t2 - 
                       9*Power(t2,2) + 4*t1*Power(t2,2) + 
                       5*Power(t1,2)*Power(t2,2) - 
                       4*Power(t1,3)*Power(t2,2) + 6*Power(t2,3) - 
                       2*t1*Power(t2,3) + 3*Power(t1,2)*Power(t2,3) - 
                       Power(t2,4) - t1*Power(t2,4) + 
                       2*Power(s1,4)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                       Power(s2,4)*t2*(2 - t1 + 2*t2) + 
                       Power(s2,3)*
                        (1 + 6*t2 + 3*Power(t2,2) + 2*Power(t2,3) + 
                        Power(t1,2)*(1 + 3*t2) - 
                        t1*(2 + 11*t2 + 6*Power(t2,2))) + 
                       Power(s2,2)*
                        (1 - 5*t2 + 8*Power(t2,2) + 5*Power(t2,3) - 
                        3*Power(t1,3)*(1 + t2) + 
                        Power(t1,2)*(8 + 19*t2 + 6*Power(t2,2)) - 
                        t1*(6 + 15*t2 + 14*Power(t2,2) + 
                        3*Power(t2,3))) + 
                       s2*(-1 - t2 + 3*Power(t2,2) - 3*Power(t2,3) + 
                        2*Power(t2,4) + Power(t1,4)*(3 + t2) - 
                        Power(t1,3)*(10 + 13*t2 + 2*Power(t2,2)) + 
                        Power(t1,2)*
                        (9 + 15*t2 + 15*Power(t2,2) + Power(t2,3)) - 
                        t1*(1 - 5*t2 + 15*Power(t2,2) + 7*Power(t2,3))\
) + Power(s1,3)*(4*Power(s2,3) - Power(s2,2)*(5 + 8*t1) + 
                        2*Power(-1 + t2,2) - Power(t1,2)*(1 + 4*t2) + 
                        t1*(2 - 6*t2 + 4*Power(t2,2)) + 
                        s2*(-2 + 4*Power(t1,2) + 6*t2 - 
                        4*Power(t2,2) + t1*(6 + 4*t2))) + 
                       Power(s1,2)*
                        (2*Power(s2,4) - Power(t1,4) - 
                        4*Power(-1 + t2,2)*t2 + 
                        2*Power(t1,3)*(2 + t2) - 
                        Power(s2,3)*(3 + 3*t1 + 5*t2) + 
                        Power(t1,2)*(5 - 3*t2 + Power(t2,2)) + 
                        t1*(-7 + 4*t2 + 5*Power(t2,2) - 
                       2*Power(t2,3)) - 
                        Power(s2,2)*
                        (5 + Power(t1,2) - 13*t2 + 5*Power(t2,2) - 
                        2*t1*(4 + 7*t2)) + 
                        s2*(9 + 3*Power(t1,3) - 8*t2 - 
                        3*Power(t2,2) + 2*Power(t2,3) + 
                        2*t1*t2*(-5 + 2*t2) - Power(t1,2)*(9 + 11*t2))\
) - s1*(Power(t1,5) - 3*Power(t1,4)*t2 + Power(s2,4)*(1 + 4*t2) - 
                        Power(-1 + t2,2)*(-4 + t2 + 2*Power(t2,2)) + 
                        Power(t1,3)*(-6 + 4*t2 + 3*Power(t2,2)) + 
                        t1*(-2 - 2*t2 + 4*Power(t2,2)) - 
                        Power(t1,2)*
                        (-3 - 7*t2 + 2*Power(t2,2) + Power(t2,3)) - 
                        Power(s2,3)*
                        (-3 + Power(t1,2) + t2 - Power(t2,2) + 
                        t1*(4 + 8*t2)) + 
                        Power(s2,2)*
                        (-4 + 3*Power(t1,3) + t2 + 13*Power(t2,2) - 
                        3*Power(t2,3) + Power(t1,2)*(5 + t2) + 
                        t1*(-8 + 3*Power(t2,2))) - 
                        s2*(1 + 3*Power(t1,4) + 
                        Power(t1,3)*(2 - 6*t2) - 11*t2 + 
                        13*Power(t2,2) - 3*Power(t2,3) + 
                        Power(t1,2)*(-11 + 3*t2 + 7*Power(t2,2)) - 
                        2*t1*(1 - 5*t2 - 5*Power(t2,2) + 2*Power(t2,3))\
)))))/(Power(-1 + s - s*s1 + s1*(s2 - t1) + t2,2)*(-1 + s2 - t1 + t2))))/
           ((-1 + t1)*(-1 + t2))) + 
       ((64*(-24 + 4*Power(s,4) + 22*s1 + 6*Power(s1,3) - 88*s2 - 
               28*s1*s2 + 4*Power(s1,3)*s2 - 16*Power(s2,2) - 
               4*s1*Power(s2,2) + 8*Power(s2,3) + 88*t1 + 18*s1*t1 - 
               4*Power(s1,3)*t1 + 32*s2*t1 + 4*s1*s2*t1 - 
               24*Power(s2,2)*t1 + 4*s1*Power(s2,2)*t1 - 
               16*Power(t1,2) + 24*s2*Power(t1,2) - 
               8*s1*s2*Power(t1,2) - 8*Power(t1,3) + 
               4*s1*Power(t1,3) - 
               4*Power(s,3)*(3*s1 + 2*s2 - t1 - 2*t2) - 54*t2 + 
               22*s1*t2 + 12*s2*t2 + 8*s1*s2*t2 + 20*Power(s2,2)*t2 - 
               8*s1*Power(s2,2)*t2 - 2*t1*t2 - 18*s1*t1*t2 - 
               36*s2*t1*t2 + 12*s1*s2*t1*t2 - 4*Power(s2,2)*t1*t2 + 
               16*Power(t1,2)*t2 - 4*s1*Power(t1,2)*t2 + 
               8*s2*Power(t1,2)*t2 - 4*Power(t1,3)*t2 - 
               8*Power(t2,2) + 18*s1*Power(t2,2) + 6*s2*Power(t2,2) + 
               4*s1*s2*Power(t2,2) + 4*Power(s2,2)*Power(t2,2) + 
               2*t1*Power(t2,2) - 4*s1*t1*Power(t2,2) - 
               8*s2*t1*Power(t2,2) + 4*Power(t1,2)*Power(t2,2) - 
               6*Power(t2,3) + 
               2*Power(s,2)*(-11 + 6*Power(s1,2) + 2*Power(s2,2) + 
                  4*t1 - 2*Power(t1,2) + 
                  s1*(1 + 10*s2 - 6*t1 - 8*t2) - 2*t2 + 6*t1*t2 + 
                  2*Power(t2,2) - 2*s2*(1 + 4*t2)) + 
               2*Power(s1,2)*
                (-11 + 2*Power(s2,2) - 9*t2 + 4*t1*(1 + t2) - 
                  s2*(3 + 2*t1 + 4*t2)) - 
               ((2*Power(s,2) + s1 + s2 + 2*s1*s2 + 2*Power(s2,2) - 
                    2*s*(1 + s1 + 2*s2 - 2*t1) - t1 - 2*s1*t1 - 
                    4*s2*t1 + 2*Power(t1,2) - t2)*(s - s1 + t2)*
                  (-s2 + 2*Power(s1,2)*(s2 - t1) + t1 + s2*t1 - 
                    Power(t1,2) + Power(s,2)*(-2 + 2*s1 + t1 - t2) + 
                    2*t2 - s2*t2 - Power(s2,2)*t2 + t1*t2 + 
                    s2*t1*t2 - 2*Power(t2,2) + 
                    s1*(-2 + Power(s2,2) - Power(t1,2) + 2*t2 - 
                       2*s2*t2 + 2*t1*t2) + 
                    s*(3 - 2*Power(s1,2) - 3*t1 + Power(t1,2) - t2 - 
                       t1*t2 + s2*(2 - t1 + 2*t2) + 
                       s1*(-3*s2 + t1 + 2*t2))))/
                (1 + s*(-1 + s1) + s1*(-s2 + t1) - t2) - 
               2*s*(-17 + 2*Power(s1,3) + 17*t1 + 4*Power(t1,2) + 
                  2*Power(t1,3) + 
                  Power(s1,2)*(4 + 8*s2 - 6*t1 - 4*t2) + 
                  2*Power(s2,2)*(1 + t1 - 2*t2) + 7*t2 - 17*t1*t2 + 
                  9*Power(t2,2) - 4*t1*Power(t2,2) + 
                  s2*(-22 - 6*t1 - 4*Power(t1,2) + 11*t2 + 4*t1*t2 + 
                     4*Power(t2,2)) + 
                  s1*(-18 + 4*Power(s2,2) + 8*t1 - 2*Power(t1,2) - 
                     9*t2 + 10*t1*t2 + 2*Power(t2,2) - 
                     s2*(5 + 2*t1 + 12*t2)))))/
           ((-1 + s1)*(-s + s1 - t2)) + 
          (64*(-(((8 + 6*s2 - Power(s2,2) - Power(s2,3) + 
                      Power(s1,2)*(1 + s2 - t1) - 6*t1 + 2*s2*t1 + 
                      3*Power(s2,2)*t1 - Power(t1,2) - 
                      3*s2*Power(t1,2) + Power(t1,3) + 
                      Power(s,2)*(-1 + s1 - s2 + t1 - t2) + 4*t2 - 
                      s2*t2 - Power(s2,2)*t2 + t1*t2 + 2*s2*t1*t2 - 
                      Power(t1,2)*t2 - Power(t2,2) + 
                      s1*(-4 + t1 + t1*t2 - s2*(1 + t2)) + 
                      s*(-2 - Power(s1,2) + s2 + 2*Power(s2,2) - 
                       t1 - 4*s2*t1 + 2*Power(t1,2) - t2 + 
                       2*s2*t2 - 2*t1*t2 + s1*(1 - s2 + t1 + t2)))*
                    (-s2 + 2*Power(s1,2)*(s2 - t1) + t1 + s2*t1 - 
                      Power(t1,2) + 
                      Power(s,2)*(-2 + 2*s1 + t1 - t2) + 2*t2 - 
                      s2*t2 - Power(s2,2)*t2 + t1*t2 + s2*t1*t2 - 
                      2*Power(t2,2) + 
                      s1*(-2 + Power(s2,2) - Power(t1,2) + 2*t2 - 
                       2*s2*t2 + 2*t1*t2) + 
                      s*(3 - 2*Power(s1,2) - 3*t1 + Power(t1,2) - 
                        t2 - t1*t2 + s2*(2 - t1 + 2*t2) + 
                        s1*(-3*s2 + t1 + 2*t2))))/
                  (1 + s*(-1 + s1) + s1*(-s2 + t1) - t2)) + 
               2*(-21 - 25*s2 + 4*Power(s2,2) + 4*Power(s2,3) + 
                  28*t1 - 6*s2*t1 - 13*Power(s2,2)*t1 + 
                  2*Power(t1,2) + 14*s2*Power(t1,2) - 5*Power(t1,3) + 
                  Power(s1,3)*(-1 - s2 + t1) - 5*t2 + 19*s2*t2 + 
                  7*Power(s2,2)*t2 - Power(s2,3)*t2 - 24*t1*t2 - 
                  18*s2*t1*t2 + 3*Power(s2,2)*t1*t2 + 
                  11*Power(t1,2)*t2 - 3*s2*Power(t1,2)*t2 + 
                  Power(t1,3)*t2 + 9*Power(t2,2) + 3*s2*Power(t2,2) - 
                  Power(s2,2)*Power(t2,2) - 2*t1*Power(t2,2) + 
                  2*s2*t1*Power(t2,2) - Power(t1,2)*Power(t2,2) - 
                  Power(t2,3) + Power(s,3)*(s1 - 2*(s2 - t1 + t2)) + 
                  Power(s1,2)*
                   (1 + Power(t1,2) + t2 - t1*(1 + 2*t2) + 
                     s2*(2 - t1 + 2*t2)) + 
                  Power(s,2)*
                   (-8 - 2*Power(s1,2) + 4*Power(s2,2) + t1 + 
                     4*Power(t1,2) + 2*t2 - 3*t1*t2 - Power(t2,2) + 
                     s2*(3 - 8*t1 + 3*t2) + s1*(s2 - 2*t1 + 4*t2)) + 
                  s1*(-13 + Power(s2,3) - 2*Power(t1,3) + 
                     Power(t1,2)*(-2 + t2) - 7*t2 + Power(t2,2) + 
                     Power(s2,2)*(1 - 4*t1 + t2) + 
                     s2*(-18 + t1 + 5*Power(t1,2) - t2 - 2*t1*t2 - 
                        Power(t2,2)) + t1*(25 - 2*t2 + Power(t2,2))) + 
                  s*(21 + Power(s1,3) - 2*Power(s2,3) - 20*t1 - 
                     5*Power(t1,2) + 2*Power(t1,3) + 
                     Power(s2,2)*(-7 + 6*t1) + 
                     Power(s1,2)*(1 + 2*s2 - t1 - 2*t2) - t2 + 
                     16*t1*t2 - 5*Power(t2,2) - 2*t1*Power(t2,2) + 
                     s2*(13 + 12*t1 - 6*Power(t1,2) - 10*t2 + 
                        2*Power(t2,2)) - 
                     s1*(-11 + 3*Power(s2,2) + 2*t1 + 5*Power(t1,2) + 
                        t2 - 5*t1*t2 - Power(t2,2) + 
                        s2*(4 - 8*t1 + 5*t2))))))/
           ((-1 + s2)*(-s + s2 - t1)) + 
          ((-64*(2*Power(s,3) - 2*Power(s1,2)*(1 + s2 - t1) + 
                  t2*(-3 + 2*Power(s2,2) + s2*(2 - 4*t1) - 2*t1 + 
                     2*Power(t1,2) + t2) + 
                  Power(s,2)*(-1 - 4*s1 - 4*s2 + 4*t1 + 2*t2) + 
                  s*(-3 + 2*Power(s1,2) + 2*Power(s2,2) - 2*t1 + 
                     2*Power(t1,2) + s2*(2 - 4*t1 - 4*t2) + 
                     s1*(3 + 6*s2 - 6*t1 - 2*t2) + 4*t1*t2) + 
                  s1*(3 - 2*Power(s2,2) - 2*Power(t1,2) - 
                     2*t1*(-1 + t2) + t2 + 2*s2*(-1 + 2*t1 + t2)))*
                (-s2 + 2*Power(s1,2)*(s2 - t1) + t1 + s2*t1 - 
                  Power(t1,2) + Power(s,2)*(-2 + 2*s1 + t1 - t2) + 
                  2*t2 - s2*t2 - Power(s2,2)*t2 + t1*t2 + s2*t1*t2 - 
                  2*Power(t2,2) + 
                  s1*(-2 + Power(s2,2) - Power(t1,2) + 2*t2 - 
                     2*s2*t2 + 2*t1*t2) + 
                  s*(3 - 2*Power(s1,2) - 3*t1 + Power(t1,2) - t2 - 
                     t1*t2 + s2*(2 - t1 + 2*t2) + 
                     s1*(-3*s2 + t1 + 2*t2))))/
              (1 + s*(-1 + s1) + s1*(-s2 + t1) - t2) + 
             128*(-20 + 2*Power(s,4) - 24*s2 + 4*Power(s2,2) + 
                4*Power(s2,3) + 2*Power(s1,3)*(1 + s2 - t1) + 24*t1 - 
                8*s2*t1 - 12*Power(s2,2)*t1 + 4*Power(t1,2) + 
                12*s2*Power(t1,2) - 4*Power(t1,3) - 8*t2 + 16*s2*t2 + 
                6*Power(s2,2)*t2 - 10*t1*t2 - 12*s2*t1*t2 - 
                2*Power(s2,2)*t1*t2 + 6*Power(t1,2)*t2 + 
                4*s2*Power(t1,2)*t2 - 2*Power(t1,3)*t2 - 
                5*Power(t2,2) + s2*Power(t2,2) + 
                2*Power(s2,2)*Power(t2,2) - t1*Power(t2,2) - 
                4*s2*t1*Power(t2,2) + 2*Power(t1,2)*Power(t2,2) + 
                Power(t2,3) + 
                Power(s,3)*(1 - 6*s1 - 4*s2 + 2*t1 + 4*t2) + 
                Power(s,2)*(-15 + 6*Power(s1,2) + s2 + 
                   2*Power(s2,2) - t1 - 2*Power(t1,2) + 
                   2*s1*(1 + 5*s2 - 3*t1 - 4*t2) + 5*t2 - 8*s2*t2 + 
                   6*t1*t2 + 2*Power(t2,2)) + 
                Power(s1,2)*(-11 + 2*Power(s2,2) + t1 - 3*t2 + 
                   4*t1*t2 - s2*(1 + 2*t1 + 4*t2)) + 
                2*s1*(-4 + Power(t1,3) + 
                   Power(s2,2)*(1 + t1 - 2*t2) - 
                   Power(t1,2)*(-1 + t2) + 10*t2 - 
                   t1*(-3 + 2*t2 + Power(t2,2)) + 
                   s2*(-6 - 2*Power(t1,2) + 2*t2 + Power(t2,2) + 
                      t1*(-2 + 3*t2))) + 
                s*(16 - 2*Power(s1,3) - 2*t1 - 6*Power(t1,2) - 
                   2*Power(t1,3) - 2*Power(s2,2)*(3 + t1 - 2*t2) - 
                   24*t2 + 10*t1*t2 + Power(t2,2) + 4*t1*Power(t2,2) + 
                   Power(s1,2)*(-5 - 8*s2 + 6*t1 + 4*t2) + 
                   2*s2*(4 + 6*t1 + 2*Power(t1,2) - 5*t2 - 2*t1*t2 - 
                      2*Power(t2,2)) - 
                   2*s1*(-15 + 2*Power(s2,2) - Power(t1,2) + 5*t1*t2 + 
                      Power(t2,2) - s2*(t1 + 6*t2)))))/
           ((s - s2 + t1)*(s - s1 + t2)) + 
          (64*(((-3 + Power(s,3) - Power(s1,2) - 2*s2 + 
                    Power(s2,2) + 2*t1 - 2*s2*t1 + Power(t1,2) + 
                    s1*(2 + s2*(-3 + t2) - t1*(-3 + t2) - 2*t2) - 
                    4*t2 + Power(s2,2)*t2 - 2*s2*t1*t2 + 
                    Power(t1,2)*t2 + 
                    Power(s,2)*(1 - s1 - 2*s2 + 2*t1 + t2) + 
                    s*(-2 + Power(s2,2) + 2*t1 + Power(t1,2) + 
                       s1*(1 + s2 - t1 - t2) + 2*t1*t2 - 
                       2*s2*(1 + t1 + t2)))*
                  (-s2 + 2*Power(s1,2)*(s2 - t1) + t1 + s2*t1 - 
                    Power(t1,2) + Power(s,2)*(-2 + 2*s1 + t1 - t2) + 
                    2*t2 - s2*t2 - Power(s2,2)*t2 + t1*t2 + 
                    s2*t1*t2 - 2*Power(t2,2) + 
                    s1*(-2 + Power(s2,2) - Power(t1,2) + 2*t2 - 
                       2*s2*t2 + 2*t1*t2) + 
                    s*(3 - 2*Power(s1,2) - 3*t1 + Power(t1,2) - t2 - 
                       t1*t2 + s2*(2 - t1 + 2*t2) + 
                       s1*(-3*s2 + t1 + 2*t2))))/
                (1 + s*(-1 + s1) + s1*(-s2 + t1) - t2) - 
               2*(13 + Power(s,4) + Power(s1,3) + 44*s2 + 
                  17*Power(s2,2) - Power(s2,3) - 45*t1 - 40*s2*t1 - 
                  Power(s2,2)*t1 + 23*Power(t1,2) + 5*s2*Power(t1,2) - 
                  3*Power(t1,3) + 19*t2 + 5*s2*t2 - 2*Power(s2,2)*t2 - 
                  3*t1*t2 + 3*s2*t1*t2 - Power(s2,2)*t1*t2 - 
                  Power(t1,2)*t2 + 2*s2*Power(t1,2)*t2 - 
                  Power(t1,3)*t2 - 8*Power(t2,2) - s2*Power(t2,2) + 
                  Power(s2,2)*Power(t2,2) + 3*t1*Power(t2,2) - 
                  2*s2*t1*Power(t2,2) + Power(t1,2)*Power(t2,2) - 
                  2*Power(t2,3) + 
                  Power(s,3)*(1 - 2*s1 - 3*s2 + 2*t1 + t2) - 
                  Power(s1,2)*
                   (-7 + s2*(-4 + t2) - t1*(-3 + t2) + t2) + 
                  Power(s,2)*
                   (3 + Power(s1,2) + 3*Power(s2,2) + 3*t1 + 
                     Power(t1,2) + s1*(4*s2 - 3*t1 - 2*t2) - 4*t2 + 
                     t1*t2 + Power(t2,2) - 2*s2*(2 + 2*t1 + t2)) + 
                  s1*(-7 - 2*Power(t1,2) - Power(s2,2)*t2 + 
                     2*Power(t2,2) - t1*(4 - 6*t2 + Power(t2,2)) + 
                     s2*(6 - 4*t2 + Power(t2,2) + t1*(2 + t2))) - 
                  s*(18 + Power(s2,3) - 22*t1 + 3*Power(t1,2) + 
                     Power(s1,2)*(2 + s2 - t1 - t2) + 7*t2 - t1*t2 + 
                     Power(t1,2)*t2 + 3*Power(t2,2) - 
                     2*t1*Power(t2,2) - Power(s2,2)*(4 + 2*t1 + t2) + 
                     s1*(8 + 2*Power(s2,2) + Power(t1,2) - 8*t2 + 
                        Power(t2,2) + 2*t1*(1 + t2) - 
                        s2*(2 + 3*t1 + 3*t2)) + 
                     s2*(t1 + Power(t1,2) + 2*(7 - 2*t2 + Power(t2,2)))\
))))/((-1 + s1)*(-1 + t2)) + (32*
             ((-2*(-s2 + 2*Power(s1,2)*(s2 - t1) + t1 + s2*t1 - 
                    Power(t1,2) + 
                    Power(s,2)*(-2 + 2*s1 + t1 - t2) + 2*t2 - 
                    s2*t2 - Power(s2,2)*t2 + t1*t2 + s2*t1*t2 - 
                    2*Power(t2,2) + 
                    s1*(-2 + Power(s2,2) - Power(t1,2) + 2*t2 - 
                       2*s2*t2 + 2*t1*t2) + 
                    s*(3 - 2*Power(s1,2) - 3*t1 + Power(t1,2) - 
                       t2 - t1*t2 + s2*(2 - t1 + 2*t2) + 
                       s1*(-3*s2 + t1 + 2*t2)))*
                  (26 + 2*Power(s1,3) + 15*s2 - Power(s2,2) - 
                    Power(s2,3) - 17*t1 + 2*s2*t1 + 
                    2*Power(s2,2)*t1 - Power(t1,2) - 
                    s2*Power(t1,2) + 
                    Power(s1,2)*(-6 + s2 + t1 - 4*t2) - 4*t2 - 
                    2*s2*t2 + Power(s2,2)*t2 + 7*t1*t2 - 
                    3*s2*t1*t2 + 2*Power(t1,2)*t2 - 3*Power(t2,2) + 
                    2*s2*Power(t2,2) - 2*t1*Power(t2,2) + 
                    s1*(9 - 2*Power(s2,2) - 5*t1 - 3*Power(t1,2) + 
                       s2*(-3 + 5*t1 - 3*t2) + 6*t2 + 2*t1*t2 + 
                       2*Power(t2,2)) + 
                    Power(s,2)*(2*s1 - 3*(s2 - t1 + t2)) + 
                    s*(-21 - 4*Power(s1,2) + 4*Power(s2,2) + 5*t1 - 
                       7*s2*t1 + 3*Power(t1,2) - 2*t2 + 2*s2*t2 - 
                       t1*t2 - 2*Power(t2,2) + 
                       s1*(6 + 2*s2 - 4*t1 + 7*t2))))/
                (1 + s*(-1 + s1) + s1*(-s2 + t1) - t2) + 
               4*(-12 - Power(s1,4) - 39*s2 - 13*Power(s2,2) + 
                  Power(s2,3) + 38*t1 + 38*s2*t1 + 3*Power(s2,2)*t1 - 
                  24*Power(t1,2) - 10*s2*Power(t1,2) + 
                  6*Power(t1,3) - Power(s1,3)*(-2 + s2 + t1 - 3*t2) - 
                  10*s2*t2 - 2*Power(s2,2)*t2 - Power(s2,3)*t2 + 
                  16*t1*t2 + 8*s2*t1*t2 + 2*Power(s2,2)*t1*t2 - 
                  10*Power(t1,2)*t2 - s2*Power(t1,2)*t2 - 
                  9*Power(t2,2) - 2*s2*Power(t2,2) + 
                  8*t1*Power(t2,2) - s2*t1*Power(t2,2) + 
                  Power(t1,2)*Power(t2,2) - Power(t2,3) + 
                  s2*Power(t2,3) - t1*Power(t2,3) + 
                  Power(s,3)*(s1 - 2*(s2 - t1 + t2)) + 
                  Power(s1,2)*
                   (-6 + Power(s2,2) + 3*t1 + 2*Power(t1,2) - 2*t2 - 
                     3*Power(t2,2) + s2*(5 - 4*t1 + 3*t2)) + 
                  Power(s,2)*
                   (-14 - 3*Power(s1,2) + 4*Power(s2,2) + 5*t1 + 
                     2*Power(t1,2) - t2 + t1*t2 - 3*Power(t2,2) + 
                     s2*(1 - 6*t1 + t2) + s1*(2 + 3*s2 - 5*t1 + 7*t2)\
) + s1*(-20 + Power(s2,3) + Power(t1,2)*(1 - 3*t2) + 7*t2 + 
                     Power(t2,2) + Power(t2,3) - 
                     Power(s2,2)*(-4 + 3*t1 + t2) + 
                     t1*(10 - 7*t2 + 2*Power(t2,2)) + 
                     s2*(-12 + 2*Power(t1,2) - 2*t2 - 
                        3*Power(t2,2) + t1*(-1 + 4*t2))) + 
                  s*(40 + 3*Power(s1,3) - 2*Power(s2,3) - 32*t1 + 
                     5*Power(t1,2) + 4*Power(s1,2)*(-1 + t1 - 2*t2) + 
                     t2 + 3*Power(t1,2)*t2 + 2*Power(t2,2) - 
                     2*t1*Power(t2,2) - Power(t2,3) + 
                     2*Power(s2,2)*(-1 + 2*t1 + t2) - 
                     s1*(-16 + 5*Power(s2,2) + 4*Power(t1,2) - 
                        3*t2 - 6*Power(t2,2) + t1*(8 + t2) + 
                        s2*(6 - 10*t1 + 4*t2)) + 
                     s2*(26 - 2*Power(t1,2) + 4*t2 + 3*Power(t2,2) - 
                        t1*(6 + 5*t2)))) + 
               (2*(-6 - Power(s1,2) + Power(s2,2) + 2*t1 + s*t1 + 
                    Power(t1,2) - s2*(2 + s + 2*t1 - t2) - 2*t2 - 
                    s*t2 - t1*t2 + s1*(4 + s + t2))*
                  (Power(s2,2) + Power(s2,3) - 2*s2*t1 - 
                    5*Power(s2,2)*t1 - 2*Power(s2,3)*t1 + 
                    Power(t1,2) + 7*s2*Power(t1,2) + 
                    7*Power(s2,2)*Power(t1,2) + 
                    Power(s2,3)*Power(t1,2) - 3*Power(t1,3) - 
                    8*s2*Power(t1,3) - 3*Power(s2,2)*Power(t1,3) + 
                    3*Power(t1,4) + 3*s2*Power(t1,4) - Power(t1,5) + 
                    4*s2*t2 - 5*Power(s2,2)*t2 + 2*Power(s2,3)*t2 + 
                    2*Power(s2,4)*t2 - 4*t1*t2 + 4*s2*t1*t2 - 
                    8*Power(s2,3)*t1*t2 - 2*Power(s2,4)*t1*t2 + 
                    Power(t1,2)*t2 - 4*s2*Power(t1,2)*t2 + 
                    11*Power(s2,2)*Power(t1,2)*t2 + 
                    6*Power(s2,3)*Power(t1,2)*t2 + 2*Power(t1,3)*t2 - 
                    6*s2*Power(t1,3)*t2 - 
                    6*Power(s2,2)*Power(t1,3)*t2 + Power(t1,4)*t2 + 
                    2*s2*Power(t1,4)*t2 - 2*Power(t2,2) - 
                    4*s2*Power(t2,2) + 7*Power(s2,2)*Power(t2,2) - 
                    3*Power(s2,3)*Power(t2,2) + 
                    Power(s2,4)*Power(t2,2) + 
                    Power(s2,5)*Power(t2,2) + 4*t1*Power(t2,2) - 
                    Power(s2,2)*t1*Power(t2,2) - 
                    4*Power(s2,3)*t1*Power(t2,2) - 
                    3*Power(s2,4)*t1*Power(t2,2) - 
                    7*Power(t1,2)*Power(t2,2) + 
                    7*s2*Power(t1,2)*Power(t2,2) + 
                    5*Power(s2,2)*Power(t1,2)*Power(t2,2) + 
                    3*Power(s2,3)*Power(t1,2)*Power(t2,2) - 
                    3*Power(t1,3)*Power(t2,2) - 
                    2*s2*Power(t1,3)*Power(t2,2) - 
                    Power(s2,2)*Power(t1,3)*Power(t2,2) + 
                    6*Power(t2,3) - 4*s2*Power(t2,3) - 
                    5*Power(s2,2)*Power(t2,3) + 
                    4*Power(s2,3)*Power(t2,3) + 
                    Power(s2,4)*Power(t2,3) + 4*t1*Power(t2,3) - 
                    6*Power(s2,2)*t1*Power(t2,3) - 
                    2*Power(s2,3)*t1*Power(t2,3) + 
                    5*Power(t1,2)*Power(t2,3) + 
                    2*s2*Power(t1,2)*Power(t2,3) + 
                    Power(s2,2)*Power(t1,2)*Power(t2,3) - 
                    6*Power(t2,4) + 4*s2*Power(t2,4) + 
                    2*Power(s2,2)*Power(t2,4) - 4*t1*Power(t2,4) - 
                    2*s2*t1*Power(t2,4) + 2*Power(t2,5) + 
                    2*Power(s1,4)*Power(s2 - t1,2)*
                     (-1 + s2 - t1 + t2) + 
                    Power(s,4)*(-1 + s2 - t1 + t2)*
                     (2 + 2*Power(s1,2) + Power(t1,2) + 
                       2*s1*(-2 + t1 - t2) + 2*t2 + Power(t2,2) - 
                       2*t1*(1 + t2)) + 
                    2*Power(s1,3)*(s2 - t1)*
                     (Power(s2,3) + Power(t1,3) + 
                       2*Power(-1 + t2,2) + 2*t1*Power(-1 + t2,2) - 
                       3*Power(t1,2)*t2 - 
                       Power(s2,2)*(2 + t1 + t2) - 
                       s2*(Power(t1,2) + 2*Power(-1 + t2,2) - 
                        2*t1*(1 + 2*t2))) + 
                    Power(s1,2)*
                     (Power(s2,5) - Power(t1,5) + 
                       Power(s2,4)*(1 - 3*t1 - 3*t2) + 
                       2*Power(-1 + t2,3) + 
                       Power(t1,4)*(-1 + 3*t2) + 
                       2*t1*Power(-1 + t2,2)*(-1 + 4*t2) + 
                       Power(t1,3)*(-6 + 6*t2 - 4*Power(t2,2)) + 
                       2*Power(t1,2)*
                       (2 + 4*t2 - 7*Power(t2,2) + Power(t2,3)) + 
                       2*Power(s2,3)*
                       (-2 + 2*Power(t1,2) + 5*t2 - Power(t2,2) + 
                       t1*(-2 + 4*t2)) - 
                       2*Power(s2,2)*
                       (-4 + 2*Power(t1,3) + 
                       2*Power(t1,2)*(-1 + t2) + 5*Power(t2,2) - 
                       Power(t2,3) + t1*(-1 + 7*t2)) + 
                       s2*(3*Power(t1,4) - 4*Power(t1,3)*t2 - 
                        2*Power(-1 + t2,2)*(-1 + 4*t2) + 
                        Power(t1,2)*(8 - 2*t2 + 6*Power(t2,2)) - 
                        4*t1*(3 + 2*t2 - 6*Power(t2,2) + Power(t2,3))\
)) - 2*s1*(-2*Power(t1,4)*(-1 + t2) + Power(s2,5)*t2 + 
                       2*Power(-1 + t2,3)*t2 + 
                       t1*Power(-1 + t2,2)*
                       (-2 - 3*t2 + 2*Power(t2,2)) + 
                       Power(t1,3)*(-2 - 3*t2 + 4*Power(t2,2)) - 
                       2*Power(t1,2)*
                       (-1 + t2 - 2*Power(t2,2) + 2*Power(t2,3)) + 
                       Power(s2,4)*(1 + t2 - t1*(1 + 3*t2)) + 
                       Power(s2,3)*
                       (1 - 4*t2 + 5*Power(t2,2) - Power(t2,3) + 
                       3*Power(t1,2)*(1 + t2) + 
                       t1*(-4 - 3*t2 + Power(t2,2))) - 
                       Power(s2,2)*
                       (1 - 7*t2 + 5*Power(t2,2) + Power(t2,3) + 
                       Power(t1,3)*(3 + t2) + 
                       Power(t1,2)*(-7 - t2 + 2*Power(t2,2)) + 
                       t1*(3 - 3*t2 + 5*Power(t2,2) - 2*Power(t2,3))\
) + s2*(Power(t1,4) + Power(t1,3)*(-6 + 3*t2 + Power(t2,2)) - 
                        Power(-1 + t2,2)*
                       (-2 - 3*t2 + 2*Power(t2,2)) - 
                        Power(t1,2)*
                       (-4 - 4*t2 + 4*Power(t2,2) + Power(t2,3)) + 
                        t1*(-1 - 5*t2 + Power(t2,2) + 5*Power(t2,3)))\
) - 2*Power(s,3)*(-3 + 2*t1 + 2*Power(t1,2) - 3*Power(t1,3) + 
                       Power(t1,4) + 3*t2 - 2*t1*t2 + 
                       7*Power(t1,2)*t2 - 3*Power(t1,3)*t2 - 
                       5*t1*Power(t2,2) + 
                       3*Power(t1,2)*Power(t2,2) + Power(t2,3) - 
                       t1*Power(t2,3) + 
                       2*Power(s1,3)*(-1 + s2 - t1 + t2) + 
                       Power(s1,2)*(-2 + 3*s2 - 3*t2)*
                       (-1 + s2 - t1 + t2) + 
                       Power(s2,2)*
                       (Power(t1,2) + 2*Power(1 + t2,2) - 
                       t1*(2 + 3*t2)) + 
                       s2*(1 - 2*Power(t1,3) - 2*t2 + 
                       3*Power(t2,2) + 2*Power(t2,3) + 
                       Power(t1,2)*(5 + 6*t2) - 
                       t1*(5 + 8*t2 + 6*Power(t2,2))) + 
                       s1*(3 + Power(t1,3) + 
                        Power(s2,2)*(-5 + 2*t1 - 4*t2) - 6*t2 + 
                        Power(t2,2) + Power(t2,3) - 
                        Power(t1,2)*(5 + t2) - 
                        t1*(2 - 4*t2 + Power(t2,2)) + 
                        s2*(2 - 3*Power(t1,2) + 2*t2 - 
                        3*Power(t2,2) + t1*(8 + 6*t2)))) + 
                    Power(s,2)*
                     (-5 + 7*t1 + 3*Power(t1,2) - 11*Power(t1,3) + 
                       7*Power(t1,4) - Power(t1,5) + 13*t2 - 
                       6*t1*t2 + 15*Power(t1,2)*t2 - 
                       18*Power(t1,3)*t2 + 3*Power(t1,4)*t2 - 
                       15*Power(t2,2) - 9*t1*Power(t2,2) + 
                       17*Power(t1,2)*Power(t2,2) - 
                       3*Power(t1,3)*Power(t2,2) + 5*Power(t2,3) - 
                       8*t1*Power(t2,3) + Power(t1,2)*Power(t2,3) + 
                       2*Power(t2,4) + 
                       2*Power(s1,4)*(-1 + s2 - t1 + t2) + 
                       Power(s2,3)*
                       (2 + Power(t1,2) + 10*t2 + 6*Power(t2,2) - 
                       2*t1*(1 + 3*t2)) + 
                       2*Power(s1,3)*
                       (5*Power(s2,2) + 3*Power(t1,2) - 
                       t1*(-3 + t2) - 2*(-1 + t2)*t2 + 
                       s2*(-5 - 8*t1 + 3*t2)) + 
                       Power(s2,2)*
                       (6 - 3*Power(t1,3) + 6*t2 + 10*Power(t2,2) + 
                       6*Power(t2,3) + Power(t1,2)*(11 + 15*t2) - 
                       2*t1*(7 + 17*t2 + 9*Power(t2,2))) + 
                       Power(s1,2)*
                       (7*Power(s2,3) + 5*Power(t1,3) - 
                       3*Power(t1,2)*(3 + 5*t2) - 
                       Power(s2,2)*(15 + 7*t1 + 9*t2) + 
                       s2*(-2 + 22*t1 - 5*Power(t1,2) + 18*t2 + 
                       26*t1*t2 - 14*Power(t2,2)) + 
                       t1*(-8 - 2*t2 + 8*Power(t2,2)) + 
                       2*(5 - 8*t2 + 2*Power(t2,2) + Power(t2,3))) + 
                       s2*(-3 + 3*Power(t1,4) - 10*t2 + 
                       21*Power(t2,2) + 8*Power(t2,3) - 
                       4*Power(t1,3)*(4 + 3*t2) + 
                       Power(t1,2)*(23 + 42*t2 + 15*Power(t2,2)) - 
                       2*t1*
                       (5 + 7*t2 + 17*Power(t2,2) + 3*Power(t2,3))) \
+ 2*s1*(Power(t1,4) + Power(s2,3)*(-4 + t1 - 6*t2) - 
                        4*Power(t1,3)*(-1 + t2) - 
                        2*Power(-1 + t2,2)*(1 + 2*t2) + 
                        Power(t1,2)*(-4 - t2 + 5*Power(t2,2)) + 
                        t1*(-2 + 10*t2 + Power(t2,2) - 
                       2*Power(t2,3)) - 
                        Power(s2,2)*
                       (2 + Power(t1,2) - 5*t2 + 3*Power(t2,2) - 
                       t1*(10 + 9*t2)) - 
                        s2*(-8 + Power(t1,3) - 
                        Power(t1,2)*(-10 + t2) + 11*t2 + 
                        7*Power(t2,2) - 3*Power(t2,3) + 
                        3*t1*(-1 + Power(t2,2))))) - 
                    2*s*(t1 - 4*Power(t1,3) + 4*Power(t1,4) - 
                       Power(t1,5) + 4*t2 - 2*t1*t2 + 
                       Power(t1,2)*t2 - 6*Power(t1,3)*t2 + 
                       3*Power(t1,4)*t2 - 9*Power(t2,2) + 
                       4*t1*Power(t2,2) + 5*Power(t1,2)*Power(t2,2) - 
                       4*Power(t1,3)*Power(t2,2) + 6*Power(t2,3) - 
                       2*t1*Power(t2,3) + 3*Power(t1,2)*Power(t2,3) - 
                       Power(t2,4) - t1*Power(t2,4) + 
                       2*Power(s1,4)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                       Power(s2,4)*t2*(2 - t1 + 2*t2) + 
                       Power(s2,3)*
                        (1 + 6*t2 + 3*Power(t2,2) + 2*Power(t2,3) + 
                        Power(t1,2)*(1 + 3*t2) - 
                        t1*(2 + 11*t2 + 6*Power(t2,2))) + 
                       Power(s2,2)*
                        (1 - 5*t2 + 8*Power(t2,2) + 5*Power(t2,3) - 
                        3*Power(t1,3)*(1 + t2) + 
                        Power(t1,2)*(8 + 19*t2 + 6*Power(t2,2)) - 
                        t1*(6 + 15*t2 + 14*Power(t2,2) + 
                        3*Power(t2,3))) + 
                       s2*(-1 - t2 + 3*Power(t2,2) - 3*Power(t2,3) + 
                        2*Power(t2,4) + Power(t1,4)*(3 + t2) - 
                        Power(t1,3)*(10 + 13*t2 + 2*Power(t2,2)) + 
                        Power(t1,2)*
                       (9 + 15*t2 + 15*Power(t2,2) + Power(t2,3)) - 
                        t1*(1 - 5*t2 + 15*Power(t2,2) + 
                       7*Power(t2,3))) + 
                       Power(s1,3)*
                        (4*Power(s2,3) - Power(s2,2)*(5 + 8*t1) + 
                        2*Power(-1 + t2,2) - 
                        Power(t1,2)*(1 + 4*t2) + 
                        t1*(2 - 6*t2 + 4*Power(t2,2)) + 
                        s2*(-2 + 4*Power(t1,2) + 6*t2 - 
                        4*Power(t2,2) + t1*(6 + 4*t2))) + 
                       Power(s1,2)*
                        (2*Power(s2,4) - Power(t1,4) - 
                        4*Power(-1 + t2,2)*t2 + 
                        2*Power(t1,3)*(2 + t2) - 
                        Power(s2,3)*(3 + 3*t1 + 5*t2) + 
                        Power(t1,2)*(5 - 3*t2 + Power(t2,2)) + 
                        t1*(-7 + 4*t2 + 5*Power(t2,2) - 
                       2*Power(t2,3)) - 
                        Power(s2,2)*
                       (5 + Power(t1,2) - 13*t2 + 5*Power(t2,2) - 
                       2*t1*(4 + 7*t2)) + 
                        s2*(9 + 3*Power(t1,3) - 8*t2 - 
                        3*Power(t2,2) + 2*Power(t2,3) + 
                        2*t1*t2*(-5 + 2*t2) - Power(t1,2)*(9 + 11*t2)\
)) - s1*(Power(t1,5) - 3*Power(t1,4)*t2 + Power(s2,4)*(1 + 4*t2) - 
                        Power(-1 + t2,2)*(-4 + t2 + 2*Power(t2,2)) + 
                        Power(t1,3)*(-6 + 4*t2 + 3*Power(t2,2)) + 
                        t1*(-2 - 2*t2 + 4*Power(t2,2)) - 
                        Power(t1,2)*
                        (-3 - 7*t2 + 2*Power(t2,2) + Power(t2,3)) - 
                        Power(s2,3)*
                        (-3 + Power(t1,2) + t2 - Power(t2,2) + 
                        t1*(4 + 8*t2)) + 
                        Power(s2,2)*
                        (-4 + 3*Power(t1,3) + t2 + 13*Power(t2,2) - 
                        3*Power(t2,3) + Power(t1,2)*(5 + t2) + 
                        t1*(-8 + 3*Power(t2,2))) - 
                        s2*(1 + 3*Power(t1,4) + 
                        Power(t1,3)*(2 - 6*t2) - 11*t2 + 
                        13*Power(t2,2) - 3*Power(t2,3) + 
                        Power(t1,2)*(-11 + 3*t2 + 7*Power(t2,2)) - 
                        2*t1*
                        (1 - 5*t2 - 5*Power(t2,2) + 2*Power(t2,3))))))\
)/(Power(-1 + s - s*s1 + s1*(s2 - t1) + t2,2)*(-1 + s2 - t1 + t2))))/
           ((-1 + s2)*(-1 + t1)) + 
          (32*((2*(-s2 + 2*Power(s1,2)*(s2 - t1) + t1 + s2*t1 - 
                    Power(t1,2) + Power(s,2)*(-2 + 2*s1 + t1 - t2) + 
                    2*t2 - s2*t2 - Power(s2,2)*t2 + t1*t2 + 
                    s2*t1*t2 - 2*Power(t2,2) + 
                    s1*(-2 + Power(s2,2) - Power(t1,2) + 2*t2 - 
                       2*s2*t2 + 2*t1*t2) + 
                    s*(3 - 2*Power(s1,2) - 3*t1 + Power(t1,2) - t2 - 
                       t1*t2 + s2*(2 - t1 + 2*t2) + 
                       s1*(-3*s2 + t1 + 2*t2)))*
                  (-17 + 2*Power(s,3) - 16*s2 + 2*Power(s2,2) + 
                    14*t1 + s2*t1 - 3*Power(t1,2) + 5*t2 + 
                    2*Power(s1,2)*t2 - 4*s2*t2 + Power(s2,2)*t2 + 
                    2*t1*t2 - Power(t1,2)*t2 - 2*s2*Power(t2,2) + 
                    2*t1*Power(t2,2) + 
                    Power(s,2)*(4 - 4*s1 - 4*s2 + 2*t1 + 3*t2) - 
                    s1*(11 + t2 + 2*Power(t2,2) + t1*(2 + t2) - 
                       s2*(4 + 3*t2)) + 
                    s*(13 + 2*Power(s1,2) + 2*Power(s2,2) + 
                       s1*(-4 + 4*s2 - 2*t1 - 5*t2) + 3*t2 + 
                       2*t1*t2 + 2*Power(t2,2) - s2*(7 + 2*t1 + 4*t2))\
))/(1 + s*(-1 + s1) + s1*(-s2 + t1) - t2) - 
               4*(12 + Power(s,4) + 29*s2 + 12*Power(s2,2) - 
                  Power(s2,3) - 26*t1 - 19*s2*t1 - 2*Power(s2,2)*t1 + 
                  6*Power(t1,2) + 3*s2*Power(t1,2) - 5*t2 - 
                  Power(s1,3)*t2 - 16*s2*t2 + Power(s2,2)*t2 + 
                  20*t1*t2 + 6*s2*t1*t2 - Power(s2,2)*t1*t2 - 
                  7*Power(t1,2)*t2 + s2*Power(t1,2)*t2 - 
                  10*Power(t2,2) - 2*s2*Power(t2,2) + 
                  Power(s2,2)*Power(t2,2) + 5*t1*Power(t2,2) - 
                  Power(t1,2)*Power(t2,2) - Power(t2,3) - 
                  s2*Power(t2,3) + t1*Power(t2,3) + 
                  Power(s,3)*(2 - 3*s1 - 3*s2 + t1 + 2*t2) + 
                  Power(s1,2)*(5 + 2*Power(t2,2) - 2*s2*(1 + t2)) + 
                  Power(s,2)*
                   (9 + 3*Power(s1,2) + 3*Power(s2,2) - 2*t1 + 
                     s1*(-4 + 6*s2 - 2*t1 - 5*t2) + t2 + t1*t2 + 
                     2*Power(t2,2) - 2*s2*(3 + t1 + 2*t2)) - 
                  s*(25 + Power(s1,3) + Power(s2,3) - 21*t1 + 
                     5*Power(t1,2) + 
                     Power(s1,2)*(-2 + 3*s2 - t1 - 4*t2) - 4*t2 + 
                     3*t1*t2 + Power(t1,2)*t2 - 2*Power(t2,2) - 
                     t1*Power(t2,2) - Power(t2,3) - 
                     Power(s2,2)*(5 + t1 + 2*t2) + 
                     s2*(24 - 4*t1 + t2 + 3*Power(t2,2)) + 
                     s1*(18 + 3*Power(s2,2) + t1*(-2 + t2) + t2 + 
                        4*Power(t2,2) - 2*s2*(4 + t1 + 3*t2))) - 
                  s1*(-17 - 5*t2 - Power(t2,2) + Power(t2,3) + 
                     Power(s2,2)*(3 + t2) - Power(t1,2)*(5 + t2) + 
                     t1*(25 + t2 + Power(t2,2)) + 
                     s2*(t1*(4 + t2) - 3*(9 + t2 + Power(t2,2))))) - 
               (2*(4 + Power(s,2) + t1 + t2 - s1*t2 + t1*t2 - 
                    s2*(1 + t2) + s*(2 - s1 - s2 + t1 + t2))*
                  (Power(s2,2) + Power(s2,3) - 2*s2*t1 - 
                    5*Power(s2,2)*t1 - 2*Power(s2,3)*t1 + 
                    Power(t1,2) + 7*s2*Power(t1,2) + 
                    7*Power(s2,2)*Power(t1,2) + 
                    Power(s2,3)*Power(t1,2) - 3*Power(t1,3) - 
                    8*s2*Power(t1,3) - 3*Power(s2,2)*Power(t1,3) + 
                    3*Power(t1,4) + 3*s2*Power(t1,4) - Power(t1,5) + 
                    4*s2*t2 - 5*Power(s2,2)*t2 + 2*Power(s2,3)*t2 + 
                    2*Power(s2,4)*t2 - 4*t1*t2 + 4*s2*t1*t2 - 
                    8*Power(s2,3)*t1*t2 - 2*Power(s2,4)*t1*t2 + 
                    Power(t1,2)*t2 - 4*s2*Power(t1,2)*t2 + 
                    11*Power(s2,2)*Power(t1,2)*t2 + 
                    6*Power(s2,3)*Power(t1,2)*t2 + 2*Power(t1,3)*t2 - 
                    6*s2*Power(t1,3)*t2 - 
                    6*Power(s2,2)*Power(t1,3)*t2 + Power(t1,4)*t2 + 
                    2*s2*Power(t1,4)*t2 - 2*Power(t2,2) - 
                    4*s2*Power(t2,2) + 7*Power(s2,2)*Power(t2,2) - 
                    3*Power(s2,3)*Power(t2,2) + 
                    Power(s2,4)*Power(t2,2) + 
                    Power(s2,5)*Power(t2,2) + 4*t1*Power(t2,2) - 
                    Power(s2,2)*t1*Power(t2,2) - 
                    4*Power(s2,3)*t1*Power(t2,2) - 
                    3*Power(s2,4)*t1*Power(t2,2) - 
                    7*Power(t1,2)*Power(t2,2) + 
                    7*s2*Power(t1,2)*Power(t2,2) + 
                    5*Power(s2,2)*Power(t1,2)*Power(t2,2) + 
                    3*Power(s2,3)*Power(t1,2)*Power(t2,2) - 
                    3*Power(t1,3)*Power(t2,2) - 
                    2*s2*Power(t1,3)*Power(t2,2) - 
                    Power(s2,2)*Power(t1,3)*Power(t2,2) + 
                    6*Power(t2,3) - 4*s2*Power(t2,3) - 
                    5*Power(s2,2)*Power(t2,3) + 
                    4*Power(s2,3)*Power(t2,3) + 
                    Power(s2,4)*Power(t2,3) + 4*t1*Power(t2,3) - 
                    6*Power(s2,2)*t1*Power(t2,3) - 
                    2*Power(s2,3)*t1*Power(t2,3) + 
                    5*Power(t1,2)*Power(t2,3) + 
                    2*s2*Power(t1,2)*Power(t2,3) + 
                    Power(s2,2)*Power(t1,2)*Power(t2,3) - 
                    6*Power(t2,4) + 4*s2*Power(t2,4) + 
                    2*Power(s2,2)*Power(t2,4) - 4*t1*Power(t2,4) - 
                    2*s2*t1*Power(t2,4) + 2*Power(t2,5) + 
                    2*Power(s1,4)*Power(s2 - t1,2)*
                     (-1 + s2 - t1 + t2) + 
                    Power(s,4)*(-1 + s2 - t1 + t2)*
                     (2 + 2*Power(s1,2) + Power(t1,2) + 
                       2*s1*(-2 + t1 - t2) + 2*t2 + Power(t2,2) - 
                       2*t1*(1 + t2)) + 
                    2*Power(s1,3)*(s2 - t1)*
                     (Power(s2,3) + Power(t1,3) + 
                       2*Power(-1 + t2,2) + 2*t1*Power(-1 + t2,2) - 
                       3*Power(t1,2)*t2 - Power(s2,2)*(2 + t1 + t2) - 
                       s2*(Power(t1,2) + 2*Power(-1 + t2,2) - 
                        2*t1*(1 + 2*t2))) + 
                    Power(s1,2)*
                     (Power(s2,5) - Power(t1,5) + 
                       Power(s2,4)*(1 - 3*t1 - 3*t2) + 
                       2*Power(-1 + t2,3) + Power(t1,4)*(-1 + 3*t2) + 
                       2*t1*Power(-1 + t2,2)*(-1 + 4*t2) + 
                       Power(t1,3)*(-6 + 6*t2 - 4*Power(t2,2)) + 
                       2*Power(t1,2)*
                        (2 + 4*t2 - 7*Power(t2,2) + Power(t2,3)) + 
                       2*Power(s2,3)*
                        (-2 + 2*Power(t1,2) + 5*t2 - Power(t2,2) + 
                        t1*(-2 + 4*t2)) - 
                       2*Power(s2,2)*
                        (-4 + 2*Power(t1,3) + 
                        2*Power(t1,2)*(-1 + t2) + 5*Power(t2,2) - 
                        Power(t2,3) + t1*(-1 + 7*t2)) + 
                       s2*(3*Power(t1,4) - 4*Power(t1,3)*t2 - 
                        2*Power(-1 + t2,2)*(-1 + 4*t2) + 
                        Power(t1,2)*(8 - 2*t2 + 6*Power(t2,2)) - 
                        4*t1*(3 + 2*t2 - 6*Power(t2,2) + Power(t2,3)))\
) - 2*s1*(-2*Power(t1,4)*(-1 + t2) + Power(s2,5)*t2 + 
                       2*Power(-1 + t2,3)*t2 + 
                       t1*Power(-1 + t2,2)*
                        (-2 - 3*t2 + 2*Power(t2,2)) + 
                       Power(t1,3)*(-2 - 3*t2 + 4*Power(t2,2)) - 
                       2*Power(t1,2)*
                        (-1 + t2 - 2*Power(t2,2) + 2*Power(t2,3)) + 
                       Power(s2,4)*(1 + t2 - t1*(1 + 3*t2)) + 
                       Power(s2,3)*
                        (1 - 4*t2 + 5*Power(t2,2) - Power(t2,3) + 
                        3*Power(t1,2)*(1 + t2) + 
                        t1*(-4 - 3*t2 + Power(t2,2))) - 
                       Power(s2,2)*
                        (1 - 7*t2 + 5*Power(t2,2) + Power(t2,3) + 
                        Power(t1,3)*(3 + t2) + 
                        Power(t1,2)*(-7 - t2 + 2*Power(t2,2)) + 
                        t1*(3 - 3*t2 + 5*Power(t2,2) - 2*Power(t2,3))\
) + s2*(Power(t1,4) + Power(t1,3)*(-6 + 3*t2 + Power(t2,2)) - 
                        Power(-1 + t2,2)*
                       (-2 - 3*t2 + 2*Power(t2,2)) - 
                        Power(t1,2)*
                        (-4 - 4*t2 + 4*Power(t2,2) + Power(t2,3)) + 
                        t1*(-1 - 5*t2 + Power(t2,2) + 5*Power(t2,3)))) \
- 2*Power(s,3)*(-3 + 2*t1 + 2*Power(t1,2) - 3*Power(t1,3) + 
                       Power(t1,4) + 3*t2 - 2*t1*t2 + 
                       7*Power(t1,2)*t2 - 3*Power(t1,3)*t2 - 
                       5*t1*Power(t2,2) + 3*Power(t1,2)*Power(t2,2) + 
                       Power(t2,3) - t1*Power(t2,3) + 
                       2*Power(s1,3)*(-1 + s2 - t1 + t2) + 
                       Power(s1,2)*(-2 + 3*s2 - 3*t2)*
                        (-1 + s2 - t1 + t2) + 
                       Power(s2,2)*
                        (Power(t1,2) + 2*Power(1 + t2,2) - 
                        t1*(2 + 3*t2)) + 
                       s2*(1 - 2*Power(t1,3) - 2*t2 + 
                        3*Power(t2,2) + 2*Power(t2,3) + 
                        Power(t1,2)*(5 + 6*t2) - 
                        t1*(5 + 8*t2 + 6*Power(t2,2))) + 
                       s1*(3 + Power(t1,3) + 
                        Power(s2,2)*(-5 + 2*t1 - 4*t2) - 6*t2 + 
                        Power(t2,2) + Power(t2,3) - 
                        Power(t1,2)*(5 + t2) - 
                        t1*(2 - 4*t2 + Power(t2,2)) + 
                        s2*(2 - 3*Power(t1,2) + 2*t2 - 
                        3*Power(t2,2) + t1*(8 + 6*t2)))) + 
                    Power(s,2)*
                     (-5 + 7*t1 + 3*Power(t1,2) - 11*Power(t1,3) + 
                       7*Power(t1,4) - Power(t1,5) + 13*t2 - 
                       6*t1*t2 + 15*Power(t1,2)*t2 - 
                       18*Power(t1,3)*t2 + 3*Power(t1,4)*t2 - 
                       15*Power(t2,2) - 9*t1*Power(t2,2) + 
                       17*Power(t1,2)*Power(t2,2) - 
                       3*Power(t1,3)*Power(t2,2) + 5*Power(t2,3) - 
                       8*t1*Power(t2,3) + Power(t1,2)*Power(t2,3) + 
                       2*Power(t2,4) + 
                       2*Power(s1,4)*(-1 + s2 - t1 + t2) + 
                       Power(s2,3)*
                        (2 + Power(t1,2) + 10*t2 + 6*Power(t2,2) - 
                        2*t1*(1 + 3*t2)) + 
                       2*Power(s1,3)*
                        (5*Power(s2,2) + 3*Power(t1,2) - 
                        t1*(-3 + t2) - 2*(-1 + t2)*t2 + 
                        s2*(-5 - 8*t1 + 3*t2)) + 
                       Power(s2,2)*
                        (6 - 3*Power(t1,3) + 6*t2 + 10*Power(t2,2) + 
                        6*Power(t2,3) + Power(t1,2)*(11 + 15*t2) - 
                        2*t1*(7 + 17*t2 + 9*Power(t2,2))) + 
                       Power(s1,2)*
                        (7*Power(s2,3) + 5*Power(t1,3) - 
                        3*Power(t1,2)*(3 + 5*t2) - 
                        Power(s2,2)*(15 + 7*t1 + 9*t2) + 
                        s2*(-2 + 22*t1 - 5*Power(t1,2) + 18*t2 + 
                       26*t1*t2 - 14*Power(t2,2)) + 
                        t1*(-8 - 2*t2 + 8*Power(t2,2)) + 
                        2*(5 - 8*t2 + 2*Power(t2,2) + Power(t2,3))) + 
                       s2*(-3 + 3*Power(t1,4) - 10*t2 + 
                        21*Power(t2,2) + 8*Power(t2,3) - 
                        4*Power(t1,3)*(4 + 3*t2) + 
                        Power(t1,2)*(23 + 42*t2 + 15*Power(t2,2)) - 
                        2*t1*
                        (5 + 7*t2 + 17*Power(t2,2) + 3*Power(t2,3))) \
+ 2*s1*(Power(t1,4) + Power(s2,3)*(-4 + t1 - 6*t2) - 
                        4*Power(t1,3)*(-1 + t2) - 
                        2*Power(-1 + t2,2)*(1 + 2*t2) + 
                        Power(t1,2)*(-4 - t2 + 5*Power(t2,2)) + 
                        t1*(-2 + 10*t2 + Power(t2,2) - 
                       2*Power(t2,3)) - 
                        Power(s2,2)*
                        (2 + Power(t1,2) - 5*t2 + 3*Power(t2,2) - 
                        t1*(10 + 9*t2)) - 
                        s2*(-8 + Power(t1,3) - 
                        Power(t1,2)*(-10 + t2) + 11*t2 + 
                        7*Power(t2,2) - 3*Power(t2,3) + 
                        3*t1*(-1 + Power(t2,2))))) - 
                    2*s*(t1 - 4*Power(t1,3) + 4*Power(t1,4) - 
                       Power(t1,5) + 4*t2 - 2*t1*t2 + Power(t1,2)*t2 - 
                       6*Power(t1,3)*t2 + 3*Power(t1,4)*t2 - 
                       9*Power(t2,2) + 4*t1*Power(t2,2) + 
                       5*Power(t1,2)*Power(t2,2) - 
                       4*Power(t1,3)*Power(t2,2) + 6*Power(t2,3) - 
                       2*t1*Power(t2,3) + 3*Power(t1,2)*Power(t2,3) - 
                       Power(t2,4) - t1*Power(t2,4) + 
                       2*Power(s1,4)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                       Power(s2,4)*t2*(2 - t1 + 2*t2) + 
                       Power(s2,3)*
                        (1 + 6*t2 + 3*Power(t2,2) + 2*Power(t2,3) + 
                        Power(t1,2)*(1 + 3*t2) - 
                        t1*(2 + 11*t2 + 6*Power(t2,2))) + 
                       Power(s2,2)*
                        (1 - 5*t2 + 8*Power(t2,2) + 5*Power(t2,3) - 
                        3*Power(t1,3)*(1 + t2) + 
                        Power(t1,2)*(8 + 19*t2 + 6*Power(t2,2)) - 
                        t1*(6 + 15*t2 + 14*Power(t2,2) + 
                        3*Power(t2,3))) + 
                       s2*(-1 - t2 + 3*Power(t2,2) - 3*Power(t2,3) + 
                        2*Power(t2,4) + Power(t1,4)*(3 + t2) - 
                        Power(t1,3)*(10 + 13*t2 + 2*Power(t2,2)) + 
                        Power(t1,2)*
                        (9 + 15*t2 + 15*Power(t2,2) + Power(t2,3)) - 
                        t1*(1 - 5*t2 + 15*Power(t2,2) + 7*Power(t2,3))\
) + Power(s1,3)*(4*Power(s2,3) - Power(s2,2)*(5 + 8*t1) + 
                        2*Power(-1 + t2,2) - Power(t1,2)*(1 + 4*t2) + 
                        t1*(2 - 6*t2 + 4*Power(t2,2)) + 
                        s2*(-2 + 4*Power(t1,2) + 6*t2 - 
                        4*Power(t2,2) + t1*(6 + 4*t2))) + 
                       Power(s1,2)*
                        (2*Power(s2,4) - Power(t1,4) - 
                        4*Power(-1 + t2,2)*t2 + 
                        2*Power(t1,3)*(2 + t2) - 
                        Power(s2,3)*(3 + 3*t1 + 5*t2) + 
                        Power(t1,2)*(5 - 3*t2 + Power(t2,2)) + 
                        t1*(-7 + 4*t2 + 5*Power(t2,2) - 
                       2*Power(t2,3)) - 
                        Power(s2,2)*
                        (5 + Power(t1,2) - 13*t2 + 5*Power(t2,2) - 
                        2*t1*(4 + 7*t2)) + 
                        s2*(9 + 3*Power(t1,3) - 8*t2 - 
                        3*Power(t2,2) + 2*Power(t2,3) + 
                        2*t1*t2*(-5 + 2*t2) - Power(t1,2)*(9 + 11*t2))\
) - s1*(Power(t1,5) - 3*Power(t1,4)*t2 + Power(s2,4)*(1 + 4*t2) - 
                        Power(-1 + t2,2)*(-4 + t2 + 2*Power(t2,2)) + 
                        Power(t1,3)*(-6 + 4*t2 + 3*Power(t2,2)) + 
                        t1*(-2 - 2*t2 + 4*Power(t2,2)) - 
                        Power(t1,2)*
                        (-3 - 7*t2 + 2*Power(t2,2) + Power(t2,3)) - 
                        Power(s2,3)*
                        (-3 + Power(t1,2) + t2 - Power(t2,2) + 
                        t1*(4 + 8*t2)) + 
                        Power(s2,2)*
                        (-4 + 3*Power(t1,3) + t2 + 13*Power(t2,2) - 
                        3*Power(t2,3) + Power(t1,2)*(5 + t2) + 
                        t1*(-8 + 3*Power(t2,2))) - 
                        s2*(1 + 3*Power(t1,4) + 
                        Power(t1,3)*(2 - 6*t2) - 11*t2 + 
                        13*Power(t2,2) - 3*Power(t2,3) + 
                        Power(t1,2)*(-11 + 3*t2 + 7*Power(t2,2)) - 
                        2*t1*(1 - 5*t2 - 5*Power(t2,2) + 2*Power(t2,3))\
)))))/(Power(-1 + s - s*s1 + s1*(s2 - t1) + t2,2)*(-1 + s2 - t1 + t2))))/
           ((-1 + t1)*(-1 + t2)))/(s - s1 + t2))*
     B1(1 - s + s2 - t1,1 - s2 + t1 - t2,1 - s + s1 - t2))/
   (128.*Power(Pi,2)) + (((-64*s*
          (-1 - 2*Power(s2,2)*(-1 + t1) - s*t1 + Power(t1,2) - 
            s*Power(t1,2) - s1*(s2 - t1)*(1 + t1) + t2 - s*t2 - 
            3*t1*t2 + s*t1*t2 + 
            s2*(2*Power(t1,2) + t1*(-2 + 2*s - t2) + 3*t2))*
          ((-3 - 4*s2 - 5*Power(s2,2) + 3*Power(s2,3) + 
               Power(s1,2)*(3 - 2*s2*(-2 + t1) - 3*t1) + 
               2*Power(s,3)*(-1 + t1) + 5*t1 + 16*s2*t1 + 
               4*Power(s2,2)*t1 - 2*Power(s2,3)*t1 + 2*Power(t1,2) - 
               2*s2*Power(t1,2) + 
               Power(s,2)*(3 - 4*s1*(-1 + t1) - 6*s2*(-1 + t1) + 
                  4*t1 - t2) - t2 + s2*t2 - Power(s2,2)*t2 + 
               s2*t1*t2 + s1*
                (2 + Power(s2,2)*(7 - 4*t1) - 2*Power(t1,2) + 
                  s2*(1 + t1 - t2) + t1*(11 + t2)) + 
               s*(2*Power(s1,2)*(-1 + t1) + Power(s2,2)*(-7 + 6*t1) + 
                  t1*(-15 + 2*t1 - t2) + 
                  s1*(-6 - t1 + 2*s2*(-5 + 4*t1) + t2) + 
                  s2*(2 - 9*t1 + 2*t2)))/((-1 + s2)*(-s + s2 - t1)) + 
            (5*s2 + 5*Power(s2,2) + 4*Power(s2,3) - 4*s2*t1 - 
               2*Power(s2,3)*t1 - 2*Power(t1,2) + 2*s2*Power(t1,2) + 
               Power(s1,2)*(-4 + s2 + t1 - 2*s2*t1) + 11*s2*t2 + 
               2*Power(s2,2)*t2 - s2*t1*t2 + 
               Power(s,2)*(3 - 2*s2*(-2 + t1) - t1 + 4*t2) + 
               s*(-5 + 4*Power(s2,2)*(-2 + t1) + 8*t1 - 
                  2*Power(t1,2) + s2*(-10 + t1 - 6*t2) + 
                  s1*(1 + s2*(-5 + 4*t1) - 4*t2) - 11*t2 + t1*t2) + 
               s1*(5 + Power(s2,2)*(5 - 4*t1) + 2*Power(t1,2) + 
                  11*t2 - t1*(8 + t2) + s2*(3 + t1 + 2*t2)))/
             ((-1 + s2)*(-1 + t1)) + 
            (-2 + s2 + 12*Power(s2,2) + 6*Power(s2,3) - 
               2*Power(s2,4) + t1 - 6*s2*t1 + Power(s2,2)*t1 - 
               Power(t1,2) - 2*t2 + 3*s2*t2 + 8*Power(s2,2)*t2 - 
               2*Power(s2,3)*t2 - 3*t1*t2 + 
               Power(s,2)*(2 + s2 - 2*Power(s2,2) + 2*t1 + t2 - 
                  2*s2*t2) + 
               Power(s1,2)*(1 + s2 - 2*Power(s2,2) + 2*t1 + 2*t2 - 
                  2*s2*t2) + 
               s1*(3 - 4*Power(s2,3) + Power(t1,2) + 
                  Power(s2,2)*(7 - 4*t2) - t2 - t1*(4 + t2) + 
                  s2*(7 + 2*t1 + 10*t2)) + 
               s*(-2 + 4*Power(s2,3) + 3*t1 - Power(t1,2) + t2 + 
                  t1*t2 + Power(s2,2)*(-7 + 4*t2) - 
                  s2*(8 + 2*t1 + 9*t2) + 
                  s1*(-3 + 4*Power(s2,2) - 4*t1 - 3*t2 + 
                     s2*(-2 + 4*t2))))/((-1 + t1)*(-1 + t2)) + 
            (-6 + 2*Power(s1,3) - 16*s2 + 5*Power(s2,2) + 
               4*Power(s2,3) - 2*Power(s2,4) + 5*t1 + 14*s2*t1 + 
               5*Power(s2,2)*t1 + Power(t1,2) - 2*t2 + 8*s2*t2 - 
               3*Power(s2,2)*t2 - 2*Power(s2,3)*t2 + 3*t1*t2 + 
               2*Power(s,3)*(-3 + s2 + t2) - 
               Power(s1,2)*(-1 + 2*Power(s2,2) + 3*t1 + 
                  2*s2*(-5 + t2) + 2*t2) + 
               s1*(-4*Power(s2,3) - Power(t1,2) + 
                  Power(s2,2)*(13 - 4*t2) + 2*s2*(3 + t1 - 2*t2) + 
                  7*(-1 + t2) + t1*(5 + t2)) - 
               2*Power(s,2)*
                (2 + 3*Power(s2,2) - 2*t1 + 3*s2*(-3 + t2) + t2 + 
                  2*s1*(-4 + s2 + t2)) + 
               s*(14 + 6*Power(s2,3) - 11*t1 - 7*t2 - t1*t2 + 
                  2*Power(s1,2)*(-6 + s2 + t2) + s2*(-10*t1 + 4*t2) + 
                  Power(s2,2)*(-17 + 6*t2) + 
                  s1*(2 + 8*Power(s2,2) + 4*t2 + s2*(-30 + 8*t2))))/
             ((s - s2 + t1)*(s - s1 + t2)) + 
            (-4 + 4*Power(s1,3) - 11*s2 + 10*Power(s2,2) + 
               10*Power(s2,3) - 2*Power(s2,4) + 3*t1 - s2*t1 - 
               2*Power(s2,2)*t1 + Power(t1,2) + 2*s2*Power(t1,2) - 
               5*t2 + 15*s2*t2 + 6*Power(s2,2)*t2 - 
               2*Power(s2,3)*t2 - 2*s2*t1*t2 + Power(t2,2) + 
               s2*Power(t2,2) + 2*Power(s,3)*(-3 + s2 + t2) - 
               2*Power(s1,2)*(4 + Power(s2,2) + s2*(-5 + t2) + t2) + 
               s1*(-3 - 4*Power(s2,3) + Power(t1,2) - 
                  s2*(4 + t1 - 5*t2) - 4*Power(s2,2)*(-4 + t2) + 
                  9*t2 + 2*Power(t2,2) - t1*(3 + 2*t2)) - 
               2*Power(s,2)*
                (-1 + 3*Power(s2,2) + 3*s2*(-3 + t2) - t2 + 
                  s1*(-7 + 2*s2 + 2*t2)) + 
               s*(6 + 6*Power(s2,3) + t1 - 2*Power(t1,2) + 
                  s2*(-7 + 2*t1 - 9*t2) - 9*t2 + 2*t1*t2 - 
                  2*Power(t2,2) + 2*Power(s1,2)*(-6 + s2 + t2) + 
                  Power(s2,2)*(-22 + 6*t2) + 
                  s1*(7 + 8*Power(s2,2) - t1 + s2*(-26 + 8*t2))))/
             ((-1 + s1)*(-s + s1 - t2)) + 
            (-2 + 2*Power(s1,3) - 2*s2 + 8*Power(s2,2) + 
               9*Power(s2,3) - 2*Power(s2,4) + t1 - s2*t1 + 
               2*Power(s2,2)*t1 - Power(t1,2) - 2*s2*Power(t1,2) - 
               5*t2 + 2*s2*t2 + 9*Power(s2,2)*t2 - 2*Power(s2,3)*t2 + 
               2*s2*t1*t2 - Power(t2,2) - s2*Power(t2,2) + 
               2*Power(s,3)*(-1 + s2 + t2) + 
               Power(s,2)*(2 - 6*Power(s2,2) + s2*(8 - 6*t2) + 
                  s1*(6 - 4*s2 - 4*t2) + 4*t2) + 
               Power(s1,2)*(-5 - 2*Power(s2,2) + s2*(7 - 2*t2) + 
                  4*t2) - s1*
                (-2 + 4*Power(s2,3) + t1 + Power(t1,2) + 
                  s2*(3 - 3*t1 - 10*t2) + 4*Power(s2,2)*(-3 + t2) - 
                  2*t1*t2 + 2*Power(t2,2)) + 
               s*(6*Power(s2,3) + 2*Power(s1,2)*(-3 + s2 + t2) + 
                  Power(s2,2)*(-13 + 6*t2) - s2*(7 + 2*t1 + 10*t2) + 
                  2*(Power(t1,2) - t1*t2 + Power(t2,2)) + 
                  s1*(3 + 8*Power(s2,2) - 8*t2 + s2*(-15 + 8*t2))))/
             ((-1 + s1)*(-1 + t2))))/
        (Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
          Power(s,2)*(-1 + s2)*(s2*t1 - Power(t1,2) - t2 + t1*t2) - 
          s*(-1 + s2)*(Power(s2,2)*(-1 + t1) + 2*t1*(1 - t1 + t2) + 
             s1*(1 + s2 - 2*t1 + s2*t1 - Power(t1,2) - t2 + t1*t2) - 
             s2*(1 + Power(t1,2) + 3*t2 - t1*(2 + t2)))) + 
       64*((-2 - 5*t1 + 8*s1*t1 + 22*Power(t1,2) + 
             4*Power(s,3)*Power(t1,2) - 7*s1*Power(t1,2) - 
             19*Power(t1,3) + 5*s1*Power(t1,3) + 4*Power(t1,4) + 
             Power(s2,3)*(-2 + 6*t1 - 4*Power(t1,2)) + 2*t2 - 
             5*t1*t2 - Power(t1,2)*t2 - 2*Power(t1,3)*t2 + 
             Power(s2,2)*(10 - 2*Power(t1,2) + 4*Power(t1,3) + 
                s1*(-7 + 9*t1 - 4*Power(t1,2)) + t1*(-12 + t2) + t2) + 
             s2*(3 - 8*Power(t1,3) + 
                2*s1*(-4 + 7*t1 - 7*Power(t1,2) + 2*Power(t1,3)) - 
                Power(t1,2)*(-35 + t2) + 7*t2 - 2*t1*(15 + t2)) + 
             Power(s,2)*((3 - 4*s1 - 12*s2)*Power(t1,2) + 
                4*Power(t1,3) + t2 + t1*(1 + 2*s1 + 6*s2 + t2)) + 
             s*(1 + 2*t1 - 26*Power(t1,2) + 9*Power(t1,3) + 
                2*Power(s2,2)*(1 - 6*t1 + 6*Power(t1,2)) + 
                s1*(t1*(-9 + 9*t1 - 4*Power(t1,2)) + 
                   s2*(3 - 11*t1 + 8*Power(t1,2))) - 3*t2 + 2*t1*t2 + 
                Power(t1,2)*t2 - 
                s2*(-13*t1 + Power(t1,2) + 8*Power(t1,3) + 2*t2 + 
                   2*t1*t2)))/
           ((-1 + s2)*(-s + s2 - t1)*(s2*(-1 + t1) - t1*(-1 + s + t1))) \
- (2 + 4*Power(s2,3)*Power(-1 + t1,2) + 4*Power(t1,4) + 
             Power(t1,3)*(-23 + 11*s + s1 - 2*t2) - 2*t2 + 2*s*t2 + 
             Power(s2,2)*(6 - 8*s*Power(-1 + t1,2) - 14*t1 + 
                12*Power(t1,2) - 4*Power(t1,3) + 
                s1*(3 - 5*t1 + 4*Power(t1,2)) + 3*t2 - 5*t1*t2) + 
             Power(t1,2)*(24 + 10*Power(s,2) - 5*s1 + 13*t2 - 
                s*(33 + 7*t2)) + 
             t1*(-7 + 4*Power(s,3) - 2*s1 - 3*t2 - 
                2*Power(s,2)*(8 + s1 + 2*t2) + s*(2*s1 + 11*(2 + t2))) \
+ s2*(9 - 32*t1 + 31*Power(t1,2) - 8*Power(t1,3) + 
                4*Power(s,2)*(1 - 3*t1 + Power(t1,2)) + 
                s1*(2 + 2*t1 + 4*Power(t1,2) - 4*Power(t1,3)) + t2 - 
                14*t1*t2 + 9*Power(t1,2)*t2 + 
                s*(-21*Power(t1,2) + 4*Power(t1,3) + 
                   s1*(-2 + 4*t1 - 4*Power(t1,2)) - 3*(4 + t2) + 
                   t1*(33 + 9*t2))))/
           ((-1 + s2)*(-1 + t1)*(s2*(-1 + t1) - t1*(-1 + s + t1))) + 
          (4*Power(s,6)*Power(t1,2) + 
             4*(s2 - t1)*(-1 + t1)*
              Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
             2*Power(s,5)*t1*
              (-3*Power(t1,2) + 2*s2*(-2 + 5*t1) + t2 + 
                t1*(-1 + 4*s1 + t2)) + 
             2*Power(s,4)*(2*Power(s2,2)*
                 (1 - 8*t1 + 10*Power(t1,2)) - 
                s2*(13*Power(t1,3) + 
                   Power(t1,2)*(1 - 15*s1 - 4*t2) + 
                   t1*(-7 + 9*s1 - 3*t2) + t2) + 
                t1*(-1 + Power(t1,3) + 
                   Power(t1,2)*(2 - 6*s1 - 3*t2) + s1*t2 + 
                   t1*(-9 + 2*Power(s1,2) - t2 + s1*t2))) - 
             Power(s,3)*(4*Power(s2,3)*(3 - 12*t1 + 10*Power(t1,2)) + 
                Power(t1,2)*(5 - 38*s1 + 6*Power(s1,2) - 12*t2) - 
                2*Power(t2,2) + Power(t1,4)*(7 + 4*s1 + 4*t2) + 
                Power(t1,3)*
                 (-18 - 6*Power(s1,2) + s1*(20 - 6*t2) + 7*t2) + 
                t1*(2 + 2*s1*(-1 + t2) - 7*t2 + 2*Power(t2,2)) + 
                2*Power(s2,2)*
                 (s1*(5 - 26*t1 + 21*Power(t1,2)) - 
                   (1 + t1)*
                    (6 - 26*t1 + 21*Power(t1,2) + 3*t2 - 6*t1*t2)) + 
                2*s2*(1 + 4*Power(t1,4) + 
                   Power(t1,3)*(7 - 18*s1 - 9*t2) + t2 - s1*t2 + 
                   Power(t1,2)*
                    (-31 + 5*Power(s1,2) + 3*s1*(3 + t2)) + 
                   t1*(7 - 5*Power(s1,2) + s1*(11 + 2*t2)))) + 
             Power(s,2)*(-9*Power(t1,5) + 
                4*Power(s2,4)*(3 - 8*t1 + 5*Power(t1,2)) - 
                4*(-1 + t2)*t2 - 
                2*Power(t1,3)*
                 (32 + Power(s1,2) + s1*(-21 + t2) - 25*t2 + 
                   Power(t2,2)) - 
                Power(t1,2)*
                 (-24 + 29*s1 + 4*Power(s1,2) + 27*t2 + 
                   2*Power(t2,2)) + 
                t1*(9 - 8*t2 + 12*Power(t2,2) - 2*s1*(2 + t2)) + 
                Power(t1,4)*
                 (40 + 2*Power(s1,2) - 7*t2 + s1*(-21 + 4*t2)) + 
                2*Power(s2,3)*
                 (-11 - 15*Power(t1,3) + 
                   s1*(10 - 25*t1 + 13*Power(t1,2)) - 3*t2 + 
                   t1*(17 + t2) + Power(t1,2)*(9 + 4*t2)) + 
                2*Power(s2,2)*
                 (4 + 5*Power(t1,4) + 
                   Power(s1,2)*(3 - 9*t1 + 4*Power(t1,2)) + 
                   Power(t1,3)*(5 - 9*t2) - t1*(-21 + t2) + 3*t2 + 
                   Power(t1,2)*(-35 + 3*t2) + 
                   s1*(-18*Power(t1,3) - 2*(6 + t2) + 
                      3*Power(t1,2)*(6 + t2) + t1*(16 + t2))) + 
                s2*(-4 - 2*Power(s1,2)*t1*
                    (1 - 10*t1 + 5*Power(t1,2)) + 
                   Power(t1,2)*(45 - 41*t2) + t2 - 10*Power(t2,2) + 
                   Power(t1,3)*(-23 + 6*t2) + 
                   Power(t1,4)*(15 + 8*t2) + 
                   3*t1*(-11 + 6*t2 + 2*Power(t2,2)) + 
                   s1*(2 + 10*Power(t1,4) + 
                      Power(t1,3)*(37 - 12*t2) + t1*(59 - 2*t2) + 
                      4*t2 + 2*Power(t1,2)*(-50 + 3*t2)))) + 
             s*(-4*Power(s2,5)*Power(-1 + t1,2) + 
                t1*(5 + 2*s1 - 14*t2)*(-1 + t2) + 2*Power(-1 + t2,2) - 
                Power(t1,5)*(-9 + s1 + 2*t2) + 
                Power(t1,4)*(-24 + 4*Power(s1,2) + s1*(9 - 2*t2) + 
                   19*t2) - 2*Power(s2,4)*(-1 + t1)*
                 (4 - 4*Power(t1,2) + s1*(-5 + 3*t1) + t2 + t1*t2) + 
                Power(t1,3)*(20 - 8*Power(s1,2) + s1*(11 - 10*t2) - 
                   49*t2 + 6*Power(t2,2)) + 
                Power(t1,2)*(-2 + 8*Power(s1,2) + 17*t2 + 
                   10*Power(t2,2) + s1*(-17 + 2*t2)) - 
                2*Power(s2,3)*(-1 + t1)*
                 (4 + Power(s1,2)*(-3 + t1) + 2*Power(t1,3) + 
                   Power(t1,2)*(2 - 3*t2) - 2*t2 - t1*(8 + t2) + 
                   s1*(9 - 6*Power(t1,2) + t2 + t1*(3 + t2))) + 
                s2*(1 + 10*Power(t1,5) - 
                   2*Power(s1,2)*t1*(8 - 5*t1 + Power(t1,3)) - 
                   11*t2 + 10*Power(t2,2) + 
                   Power(t1,4)*(-31 + 5*t2) - 
                   2*Power(t1,2)*(25 - 51*t2 + Power(t2,2)) + 
                   2*Power(t1,3)*(26 - 23*t2 + Power(t2,2)) - 
                   2*t1*(-9 + 25*t2 + 9*Power(t2,2)) - 
                   2*s1*(-1 + Power(t1,3)*(27 - 4*t2) + 
                      Power(t1,2)*(-10 + t2) + t2 + 
                      Power(t1,4)*(-9 + 2*t2) - t1*(7 + 8*t2))) + 
                Power(s2,2)*(-12 + 
                   4*Power(s1,2)*
                    (2 + t1 - 3*Power(t1,2) + Power(t1,3)) + 25*t2 + 
                   12*Power(t2,2) + Power(t1,3)*(14 + t2) - 
                   2*Power(t1,4)*(5 + 2*t2) + 
                   Power(t1,2)*(-10 + 19*t2) + 
                   t1*(18 - 41*t2 - 8*Power(t2,2)) + 
                   s1*(3 - 6*Power(t1,4) + Power(t1,2)*(67 - 6*t2) - 
                      18*t2 + 3*Power(t1,3)*(-5 + 2*t2) + 
                      t1*(-49 + 10*t2)))))/
           (s*(s - s2 + t1)*Power(s2 - s2*t1 + t1*(-1 + s + t1),2)*
             (s - s1 + t2)) - 
          (2*Power(s,5)*(1 + 2*s2)*Power(t1,2) - 
             4*(s2 - t1)*(-1 + t1)*
              Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
             2*Power(s,4)*t1*
              (Power(s2,2)*(4 - 8*t1) + 4*Power(t1,2) + t2 + 
                t1*(2 + 3*t2) + 
                s2*(2 + 3*Power(t1,2) - t2 - t1*(2 + 4*s1 + t2))) + 
             Power(s,3)*(7*Power(t1,4) + 
                4*Power(s2,3)*(1 - 6*t1 + 6*Power(t1,2)) - 
                2*Power(t2,2) + Power(t1,3)*(15 - 4*s1 + 13*t2) + 
                2*t1*(1 - (3 + s1)*t2 + Power(t2,2)) - 
                Power(t1,2)*(20 + 9*t2 + 2*s1*(5 + 3*t2)) + 
                Power(s2,2)*
                 (2 - 20*Power(t1,3) - 2*t2 + 
                   t1*(6 - 18*s1 + 4*t2) + 
                   Power(t1,2)*(4 + 22*s1 + 6*t2)) + 
                s2*(2*Power(t1,4) + 2*t2 - 
                   Power(t1,3)*(1 + 12*s1 + 6*t2) + 
                   t1*(8 + 7*t2 + 2*s1*(1 + t2)) + 
                   Power(t1,2)*
                    (-37 + 4*Power(s1,2) - 21*t2 + 2*s1*(6 + t2)))) + 
             Power(s,2)*(3*Power(t1,5) - 
                8*Power(s2,4)*(1 - 3*t1 + 2*Power(t1,2)) + 
                4*(-1 + t2)*t2 + Power(t1,4)*(33 - s1 + t2) + 
                2*t1*(-4 + s1 + 4*t2 + 2*s1*t2 - 6*Power(t2,2)) - 
                2*Power(s2,3)*
                 (-5 - 11*Power(t1,3) + 
                   s1*(5 - 17*t1 + 10*Power(t1,2)) - 2*t2 + 
                   3*Power(t1,2)*(1 + t2) + t1*(13 + t2)) + 
                Power(t1,2)*
                 (23 + 2*Power(s1,2) + 2*t2 + 2*Power(t2,2) + 
                   2*s1*(7 + 6*t2)) + 
                Power(t1,3)*
                 (-51 + 2*Power(s1,2) - 17*t2 + 2*Power(t2,2) - 
                   s1*(5 + 16*t2)) + 
                Power(s2,2)*
                 (4 + 2*Power(s1,2)*(5 - 3*t1)*t1 - 6*Power(t1,4) + 
                   t2 - 23*t1*(3 + t2) + Power(t1,3)*(-23 + 12*t2) + 
                   Power(t1,2)*(94 + 20*t2) + 
                   s1*(24*Power(t1,3) + 2*(1 + t2) - t1*(1 + 2*t2) - 
                      Power(t1,2)*(35 + 4*t2))) - 
                s2*(-2 + 6*s1*t2 - 10*Power(t2,2) + 
                   4*Power(t1,4)*(s1 + t2) + 
                   Power(t1,3)*
                    (87 - 6*Power(s1,2) + 25*t2 - 6*s1*t2) + 
                   2*Power(t1,2)*
                    (-55 + 6*Power(s1,2) - 17*t2 - s1*(11 + 6*t2)) + 
                   t1*(25 + 2*Power(s1,2) + t2 + 6*Power(t2,2) + 
                      s1*(22 + 8*t2)))) + 
             s*(4*Power(s2,5)*Power(-1 + t1,2) - 
                2*t1*(3 + s1 - 7*t2)*(-1 + t2) - 2*Power(-1 + t2,2) + 
                2*Power(s2,4)*(-1 + t1)*
                 (6 - 4*Power(t1,2) + s1*(-5 + 3*t1) + t1*(-2 + t2) + 
                   t2) - Power(t1,5)*(-13 + s1 + 2*t2) - 
                Power(t1,4)*(39 - 9*t2 + 2*s1*(5 + 3*t2)) - 
                Power(t1,2)*(19 + 4*Power(s1,2) - 3*t2 + 
                   10*Power(t2,2) + 2*s1*(-7 + 5*t2)) + 
                Power(t1,3)*(41 + 6*t2 - 6*Power(t2,2) + 
                   s1*(-5 + 26*t2)) + 
                Power(s2,3)*(-1 + t1)*
                 (2*Power(s1,2)*(-3 + t1) + 4*Power(t1,3) + 
                   Power(t1,2)*(20 - 6*t2) + 3*(6 + t2) - 
                   t1*(42 + 11*t2) + 
                   s1*(13 - 12*Power(t1,2) + 2*t2 + t1*(13 + 2*t2))) - 
                Power(s2,2)*(7 + 
                   4*Power(s1,2)*
                    (1 + 3*t1 - 4*Power(t1,2) + Power(t1,3)) + 
                   Power(t1,4)*(10 - 4*t2) + 7*t2 + 12*Power(t2,2) - 
                   3*Power(t1,3)*(23 + 4*t2) + 
                   Power(t1,2)*(115 + 21*t2) - 
                   t1*(63 + 12*t2 + 8*Power(t2,2)) + 
                   s1*(6 - 6*Power(t1,4) - 10*t2 + 
                      2*Power(t1,2)*(22 + t2) - 3*t1*(15 + 2*t2) + 
                      Power(t1,3)*(1 + 6*t2))) + 
                s2*(-2*Power(t1,5) + 
                   2*Power(s1,2)*t1*
                    (4 + 3*t1 - 4*Power(t1,2) + Power(t1,3)) - 
                   5*Power(t1,4)*(8 + t2) - 
                   2*Power(t1,3)*(-52 - 2*t2 + Power(t2,2)) + 
                   Power(t1,2)*(-82 - 23*t2 + 2*Power(t2,2)) - 
                   2*(1 - 6*t2 + 5*Power(t2,2)) + 
                   2*t1*(11 + 6*t2 + 9*Power(t2,2)) + 
                   s1*(-8*t1 + 2*(-1 + t2) + Power(t1,4)*(-7 + 4*t2) + 
                      Power(t1,3)*(44 + 8*t2) - 
                      3*Power(t1,2)*(9 + 10*t2)))))/
           (s*(-1 + t1)*Power(s2 - s2*t1 + t1*(-1 + s + t1),2)*
             (-1 + t2)) + (32*s*(-1 + t1)*(s1 - s2 + t1 - t2) - 
             16*(1 - 2*Power(s2,3) - t1 - s1*t1 + 2*Power(t1,2) - 
                2*s1*Power(t1,2) + Power(s,2)*(-2 + s1 + 2*t1 - t2) - 
                t2 + 2*t1*t2 - Power(s2,2)*(-2 + s1 - 4*t1 + t2) + 
                s2*(-3 + t1 - 4*Power(t1,2) + s1*(3 + t1) - 2*t2 + 
                   2*t1*t2) - 
                s*(2 + Power(s1,2) - 2*Power(s2,2) + 13*t1 - 
                   18*Power(t1,2) + s2*(-15 + 20*t1 - 2*t2) + 
                   s1*(11 + s2 - 13*t1 - t2) - 13*t2 + 16*t1*t2)) + 
             (2*(-2 + 2*Power(s,2) + Power(s1,2) - 10*s2 + 
                  3*Power(s2,2) - 3*t2 + s2*t2 + 
                  s1*(-3 + 4*s2 + t2) - s*(-6 + 3*s1 + 5*s2 + t2))*
                (1 + 2*Power(s2,2)*(-1 + t1) + s*t1 - Power(t1,2) + 
                  s*Power(t1,2) + s1*(s2 - t1)*(1 + t1) - t2 + 
                  s*t2 + 3*t1*t2 - s*t1*t2 + 
                  s2*(-2*Power(t1,2) - 3*t2 + t1*(2 - 2*s + t2))))/
              (s2*(-1 + t1) - t1*(-1 + s + t1)) + 
             8*(6 + Power(s,3) - 25*s2 + 15*Power(s2,2) - 
                19*Power(s2,3) - 2*t1 + 17*s2*t1 + 
                36*Power(s2,2)*t1 + 8*Power(t1,2) - 
                36*s2*Power(t1,2) + Power(s1,2)*(-s2 + t1) + 
                Power(s,2)*(-15 + 6*s1 - s2 + 18*t1 - 8*t2) - 12*t2 - 
                20*s2*t2 - 10*Power(s2,2)*t2 + 24*t1*t2 + 
                18*s2*t1*t2 - 
                s1*(1 + 12*Power(s2,2) + 13*t1 + 18*Power(t1,2) - 
                   2*t2 + s2*(-33 - 11*t1 + 2*t2)) + 
                s*(-14 - 7*Power(s1,2) + 19*Power(s2,2) - 60*t1 + 
                   103*Power(t1,2) + 61*t2 - 85*t1*t2 + 
                   s1*(-47 - 5*s2 + 56*t1 + 9*t2) + 
                   s2*(73 - 121*t1 + 18*t2))) + 
             2*(-6 + 4*Power(s,3) - 22*s2 + 18*Power(s2,2) - 
                40*Power(s2,3) + 14*t1 + 62*s2*t1 + 
                80*Power(s2,2)*t1 - 10*Power(t1,2) - 
                80*s2*Power(t1,2) + 
                2*Power(s,2)*(-15 + 3*s1 + 2*s2 + 20*t1 - 5*t2) - 
                26*t2 - 57*s2*t2 - 22*Power(s2,2)*t2 + 75*t1*t2 + 
                40*s2*t1*t2 + 2*Power(s1,2)*(1 - s2 + t1 + t2) + 
                s1*(-28*Power(s2,2) - 41*t1 - 40*Power(t1,2) + 
                   s2*(89 + 28*t1 - 6*t2) + 8*t2) + 
                s*(-21 - 10*Power(s1,2) + 32*Power(s2,2) - 83*t1 + 
                   176*Power(t1,2) + 90*t2 - 136*t1*t2 + 
                   s1*(-68 - 6*s2 + 68*t1 + 18*t2) + 
                   s2*(110 - 216*t1 + 32*t2))) - 
             4*(11 + 3*Power(s,3) - 56*s2 + 19*Power(s2,2) - 
                47*Power(s2,3) + 7*t1 + 62*s2*t1 + 
                96*Power(s2,2)*t1 + 2*Power(t1,2) - 
                96*s2*Power(t1,2) + 
                Power(s,2)*(-43 + 13*s1 + 5*s2 + 48*t1 - 15*t2) - 
                Power(s1,2)*(2 + s2 - 3*t1 - 2*t2) - 39*t2 - 
                66*s2*t2 - 25*Power(s2,2)*t2 + 78*t1*t2 + 
                48*s2*t1*t2 + 
                s1*(1 - 29*Power(s2,2) - 42*t1 - 48*Power(t1,2) + 
                   4*s2*(22 + 8*t1 - t2) + 6*t2) + 
                s*(-31 - 16*Power(s1,2) + 39*Power(s2,2) - 118*t1 + 
                   230*Power(t1,2) - 
                   2*s1*(42 + 8*s2 - 51*t1 - 10*t2) + 128*t2 - 
                   182*t1*t2 + s2*(166 - 278*t1 + 40*t2))) + 
             ((3 + s - s1 - s2)*
                (4*(s2 - t1)*(-1 + t1)*
                   Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                  2*Power(s,3)*
                   (Power(t1,4) + Power(t2,2) - 
                     2*Power(t1,3)*(s2 + t2) - 2*t1*t2*(s2 + t2) + 
                     Power(t1,2)*
                      (2*Power(s2,2) + 2*s2*(-1 + t2) + 
                        Power(1 + t2,2))) - 
                  4*Power(s,2)*
                   ((2 + s1 + s2)*Power(t1,4) + 
                     t2*(-1 + s2 - s1*s2 + Power(s2,2) + t2 + 
                       3*s2*t2) + 
                     t1*(-1 + s2 + Power(s2,2) - 2*Power(s2,3) + 
                       s1*(1 + Power(s2,2) + s2*(-1 + t2)) + t2 - 
                       5*s2*t2 - 5*Power(s2,2)*t2 - 4*Power(t2,2) - 
                       4*s2*Power(t2,2)) - 
                     Power(t1,3)*
                      (2 + 3*Power(s2,2) + 5*t2 + 2*s2*(1 + t2) + 
                       s1*(-1 + 2*s2 + t2)) + 
                     Power(t1,2)*
                      (1 + 2*Power(s2,3) + 3*t2 + 3*Power(t2,2) + 
                        2*Power(s2,2)*(1 + t2) + s2*t2*(10 + t2) + 
                        s1*(-1 - s2 + Power(s2,2) + t2))) + 
                  2*s*(2*Power(s2,4)*Power(-1 + t1,2) + 
                     (5 + 4*s1 + Power(s1,2))*Power(t1,4) + 
                     Power(-1 + t2,2) - 
                     2*t1*(1 - 5*t2 + 4*Power(t2,2)) + 
                     Power(t1,2)*
                      (6 + Power(s1,2) + 4*s1*(-3 + t2) + 2*t2 + 
                        9*Power(t2,2)) - 
                     2*Power(t1,3)*(4*s1*(-1 + t2) + 5*(1 + t2)) + 
                     2*Power(s2,3)*(-1 + t1)*
                      (-2*Power(t1,2) + s1*(1 + t1) - 3*t2 + 
                        t1*(2 + t2)) + 
                     Power(s2,2)*
                      (-2 + 2*Power(t1,4) + 
                        Power(s1,2)*(1 + Power(t1,2)) + 10*t2 + 
                        9*Power(t2,2) - 2*Power(t1,3)*(3 + t2) - 
                        4*s1*
                       (-Power(t1,2) + Power(t1,3) + 2*t2 - t1*t2) + 
                        t1*(2 - 26*t2 - 8*Power(t2,2)) + 
                        Power(t1,2)*(4 + 18*t2 + Power(t2,2))) + 
                     2*s2*((1 + s1)*Power(t1,4) + 3*(-1 + t2)*t2 - 
                        Power(t1,3)*(3 + 4*s1 + Power(s1,2) + 6*t2) - 
                        t1*(1 + Power(s1,2) + 8*t2 + 8*Power(t2,2) - 
                        2*s1*(3 + t2)) + 
                        Power(t1,2)*
                        (3 + 17*t2 + 3*Power(t2,2) + s1*(-3 + 2*t2))))\
))/(s*Power(s2 - s2*t1 + t1*(-1 + s + t1),2)))/
           ((-1 + s1)*(-s + s1 - t2)) + 
          (32*s*(-1 + t1)*(-1 + s2 - t1 + t2) + 
             (2*(-2 + 2*Power(s,2) + Power(s1,2) - 10*s2 + 
                  3*Power(s2,2) - 3*t2 + s2*t2 + 
                  s1*(-3 + 4*s2 + t2) - s*(-5 + 3*s1 + 5*s2 + t2))*
                (1 + 2*Power(s2,2)*(-1 + t1) + s*t1 - Power(t1,2) + 
                  s*Power(t1,2) + s1*(s2 - t1)*(1 + t1) - t2 + s*t2 + 
                  3*t1*t2 - s*t1*t2 + 
                  s2*(-2*Power(t1,2) - 3*t2 + t1*(2 - 2*s + t2))))/
              (s2*(-1 + t1) - t1*(-1 + s + t1)) - 
             8*(Power(s,3) + 46*s2 - 38*Power(s2,2) + Power(s2,3) + 
                Power(s1,2)*(s2 - t1) - 29*t1 + 12*s2*t1 + 
                36*Power(s2,2)*t1 + 8*Power(t1,2) - 
                36*s2*Power(t1,2) - t2 - 42*s2*t2 + Power(s2,2)*t2 + 
                24*t1*t2 + 18*s2*t1*t2 - 
                Power(s,2)*(8 + s1 + s2 - 18*t1 + 8*t2) + 
                s1*(3 + 3*Power(s2,2) + 16*t1 - 18*Power(t1,2) - 
                   2*t2 + 2*s2*(-10 + 8*t1 + t2)) + 
                s*(-69 - Power(s2,2) - 5*t1 + 103*Power(t1,2) + 
                   72*t2 - 85*t1*t2 + s2*(97 - 121*t1 + 7*t2) + 
                   s1*(-8 + t1 + 9*t2))) - 
             2*(11 + 4*Power(s,3) + 50*s2 - 83*Power(s2,2) - 42*t1 + 
                58*s2*t1 + 80*Power(s2,2)*t1 - 10*Power(t1,2) - 
                80*s2*Power(t1,2) + 
                2*Power(s1,2)*(-2 + s2 - t1 - t2) - 21*t2 - 
                106*s2*t2 + 2*Power(s2,2)*t2 + 75*t1*t2 + 
                40*s2*t1*t2 - 
                2*Power(s,2)*(13 + 2*s1 - 2*s2 - 20*t1 + 5*t2) + 
                s1*(15 + 8*Power(s2,2) + 19*t1 - 40*Power(t1,2) - 
                   4*t2 + s2*(-44 + 32*t1 + 6*t2)) + 
                s*(-97 - 8*Power(s2,2) - 18*t1 + 176*Power(t1,2) + 
                   121*t2 - 136*t1*t2 + s2*(183 - 216*t1 + 8*t2) + 
                   s1*(-7 - 4*s2 + 6*t1 + 18*t2))) + 
             4*(3 + 5*Power(s,3) + 99*s2 - 90*Power(s2,2) - 
                Power(s2,3) - 63*t1 + 54*s2*t1 + 96*Power(s2,2)*t1 + 
                2*Power(t1,2) - 96*s2*Power(t1,2) + 
                Power(s1,2)*(s2 - 3*t1 - 2*t2) - 13*t2 - 114*s2*t2 + 
                Power(s2,2)*t2 + 78*t1*t2 + 48*s2*t1*t2 - 
                Power(s,2)*(22 + 6*s1 + 5*s2 - 48*t1 + 19*t2) + 
                s1*(5*Power(s2,2) - 
                   2*(-6 - 17*t1 + 24*Power(t1,2) + t2) + 
                   s2*(-45 + 40*t1 + 4*t2)) + 
                s*(-141 + Power(s1,2) + Power(s2,2) - 19*t1 + 
                   230*Power(t1,2) + 152*t2 - 182*t1*t2 + 
                   s2*(214 - 278*t1 + 18*t2) + 
                   s1*(-20 + 5*s2 + 5*t1 + 26*t2))) + 
             16*(Power(s,2)*(-1 + 2*t1 - t2) + 
                s*(s1*(-1 + t2) + s2*(17 - 20*t1 + t2) + 
                   2*(-7 + 9*Power(t1,2) + 7*t2 - 8*t1*t2)) + 
                2*(2*Power(s2,2)*(-1 + t1) + 
                   t1*(-2 + s1 + t1 - s1*t1 + t2) + 
                   s2*(3 + s1*(-1 + t1) - 2*Power(t1,2) - 2*t2 + t1*t2)\
)) + ((3 + s - s1 - s2)*(4*(s2 - t1)*(-1 + t1)*
                   Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                  2*Power(s,3)*
                   (Power(t1,4) + Power(t2,2) - 
                     2*Power(t1,3)*(s2 + t2) - 2*t1*t2*(s2 + t2) + 
                     Power(t1,2)*
                      (2*Power(s2,2) + 2*s2*(-1 + t2) + 
                        Power(1 + t2,2))) - 
                  4*Power(s,2)*
                   ((2 + s1 + s2)*Power(t1,4) + 
                     t2*(-1 + s2 - s1*s2 + Power(s2,2) + t2 + 
                       3*s2*t2) + 
                     t1*(-1 + s2 + Power(s2,2) - 2*Power(s2,3) + 
                        s1*(1 + Power(s2,2) + s2*(-1 + t2)) + t2 - 
                        5*s2*t2 - 5*Power(s2,2)*t2 - 4*Power(t2,2) - 
                        4*s2*Power(t2,2)) - 
                     Power(t1,3)*
                      (2 + 3*Power(s2,2) + 5*t2 + 2*s2*(1 + t2) + 
                        s1*(-1 + 2*s2 + t2)) + 
                     Power(t1,2)*
                      (1 + 2*Power(s2,3) + 3*t2 + 3*Power(t2,2) + 
                        2*Power(s2,2)*(1 + t2) + s2*t2*(10 + t2) + 
                        s1*(-1 - s2 + Power(s2,2) + t2))) + 
                  2*s*(2*Power(s2,4)*Power(-1 + t1,2) + 
                     (5 + 4*s1 + Power(s1,2))*Power(t1,4) + 
                     Power(-1 + t2,2) - 
                     2*t1*(1 - 5*t2 + 4*Power(t2,2)) + 
                     Power(t1,2)*
                      (6 + Power(s1,2) + 4*s1*(-3 + t2) + 2*t2 + 
                        9*Power(t2,2)) - 
                     2*Power(t1,3)*(4*s1*(-1 + t2) + 5*(1 + t2)) + 
                     2*Power(s2,3)*(-1 + t1)*
                      (-2*Power(t1,2) + s1*(1 + t1) - 3*t2 + 
                        t1*(2 + t2)) + 
                     Power(s2,2)*
                      (-2 + 2*Power(t1,4) + 
                        Power(s1,2)*(1 + Power(t1,2)) + 10*t2 + 
                        9*Power(t2,2) - 2*Power(t1,3)*(3 + t2) - 
                        4*s1*
                        (-Power(t1,2) + Power(t1,3) + 2*t2 - t1*t2) + 
                        t1*(2 - 26*t2 - 8*Power(t2,2)) + 
                        Power(t1,2)*(4 + 18*t2 + Power(t2,2))) + 
                     2*s2*((1 + s1)*Power(t1,4) + 3*(-1 + t2)*t2 - 
                        Power(t1,3)*(3 + 4*s1 + Power(s1,2) + 6*t2) - 
                        t1*(1 + Power(s1,2) + 8*t2 + 8*Power(t2,2) - 
                        2*s1*(3 + t2)) + 
                        Power(t1,2)*
                        (3 + 17*t2 + 3*Power(t2,2) + s1*(-3 + 2*t2)))))\
)/(s*Power(s2 - s2*t1 + t1*(-1 + s + t1),2)))/((-1 + s1)*(-1 + t2))) + 
       (((-64*(11 + Power(s,2) + s*(-5 + s2*(-2 + t1)) - 
                  Power(s2,2)*(-1 + t1) - t1 - 2*s2*t1)*
                (1 + 2*Power(s2,2)*(-1 + t1) + s*t1 - Power(t1,2) + 
                  s*Power(t1,2) + s1*(s2 - t1)*(1 + t1) - t2 + 
                  s*t2 + 3*t1*t2 - s*t1*t2 + 
                  s2*(-2*Power(t1,2) - 3*t2 + t1*(2 - 2*s + t2))))/
              (s2*(-1 + t1) - t1*(-1 + s + t1)) - 
             128*(-6 - 11*s2 - 3*Power(s2,3) + 8*t1 - 5*s2*t1 + 
                Power(s2,3)*t1 - 2*Power(t1,2) + 2*s2*Power(t1,2) + 
                s1*(7 + Power(s2,2)*(-1 + t1) - 3*t1 + 2*s2*t1 + 
                   2*Power(t1,2)) + 
                Power(s,2)*(-1 - 3*t1 + Power(t1,2) + 
                   s2*(-3 + 2*t1) - 4*t2) - 20*t2 - 2*s2*t2 - 
                2*Power(s2,2)*t2 + 
                s*(7 - 3*Power(s2,2)*(-2 + t1) - 3*t1 + 
                   s1*(-5 + s2 + 3*t1 - 2*s2*t1 - Power(t1,2)) + 
                   17*t2 - t1*t2 + s2*(6 + 4*t1 - Power(t1,2) + 4*t2)))\
)/(-1 + t1) - (64*(7 - 4*Power(s2,4)*(-1 + t1) - t1 - 17*s1*t1 + 
               4*Power(s,4)*(-1 + t1)*t1 - 35*Power(t1,2) + 
               2*s1*Power(t1,2) + 33*Power(t1,3) - 17*s1*Power(t1,3) - 
               4*Power(t1,4) + 4*s1*Power(t1,4) - 7*t2 + 14*t1*t2 + 
               21*Power(t1,2)*t2 + 
               Power(s2,3)*(-1 + t1)*
                (-8 + s1*(-7 + t1) - t1*(-4 + t2) + 3*t2) + 
               s2*(4 + 2*Power(t1,3) + 4*Power(t1,4) + 
                  s1*(17 - 6*t1 + 43*Power(t1,2) - 14*Power(t1,3)) + 
                  Power(t1,2)*(-52 + t2) - 17*t2 - 6*t1*(-7 + 4*t2)) + 
               Power(s2,2)*(-15 - 3*Power(t1,3) - 
                  s1*(-4 + 33*t1 - 18*Power(t1,2) + Power(t1,3)) + 
                  Power(t1,2)*(11 - 5*t2) + 11*t2 + t1*(7 + 6*t2)) + 
               Power(s,3)*(-4*s2*(1 - 5*t1 + 3*Power(t1,2)) + 
                  t1*(18 - 4*s1*(-1 + t1) - 13*t1 + 5*Power(t1,2) - 
                     3*t2 + t1*t2)) + 
               Power(s,2)*(-5*(1 + s1)*Power(t1,3) + 2*Power(t1,4) + 
                  4*Power(s2,2)*(3 - 8*t1 + 3*Power(t1,2)) + 
                  Power(t1,2)*(-9 + 23*s1 - 7*t2) - 
                  2*t1*(17 + 8*s1 - 4*t2) - t2 + 
                  s2*(14 - 10*Power(t1,3) + 
                     s1*(4 - 19*t1 + 9*Power(t1,2)) + 
                     Power(t1,2)*(35 - 3*t2) - 3*t2 + t1*(-41 + 10*t2)\
)) - s*(1 + 2*s1*Power(t1,4) + 
                  4*Power(s2,3)*(3 - 5*t1 + Power(t1,2)) - 8*t2 + 
                  Power(t1,2)*(-59 + 31*s1 + t2) + 
                  Power(t1,3)*(23 - 20*s1 + 2*t2) + 
                  t1*(-33 - 23*s1 + 7*t2) + 
                  s2*(19 - 10*Power(t1,3) + 2*Power(t1,4) + 
                     s1*(13 - 44*t1 + 45*Power(t1,2) - 
                       6*Power(t1,3)) - 4*t2 + 4*t1*(4 + 3*t2) - 
                     3*Power(t1,2)*(1 + 4*t2)) + 
                  Power(s2,2)*
                   (10 - 5*Power(t1,3) + 
                     s1*(11 - 23*t1 + 6*Power(t1,2)) + 
                     Power(t1,2)*(26 - 3*t2) - 6*t2 + t1*(-39 + 11*t2))\
)))/((s - s2 + t1)*(s2 - s2*t1 + t1*(-1 + s + t1))) + 
          (32*(-1 + s2)*((2*
                  (-5 - 2*Power(s,3) - 12*s2 + 2*Power(s2,2) + 
                    2*Power(s2,3) + 
                    s1*(-2 + 2*s2 + Power(s2,2) - 3*t1) - 5*t1 - 
                    6*s2*t1 + 
                    s*(10 - 6*Power(s2,2) + 5*t1 + 
                       s1*(-2 - 3*s2 + t1) + s2*(2 + t1 - 2*t2) - 
                       2*t2) - 3*t2 + 4*s2*t2 + Power(s2,2)*t2 + 
                    Power(s,2)*(-1 + 2*s1 + 6*s2 - t1 + t2))*
                  (1 + 2*Power(s2,2)*(-1 + t1) + s*t1 - 
                    Power(t1,2) + s*Power(t1,2) + 
                    s1*(s2 - t1)*(1 + t1) - t2 + s*t2 + 3*t1*t2 - 
                    s*t1*t2 + 
                    s2*(-2*Power(t1,2) - 3*t2 + t1*(2 - 2*s + t2))))/
                (s2*(-1 + t1) - t1*(-1 + s + t1)) + 
               4*(-8 - s2 + 7*Power(s2,2) + Power(s2,3) - 
                  Power(s2,4) + Power(s1,2)*(3 + 2*s2 - t1) + 9*t1 + 
                  11*s2*t1 + 3*Power(s2,2)*t1 - Power(t1,2) + 10*t2 - 
                  s2*t2 - 4*Power(s2,2)*t2 - Power(s2,3)*t2 + t1*t2 + 
                  3*s2*t1*t2 + 2*Power(s,3)*(-3 + s2 + t2) + 
                  Power(s,2)*
                   (7 - 5*Power(s2,2) + s2*(13 + t1 - 5*t2) - 6*t2 + 
                     t1*t2 - 2*s1*(-5 + s2 + t2)) - 
                  s1*(Power(s2,3) + s2*(2*t1 + 5*(-2 + t2)) + 
                     Power(s2,2)*(-2 + t2) - (-2 + t1)*(t1 + 2*t2)) - 
                  s*(-10 + 2*Power(s1,2) - 4*Power(s2,3) + 9*t1 + 
                     s1*(16 - 3*Power(s2,2) + s2*(11 + t1 - 3*t2) + 
                        t1*(-5 + t2) - 6*t2) + 
                     Power(s2,2)*(9 + t1 - 4*t2) + t2 + 3*t1*t2 + 
                     s2*(-9*(-2 + t2) + t1*(3 + t2)))) - 
               ((-5 + Power(s,2) + 3*s2 + Power(s2,2) - 
                    s*(1 + 2*s2) - 3*t1)*
                  (4*(s2 - t1)*(-1 + t1)*
                     Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                    2*Power(s,3)*
                     (Power(t1,4) + Power(t2,2) - 
                       2*Power(t1,3)*(s2 + t2) - 2*t1*t2*(s2 + t2) + 
                       Power(t1,2)*
                        (2*Power(s2,2) + 2*s2*(-1 + t2) + 
                        Power(1 + t2,2))) - 
                    4*Power(s,2)*
                     ((2 + s1 + s2)*Power(t1,4) + 
                       t2*(-1 + s2 - s1*s2 + Power(s2,2) + t2 + 
                       3*s2*t2) + 
                       t1*(-1 + s2 + Power(s2,2) - 2*Power(s2,3) + 
                       s1*(1 + Power(s2,2) + s2*(-1 + t2)) + t2 - 
                       5*s2*t2 - 5*Power(s2,2)*t2 - 4*Power(t2,2) - 
                       4*s2*Power(t2,2)) - 
                       Power(t1,3)*
                       (2 + 3*Power(s2,2) + 5*t2 + 2*s2*(1 + t2) + 
                       s1*(-1 + 2*s2 + t2)) + 
                       Power(t1,2)*
                        (1 + 2*Power(s2,3) + 3*t2 + 3*Power(t2,2) + 
                        2*Power(s2,2)*(1 + t2) + s2*t2*(10 + t2) + 
                        s1*(-1 - s2 + Power(s2,2) + t2))) + 
                    2*s*(2*Power(s2,4)*Power(-1 + t1,2) + 
                       (5 + 4*s1 + Power(s1,2))*Power(t1,4) + 
                       Power(-1 + t2,2) - 
                       2*t1*(1 - 5*t2 + 4*Power(t2,2)) + 
                       Power(t1,2)*
                        (6 + Power(s1,2) + 4*s1*(-3 + t2) + 2*t2 + 
                        9*Power(t2,2)) - 
                       2*Power(t1,3)*(4*s1*(-1 + t2) + 5*(1 + t2)) + 
                       2*Power(s2,3)*(-1 + t1)*
                        (-2*Power(t1,2) + s1*(1 + t1) - 3*t2 + 
                        t1*(2 + t2)) + 
                       Power(s2,2)*
                        (-2 + 2*Power(t1,4) + 
                        Power(s1,2)*(1 + Power(t1,2)) + 10*t2 + 
                        9*Power(t2,2) - 2*Power(t1,3)*(3 + t2) - 
                        4*s1*
                       (-Power(t1,2) + Power(t1,3) + 2*t2 - t1*t2) + 
                        t1*(2 - 26*t2 - 8*Power(t2,2)) + 
                        Power(t1,2)*(4 + 18*t2 + Power(t2,2))) + 
                       2*s2*((1 + s1)*Power(t1,4) + 3*(-1 + t2)*t2 - 
                        Power(t1,3)*(3 + 4*s1 + Power(s1,2) + 6*t2) - 
                        t1*(1 + Power(s1,2) + 8*t2 + 8*Power(t2,2) - 
                        2*s1*(3 + t2)) + 
                        Power(t1,2)*
                        (3 + 17*t2 + 3*Power(t2,2) + s1*(-3 + 2*t2))))\
))/(s*Power(s2 - s2*t1 + t1*(-1 + s + t1),2))))/
           ((s - s2 + t1)*(s - s1 + t2)) + 
          (32*(-1 + s2)*((2*
                  (1 + 2*Power(s2,2)*(-1 + t1) + s*t1 - 
                    Power(t1,2) + s*Power(t1,2) + 
                    s1*(s2 - t1)*(1 + t1) - t2 + s*t2 + 3*t1*t2 - 
                    s*t1*t2 + 
                    s2*(-2*Power(t1,2) - 3*t2 + t1*(2 - 2*s + t2)))*
                  (-6 - 16*s2 - Power(s2,2) + 2*Power(s2,3) + t1 + 
                    Power(s,2)*(-1 + 2*s2 + t1) + 
                    s1*(-2 + 2*s2 + Power(s2,2) + t1) - 7*t2 + 
                    Power(s2,2)*t2 + 
                    s*(2 - 4*Power(s2,2) + 3*t1 - 
                       s1*(-1 + 2*s2 + t1) + 2*t2 - s2*(-5 + t1 + t2)\
)))/(s2*(-1 + t1) - t1*(-1 + s + t1)) - 
               4*(-3 - 4*s2 - 11*Power(s2,2) - 2*Power(s2,3) + 
                  Power(s2,4) - Power(s1,2)*(-1 + t1) + t1 - s2*t1 - 
                  Power(s2,2)*t1 + t2 - 12*s2*t2 - 2*Power(s2,2)*t2 + 
                  Power(s2,3)*t2 + 2*t1*t2 + 2*s2*t1*t2 + 
                  s1*(-1 + Power(s2,3) - 4*t2 + 
                     Power(s2,2)*(1 + t2) + s2*(-1 - t1 + t2) + 
                     t1*(-3 + 2*t2)) + 
                  Power(s,2)*
                   (2*Power(s2,2) + t1*(-3 + t2) - 2*t2 + 
                     s2*(-3 + t1 + 2*t2)) + 
                  s*(3 - 3*Power(s2,3) + Power(t1,2) + t2 - 2*t1*t2 - 
                     Power(s2,2)*(-6 + t1 + 3*t2) + 
                     s2*(5 - t1*(-5 + t2) + 7*t2) - 
                     s1*(2 + 2*Power(s2,2) + t1*(-5 + t2) - 3*t2 + 
                        s2*(-2 + t1 + 2*t2)))) + 
               ((7 + s*(-1 + s2) - s2 - Power(s2,2) - t1)*
                  (4*(s2 - t1)*(-1 + t1)*
                     Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                    2*Power(s,3)*
                     (Power(t1,4) + Power(t2,2) - 
                       2*Power(t1,3)*(s2 + t2) - 2*t1*t2*(s2 + t2) + 
                       Power(t1,2)*
                        (2*Power(s2,2) + 2*s2*(-1 + t2) + 
                        Power(1 + t2,2))) - 
                    4*Power(s,2)*
                     ((2 + s1 + s2)*Power(t1,4) + 
                       t2*(-1 + s2 - s1*s2 + Power(s2,2) + t2 + 
                       3*s2*t2) + 
                       t1*(-1 + s2 + Power(s2,2) - 2*Power(s2,3) + 
                       s1*(1 + Power(s2,2) + s2*(-1 + t2)) + t2 - 
                       5*s2*t2 - 5*Power(s2,2)*t2 - 4*Power(t2,2) - 
                       4*s2*Power(t2,2)) - 
                       Power(t1,3)*
                       (2 + 3*Power(s2,2) + 5*t2 + 2*s2*(1 + t2) + 
                       s1*(-1 + 2*s2 + t2)) + 
                       Power(t1,2)*
                        (1 + 2*Power(s2,3) + 3*t2 + 3*Power(t2,2) + 
                        2*Power(s2,2)*(1 + t2) + s2*t2*(10 + t2) + 
                        s1*(-1 - s2 + Power(s2,2) + t2))) + 
                    2*s*(2*Power(s2,4)*Power(-1 + t1,2) + 
                       (5 + 4*s1 + Power(s1,2))*Power(t1,4) + 
                       Power(-1 + t2,2) - 
                       2*t1*(1 - 5*t2 + 4*Power(t2,2)) + 
                       Power(t1,2)*
                        (6 + Power(s1,2) + 4*s1*(-3 + t2) + 2*t2 + 
                        9*Power(t2,2)) - 
                       2*Power(t1,3)*(4*s1*(-1 + t2) + 5*(1 + t2)) + 
                       2*Power(s2,3)*(-1 + t1)*
                        (-2*Power(t1,2) + s1*(1 + t1) - 3*t2 + 
                        t1*(2 + t2)) + 
                       Power(s2,2)*
                        (-2 + 2*Power(t1,4) + 
                        Power(s1,2)*(1 + Power(t1,2)) + 10*t2 + 
                        9*Power(t2,2) - 2*Power(t1,3)*(3 + t2) - 
                        4*s1*
                       (-Power(t1,2) + Power(t1,3) + 2*t2 - t1*t2) + 
                        t1*(2 - 26*t2 - 8*Power(t2,2)) + 
                        Power(t1,2)*(4 + 18*t2 + Power(t2,2))) + 
                       2*s2*((1 + s1)*Power(t1,4) + 3*(-1 + t2)*t2 - 
                        Power(t1,3)*(3 + 4*s1 + Power(s1,2) + 6*t2) - 
                        t1*(1 + Power(s1,2) + 8*t2 + 8*Power(t2,2) - 
                        2*s1*(3 + t2)) + 
                        Power(t1,2)*
                        (3 + 17*t2 + 3*Power(t2,2) + s1*(-3 + 2*t2))))\
))/(s*Power(s2 - s2*t1 + t1*(-1 + s + t1),2))))/((-1 + t1)*(-1 + t2)) + 
          (16*(-1 + s2)*((4*
                  (-1 - 2*Power(s2,2)*(-1 + t1) - s*t1 + 
                    Power(t1,2) - s*Power(t1,2) - 
                    s1*(s2 - t1)*(1 + t1) + t2 - s*t2 - 3*t1*t2 + 
                    s*t1*t2 + 
                    s2*(2*Power(t1,2) + t1*(-2 + 2*s - t2) + 3*t2))*
                  (2 + 2*Power(s,3) + Power(s1,2)*(-1 + s2) + 
                    20*s2 + 3*Power(s2,2) - Power(s2,3) - 2*t1 + 
                    s1*(4 - 5*s2 + Power(s2,2) + 2*t1 - 2*t2) + 
                    10*t2 + s2*t2 - Power(s2,2)*t2 - 
                    2*Power(s,2)*(-1 + s1 + 3*s2 + t2) + 
                    s*(-15 + 5*Power(s2,2) + 3*s2*(-1 + t2) + 
                       s1*(3 + 2*s2 + t2))))/
                (s2*(-1 + t1) - t1*(-1 + s + t1)) + 
               8*(1 + 4*s2 + 11*Power(s2,2) + Power(s2,3) - 4*t1 - 
                  4*s2*t1 + 3*Power(t1,2) + 
                  s1*(-13 + Power(s2,3) + Power(t1,2) + 
                     s2*(5 + t1 - 5*t2) + Power(s2,2)*(-5 + t2) - 
                     t1*(-2 + t2) - 3*t2) + 
                  Power(s1,2)*
                   (7 + Power(s2,2) + t1 + s2*(-1 + t2) - t2) + 
                  12*t2 + 13*s2*t2 + 2*Power(s2,2)*t2 - 2*t1*t2 - 
                  s2*t1*t2 + 3*Power(t2,2) + s2*Power(t2,2) + 
                  2*Power(s,3)*(-3 + s2 + t2) - 
                  2*Power(s,2)*
                   (-7 + 2*Power(s2,2) + 2*s2*(-2 + t2) + t2 + 
                     s1*(-4 + s2 + t2)) + 
                  s*(-5 - 4*Power(s1,2) + 2*Power(s2,3) + t1 - 
                     2*Power(t1,2) + s2*(-20 + t1 - t2) - 11*t2 + 
                     2*t1*t2 - 2*Power(t2,2) + 
                     Power(s2,2)*(-3 + 2*t2) + 
                     s1*(-8 + Power(s2,2) + s2*(-1 + t2) + 7*t2))) + 
               ((1 + s - s2)*
                  (1 + 2*Power(s2,2)*(-1 + t1) + s*t1 - 
                    Power(t1,2) + s*Power(t1,2) + 
                    s1*(s2 - t1)*(1 + t1) - t2 + s*t2 + 3*t1*t2 - 
                    s*t1*t2 + 
                    s2*(-2*Power(t1,2) - 3*t2 + t1*(2 - 2*s + t2)))*
                  (4*s*Power(1 + 2*Power(s2,2)*(-1 + t1) + s*t1 - 
                       Power(t1,2) + s*Power(t1,2) + 
                       s1*(s2 - t1)*(1 + t1) - t2 + s*t2 + 
                       3*t1*t2 - s*t1*t2 + 
                       s2*(-2*Power(t1,2) - 3*t2 + 
                       t1*(2 - 2*s + t2)),2) - 
                    12*(s2*(-1 + t1) - t1*(-1 + s + t1))*
                     (-Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,
                       2) - 
                       Power(s,2)*(-1 + s2)*
                       (s2*t1 - Power(t1,2) - t2 + t1*t2) + 
                       s*(-1 + s2)*
                        (Power(s2,2)*(-1 + t1) + 
                        2*t1*(1 - t1 + t2) + 
                        s1*(1 + s2 - 2*t1 + s2*t1 - Power(t1,2) - 
                       t2 + t1*t2) - 
                        s2*(1 + Power(t1,2) + 3*t2 - t1*(2 + t2))))))/
                (s*Power(s2 - s2*t1 + t1*(-1 + s + t1),3)) - 
               (2*(-9 + 2*Power(s,2) + 2*s1 - 3*s2 + 2*Power(s2,2) - 
                    t2 + s2*t2 - s*(s1 + 4*s2 + t2))*
                  (4*(s2 - t1)*(-1 + t1)*
                     Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                    2*Power(s,3)*
                     (Power(t1,4) + Power(t2,2) - 
                       2*Power(t1,3)*(s2 + t2) - 2*t1*t2*(s2 + t2) + 
                       Power(t1,2)*
                        (2*Power(s2,2) + 2*s2*(-1 + t2) + 
                        Power(1 + t2,2))) - 
                    4*Power(s,2)*
                     ((2 + s1 + s2)*Power(t1,4) + 
                       t2*(-1 + s2 - s1*s2 + Power(s2,2) + t2 + 
                       3*s2*t2) + 
                       t1*(-1 + s2 + Power(s2,2) - 2*Power(s2,3) + 
                       s1*(1 + Power(s2,2) + s2*(-1 + t2)) + t2 - 
                       5*s2*t2 - 5*Power(s2,2)*t2 - 4*Power(t2,2) - 
                       4*s2*Power(t2,2)) - 
                       Power(t1,3)*
                       (2 + 3*Power(s2,2) + 5*t2 + 2*s2*(1 + t2) + 
                       s1*(-1 + 2*s2 + t2)) + 
                       Power(t1,2)*
                        (1 + 2*Power(s2,3) + 3*t2 + 3*Power(t2,2) + 
                        2*Power(s2,2)*(1 + t2) + s2*t2*(10 + t2) + 
                        s1*(-1 - s2 + Power(s2,2) + t2))) + 
                    2*s*(2*Power(s2,4)*Power(-1 + t1,2) + 
                       (5 + 4*s1 + Power(s1,2))*Power(t1,4) + 
                       Power(-1 + t2,2) - 
                       2*t1*(1 - 5*t2 + 4*Power(t2,2)) + 
                       Power(t1,2)*
                        (6 + Power(s1,2) + 4*s1*(-3 + t2) + 2*t2 + 
                        9*Power(t2,2)) - 
                       2*Power(t1,3)*(4*s1*(-1 + t2) + 5*(1 + t2)) + 
                       2*Power(s2,3)*(-1 + t1)*
                        (-2*Power(t1,2) + s1*(1 + t1) - 3*t2 + 
                        t1*(2 + t2)) + 
                       Power(s2,2)*
                        (-2 + 2*Power(t1,4) + 
                        Power(s1,2)*(1 + Power(t1,2)) + 10*t2 + 
                        9*Power(t2,2) - 2*Power(t1,3)*(3 + t2) - 
                        4*s1*
                       (-Power(t1,2) + Power(t1,3) + 2*t2 - t1*t2) + 
                        t1*(2 - 26*t2 - 8*Power(t2,2)) + 
                        Power(t1,2)*(4 + 18*t2 + Power(t2,2))) + 
                       2*s2*((1 + s1)*Power(t1,4) + 3*(-1 + t2)*t2 - 
                        Power(t1,3)*(3 + 4*s1 + Power(s1,2) + 6*t2) - 
                        t1*(1 + Power(s1,2) + 8*t2 + 8*Power(t2,2) - 
                        2*s1*(3 + t2)) + 
                        Power(t1,2)*
                        (3 + 17*t2 + 3*Power(t2,2) + s1*(-3 + 2*t2))))\
))/(s*Power(s2 - s2*t1 + t1*(-1 + s + t1),2))))/
           ((-1 + s1)*(-s + s1 - t2)) + 
          (8*(-1 + s2)*((8*(-1 - 2*Power(s2,2)*(-1 + t1) - s*t1 + 
                    Power(t1,2) - s*Power(t1,2) - 
                    s1*(s2 - t1)*(1 + t1) + t2 - s*t2 - 3*t1*t2 + 
                    s*t1*t2 + 
                    s2*(2*Power(t1,2) + t1*(-2 + 2*s - t2) + 3*t2))*
                  (3 + 2*Power(s,3) + Power(s1,2)*(-1 + s2) + 19*s2 + 
                    3*Power(s2,2) - Power(s2,3) - 2*t1 + 
                    s1*(5 - 6*s2 + Power(s2,2) + 2*t1 - 2*t2) + 
                    8*t2 + 3*s2*t2 - Power(s2,2)*t2 - 
                    2*Power(s,2)*(s1 + 3*s2 + t2) + 
                    s*(5*Power(s2,2) - 2*(6 + t1 - t2) + 
                       3*s2*(1 + t2) + s1*(2 + 2*s2 + t2))))/
                (s2*(-1 + t1) - t1*(-1 + s + t1)) + 
               16*(6 + 2*s2 + 10*Power(s2,2) + Power(s2,3) + 3*t1 + 
                  s2*t1 - 3*Power(t1,2) - t2 + 9*s2*t2 + 
                  Power(s2,2)*t2 + 2*t1*t2 + s2*t1*t2 - 
                  3*Power(t2,2) - s2*Power(t2,2) + 
                  2*Power(s,3)*(-1 + s2 + t2) + 
                  s1*(-8 + Power(s2,3) - Power(t1,2) + 
                     s2*(8 + 2*t1 - 5*t2) + Power(s2,2)*(-5 + t2) + 
                     t1*(-1 + t2) + 7*t2) + 
                  Power(s1,2)*(4 + Power(s2,2) - t1 - t2 + s2*t2) - 
                  2*Power(s,2)*
                   (-3 + 2*Power(s2,2) + s1*(-2 + s2 + t2) + 
                     s2*(-1 + 2*t2)) + 
                  s*(-2*Power(s1,2) + 2*Power(s2,3) + 
                     s1*(-3 + Power(s2,2) + 2*t1 + s2*(-2 + t2) + 
                       t2) + Power(s2,2)*(1 + 2*t2) + 
                     s2*(-13 - 4*t1 + 2*t2) + 
                     2*(-3 + Power(t1,2) - 4*t2 - t1*t2 + Power(t2,2)))\
) + (2*(1 + s - s2)*(1 + 2*Power(s2,2)*(-1 + t1) + s*t1 - 
                    Power(t1,2) + s*Power(t1,2) + 
                    s1*(s2 - t1)*(1 + t1) - t2 + s*t2 + 3*t1*t2 - 
                    s*t1*t2 + 
                    s2*(-2*Power(t1,2) - 3*t2 + t1*(2 - 2*s + t2)))*
                  (4*s*Power(1 + 2*Power(s2,2)*(-1 + t1) + s*t1 - 
                       Power(t1,2) + s*Power(t1,2) + 
                       s1*(s2 - t1)*(1 + t1) - t2 + s*t2 + 
                       3*t1*t2 - s*t1*t2 + 
                       s2*(-2*Power(t1,2) - 3*t2 + 
                       t1*(2 - 2*s + t2)),2) - 
                    12*(s2*(-1 + t1) - t1*(-1 + s + t1))*
                     (-Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                       Power(s,2)*(-1 + s2)*
                        (s2*t1 - Power(t1,2) - t2 + t1*t2) + 
                       s*(-1 + s2)*
                        (Power(s2,2)*(-1 + t1) + 2*t1*(1 - t1 + t2) + 
                        s1*(1 + s2 - 2*t1 + s2*t1 - Power(t1,2) - 
                        t2 + t1*t2) - 
                        s2*(1 + Power(t1,2) + 3*t2 - t1*(2 + t2))))))/
                (s*Power(s2 - s2*t1 + t1*(-1 + s + t1),3)) - 
               (4*(-9 + 2*Power(s,2) + 2*s1 - 3*s2 + 2*Power(s2,2) - 
                    t2 + s2*t2 - s*(2 + s1 + 4*s2 + t2))*
                  (4*(s2 - t1)*(-1 + t1)*
                     Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                    2*Power(s,3)*
                     (Power(t1,4) + Power(t2,2) - 
                       2*Power(t1,3)*(s2 + t2) - 2*t1*t2*(s2 + t2) + 
                       Power(t1,2)*
                        (2*Power(s2,2) + 2*s2*(-1 + t2) + 
                        Power(1 + t2,2))) - 
                    4*Power(s,2)*
                     ((2 + s1 + s2)*Power(t1,4) + 
                       t2*(-1 + s2 - s1*s2 + Power(s2,2) + t2 + 
                        3*s2*t2) + 
                       t1*(-1 + s2 + Power(s2,2) - 2*Power(s2,3) + 
                        s1*(1 + Power(s2,2) + s2*(-1 + t2)) + t2 - 
                        5*s2*t2 - 5*Power(s2,2)*t2 - 4*Power(t2,2) - 
                        4*s2*Power(t2,2)) - 
                       Power(t1,3)*
                        (2 + 3*Power(s2,2) + 5*t2 + 2*s2*(1 + t2) + 
                        s1*(-1 + 2*s2 + t2)) + 
                       Power(t1,2)*
                        (1 + 2*Power(s2,3) + 3*t2 + 3*Power(t2,2) + 
                        2*Power(s2,2)*(1 + t2) + s2*t2*(10 + t2) + 
                        s1*(-1 - s2 + Power(s2,2) + t2))) + 
                    2*s*(2*Power(s2,4)*Power(-1 + t1,2) + 
                       (5 + 4*s1 + Power(s1,2))*Power(t1,4) + 
                       Power(-1 + t2,2) - 
                       2*t1*(1 - 5*t2 + 4*Power(t2,2)) + 
                       Power(t1,2)*
                        (6 + Power(s1,2) + 4*s1*(-3 + t2) + 2*t2 + 
                        9*Power(t2,2)) - 
                       2*Power(t1,3)*(4*s1*(-1 + t2) + 5*(1 + t2)) + 
                       2*Power(s2,3)*(-1 + t1)*
                        (-2*Power(t1,2) + s1*(1 + t1) - 3*t2 + 
                        t1*(2 + t2)) + 
                       Power(s2,2)*
                        (-2 + 2*Power(t1,4) + 
                        Power(s1,2)*(1 + Power(t1,2)) + 10*t2 + 
                        9*Power(t2,2) - 2*Power(t1,3)*(3 + t2) - 
                        4*s1*
                        (-Power(t1,2) + Power(t1,3) + 2*t2 - t1*t2) + 
                        t1*(2 - 26*t2 - 8*Power(t2,2)) + 
                        Power(t1,2)*(4 + 18*t2 + Power(t2,2))) + 
                       2*s2*((1 + s1)*Power(t1,4) + 3*(-1 + t2)*t2 - 
                        Power(t1,3)*(3 + 4*s1 + Power(s1,2) + 6*t2) - 
                        t1*(1 + Power(s1,2) + 8*t2 + 8*Power(t2,2) - 
                        2*s1*(3 + t2)) + 
                        Power(t1,2)*
                        (3 + 17*t2 + 3*Power(t2,2) + s1*(-3 + 2*t2)))))\
)/(s*Power(s2 - s2*t1 + t1*(-1 + s + t1),2))))/((-1 + s1)*(-1 + t2)))/
        Power(-1 + s2,2))*B1(t1,s,s2))/(128.*Power(Pi,2)) + 
  (((8*(-3 + 2*Power(s1,3)*Power(t1,2) - (-1 + s)*Power(t1,4) + t2 + 
            11*s*t2 - 13*s2*t2 - 2*Power(t2,2) - 2*s2*Power(t2,2) + 
            2*Power(t1,2)*(-1 + t2)*(-1 + s + t2) + 
            Power(t1,3)*(-1 + (-1 + s - s2)*t2) + 
            Power(s1,2)*t1*(6 + s2*(-1 + t1) + 3*Power(t1,2) - 2*t2 - 
               t1*(3 + 2*s + 4*t2)) + 
            t1*(1 + (-2 + 6*s2)*t2 + 2*(1 + s2)*Power(t2,2) - 
               s*(3 + 4*t2 + 2*Power(t2,2))) + 
            s1*(Power(t1,4) + Power(t1,3)*(5 - 3*s - 2*t2) - 6*t2 + 
               Power(t1,2)*(2 + 5*s - 5*t2 + 3*s*t2) + 
               t1*(-2 - 8*s + 5*t2 + 3*s*t2 + 4*Power(t2,2)) + 
               s2*(3 + Power(t1,3) + 2*t2 + t1*(12 + t2) - 
                  Power(t1,2)*(8 + 3*t2)))))/
        ((-1 + t1)*(s1*t1 - t2)*(-1 + t2)) + 
       (8*(-1 - 3*(-1 + s)*Power(t1,3) + 5*t2 + s*t2 - 5*s2*t2 - 
            2*Power(t2,2) - 8*s*Power(t2,2) + 8*s2*Power(t2,2) + 
            2*Power(t2,3) + 2*Power(s1,3)*t1*(3 - t1 + t2) + 
            Power(t1,2)*(1 - (3 + s2)*t2 + 2*Power(t2,2) + 
               s*(-4 + 3*t2)) - 
            t1*(3 + (2 + 8*s2)*t2 - 2*(-1 + s2)*Power(t2,2) + 
               2*Power(t2,3) + s*(1 - 12*t2 + 4*Power(t2,2))) - 
            Power(s1,2)*(2*Power(t1,3) + 2*t2*(3 + t2) + 
               Power(t1,2)*(-9 - 2*s + 3*t2) + 
               t1*(1 + 5*t2 + 2*s*(2 + t2)) + 
               s2*(-1 + 2*Power(t1,2) - t2 - t1*(5 + 3*t2))) + 
            s1*(-1 + Power(t1,3)*(1 + 2*s - 2*t2) + 5*s*t2 + 
               7*Power(t2,2) + 3*s*Power(t2,2) + 
               Power(t1,2)*(3 + s*(-7 + t2) + 5*t2 + 2*Power(t2,2)) + 
               t1*(-3 - 7*t2 + Power(t2,2) + 
                  s*(-1 + 4*t2 + Power(t2,2))) - 
               s2*(-1 - 7*Power(t1,2) + 2*Power(t1,3) + 7*t2 + 
                  3*Power(t2,2) + t1*(-8 + 7*t2 + Power(t2,2))))))/
        ((-1 + s1)*(s1*t1 - t2)*(-1 + t2)) + 
       (8*(-1 + s1 + t1 - t2)*
          (1 - Power(t1,2) + s1*(s2*(-1 + t1) - t1*(1 + t1 - 2*t2)) + 
            s*(1 + t1)*(t1 - t2) + t2 - s2*t2 + t1*t2 + s2*t1*t2 - 
            2*Power(t2,2))*((-2 + 3*Power(s2,2) - 
               Power(s1,3)*(-3 + t1) + t1 + 3*s2*t1 + 
               3*Power(s2,2)*t1 + 4*Power(t1,2) + s2*Power(t1,2) - 
               Power(s2,2)*Power(t1,2) + Power(t1,3) - 2*t2 + 
               4*s2*t2 - 4*Power(s2,2)*t2 + t1*t2 + 6*s2*t1*t2 - 
               Power(t1,2)*t2 - s2*Power(t1,2)*t2 - 
               5*s2*Power(t2,2) + s2*t1*Power(t2,2) - 2*Power(t2,3) - 
               Power(s,2)*(1 + s1*(-1 + t1) + t1 + Power(t1,2) + 
                  t1*t2) + Power(s1,2)*
                (-3 - Power(t1,2) - 2*t1*(-2 + t2) + 2*t2 + 
                  Power(t2,2) + s2*(6 - 2*t1 + t2)) + 
               s*(2*Power(s1,2)*(-2 + t1) - t1 - Power(t1,2) - 
                  3*t2 - 6*t1*t2 + Power(t1,2)*t2 + 2*Power(t2,2) - 
                  t1*Power(t2,2) + 
                  s2*(-1 + 2*Power(t1,2) + t1*(-2 + t2) + 3*t2) + 
                  s1*(3 - 2*t1 + 2*Power(t1,2) + 
                     s2*(-4 + 2*t1 - t2) - 2*t2 + 3*t1*t2 - 
                     Power(t2,2))) + 
               s1*(-1 + s2*(-2*Power(t1,2) - 2*t1*(-3 + t2) + 
                     Power(-1 + t2,2)) - Power(t1,2)*(-3 + t2) + t2 - 
                  3*Power(t2,2) + Power(s2,2)*(3 - t1 + t2) + 
                  t1*(4 + 5*t2 + Power(t2,2))))/((-1 + s1)*(-1 + t2)) + 
            (-9 - 2*s2 + 14*Power(s2,2) - Power(s2,3) + 6*t1 - 
               7*s2*t1 - 3*Power(s2,2)*t1 - 6*Power(t1,2) + 
               s2*Power(t1,2) + Power(t1,3) - t2 + 11*s2*t2 + 
               5*Power(s2,2)*t2 - Power(s2,3)*t2 + 7*t1*t2 - 
               2*s2*t1*t2 - 5*Power(t1,2)*t2 + 3*Power(t2,2) + 
               6*s2*Power(t2,2) - Power(s2,2)*Power(t2,2) + 
               3*t1*Power(t2,2) + Power(t2,3) + 
               Power(s,2)*(4 + s2 - 2*t1 + 3*t2 - s2*t2 - 
                  Power(t2,2)) - 
               Power(s1,2)*(-1 - t1 + 2*t2 + Power(t2,2) + 
                  s2*(2 + t2)) + 
               s1*(2 + 3*Power(t1,2) + 12*t2 + 5*Power(t2,2) - 
                  5*t1*(2 + t2) - Power(s2,2)*(3 + 2*t2) + 
                  s2*(16 - 3*t1 + 3*t2 - 2*Power(t2,2))) + 
               s*(3 + 5*t1 - Power(t1,2) - 11*t2 + 2*Power(s2,2)*t2 + 
                  2*t1*t2 - 5*Power(t2,2) + 
                  s2*(-18 + 5*t1 - 8*t2 + 2*Power(t2,2)) + 
                  s1*(-6 + s2 + 2*t1 - t2 + 2*s2*t2 + 2*Power(t2,2))))/
             ((-1 + s2)*(-1 + t1)) - 
            (6 + 5*s2 - 6*Power(s2,2) - 2*Power(s2,3) + 
               2*Power(s1,3)*(-3 + t1) + 2*t1 - s2*t1 + 
               Power(s2,2)*t1 + Power(s2,3)*t1 - 4*Power(t1,2) - 
               5*s2*Power(t1,2) - 4*Power(t1,3) + s2*Power(t1,3) + 
               t2 - 4*s2*t2 - 4*Power(s2,2)*t2 + t1*t2 + 6*s2*t1*t2 + 
               Power(s2,2)*t1*t2 + 2*Power(t1,2)*t2 - 
               s2*Power(t1,2)*t2 - 4*Power(t2,2) - 4*s2*Power(t2,2) + 
               t1*Power(t2,2) - Power(t2,3) + 
               Power(s,2)*(s1 + t1)*(-5 + s2 + t1 + t2) + 
               Power(s1,2)*(Power(s2,2) + 3*Power(t1,2) + 
                  3*(2 + t2) - t1*(12 + t2) + s2*(-7 + 2*t1 + t2)) + 
               s1*(4 + Power(s2,3) + Power(t1,3) - 3*t2 - 
                  Power(t1,2)*(11 + t2) + 
                  Power(s2,2)*(-3 + t1 + t2) + 
                  s2*(-5 - 9*t1 + 3*Power(t1,2) + 2*t2) + 
                  t1*(2 + 5*t2)) - 
               s*(4 + 2*Power(s2,2)*(-1 + t1) - 5*Power(t1,2) + 
                  Power(t1,3) - 4*t2 + 3*t1*t2 - Power(t1,2)*t2 - 
                  Power(t2,2) + Power(s1,2)*(-11 + s2 + 3*t1 + t2) + 
                  s2*(Power(t1,2) + 2*t1*(-2 + t2) - 2*(4 + t2)) + 
                  s1*(4 + 2*Power(s2,2) - 15*t1 + 4*Power(t1,2) + 
                     3*t2 + 2*s2*(-3 + t1 + t2))))/
             ((-1 + s1)*(-s + s1 - t2)) + 
            (2 - 7*s2 + 2*Power(s2,2) - 2*Power(s2,3) - 10*t1 - 
               2*s2*t1 + 9*Power(s2,2)*t1 - Power(s2,3)*t1 + 
               8*Power(t1,2) + 4*s2*Power(t1,2) - 
               2*Power(s2,2)*Power(t1,2) + s2*Power(t1,3) - 4*t2 + 
               7*s2*t2 - 2*Power(s2,2)*t2 + 2*t1*t2 + 2*s2*t1*t2 - 
               Power(s2,2)*t1*t2 + 2*Power(t1,2)*t2 - 
               s2*Power(t1,2)*t2 - s2*Power(t2,2) - 
               2*Power(s1,3)*(-3 + s2 + t2) + 
               Power(s1,2)*(-1 - 3*Power(s2,2) + Power(t1,2) + 
                  t1*(7 - 3*t2) + s2*(15 - 4*t1 - 3*t2) + t2) - 
               Power(s,2)*(4 - 7*t1 + Power(t1,2) + s2*(2 + t1) + 
                  2*t2 + t1*t2 + s1*(-7 + s2 + t1 + t2)) + 
               s1*(-15 - Power(s2,3) + Power(t1,3) - 
                  Power(t1,2)*(-2 + t2) + 5*t2 + Power(t2,2) + 
                  2*t1*(4 + t2) - Power(s2,2)*(-7 + 5*t1 + t2) + 
                  s2*(-3 - Power(t1,2) + t1*(18 - 4*t2) + t2)) + 
               s*(3 + 2*t1 - 4*Power(t1,2) - Power(t1,3) + 
                  2*Power(s2,2)*(2 + t1) - 4*t2 + Power(t1,2)*t2 + 
                  Power(s1,2)*(-13 + 3*s2 + t1 + 3*t2) + 
                  s2*(3*Power(t1,2) + 2*t1*(-8 + t2) + 3*(1 + t2)) + 
                  s1*(7 + 2*Power(s2,2) + 4*t1*(-4 + t2) + t2 + 
                     2*s2*(-7 + 3*t1 + t2))))/
             ((s - s2 + t1)*(s - s1 + t2)) + 
            (-3 + 3*s2 + 5*Power(s2,2) + 9*t1 + s2*t1 - 
               4*Power(s2,2)*t1 - Power(t1,2) + 2*s2*Power(t1,2) - 
               Power(t1,3) - 5*t2 + 4*s2*t2 + 2*Power(s2,2)*t2 + 
               4*t1*t2 + s2*t1*t2 - Power(s2,2)*t1*t2 + 
               3*Power(t1,2)*t2 - 2*Power(t2,2) + s2*Power(t2,2) - 
               4*t1*Power(t2,2) + 
               Power(s1,2)*(3 + 2*t2 - t1*(2 + t2)) - 
               Power(s,2)*(-1 + t1*(2 + t2)) + 
               s*(s1*(-1 + t1)*(3 + 2*t2) - 
                  2*(1 + t1 + Power(t1,2) + 2*t2 + t1*t2) + 
                  s2*(-7 - t2 + 2*t1*(3 + t2))) + 
               s1*(2 + 3*t2 + t1*(2 + 5*t2) + 
                  s2*(8 + 3*t2 - t1*(5 + 2*t2))))/((-1 + t1)*(-1 + t2)) \
+ (-3 - 4*s2 - 2*Power(s2,3) - t1 + 10*Power(s2,2)*t1 - 
               Power(s2,3)*t1 + Power(t1,2) - s2*Power(t1,2) - 
               3*Power(t1,3) - 4*t2 + 3*s2*t2 + 3*t1*t2 + 7*s2*t1*t2 - 
               3*Power(s2,2)*t1*t2 + Power(t1,2)*t2 + 
               s2*Power(t1,2)*t2 + 4*Power(t2,2) + s2*Power(t2,2) + 
               t1*Power(t2,2) - s2*t1*Power(t2,2) + Power(t2,3) - 
               Power(s1,3)*(-2 + s2 + t2) - 
               Power(s,2)*(1 - 3*t1 + s2*(2 + t1) + t2 + 2*t1*t2 + 
                  s1*(-4 + s2 + t2)) - 
               Power(s1,2)*(-1 + 2*Power(s2,2) + (-4 + t1)*t2 + 
                  Power(t2,2) + s2*(-9 + t1 + 3*t2)) + 
               s*(2 - 2*t1 + Power(t1,2) + 2*Power(s2,2)*(2 + t1) - 
                  t2 - 7*t1*t2 - Power(t1,2)*t2 + Power(t2,2) + 
                  t1*Power(t2,2) + 2*Power(s1,2)*(-3 + s2 + t2) + 
                  s2*(-13*t1 + 2*t2 + 5*t1*t2) + 
                  s1*(1 + 2*Power(s2,2) - 4*t1 - 3*t2 + 3*t1*t2 + 
                     Power(t2,2) + s2*(-9 + 2*t1 + 3*t2))) - 
               s1*(4 + Power(s2,3) - Power(t1,2)*(-5 + t2) - 3*t2 + 
                  Power(s2,2)*(-5 + 2*t1 + 2*t2) + 
                  t1*(-4 - 8*t2 + Power(t2,2)) + 
                  s2*((-4 + t2)*t2 + t1*(-11 + 4*t2))))/
             ((-1 + s2)*(-s + s2 - t1))))/
        (-Power(s*(-1 + s2) - (s2 - t1)*(-1 + s1 + t1 - t2),2) + 
          (-1 + t1)*(-1 + s1 + t1 - t2)*Power(1 - s2 + t1 - t2,3) - 
          Power(-1 + s2 - t1 + t2,2)*
           (3 + Power(s2,2) - s2*(4 - s + s2 - t1) - 
             2*(-1 + s - s2 + t1) + 2*Power(t1,2)*(-1 + s1 + t1 - t2) + 
             (2 + s - s2 + t1)*(-1 + s - s1 + t2) - 
             t1*(-1 + s1 + t1 - t2)*(-1 + s - s1 + s2 + t2)) + 
          (1 - s2 + t1 - t2)*(-2*Power(s,2)*(-1 + s2) + 
             s*(t1 - Power(t1,3) - 
                s1*(2 + t1 + Power(t1,2) - s2*(3 + t1)) + 2*t2 + 
                t1*t2 + Power(t1,2)*t2 + 
                s2*(-1 + Power(t1,2) - 3*t2 - t1*t2)) - 
             (s2 - t1)*(3 + Power(t1,3) + Power(s1,2)*(1 + t1) + 
                Power(t1,2)*(1 - 2*t2) + 4*t2 + Power(t2,2) + 
                2*s1*(-2 + t1 + Power(t1,2) - t2 - t1*t2) + 
                t1*(-5 - 2*t2 + Power(t2,2))))) - 
       (8*(-1 + Power(s1,5)*t1*(s2*(-1 + t1) + t1*(-5 - 2*s + 7*t1)) + 
            5*t2 + 2*s*t2 - 2*s2*t2 + 2*s*s2*t2 - 5*Power(t2,2) + 
            2*s*Power(t2,2) - 5*Power(s,2)*Power(t2,2) + 
            18*s2*Power(t2,2) - 2*s*s2*Power(t2,2) + 
            3*Power(s2,2)*Power(t2,2) + Power(t2,3) - 
            14*s*Power(t2,3) - 3*Power(s,2)*Power(t2,3) + 
            18*s2*Power(t2,3) + 5*Power(s2,2)*Power(t2,3) + 
            4*Power(t2,4) - 6*s*Power(t2,4) + 6*s2*Power(t2,4) + 
            (-1 + s)*Power(t1,5)*(-1 + s + t2) + 
            Power(t1,4)*(-1 + Power(s,2)*(1 - 5*t2) + (4 - 3*s2)*t2 + 
               s*(4 + 3*s2 - 3*t2)*t2 + (5 + s2)*Power(t2,2)) + 
            Power(t1,3)*(-2 + (1 + s2)*t2 + 
               2*(-9 + Power(s2,2))*Power(t2,2) - 
               (9 + 2*s2)*Power(t2,3) + 
               Power(s,2)*(-1 - 6*t2 + 7*Power(t2,2)) + 
               s*(4 + t2 - (9 + 8*s2)*Power(t2,2) + 3*Power(t2,3))) - 
            t1*(-1 + (4 + s2)*t2 + 
               2*(-11 + 9*s2 + 4*Power(s2,2))*Power(t2,2) + 
               (-2 - 2*s2 + 5*Power(s2,2))*Power(t2,3) + 
               (5 - 3*s2)*Power(t2,4) + 2*Power(t2,5) + 
               Power(s,2)*t2*(-4 + 5*t2 + 10*Power(t2,2)) + 
               s*(2 - 2*t2 - 3*(2 + 5*s2)*Power(t2,2) + 
                  (9 - 13*s2)*Power(t2,3) + 7*Power(t2,4))) - 
            Power(t1,2)*(-2 - 5*(-1 + s2)*t2 + 
               (4 + s2 - 3*Power(s2,2))*Power(t2,2) + 
               2*(-9 + Power(s2,2))*Power(t2,3) - 
               (7 + s2)*Power(t2,4) + 
               Power(s,2)*(1 - 5*t2 - 15*Power(t2,2) + 
                  3*Power(t2,3)) + 
               s*t2*(10 - 14*t2 - 14*Power(t2,2) + Power(t2,3) + 
                  s2*(5 + 11*t2 - 5*Power(t2,2)))) + 
            Power(s1,4)*(Power(s2,2)*(-1 + t1 + 2*Power(t1,2)) + 
               s2*(3*t1 + 6*Power(t1,3) + t2 + 
                  Power(t1,2)*(-7 - 4*s + t2)) + 
               t1*(1 + 16*Power(t1,3) + 10*t2 + 3*s*t2 - 
                  Power(t1,2)*(27 + 11*s + 14*t2) + 
                  t1*(6 + 2*Power(s,2) - 12*t2 + s*(5 + t2)))) + 
            Power(s1,3)*(12*Power(t1,5) - 
               3*Power(t1,4)*(12 + 5*s + 7*t2) - t2*(1 + (5 + s)*t2) + 
               Power(s2,2)*(2 + 5*Power(t1,3) - 5*t1*t2 - 
                  Power(t1,2)*(5 + t2)) + 
               Power(t1,3)*(6 + 6*Power(s,2) - 7*t2 + 9*Power(t2,2) + 
                  s*(16 + 11*t2)) - 
               t1*(3 + 14*t2 + 4*Power(s,2)*t2 - 3*Power(t2,2) + 
                  s*(-2 + 8*t2 + Power(t2,2))) + 
               Power(t1,2)*(9 + 49*t2 - 2*Power(s,2)*t2 + 
                  31*Power(t2,2) + s*(-5 + 27*t2 + 2*Power(t2,2))) + 
               s2*(2 + 7*Power(t1,4) - 3*t2 - 2*s*t2 - Power(t2,2) - 
                  Power(t1,3)*(7 + 11*s + 6*t2) + 
                  Power(t1,2)*
                   (-11 - 17*t2 - 3*Power(t2,2) + 3*s*(1 + t2)) + 
                  t1*(-9 + 10*t2 - 2*Power(t2,2) + s*(-2 + 9*t2)))) + 
            Power(s1,2)*(-1 + 3*Power(t1,6) + 3*t2 + 8*Power(t2,2) + 
               3*s*Power(t2,2) + Power(s,2)*Power(t2,2) + 
               2*Power(t2,3) - Power(t1,5)*(12 + 7*s + 8*t2) + 
               Power(t1,3)*(12 + Power(s,2)*(2 - 4*t2) + 69*t2 + 
                  s*(23 - 3*t2)*t2 + 36*Power(t2,2) - 2*Power(t2,3)) + 
               Power(t1,4)*(-1 + 5*Power(s,2) - 6*t2 + 
                  7*Power(t2,2) + s*(15 + 11*t2)) + 
               Power(s2,2)*(-1 + 3*Power(t1,4) - 2*t2 + 
                  5*Power(t2,2) - Power(t1,3)*(5 + 3*t2) + 
                  Power(t1,2)*(4 - 9*t2 - 2*Power(t2,2)) + 
                  t1*(-1 + 8*t2 + 3*Power(t2,2))) - 
               Power(t1,2)*(9 + 8*t2 + 28*Power(t2,2) + 
                  18*Power(t2,3) + 
                  Power(s,2)*(1 + 18*t2 + Power(t2,2)) + 
                  s*(-6 + 43*t2 + 31*Power(t2,2) + Power(t2,3))) + 
               t1*(8 - 14*t2 - 17*Power(t2,2) + 
                  4*Power(s,2)*Power(t2,2) - 20*Power(t2,3) - 
                  s*(4 - 7*t2 + 25*Power(t2,2) + 3*Power(t2,3))) + 
               s2*(-4 + 2*Power(t1,5) + (7 + 8*s)*t2 - 
                  (3 + 5*s)*Power(t2,2) + Power(t2,3) - 
                  2*Power(t1,4)*(1 + 4*s + 3*t2) + 
                  Power(t1,3)*
                   (-19 + s - 13*t2 + 7*s*t2 + 3*Power(t2,2)) + 
                  t1*(7 + 25*t2 + 13*Power(t2,2) + 4*Power(t2,3) + 
                     s*(2 - 8*t2 - 6*Power(t2,2))) + 
                  Power(t1,2)*
                   (16 + 41*t2 + 23*Power(t2,2) + Power(t2,3) + 
                     s*(-1 + 31*t2 + 3*Power(t2,2))))) + 
            s1*(2 - (-1 + s)*Power(t1,6) - 7*t2 - 2*s*t2 + 
               5*Power(t2,2) - 2*s*Power(t2,2) + 
               4*Power(s,2)*Power(t2,2) - 5*Power(t2,3) + 
               9*s*Power(t2,3) - 3*Power(s,2)*Power(t2,3) + 
               3*Power(t2,4) + s*Power(t2,4) + 
               Power(t1,5)*(-6 + Power(s,2) - 8*t2 + s*(2 + 3*t2)) + 
               Power(t1,4)*(1 - Power(s,2)*(-2 + t2) + 31*t2 + 
                  17*Power(t2,2) + s*(-1 + 10*t2 - 3*Power(t2,2))) + 
               Power(s2,2)*t2*
                (2 - 7*Power(t1,3) - Power(t1,4) - 7*t2 - 
                  4*Power(t2,2) + t1*(-5 + 3*t2 + Power(t2,2)) + 
                  Power(t1,2)*(11 + 10*t2 + Power(t2,2))) - 
               Power(t1,3)*(-9 - 4*t2 + 24*Power(t2,2) + 
                  14*Power(t2,3) + 
                  Power(s,2)*(1 + 14*t2 + Power(t2,2)) + 
                  s*(2 + 29*t2 + 23*Power(t2,2) - Power(t2,3))) + 
               Power(t1,2)*(Power(s,2)*t2*(-3 + 10*t2 + Power(t2,2)) + 
                  t2*(-36 - 35*t2 - 10*Power(t2,2) + 4*Power(t2,3)) + 
                  s*(-2 + 6*t2 + Power(t2,2) + 10*Power(t2,3))) + 
               t1*(-7 + 16*t2 + Power(t2,2) + 15*Power(t2,3) + 
                  9*Power(t2,4) + Power(s,2)*Power(t2,2)*(19 + 2*t2) + 
                  s*(4 - 8*t2 + 37*Power(t2,2) + 26*Power(t2,3) + 
                     Power(t2,4))) - 
               s2*(-2 + (3 + 8*s)*t2 + (18 - 5*s)*Power(t2,2) + 
                  (2 - 7*s)*Power(t2,3) + Power(t2,4) + 
                  Power(t1,5)*(-1 + s + t2) - 
                  Power(t1,4)*(1 + 2*(1 + s)*t2 + 2*Power(t2,2)) + 
                  Power(t1,2)*
                   (3 + (-29 + 4*s)*t2 + (-6 + 20*s)*Power(t2,2) + 
                     2*(3 + s)*Power(t2,3)) + 
                  Power(t1,3)*
                   (1 + t2 - 4*Power(t2,2) + Power(t2,3) - 
                     s*(1 + 19*t2 + Power(t2,2))) + 
                  t1*t2*(26 + 48*t2 + 23*Power(t2,2) + Power(t2,3) + 
                     s*(-3 + 32*t2 + 3*Power(t2,2)))))))/
        ((-1 + s1)*(-s + s1 - t2)*(-1 + s1 + t1 - t2)*
          Power(-(s1*t1) + t2,2)) - 
       (8*(6 + 2*Power(s1,6)*Power(t1,2) + 
            2*Power(-1 + s,2)*Power(t1,6) + 
            Power(s1,5)*t1*(-((15 + 4*s)*t1) + 7*Power(t1,2) + 
               s2*(3 + 7*t1) - 4*t2) - 5*t2 - 5*s*t2 + 11*s2*t2 - 
            12*s*s2*t2 + 9*Power(t2,2) - 45*s*Power(t2,2) + 
            13*Power(s,2)*Power(t2,2) - 11*s2*Power(t2,2) + 
            22*s*s2*Power(t2,2) + Power(s2,2)*Power(t2,2) - 
            21*Power(t2,3) + 17*s*Power(t2,3) + 
            Power(s,2)*Power(t2,3) + 23*s2*Power(t2,3) + 
            10*s*s2*Power(t2,3) - 11*Power(s2,2)*Power(t2,3) + 
            9*Power(t2,4) + 9*s*Power(t2,4) - 3*s2*Power(t2,4) + 
            2*Power(t2,5) - Power(t1,5)*
             (8 + t2 + 4*s2*t2 + Power(s,2)*(4 + 6*t2) - 
               s*(12 + (11 + 4*s2)*t2)) + 
            Power(t1,4)*(2 + 2*(6 + 5*s2)*t2 + 
               (1 + 3*s2 + 2*Power(s2,2))*Power(t2,2) + 
               Power(s,2)*(-8 + 6*t2 + 6*Power(t2,2)) - 
               s*(-8 + (27 + 10*s2)*t2 + (11 + 8*s2)*Power(t2,2))) + 
            Power(t1,3)*(16 + t2 + 9*s2*t2 - 
               (14 + 11*s2 + 2*Power(s2,2))*Power(t2,2) - 
               (3 + 2*Power(s2,2))*Power(t2,3) + 
               Power(s,2)*(4 + 35*t2 - 2*Power(t2,3)) + 
               s*(-24 - (65 + 7*s2)*t2 + 4*(7 + 2*s2)*Power(t2,2) + 
                  (5 + 4*s2)*Power(t2,3))) + 
            Power(t1,2)*(-10 - 3*(13 + 7*s2)*t2 - 
               (12 + 40*s2 + 5*Power(s2,2))*Power(t2,2) + 
               (4 + 11*s2)*Power(t2,3) + (1 + s2)*Power(t2,4) - 
               2*Power(s,2)*
                (-3 + 5*t2 + 19*Power(t2,2) + Power(t2,3)) + 
               s*(-4 + (64 + 22*s2)*t2 + (109 + 35*s2)*Power(t2,2) + 
                  (-17 + 2*s2)*Power(t2,3) - Power(t2,4))) + 
            t1*(-8 + (32 - 5*s2)*t2 + 
               (16 + 59*s2 + 4*Power(s2,2))*Power(t2,2) + 
               (24 + 14*s2 + 9*Power(s2,2))*Power(t2,3) - 
               4*s2*Power(t2,4) + 
               Power(s,2)*t2*(-17 + t2 + 11*Power(t2,2)) + 
               s*(12 + (22 + 3*s2)*t2 - 47*(1 + s2)*Power(t2,2) - 
                  5*(13 + 4*s2)*Power(t2,3) + 4*Power(t2,4))) + 
            Power(s1,4)*(8*Power(t1,4) + 
               Power(s2,2)*(2 + t1 + 5*Power(t1,2)) + 2*Power(t2,2) - 
               2*Power(t1,3)*(19 + 8*s + t2) + 
               t1*(-3 + 30*t2 + 11*s*t2) + 
               Power(t1,2)*(1 + 18*s + 2*Power(s,2) + 5*s*t2 - 
                  2*Power(t2,2)) + 
               s2*(17*Power(t1,3) - 3*t2 - 
                  Power(t1,2)*(16 + 7*s + 9*t2) - t1*(13 + 3*s + 18*t2)\
)) + Power(s1,3)*(4*Power(t1,5) - Power(t1,4)*(22 + 21*s + 3*t2) + 
               t2*(3 - (15 + 7*s)*t2) + 
               Power(s2,2)*(-8 + 11*Power(t1,3) + 
                  Power(t1,2)*(1 - 10*t2) - 3*t2 - 11*t1*t2) - 
               t1*(-13 + s - 2*t2 + 48*s*t2 + 7*Power(s,2)*t2 + 
                  21*Power(t2,2) + 11*s*Power(t2,2) - 4*Power(t2,3)) + 
               Power(t1,3)*(-10 + 9*Power(s,2) - t2 - Power(t2,2) + 
                  s*(47 + 19*t2)) - 
               Power(t1,2)*(-11 - 85*t2 + Power(t2,2) + 
                  Power(s,2)*(1 + 5*t2) + 
                  s*(-3 - 5*t2 + 2*Power(t2,2))) + 
               s2*(-4 + 11*Power(t1,4) + 13*t2 + 7*s*t2 + 
                  11*Power(t2,2) - 2*Power(t1,3)*(20 + 9*s + 7*t2) + 
                  3*Power(t1,2)*
                   (-9 - 2*s - 3*t2 + 5*s*t2 + Power(t2,2)) + 
                  2*t1*(6 + 4*s + 26*t2 + 9*s*t2 + 8*Power(t2,2)))) + 
            Power(s1,2)*(2 + Power(t1,6) - 13*t2 - 3*s*t2 - 
               3*Power(t2,2) + 30*s*Power(t2,2) + 
               7*Power(s,2)*Power(t2,2) + 14*Power(t2,3) + 
               6*s*Power(t2,3) - 2*Power(t2,4) - 
               Power(t1,5)*(1 + 12*s + 2*t2) + 
               Power(t1,4)*(-32 + 13*Power(s,2) - 12*t2 + 
                  Power(t2,2) + s*(30 + 19*t2)) + 
               Power(t1,3)*(8 + 87*t2 + 7*Power(t2,2) - 
                  Power(s,2)*(7 + 18*t2) + 
                  s*(21 - 9*t2 - 8*Power(t2,2))) + 
               Power(s2,2)*(6*Power(t1,4) - Power(t1,3)*(1 + 17*t2) + 
                  2*(3 + 5*t2 + Power(t2,2)) + 
                  3*t1*(1 + 2*t2 + 5*Power(t2,2)) + 
                  Power(t1,2)*(-14 - 11*t2 + 7*Power(t2,2))) + 
               Power(t1,2)*(35 - 6*t2 - 17*Power(t2,2) + 
                  4*Power(t2,3) + 
                  Power(s,2)*(-4 - 15*t2 + 5*Power(t2,2)) + 
                  s*(-17 - 116*t2 - 18*Power(t2,2) + Power(t2,3))) + 
               t1*(-13 - 42*t2 - 54*Power(t2,2) + 8*Power(t2,3) + 
                  Power(s,2)*t2*(9 + 10*t2) + 
                  s*(12 - 7*t2 + 50*Power(t2,2) + Power(t2,3))) + 
               s2*(16 + Power(t1,5) - 4*(2 + 7*s)*t2 - 
                  (36 + 11*s)*Power(t2,2) - 7*Power(t2,3) - 
                  Power(t1,4)*(13 + 13*s + 2*t2) + 
                  Power(t1,3)*
                   (-2 + 5*t2 + 4*Power(t2,2) + s*(-2 + 33*t2)) + 
                  Power(t1,2)*
                   (13 + 116*t2 + 6*Power(t2,2) - Power(t2,3) + 
                     s*(30 + 20*t2 - 12*Power(t2,2))) - 
                  t1*(15 - 33*t2 + 28*Power(t2,2) + 2*Power(t2,3) + 
                     s*(5 - 3*t2 + 27*Power(t2,2))))) - 
            s1*(8 + (1 + 3*s)*Power(t1,6) - 11*t2 - 8*s*t2 - 
               31*Power(t2,2) - 4*s*Power(t2,2) + 
               20*Power(s,2)*Power(t2,2) - 7*Power(t2,3) + 
               39*s*Power(t2,3) + 3*Power(s,2)*Power(t2,3) + 
               5*Power(t2,4) - s*Power(t2,4) + 
               Power(t1,5)*(2 + s - 8*Power(s,2) + 2*t2 - 7*s*t2) + 
               Power(t1,4)*(5 - 13*t2 - 5*Power(t2,2) + 
                  9*Power(s,2)*(1 + 2*t2) + 
                  s*(-49 - 4*t2 + 5*Power(t2,2))) + 
               Power(s2,2)*t2*
                (11 + 6*Power(t1,4) - 7*t2 + Power(t2,2) - 
                  Power(t1,3)*(1 + 10*t2) + 
                  Power(t1,2)*(-13 - 4*t2 + 2*Power(t2,2)) + 
                  t1*(-3 + 9*t2 + 5*Power(t2,2))) + 
               Power(t1,3)*(-19 - 46*t2 - 4*Power(t2,2) + 
                  2*Power(t2,3) + 
                  Power(s,2)*(4 + 5*t2 - 12*Power(t2,2)) + 
                  s*(8 + 131*t2 + 6*Power(t2,2) - Power(t2,3))) + 
               Power(t1,2)*(18 + 20*t2 + 89*Power(t2,2) + 
                  4*Power(t2,3) + 
                  s*(14 + 14*t2 - 95*Power(t2,2) - 4*Power(t2,3)) + 
                  Power(s,2)*
                   (3 - 24*t2 - 21*Power(t2,2) + 2*Power(t2,3))) + 
               t1*(-15 + 48*t2 - 37*Power(t2,2) - Power(t2,3) + 
                  5*Power(t2,4) + 
                  Power(s,2)*t2*(-15 + 5*t2 + 7*Power(t2,2)) + 
                  s*(23 - 58*t2 - 60*Power(t2,2) + 10*Power(t2,3) + 
                     Power(t2,4))) + 
               s2*(12 + Power(t1,5)*(-6 + 2*s - t2) + (5 - 33*s)*t2 - 
                  (2 + 7*s)*Power(t2,2) - 4*(5 + s)*Power(t2,3) + 
                  Power(t2,4) + 
                  Power(t1,4)*(13 + s + 10*t2 - 20*s*t2) - 
                  t1*(5 - 30*(1 + 2*s)*t2 + 
                     (-107 + 28*s)*Power(t2,2) + 
                     (11 + 12*s)*Power(t2,3) + Power(t2,4)) + 
                  Power(t1,3)*
                   (11 - 61*t2 + 5*Power(t2,2) + Power(t2,3) + 
                     s*(-2 + 6*t2 + 22*Power(t2,2))) + 
                  Power(t1,2)*
                   (-25 + 17*t2 + 34*Power(t2,2) + 
                     s*(-1 + 7*t2 + 21*Power(t2,2) - 4*Power(t2,3)))))))/
        ((-1 + s2)*(-s + s2 - t1)*(-1 + s1 + t1 - t2)*
          Power(-(s1*t1) + t2,2)) - 
       (8*(6 + 4*Power(s1,6)*Power(t1,2) + 
            2*Power(-1 + s,2)*Power(t1,6) - 7*t2 - 5*s*t2 + 11*s2*t2 - 
            12*s*s2*t2 - 7*Power(t2,2) - 37*s*Power(t2,2) + 
            13*Power(s,2)*Power(t2,2) - 11*s2*Power(t2,2) + 
            22*s*s2*Power(t2,2) + Power(s2,2)*Power(t2,2) - 
            37*Power(t2,3) + 27*s*Power(t2,3) + 
            Power(s,2)*Power(t2,3) + 17*s2*Power(t2,3) + 
            10*s*s2*Power(t2,3) - 11*Power(s2,2)*Power(t2,3) + 
            5*Power(t2,4) + 11*s*Power(t2,4) - 9*s2*Power(t2,4) + 
            2*Power(s1,5)*t1*
             (7*Power(t1,2) + s2*(2 + 5*t1) - 4*t2 - 
               t1*(13 + 3*s + 2*t2)) - 
            2*Power(t1,5)*(4 + t2 + 2*s2*t2 + Power(s,2)*(2 + 3*t2) - 
               2*s*(3 + (3 + s2)*t2)) + 
            2*Power(t1,4)*(1 + (8 + 5*s2)*t2 + 
               Power(1 + s2,2)*Power(t2,2) + 
               Power(s,2)*(-4 + 3*t2 + 3*Power(t2,2)) - 
               s*(-4 + 5*(3 + s2)*t2 + (7 + 4*s2)*Power(t2,2))) - 
            Power(t1,3)*(-16 + (3 - 9*s2)*t2 + 
               (39 + 10*s2 + 2*Power(s2,2))*Power(t2,2) + 
               2*(1 + s2 + Power(s2,2))*Power(t2,3) + 
               Power(s,2)*(-4 - 35*t2 + 2*Power(t2,3)) + 
               s*(24 + (64 + 7*s2)*t2 - (31 + 8*s2)*Power(t2,2) - 
                  4*(2 + s2)*Power(t2,3))) + 
            Power(t1,2)*(-10 - (41 + 21*s2)*t2 + 
               (5 - 43*s2 - 5*Power(s2,2))*Power(t2,2) + 
               7*(6 + s2)*Power(t2,3) + 2*s2*Power(t2,4) - 
               2*Power(s,2)*
                (-3 + 5*t2 + 19*Power(t2,2) + Power(t2,3)) + 
               s*(-4 + (67 + 22*s2)*t2 + (122 + 35*s2)*Power(t2,2) + 
                  2*(-7 + s2)*Power(t2,3) - 2*Power(t2,4))) + 
            t1*(-8 + (37 - 5*s2)*t2 + 
               (39 + 60*s2 + 4*Power(s2,2))*Power(t2,2) + 
               3*(5 + 8*s2 + 3*Power(s2,2))*Power(t2,3) - 
               (15 + s2)*Power(t2,4) + 
               Power(s,2)*t2*(-17 + t2 + 11*Power(t2,2)) + 
               s*(12 + (20 + 3*s2)*t2 - (68 + 47*s2)*Power(t2,2) - 
                  (81 + 20*s2)*Power(t2,3) + Power(t2,4))) + 
            Power(s1,4)*(15*Power(t1,4) + 
               Power(s2,2)*(2 + t1 + 5*Power(t1,2)) + 4*Power(t2,2) - 
               Power(t1,3)*(68 + 25*s + 12*t2) + 
               Power(t1,2)*(15 + 27*s + 2*Power(s,2) + 2*t2 + 
                  10*s*t2) + 
               4*t1*(-1 + 13*t2 + 4*s*t2 + 2*Power(t2,2)) + 
               s2*(28*Power(t1,3) - 4*t2 - 
                  Power(t1,2)*(23 + 7*s + 16*t2) - 
                  t1*(17 + 3*s + 26*t2))) + 
            Power(s1,3)*(5*Power(t1,5) - 
               Power(t1,4)*(58 + 34*s + 7*t2) + 
               Power(s2,2)*(-8 + 11*Power(t1,3) + 
                  Power(t1,2)*(1 - 10*t2) - 3*t2 - 11*t1*t2) - 
               2*t2*(-2 + (13 + 5*s)*t2 + 2*Power(t2,2)) - 
               t1*(-17 + s + 24*t2 + 71*s*t2 + 7*Power(s,2)*t2 + 
                  46*Power(t2,2) + 22*s*Power(t2,2)) + 
               Power(t1,3)*(3 + 9*Power(s,2) + 25*t2 + 
                  2*Power(t2,2) + s*(70 + 36*t2)) - 
               Power(t1,2)*(-15 - 138*t2 - 16*Power(t2,2) + 
                  Power(s,2)*(1 + 5*t2) + 
                  s*(7 - 13*t2 + 6*Power(t2,2))) + 
               s2*(-4 + 24*Power(t1,4) + 17*t2 + 7*s*t2 + 
                  16*Power(t2,2) - Power(t1,3)*(55 + 18*s + 33*t2) + 
                  Power(t1,2)*
                   (-28 - 6*s - 15*t2 + 15*s*t2 + 8*Power(t2,2)) + 
                  t1*(17 + 8*s + 75*t2 + 18*s*t2 + 30*Power(t2,2)))) + 
            Power(s1,2)*(2 - (17 + 19*s)*Power(t1,5) - 17*t2 - 
               3*s*t2 + 9*Power(t2,2) + 44*s*Power(t2,2) + 
               7*Power(s,2)*Power(t2,2) + 30*Power(t2,3) + 
               12*s*Power(t2,3) + 
               Power(t1,4)*(-25 + 13*Power(s,2) + 18*t2 + 
                  s*(47 + 34*t2)) + 
               Power(t1,3)*(27 + 158*t2 + Power(t2,2) - 
                  Power(s,2)*(7 + 18*t2) + 
                  s*(10 - 5*t2 - 17*Power(t2,2))) + 
               Power(s2,2)*(6*Power(t1,4) - Power(t1,3)*(1 + 17*t2) + 
                  2*(3 + 5*t2 + Power(t2,2)) + 
                  3*t1*(1 + 2*t2 + 5*Power(t2,2)) + 
                  Power(t1,2)*(-14 - 11*t2 + 7*Power(t2,2))) + 
               Power(t1,2)*(31 - 46*t2 - 96*Power(t2,2) - 
                  4*Power(t2,3) + 
                  Power(s,2)*(-4 - 15*t2 + 5*Power(t2,2)) + 
                  2*s*(-8 - 82*t2 - 23*Power(t2,2) + Power(t2,3))) + 
               t1*(-18 - 59*t2 - 70*Power(t2,2) + 4*Power(t2,3) + 
                  Power(s,2)*t2*(9 + 10*t2) + 
                  s*(12 + 22*t2 + 67*Power(t2,2) + 8*Power(t2,3))) + 
               s2*(16 + 6*Power(t1,5) - (13 + 28*s)*t2 - 
                  (52 + 11*s)*Power(t2,2) - 14*Power(t2,3) - 
                  Power(t1,4)*(20 + 13*s + 15*t2) + 
                  Power(t1,3)*
                   (-3 - t2 + 13*Power(t2,2) + s*(-2 + 33*t2)) - 
                  2*Power(t1,2)*
                   (-9 - 72*t2 - 15*Power(t2,2) + Power(t2,3) + 
                     s*(-15 - 10*t2 + 6*Power(t2,2))) - 
                  t1*(17 - 23*t2 + 51*Power(t2,2) + 10*Power(t2,3) + 
                     s*(5 - 3*t2 + 27*Power(t2,2))))) - 
            s1*(8 + 4*s*Power(t1,6) - 16*t2 - 8*s*t2 - 44*Power(t2,2) + 
               15*s*Power(t2,2) + 20*Power(s,2)*Power(t2,2) + 
               55*s*Power(t2,3) + 3*Power(s,2)*Power(t2,3) + 
               8*Power(t2,4) + 2*s*Power(t2,4) - 
               2*Power(t1,5)*(-3 + s + 4*Power(s,2) - t2 + 5*s*t2) + 
               Power(t1,4)*(1 - 54*t2 - 2*Power(t2,2) + 
                  9*Power(s,2)*(1 + 2*t2) + 8*s*(-6 - t2 + Power(t2,2))\
) + Power(s2,2)*t2*(11 + 6*Power(t1,4) - 7*t2 + Power(t2,2) - 
                  Power(t1,3)*(1 + 10*t2) + 
                  Power(t1,2)*(-13 - 4*t2 + 2*Power(t2,2)) + 
                  t1*(-3 + 9*t2 + 5*Power(t2,2))) + 
               Power(t1,3)*(-21 - 22*t2 + 65*Power(t2,2) + 
                  Power(s,2)*(4 + 5*t2 - 12*Power(t2,2)) + 
                  s*(11 + 161*t2 + 24*Power(t2,2) - 2*Power(t2,3))) + 
               Power(t1,2)*(23 + 62*t2 + 115*Power(t2,2) - 
                  21*Power(t2,3) + 
                  Power(s,2)*
                   (3 - 24*t2 - 21*Power(t2,2) + 2*Power(t2,3)) - 
                  2*s*(-6 + 9*t2 + 60*Power(t2,2) + 8*Power(t2,3))) + 
               t1*(-17 + 28*t2 - 80*Power(t2,2) - 51*Power(t2,3) - 
                  2*Power(t2,4) + 
                  Power(s,2)*t2*(-15 + 5*t2 + 7*Power(t2,2)) + 
                  s*(23 - 49*t2 - 75*Power(t2,2) + Power(t2,3) + 
                     2*Power(t2,4))) + 
               s2*(12 + 2*(-3 + s)*Power(t1,5) + (3 - 33*s)*t2 - 
                  (13 + 7*s)*Power(t2,2) - 2*(19 + 2*s)*Power(t2,3) - 
                  2*Power(t2,4) + 
                  Power(t1,4)*
                   (13 + s + 16*t2 - 20*s*t2 - 2*Power(t2,2)) - 
                  t1*(5 - 5*(7 + 12*s)*t2 + 
                     2*(-57 + 14*s)*Power(t2,2) + 
                     12*(1 + s)*Power(t2,3) + 2*Power(t2,4)) + 
                  Power(t1,3)*
                   (11 - 71*t2 - 12*Power(t2,2) + 2*Power(t2,3) + 
                     s*(-2 + 6*t2 + 22*Power(t2,2))) + 
                  Power(t1,2)*
                   (-25 + 17*t2 + 51*Power(t2,2) + 12*Power(t2,3) + 
                     s*(-1 + 7*t2 + 21*Power(t2,2) - 4*Power(t2,3)))))))/
        ((s - s2 + t1)*(-1 + s1 + t1 - t2)*(s - s1 + t2)*
          Power(-(s1*t1) + t2,2)) - 
       (8*(-1 - Power(-1 + s,2)*Power(t1,6) + 10*t2 + 2*s*t2 - 
            4*s2*t2 + 2*s*s2*t2 + 12*Power(t2,2) - 20*s*Power(t2,2) - 
            5*Power(s,2)*Power(t2,2) + 30*s2*Power(t2,2) + 
            2*s*s2*Power(t2,2) - 5*Power(s2,2)*Power(t2,2) + 
            14*Power(t2,3) - 48*s*Power(t2,3) - 
            6*Power(s,2)*Power(t2,3) + 42*s2*Power(t2,3) + 
            10*s*s2*Power(t2,3) - 6*Power(s2,2)*Power(t2,3) + 
            7*Power(t2,4) - 18*s*Power(t2,4) - 
            3*Power(s,2)*Power(t2,4) + 16*s2*Power(t2,4) + 
            6*s*s2*Power(t2,4) - 3*Power(s2,2)*Power(t2,4) + 
            2*Power(t2,5) + 2*Power(s1,5)*Power(t1,2)*(4 + t2) + 
            Power(t1,5)*(2 + t2 + 3*Power(s,2)*t2 + 2*s2*t2 - 
               2*s*(1 + (3 + s2)*t2)) - 
            Power(t1,4)*(-1 + t2 + 3*s2*t2 + 
               (2 + 2*s2 + Power(s2,2))*Power(t2,2) + 
               Power(s,2)*(-2 - 2*t2 + 3*Power(t2,2)) + 
               s*(4 + t2 - 3*s2*t2 - (7 + 4*s2)*Power(t2,2))) + 
            Power(t1,3)*(-4 - 5*(-1 + s2)*t2 + 
               (4 + s2 + Power(s2,2))*Power(t2,2) + 
               (3 + s2 + Power(s2,2))*Power(t2,3) + 
               Power(s,2)*t2*(-9 - 5*t2 + Power(t2,2)) + 
               2*s*(2 + (5 + 3*s2)*t2 + Power(t2,2) - 
                  (2 + s2)*Power(t2,3))) + 
            Power(t1,2)*(1 + (-1 + 7*s2)*t2 + 
               2*(-9 + 7*s2 + 2*Power(s2,2))*Power(t2,2) + 
               (-7 - 6*s2 + Power(s2,2))*Power(t2,3) - 
               (1 + s2)*Power(t2,4) + 
               Power(s,2)*(-1 + 11*Power(t2,2) + 4*Power(t2,3)) + 
               s*(2 - (9 + 5*s2)*t2 - (10 + 21*s2)*Power(t2,2) + 
                  (4 - 5*s2)*Power(t2,3) + Power(t2,4))) + 
            t1*(2 + (-14 + 3*s2)*t2 + 
               (4 - 43*s2 + Power(s2,2))*Power(t2,2) - 
               s2*(11 + 4*s2)*Power(t2,3) + 
               (2 + 5*s2 - Power(s2,2))*Power(t2,4) - 
               Power(s,2)*t2*(-4 - 4*t2 + Power(t2,2) + Power(t2,3)) + 
               s*(-2 - 4*(-1 + s2)*t2 + (39 + 9*s2)*Power(t2,2) + 
                  (24 + 7*s2)*Power(t2,3) + (-3 + 2*s2)*Power(t2,4))) + 
            Power(s1,4)*t1*(-Power(t1,3) - 4*t2*(4 + t2) + 
               Power(t1,2)*(14 + 3*t2) - 
               t1*(27 + 10*s + 23*t2 + 4*s*t2 + 2*Power(t2,2)) + 
               s2*(5 + Power(t1,2) + t2 + 5*t1*(2 + t2))) + 
            Power(s1,3)*(-2*Power(t1,5) + 2*Power(t2,2)*(4 + t2) + 
               Power(t1,4)*(5 + 2*s + 3*t2) + 
               t1*(-5 + (53 + 25*s)*t2 + (46 + 9*s)*Power(t2,2) + 
                  4*Power(t2,3)) + 
               Power(s2,2)*(1 + t1*(1 + t2) + 
                  3*Power(t1,2)*(2 + t2)) - 
               Power(t1,3)*(26 + 21*t2 + Power(t2,2) + 7*s*(3 + t2)) + 
               Power(t1,2)*(18 - t2 + 5*Power(t2,2) + 
                  2*Power(s,2)*(1 + t2) + 5*s*(5 + 6*t2 + Power(t2,2))\
) + s2*(2*Power(t1,4) - t2*(5 + t2) + Power(t1,3)*(18 - s + 5*t2) - 
                  t1*(14 + 32*t2 + 11*Power(t2,2) + s*(2 + t2)) - 
                  Power(t1,2)*
                   (32 + 36*t2 + 6*Power(t2,2) + s*(3 + 5*t2)))) + 
            Power(s1,2)*(-Power(t1,6) + 
               Power(t1,5)*(-3 + 4*s + 2*t2) - 
               t2*(-5 + (26 + 15*s)*t2 + (23 + 5*s)*Power(t2,2) + 
                  2*Power(t2,3)) - 
               t1*(-14 + 6*(4 + 11*s + Power(s,2))*t2 + 
                  (39 + 68*s + 5*Power(s,2))*Power(t2,2) + 
                  (19 + 10*s)*Power(t2,3)) - 
               Power(t1,4)*(-1 + Power(s,2) - 5*t2 + Power(t2,2) + 
                  4*s*(3 + 2*t2)) + 
               Power(t1,3)*(-3 - 21*t2 - 5*Power(t2,2) + 
                  Power(s,2)*(3 + 4*t2) + 
                  s*(29 + 32*t2 + 5*Power(t2,2))) + 
               Power(s2,2)*(-1 - 3*t2 - Power(t2,2) + 
                  4*Power(t1,3)*(1 + t2) - 
                  Power(t1,2)*(1 + 12*t2 + 4*Power(t2,2)) - 
                  t1*(2 + 13*t2 + 7*Power(t2,2))) - 
               Power(t1,2)*(8 - 63*t2 - 53*Power(t2,2) - 
                  4*Power(t2,3) + Power(s,2)*t2*(4 + 3*t2) + 
                  s*(3 - 6*t2 + 8*Power(t2,2) + Power(t2,3))) + 
               s2*(-2 + (11 - 2*s)*Power(t1,4) + Power(t1,5) + 
                  2*(7 + 2*s)*t2 + (22 + s)*Power(t2,2) + 
                  6*Power(t2,3) - 
                  Power(t1,3)*
                   (38 + s + 33*t2 + 7*s*t2 + 3*Power(t2,2)) + 
                  t1*(11 + 83*t2 + 72*Power(t2,2) + 11*Power(t2,3) + 
                     2*s*(1 + 5*t2 + 6*Power(t2,2))) + 
                  Power(t1,2)*
                   (17 + 14*t2 + 17*Power(t2,2) + Power(t2,3) + 
                     s*(-5 + 15*t2 + 7*Power(t2,2))))) + 
            s1*(1 + 2*s*Power(t1,6) - 14*t2 - 2*s*t2 + 6*Power(t2,2) + 
               41*s*Power(t2,2) + 5*Power(s,2)*Power(t2,2) + 
               26*Power(t2,3) + 38*s*Power(t2,3) + 
               3*Power(s,2)*Power(t2,3) + 11*Power(t2,4) + 
               5*s*Power(t2,4) + 
               Power(t1,5)*(-1 - 2*Power(s,2) - 5*s*(-1 + t2) + 3*t2) + 
               Power(t1,4)*(-6 - t2 - 5*Power(t2,2) + 
                  Power(s,2)*(3 + 5*t2) + s*(-2 - 8*t2 + 4*Power(t2,2))\
) + Power(s2,2)*t2*(4 + Power(t1,4) + 6*t2 + 4*Power(t2,2) - 
                  3*Power(t1,3)*(1 + t2) + 
                  Power(t1,2)*(-9 - t2 + Power(t2,2)) + 
                  t1*(7 + 22*t2 + 7*Power(t2,2))) + 
               Power(t1,3)*(5 + 17*t2 + 4*Power(t2,2) + 
                  2*Power(t2,3) - 4*Power(s,2)*t2*(1 + t2) + 
                  s*(1 + 18*t2 + 4*Power(t2,2) - Power(t2,3))) + 
               Power(t1,2)*(13 - t2 + 16*Power(t2,2) + 
                  Power(s,2)*
                   (-1 - 7*t2 - 5*Power(t2,2) + Power(t2,3)) - 
                  2*s*(4 + 32*t2 + 29*Power(t2,2) + Power(t2,3))) + 
               t1*(-12 - 4*t2 - 51*Power(t2,2) - 38*Power(t2,3) - 
                  5*Power(t2,4) + 
                  2*Power(s,2)*t2*(1 + 5*t2 + 3*Power(t2,2)) + 
                  s*(2 + 25*t2 + 63*Power(t2,2) + 33*Power(t2,3) + 
                     Power(t2,4))) - 
               s2*(-2 + (1 + s)*Power(t1,5) + (7 + 6*s)*t2 + 
                  (51 + 9*s)*Power(t2,2) + (37 + 7*s)*Power(t2,3) + 
                  5*Power(t2,4) - 
                  Power(t1,3)*
                   (1 + s - 23*t2 + 15*s*t2 + 6*Power(t2,2) + 
                     7*s*Power(t2,2) + Power(t2,3)) + 
                  Power(t1,4)*
                   (-3 - 2*t2 + Power(t2,2) + s*(2 + 4*t2)) + 
                  t1*t2*(45 + 74*t2 + 38*Power(t2,2) + Power(t2,3) + 
                     s*(-7 + 28*t2 + 13*Power(t2,2))) + 
                  Power(t1,2)*
                   (5 - 73*t2 - 42*Power(t2,2) + 2*Power(t2,3) + 
                     2*s*(-1 - 2*Power(t2,2) + Power(t2,3)))))))/
        ((-1 + s2)*(-1 + t1)*(-1 + s1 + t1 - t2)*Power(-(s1*t1) + t2,2)) \
+ (4*((-2*(3 - (-1 + s)*Power(t1,4)*(-2 + t2) - s*t2 - s2*t2 - 
                 11*Power(t2,2) - 3*s*Power(t2,2) + 5*s2*Power(t2,2) - 
                 12*Power(t2,3) + 4*s*Power(t2,3) + 4*Power(t2,4) + 
                 Power(s1,3)*
                  (s2*(-1 + t1) + t1*(-5 + t1 - 2*t2 + 2*t1*t2)) + 
                 Power(t1,3)*
                  (-8 + 2*(5 + s2)*t2 - (1 + s2)*Power(t2,2) + 
                    s*(10 - 9*t2 + Power(t2,2))) + 
                 t1*(8 + (19 + s2)*t2 + (-7 + s2)*Power(t2,2) + 
                    2*Power(t2,3) + s*(3 - 12*t2 + 7*Power(t2,2))) + 
                 Power(t1,2)*
                  (-1 + 4*(5 + s2)*t2 - (3 + 7*s2)*Power(t2,2) + 
                    2*Power(t2,3) + s*(11 - 17*t2 + 9*Power(t2,2))) + 
                 Power(s1,2)*
                  (1 + 3*Power(t1,3)*(-1 + t2) + 5*t2 - s*t2 + 
                    2*Power(t2,2) + t1*(-5 + s - 6*t2 - s*t2) - 
                    Power(t1,2)*
                     (11 + s + 8*t2 + 2*s*t2 + 2*Power(t2,2)) + 
                    s2*(-3 + Power(t1,2)*(5 + t2) - t1*(4 + 3*t2))) + 
                 s1*(3 - 28*Power(t1,2) - 21*Power(t1,3) - 
                    4*Power(t1,4) + 4*t2 + 20*t1*t2 + 
                    11*Power(t1,2)*t2 + 4*Power(t1,3)*t2 + 
                    Power(t1,4)*t2 + 5*Power(t2,2) + 
                    20*t1*Power(t2,2) - 5*Power(t1,2)*Power(t2,2) - 
                    2*Power(t1,3)*Power(t2,2) - 2*Power(t2,3) - 
                    2*t1*Power(t2,3) + 
                    s*(Power(t1,3)*(1 - 3*t2) + (-3 + t2)*t2 + 
                       t1*(1 + t2) + 
                       Power(t1,2)*(4 - 3*t2 + Power(t2,2))) + 
                    s2*(-3 + 4*t2 + 5*Power(t2,2) - 
                       Power(t1,2)*Power(t2,2) + 
                       Power(t1,3)*(4 + t2) - 
                       t1*(7 + t2 + 2*Power(t2,2))))))/
             ((-1 + t1)*(s1*t1 - t2)) + 
            (2*(-1 + 4*Power(s1,4)*t1 - (-1 + s)*Power(t1,4) - 
                 3*s*t2 + 5*s2*t2 + 13*Power(t2,2) - s*Power(t2,2) - 
                 9*s2*Power(t2,2) + 4*Power(t2,3) + 4*s*Power(t2,3) + 
                 Power(s1,3)*(1 + t1)*
                  (-(t1*(-6 + t2)) - 4*t2 + s2*(2 + t2)) + 
                 Power(t1,3)*
                  (10 - (8 + s2)*t2 - Power(t2,2) + 
                    s*(-11 + 6*t2 + Power(t2,2))) - 
                 t1*(10 - (-12 + s2)*t2 - 13*Power(t2,2) - 
                    (-1 + s2)*Power(t2,3) + 
                    s*(1 - 16*t2 + 12*Power(t2,2) + Power(t2,3))) + 
                 Power(s1,2)*
                  (-2 - Power(t1,3)*(-5 + t2) - 7*t2 + 2*s*t2 + 
                    5*Power(t2,2) + s*Power(t2,2) + 
                    Power(t1,2)*(20 + Power(t2,2) + s*(6 + t2)) + 
                    s2*(5 + Power(t1,2)*(-3 + t2) - 6*t2 - 
                       2*Power(t2,2) + t1*(14 - 5*t2 - 2*Power(t2,2))\
) + t1*(-3 - 10*t2 + 6*Power(t2,2) + s*(-6 + 3*t2 + Power(t2,2)))) - 
                 Power(t1,2)*
                  (s*(11 - 13*t2 + 6*Power(t2,2) + Power(t2,3)) + 
                    t2*(20 + 5*t2 + Power(t2,2) - 
                       s2*(-5 + 5*t2 + Power(t2,2)))) + 
                 s1*(-5 + 3*Power(t1,4) + 7*t2 + 9*s*t2 + 
                    5*Power(t2,2) - 6*s*Power(t2,2) - 5*Power(t2,3) - 
                    s*Power(t2,3) + 
                    Power(t1,3)*
                     (18 + 6*t2 + Power(t2,2) + s*(5 + t2)) + 
                    Power(t1,2)*
                     (22 - 21*t2 + Power(t2,2) + 
                       2*s*(-4 + 4*t2 + Power(t2,2))) + 
                    s2*(1 - 5*Power(t1,3) - 16*t2 + 4*Power(t2,2) + 
                       Power(t2,3) + 
                       Power(t1,2)*(-1 + 2*t2 - 2*Power(t2,2)) + 
                       t1*(5 + 2*t2 + 2*Power(t2,2) + Power(t2,3))) - 
                    t1*(-2 + 32*t2 + 3*Power(t2,2) + Power(t2,3) + 
                       s*(1 + 8*Power(t2,2) + Power(t2,3))))))/
             ((-1 + s1)*(s1*t1 - t2)) + 
            ((-1 + t2)*((2*(1 - Power(t1,2) + 
                      s1*(s2*(-1 + t1) - t1*(1 + t1 - 2*t2)) + 
                      s*(1 + t1)*(t1 - t2) + t2 - s2*t2 + t1*t2 + 
                      s2*t1*t2 - 2*Power(t2,2))*
                    (Power(s1,3) - 13*s2 - 8*t1 + 6*s2*t1 - 
                      9*Power(t1,2) + s2*Power(t1,2) + Power(t1,3) + 
                      Power(s1,2)*(-3 + 2*s2 + 2*t1) - 8*t2 + 
                      9*s2*t2 + 3*t1*t2 - 4*s2*t1*t2 + 
                      2*Power(t1,2)*t2 + 4*Power(t2,2) - 
                      3*t1*Power(t2,2) - 
                      s*(-5 + Power(s1,2) + 5*t1 + 
                       s1*(3 + t1 - 2*t2) + 5*t2 - 2*t1*t2) + 
                      s1*(-9 + 2*Power(t1,2) + 
                        s2*(4 + 3*t1 - 4*t2) + 2*t1*(-5 + t2) + 
                        8*t2 - 3*Power(t2,2))))/(s1*t1 - t2) - 
                 4*(5 + s2 - 4*Power(s2,2) - 8*s2*t1 + 
                    6*Power(s2,2)*t1 + 2*Power(t1,2) - 
                    12*s2*Power(t1,2) + Power(t1,3) + 
                    s2*Power(t1,3) - 4*t2 - 4*s2*t2 + 
                    3*Power(s2,2)*t2 - 3*t1*t2 + 12*s2*t1*t2 - 
                    2*Power(s2,2)*t1*t2 - 10*Power(t1,2)*t2 + 
                    2*s2*Power(t1,2)*t2 - Power(t2,2) + 
                    4*s2*Power(t2,2) + 4*t1*Power(t2,2) - 
                    3*s2*t1*Power(t2,2) + Power(t1,2)*Power(t2,2) - 
                    t1*Power(t2,3) + Power(s1,3)*(-2 + s2 + t2) + 
                    Power(s,2)*
                     (-3 + s1 + 2*t1 - s1*t1 - Power(t1,2) + 2*t2) + 
                    Power(s1,2)*
                     (-3 + Power(s2,2) + s2*(-4 + t1) + 
                       Power(t1,2) + 2*t1*(-2 + t2) + t2 - 
                       Power(t2,2)) + 
                    s*(3*t1 + 5*Power(t1,2) - Power(t1,3) + 
                       s2*(7 + Power(t1,2) + 2*t1*(-4 + t2) - 
                      5*t2) + t2 - 8*t1*t2 - Power(t1,2)*t2 - 
                       2*Power(t2,2) + 2*t1*Power(t2,2) - 
                       Power(s1,2)*(-3 + s2 - t1 + t2) + 
                       s1*(5 + 7*t1 + 2*s2*(-3 + t2) - 10*t2 - 
                        2*t1*t2 + 2*Power(t2,2))) + 
                    s1*(Power(t1,3) + Power(s2,2)*(5 + t1 - 2*t2) + 
                       Power(t1,2)*(-1 + t2) - t1*(3 + 8*t2) - 
                       t2*(3 - 5*t2 + Power(t2,2)) + 
                       s2*(-10 + Power(t1,2) + 14*t2 - 
                        3*Power(t2,2) + t1*(-15 + 2*t2)))) - 
                 (2*(-9 + s1 + Power(s1,2) + 2*s1*t1 + Power(t1,2) + 
                      4*t2 - 2*s1*t2 - 2*t1*t2)*
                    (-1 + Power(-1 + s,2)*Power(t1,5) + t2 + 
                      2*s*s2*t2 - 3*Power(t2,2) + 6*s*Power(t2,2) - 
                      Power(s,2)*Power(t2,2) + 6*s2*Power(t2,2) - 
                      6*s*s2*Power(t2,2) - Power(s2,2)*Power(t2,2) + 
                      5*Power(t2,3) - 4*s*Power(t2,3) + 
                      Power(s,2)*Power(t2,3) - 4*s2*Power(t2,3) - 
                      4*s*s2*Power(t2,3) + Power(s2,2)*Power(t2,3) - 
                      2*s*Power(t2,4) - 2*s2*Power(t2,4) - 
                      2*Power(t2,5) + 
                      Power(t1,4)*
                       (-1 + Power(s,2)*(1 - 3*t2) + 
                        2*s*(3 + s2)*t2 - (1 + 2*s2)*t2) + 
                      Power(s1,3)*
                       (Power(s2,2)*Power(-1 + t1,2) - 
                        2*s2*(-1 + t1)*t1*(t1 - t2) + 
                        Power(t1,2)*
                        (1 + Power(t1,2) - 2*t2 - 2*t1*t2 + 
                        2*Power(t2,2))) + 
                      Power(t1,3)*
                       (-2 + 3*Power(t2,2) + 
                        Power(s2,2)*Power(t2,2) + 2*s2*t2*(1 + t2) + 
                        Power(s,2)*(-1 - 4*t2 + 3*Power(t2,2)) - 
                        2*s*
                        (-2 + (-3 + s2)*t2 + 2*(2 + s2)*Power(t2,2))) \
+ t1*(1 - 2*(2 + s2)*t2 + (1 - 6*s2 + Power(s2,2))*Power(t2,2) + 
                        (-2 + 6*s2)*Power(t2,3) + 
                        2*(2 + s2)*Power(t2,4) + 
                        Power(s,2)*t2*(2 - 3*t2 - 2*Power(t2,2)) - 
                        2*s*
                        (1 + t2 - s2*t2 - (4 + 5*s2)*Power(t2,2) - 
                        (3 + s2)*Power(t2,3) + Power(t2,4))) - 
                      Power(t1,2)*
                       (-2 - 2*(2 + s2)*t2 + 
                        Power(1 + s2,2)*Power(t2,2) + 
                        (5 + 2*s2 + Power(s2,2))*Power(t2,3) + 
                        Power(s,2)*
                       (1 - 3*t2 - 5*Power(t2,2) + Power(t2,3)) - 
                        2*s*t2*
                        (-5 - 5*t2 + 3*Power(t2,2) + 
                        s2*(-2 + Power(t2,2)))) + 
                      s1*(1 - 2*s*Power(t1,5) + 
                        Power(t1,4)*
                       (1 + Power(s,2) + 6*s*(-1 + t2) - 4*t2) - 
                        Power(t2,2) - 2*s*Power(t2,2) + 
                        Power(s,2)*Power(t2,2) - 2*Power(t2,3) + 
                        2*s*Power(t2,3) + 2*Power(t2,4) + 
                        Power(s2,2)*(-1 + t1)*t2*
                       (t1*(-4 + t2) + 3*t2) - 
                        2*Power(t1,3)*
                       (1 + Power(s,2)*t2 - 4*Power(t2,2) + 
                       s*(-1 - 6*t2 + 3*Power(t2,2))) + 
                        Power(t1,2)*
                       (2 - 4*t2 + 7*Power(t2,2) - 8*Power(t2,3) + 
                       Power(s,2)*(1 + Power(t2,2)) + 
                       2*s*(1 - t2 - 5*Power(t2,2) + Power(t2,3))) - 
                        2*t1*
                       (1 - 4*t2 + Power(s,2)*t2 + 4*Power(t2,2) + 
                       Power(t2,3) - 2*Power(t2,4) - 
                       s*(2 - 4*t2 + Power(t2,2) + 2*Power(t2,3))) + 
                        2*s2*(-1 + t1)*
                        (-1 + s*Power(t1,3) + t2 + 2*s*t2 - 
                        s*Power(t2,2) + 
                        Power(t1,2)*(1 + s + t2 + Power(t2,2)) - 
                        t1*t2*(-2 + 4*t2 + Power(t2,2) + s*(2 + t2)))\
) + Power(s1,2)*(Power(s2,2)*(-1 + t1)*
                       (1 + Power(t1,2) - t2 - t1*t2) - 
                        2*s2*(-1 + t1)*
                        (-1 + t1 + Power(t1,3) + s*t2 - t1*t2 + 
                        Power(t2,2) - Power(t1,2)*(-1 + s + 2*t2)) + 
                        t1*(Power(t1,4) + Power(t1,3)*(1 - 3*t2) - 
                        4*(-1 + t2)*Power(t2,2) + 
                        Power(t1,2)*(3 - 6*t2 + 4*Power(t2,2)) + 
                        t1*(-5 + 3*t2 + 4*Power(t2,2) - 
                       2*Power(t2,3)) - 
                        2*s*(1 + Power(t1,3) + 
                        Power(t1,2)*(1 - 2*t2) - t2 + Power(t2,2) + 
                        t1*(-1 - t2 + Power(t2,2)))))))/
                  ((-1 + s1 + t1 - t2)*Power(-(s1*t1) + t2,2))))/
             ((-1 + s2)*(-s + s2 - t1)) + 
            ((-1 + t2)*((2*(1 - Power(t1,2) + 
                      s1*(s2*(-1 + t1) - t1*(1 + t1 - 2*t2)) + 
                      s*(1 + t1)*(t1 - t2) + t2 - s2*t2 + t1*t2 + 
                      s2*t1*t2 - 2*Power(t2,2))*
                    (2 + 2*Power(s1,3) - 13*s2 - 10*t1 + 6*s2*t1 - 
                      9*Power(t1,2) + s2*Power(t1,2) + Power(t1,3) + 
                      2*Power(s1,2)*(-1 + s2 + 2*t1 - t2) - 6*t2 + 
                      9*s2*t2 + 4*t1*t2 - 4*s2*t1*t2 + 
                      Power(t1,2)*t2 + 4*Power(t2,2) - 
                      2*t1*Power(t2,2) - 
                      s*(-5 + Power(s1,2) + 5*t1 + 
                       s1*(3 + t1 - 2*t2) + 5*t2 - 2*t1*t2) + 
                      s1*(-13 + 3*Power(t1,2) + 
                        s2*(4 + 3*t1 - 4*t2) + 7*t2 - 
                        2*Power(t2,2) - t1*(9 + t2))))/(s1*t1 - t2) + 
                 4*(-11 + 4*Power(s2,2) + 3*t1 + 8*s2*t1 - 
                    6*Power(s2,2)*t1 + 7*Power(t1,2) + 
                    9*s2*Power(t1,2) + Power(t1,3) - s2*Power(t1,3) + 
                    Power(s1,3)*(5 - 2*s2 + t1 - 2*t2) + 
                    Power(s,2)*
                     (3 + s1*(-1 + t1) - 2*t1 + Power(t1,2) - 2*t2) + 
                    10*t2 + 2*s2*t2 - 3*Power(s2,2)*t2 - 10*t1*t2 - 
                    12*s2*t1*t2 + 2*Power(s2,2)*t1*t2 + 
                    3*Power(t1,2)*t2 - s2*Power(t1,2)*t2 - 
                    Power(t1,3)*t2 - Power(t2,2) - s2*Power(t2,2) + 
                    3*t1*Power(t2,2) + 2*s2*t1*Power(t2,2) + 
                    Power(t1,2)*Power(t2,2) - 2*Power(t2,3) + 
                    Power(s1,2)*
                     (9 - Power(s2,2) + Power(t1,2) - 
                       6*t1*(-2 + t2) - 10*t2 + 3*Power(t2,2) + 
                       s2*(3 - 3*t1 + 2*t2)) + 
                    s*(5 - 5*t1 - 8*Power(t1,2) - 
                       s2*(7 + Power(t1,2) + 2*t1*(-4 + t2) - 
                      5*t2) - 5*t2 + 15*t1*t2 + 2*Power(t1,2)*t2 - 
                       Power(t2,2) - 2*t1*Power(t2,2) + 
                       Power(s1,2)*(-4 + s2 - 2*t1 + t2) - 
                       s1*(9 + 9*t1 + 2*Power(t1,2) + 
                        2*s2*(-3 + t2) - 15*t2 - 3*t1*t2 + 
                        2*Power(t2,2))) + 
                    s1*(-10 - Power(s2,2)*(5 + t1 - 2*t2) - 
                       5*Power(t1,2)*(-2 + t2) + 3*Power(t2,2) + 
                       2*t1*(7 - 4*t2 + 2*Power(t2,2)) + 
                       s2*(-2*Power(t1,2) + t1*(9 + t2) + 
                        2*(7 - 7*t2 + Power(t2,2))))) - 
                 (2*(-9 + s1 + Power(s1,2) + 2*s1*t1 + Power(t1,2) + 
                      4*t2 - 2*s1*t2 - 2*t1*t2)*
                    (-1 + Power(-1 + s,2)*Power(t1,5) + t2 + 
                      2*s*s2*t2 - 3*Power(t2,2) + 6*s*Power(t2,2) - 
                      Power(s,2)*Power(t2,2) + 6*s2*Power(t2,2) - 
                      6*s*s2*Power(t2,2) - Power(s2,2)*Power(t2,2) + 
                      5*Power(t2,3) - 4*s*Power(t2,3) + 
                      Power(s,2)*Power(t2,3) - 4*s2*Power(t2,3) - 
                      4*s*s2*Power(t2,3) + Power(s2,2)*Power(t2,3) - 
                      2*s*Power(t2,4) - 2*s2*Power(t2,4) - 
                      2*Power(t2,5) + 
                      Power(t1,4)*
                       (-1 + Power(s,2)*(1 - 3*t2) + 
                        2*s*(3 + s2)*t2 - (1 + 2*s2)*t2) + 
                      Power(s1,3)*
                       (Power(s2,2)*Power(-1 + t1,2) - 
                        2*s2*(-1 + t1)*t1*(t1 - t2) + 
                        Power(t1,2)*
                        (1 + Power(t1,2) - 2*t2 - 2*t1*t2 + 
                        2*Power(t2,2))) + 
                      Power(t1,3)*
                       (-2 + 3*Power(t2,2) + 
                        Power(s2,2)*Power(t2,2) + 2*s2*t2*(1 + t2) + 
                        Power(s,2)*(-1 - 4*t2 + 3*Power(t2,2)) - 
                        2*s*
                        (-2 + (-3 + s2)*t2 + 2*(2 + s2)*Power(t2,2))) \
+ t1*(1 - 2*(2 + s2)*t2 + (1 - 6*s2 + Power(s2,2))*Power(t2,2) + 
                        (-2 + 6*s2)*Power(t2,3) + 
                        2*(2 + s2)*Power(t2,4) + 
                        Power(s,2)*t2*(2 - 3*t2 - 2*Power(t2,2)) - 
                        2*s*
                        (1 + t2 - s2*t2 - (4 + 5*s2)*Power(t2,2) - 
                        (3 + s2)*Power(t2,3) + Power(t2,4))) - 
                      Power(t1,2)*
                       (-2 - 2*(2 + s2)*t2 + 
                        Power(1 + s2,2)*Power(t2,2) + 
                        (5 + 2*s2 + Power(s2,2))*Power(t2,3) + 
                        Power(s,2)*
                       (1 - 3*t2 - 5*Power(t2,2) + Power(t2,3)) - 
                        2*s*t2*
                        (-5 - 5*t2 + 3*Power(t2,2) + 
                        s2*(-2 + Power(t2,2)))) + 
                      s1*(1 - 2*s*Power(t1,5) + 
                        Power(t1,4)*
                       (1 + Power(s,2) + 6*s*(-1 + t2) - 4*t2) - 
                        Power(t2,2) - 2*s*Power(t2,2) + 
                        Power(s,2)*Power(t2,2) - 2*Power(t2,3) + 
                        2*s*Power(t2,3) + 2*Power(t2,4) + 
                        Power(s2,2)*(-1 + t1)*t2*
                       (t1*(-4 + t2) + 3*t2) - 
                        2*Power(t1,3)*
                       (1 + Power(s,2)*t2 - 4*Power(t2,2) + 
                       s*(-1 - 6*t2 + 3*Power(t2,2))) + 
                        Power(t1,2)*
                       (2 - 4*t2 + 7*Power(t2,2) - 8*Power(t2,3) + 
                       Power(s,2)*(1 + Power(t2,2)) + 
                       2*s*(1 - t2 - 5*Power(t2,2) + Power(t2,3))) - 
                        2*t1*
                       (1 - 4*t2 + Power(s,2)*t2 + 4*Power(t2,2) + 
                       Power(t2,3) - 2*Power(t2,4) - 
                       s*(2 - 4*t2 + Power(t2,2) + 2*Power(t2,3))) + 
                        2*s2*(-1 + t1)*
                        (-1 + s*Power(t1,3) + t2 + 2*s*t2 - 
                        s*Power(t2,2) + 
                        Power(t1,2)*(1 + s + t2 + Power(t2,2)) - 
                        t1*t2*(-2 + 4*t2 + Power(t2,2) + s*(2 + t2)))\
) + Power(s1,2)*(Power(s2,2)*(-1 + t1)*
                       (1 + Power(t1,2) - t2 - t1*t2) - 
                        2*s2*(-1 + t1)*
                        (-1 + t1 + Power(t1,3) + s*t2 - t1*t2 + 
                        Power(t2,2) - Power(t1,2)*(-1 + s + 2*t2)) + 
                        t1*(Power(t1,4) + Power(t1,3)*(1 - 3*t2) - 
                        4*(-1 + t2)*Power(t2,2) + 
                        Power(t1,2)*(3 - 6*t2 + 4*Power(t2,2)) + 
                        t1*(-5 + 3*t2 + 4*Power(t2,2) - 
                       2*Power(t2,3)) - 
                        2*s*(1 + Power(t1,3) + 
                        Power(t1,2)*(1 - 2*t2) - t2 + Power(t2,2) + 
                        t1*(-1 - t2 + Power(t2,2)))))))/
                  ((-1 + s1 + t1 - t2)*Power(-(s1*t1) + t2,2))))/
             ((s - s2 + t1)*(s - s1 + t2)) + 
            ((-1 + t2)*((2*(1 - Power(t1,2) + 
                      s1*(s2*(-1 + t1) - t1*(1 + t1 - 2*t2)) + 
                      s*(1 + t1)*(t1 - t2) + t2 - s2*t2 + t1*t2 + 
                      s2*t1*t2 - 2*Power(t2,2))*
                    (-9 + 6*s2 - 12*t1 - 5*s2*t1 + 10*Power(t1,2) + 
                      s2*Power(t1,2) - Power(t1,3) + t2 - 3*s2*t2 - 
                      14*t1*t2 + 2*s2*t1*t2 + Power(t1,2)*t2 + 
                      t1*Power(t2,2) + Power(s1,2)*(3 + t2) - 
                      s*(3 + Power(t1,2) + t1*(-4 + t2) - 2*t2 + 
                       s1*(-1 + t1 + t2)) + 
                      s1*(-1 - Power(t1,2) + 2*(-5 + s2)*t2 + 
                        Power(t2,2) + t1*(8 + s2 + 2*t2))))/
                  (s1*t1 - t2) + 
                 4*(-5 + 11*s2 - 4*Power(s2,2) - 4*t1 + 16*s2*t1 + 
                    Power(s2,2)*t1 - 3*Power(t1,2) - 
                    5*s2*Power(t1,2) + 2*Power(t1,3) + 5*t2 - 
                    14*s2*t2 + 2*Power(s2,2)*t2 + 14*t1*t2 + 
                    9*s2*t1*t2 - Power(s2,2)*t1*t2 - 
                    9*Power(t1,2)*t2 - 2*Power(t2,2) - 
                    s2*Power(t2,2) + 11*t1*Power(t2,2) - 
                    s2*t1*Power(t2,2) - 2*Power(t2,3) + 
                    Power(s,2)*(-2 + t1 + t2) - 
                    Power(s1,2)*
                     (-2 + t1 + t2 + Power(t2,2) + s2*(2 + t2)) + 
                    s*(-3 - 8*t1 + 3*Power(t1,2) + 
                       s2*(5 + t1*(-4 + t2) - 2*t2) + 7*t2 - 
                       4*t1*t2 + t1*Power(t2,2) + 
                       s1*(-2 + 2*t1 + s2*(-1 + t2) - 5*t2 + 
                        Power(t2,2))) - 
                    s1*(-1 - Power(t1,2) - 2*t2 - 8*Power(t2,2) + 
                       Power(s2,2)*(1 + t2) + 
                       t1*(-2 + 7*t2 + Power(t2,2)) + 
                       s2*(-8 - 9*t2 + Power(t2,2) + t1*(3 + t2)))) - 
                 (2*(2 + Power(t1,2) + (-2 + s1)*t2 + 
                      t1*(-4 + s1 + t2))*
                    (-1 + Power(-1 + s,2)*Power(t1,5) + t2 + 
                      2*s*s2*t2 - 3*Power(t2,2) + 6*s*Power(t2,2) - 
                      Power(s,2)*Power(t2,2) + 6*s2*Power(t2,2) - 
                      6*s*s2*Power(t2,2) - Power(s2,2)*Power(t2,2) + 
                      5*Power(t2,3) - 4*s*Power(t2,3) + 
                      Power(s,2)*Power(t2,3) - 4*s2*Power(t2,3) - 
                      4*s*s2*Power(t2,3) + Power(s2,2)*Power(t2,3) - 
                      2*s*Power(t2,4) - 2*s2*Power(t2,4) - 
                      2*Power(t2,5) + 
                      Power(t1,4)*
                       (-1 + Power(s,2)*(1 - 3*t2) + 
                        2*s*(3 + s2)*t2 - (1 + 2*s2)*t2) + 
                      Power(s1,3)*
                       (Power(s2,2)*Power(-1 + t1,2) - 
                        2*s2*(-1 + t1)*t1*(t1 - t2) + 
                        Power(t1,2)*
                        (1 + Power(t1,2) - 2*t2 - 2*t1*t2 + 
                        2*Power(t2,2))) + 
                      Power(t1,3)*
                       (-2 + 3*Power(t2,2) + 
                        Power(s2,2)*Power(t2,2) + 2*s2*t2*(1 + t2) + 
                        Power(s,2)*(-1 - 4*t2 + 3*Power(t2,2)) - 
                        2*s*
                        (-2 + (-3 + s2)*t2 + 2*(2 + s2)*Power(t2,2))) \
+ t1*(1 - 2*(2 + s2)*t2 + (1 - 6*s2 + Power(s2,2))*Power(t2,2) + 
                        (-2 + 6*s2)*Power(t2,3) + 
                        2*(2 + s2)*Power(t2,4) + 
                        Power(s,2)*t2*(2 - 3*t2 - 2*Power(t2,2)) - 
                        2*s*
                        (1 + t2 - s2*t2 - (4 + 5*s2)*Power(t2,2) - 
                        (3 + s2)*Power(t2,3) + Power(t2,4))) - 
                      Power(t1,2)*
                       (-2 - 2*(2 + s2)*t2 + 
                        Power(1 + s2,2)*Power(t2,2) + 
                        (5 + 2*s2 + Power(s2,2))*Power(t2,3) + 
                        Power(s,2)*
                       (1 - 3*t2 - 5*Power(t2,2) + Power(t2,3)) - 
                        2*s*t2*
                        (-5 - 5*t2 + 3*Power(t2,2) + 
                        s2*(-2 + Power(t2,2)))) + 
                      s1*(1 - 2*s*Power(t1,5) + 
                        Power(t1,4)*
                       (1 + Power(s,2) + 6*s*(-1 + t2) - 4*t2) - 
                        Power(t2,2) - 2*s*Power(t2,2) + 
                        Power(s,2)*Power(t2,2) - 2*Power(t2,3) + 
                        2*s*Power(t2,3) + 2*Power(t2,4) + 
                        Power(s2,2)*(-1 + t1)*t2*
                       (t1*(-4 + t2) + 3*t2) - 
                        2*Power(t1,3)*
                       (1 + Power(s,2)*t2 - 4*Power(t2,2) + 
                       s*(-1 - 6*t2 + 3*Power(t2,2))) + 
                        Power(t1,2)*
                       (2 - 4*t2 + 7*Power(t2,2) - 8*Power(t2,3) + 
                       Power(s,2)*(1 + Power(t2,2)) + 
                       2*s*(1 - t2 - 5*Power(t2,2) + Power(t2,3))) - 
                        2*t1*
                       (1 - 4*t2 + Power(s,2)*t2 + 4*Power(t2,2) + 
                       Power(t2,3) - 2*Power(t2,4) - 
                       s*(2 - 4*t2 + Power(t2,2) + 2*Power(t2,3))) + 
                        2*s2*(-1 + t1)*
                        (-1 + s*Power(t1,3) + t2 + 2*s*t2 - 
                        s*Power(t2,2) + 
                        Power(t1,2)*(1 + s + t2 + Power(t2,2)) - 
                        t1*t2*(-2 + 4*t2 + Power(t2,2) + s*(2 + t2)))\
) + Power(s1,2)*(Power(s2,2)*(-1 + t1)*
                       (1 + Power(t1,2) - t2 - t1*t2) - 
                        2*s2*(-1 + t1)*
                        (-1 + t1 + Power(t1,3) + s*t2 - t1*t2 + 
                        Power(t2,2) - Power(t1,2)*(-1 + s + 2*t2)) + 
                        t1*(Power(t1,4) + Power(t1,3)*(1 - 3*t2) - 
                        4*(-1 + t2)*Power(t2,2) + 
                        Power(t1,2)*(3 - 6*t2 + 4*Power(t2,2)) + 
                        t1*(-5 + 3*t2 + 4*Power(t2,2) - 
                       2*Power(t2,3)) - 
                        2*s*(1 + Power(t1,3) + 
                        Power(t1,2)*(1 - 2*t2) - t2 + Power(t2,2) + 
                        t1*(-1 - t2 + Power(t2,2)))))))/
                  ((-1 + s1 + t1 - t2)*Power(-(s1*t1) + t2,2))))/
             ((-1 + s2)*(-1 + t1)) + 
            ((-1 + t2)*((2*(1 - Power(t1,2) + 
                      s1*(s2*(-1 + t1) - t1*(1 + t1 - 2*t2)) + 
                      s*(1 + t1)*(t1 - t2) + t2 - s2*t2 + t1*t2 + 
                      s2*t1*t2 - 2*Power(t2,2))*
                    (-5 + 2*s2 - 10*t1 + 2*s2*t1 + 3*Power(t1,2) + 
                      t2 - 7*s2*t2 - 6*t1*t2 + 2*s2*t1*t2 - 
                      4*Power(t2,2) + t1*Power(t2,2) + 
                      Power(s1,2)*(-2 + 2*s2 - 2*t1 + 3*t2) - 
                      s*(2 + Power(s1,2) - Power(t1,2) + 
                        s1*(-1 + t2) - 3*t2 + t1*(4 + t2)) + 
                      s1*(3 - 2*Power(t1,2) - 8*t2 + Power(t2,2) + 
                        s2*(-5 + 2*t1 + 2*t2) + t1*(8 + 3*t2))))/
                  (s1*t1 - t2) + 
                 4*(-8 + 4*Power(s1,3) + 7*t1 - 14*s*t1 - 
                    3*Power(s,2)*t1 + Power(t1,2) - 4*s*Power(t1,2) + 
                    Power(s,2)*Power(t1,2) + s*Power(t1,3) - 5*t2 + 
                    12*s*t2 + 2*Power(s,2)*t2 + 5*t1*t2 + 5*s*t1*t2 - 
                    4*Power(t1,2)*t2 - s*Power(t1,2)*t2 + 
                    Power(t1,3)*t2 - 3*s*Power(t2,2) + 
                    6*t1*Power(t2,2) + s*t1*Power(t2,2) - 
                    Power(t1,2)*Power(t2,2) + Power(t2,3) - 
                    Power(s2,2)*(1 + t1 - 4*t2 + t1*t2) + 
                    Power(s1,2)*
                     (4 + s2 - Power(s2,2) + 4*t1 + 2*s2*t1 - 
                       Power(t1,2) - 3*t2 - 3*s2*t2 + 3*t1*t2 - 
                       2*Power(t2,2) + s*(-3 + s2 - t1 + t2)) + 
                    s1*(-3 + Power(s,2)*(-1 + t1) + 15*t1 + 
                       4*Power(t1,2) - Power(t1,3) - 6*t2 - 
                       12*t1*t2 + 4*Power(t1,2)*t2 + 6*Power(t2,2) - 
                       3*t1*Power(t2,2) - 
                       Power(s2,2)*(-2 + t1 + t2) + 
                       s2*(7 - 7*t1 + 2*Power(t1,2) + 10*t2 - 
                        3*t1*t2 - Power(t2,2)) + 
                       s*(-4 + s2 - t1 - 2*t2 + s2*t2 + Power(t2,2))) \
+ s2*(5 - 3*Power(t1,2) - 11*t2 + 3*Power(t2,2) + 
                       t1*(10 + 5*t2 - Power(t2,2)) + 
                       s*(2 - Power(t1,2) - 7*t2 + t1*(6 + t2)))) - 
                 (2*(1 + Power(s1,2) + t1 - 3*t2 + t1*t2 + 
                      s1*(-3 + t1 + t2))*
                    (-1 + Power(-1 + s,2)*Power(t1,5) + t2 + 
                      2*s*s2*t2 - 3*Power(t2,2) + 6*s*Power(t2,2) - 
                      Power(s,2)*Power(t2,2) + 6*s2*Power(t2,2) - 
                      6*s*s2*Power(t2,2) - Power(s2,2)*Power(t2,2) + 
                      5*Power(t2,3) - 4*s*Power(t2,3) + 
                      Power(s,2)*Power(t2,3) - 4*s2*Power(t2,3) - 
                      4*s*s2*Power(t2,3) + Power(s2,2)*Power(t2,3) - 
                      2*s*Power(t2,4) - 2*s2*Power(t2,4) - 
                      2*Power(t2,5) + 
                      Power(t1,4)*
                       (-1 + Power(s,2)*(1 - 3*t2) + 
                        2*s*(3 + s2)*t2 - (1 + 2*s2)*t2) + 
                      Power(s1,3)*
                       (Power(s2,2)*Power(-1 + t1,2) - 
                        2*s2*(-1 + t1)*t1*(t1 - t2) + 
                        Power(t1,2)*
                        (1 + Power(t1,2) - 2*t2 - 2*t1*t2 + 
                        2*Power(t2,2))) + 
                      Power(t1,3)*
                       (-2 + 3*Power(t2,2) + 
                        Power(s2,2)*Power(t2,2) + 2*s2*t2*(1 + t2) + 
                        Power(s,2)*(-1 - 4*t2 + 3*Power(t2,2)) - 
                        2*s*(-2 + (-3 + s2)*t2 + 
                        2*(2 + s2)*Power(t2,2))) + 
                      t1*(1 - 2*(2 + s2)*t2 + 
                        (1 - 6*s2 + Power(s2,2))*Power(t2,2) + 
                        (-2 + 6*s2)*Power(t2,3) + 
                        2*(2 + s2)*Power(t2,4) + 
                        Power(s,2)*t2*(2 - 3*t2 - 2*Power(t2,2)) - 
                        2*s*(1 + t2 - s2*t2 - 
                        (4 + 5*s2)*Power(t2,2) - 
                        (3 + s2)*Power(t2,3) + Power(t2,4))) - 
                      Power(t1,2)*
                       (-2 - 2*(2 + s2)*t2 + 
                        Power(1 + s2,2)*Power(t2,2) + 
                        (5 + 2*s2 + Power(s2,2))*Power(t2,3) + 
                        Power(s,2)*
                        (1 - 3*t2 - 5*Power(t2,2) + Power(t2,3)) - 
                        2*s*t2*
                        (-5 - 5*t2 + 3*Power(t2,2) + 
                        s2*(-2 + Power(t2,2)))) + 
                      s1*(1 - 2*s*Power(t1,5) + 
                        Power(t1,4)*
                        (1 + Power(s,2) + 6*s*(-1 + t2) - 4*t2) - 
                        Power(t2,2) - 2*s*Power(t2,2) + 
                        Power(s,2)*Power(t2,2) - 2*Power(t2,3) + 
                        2*s*Power(t2,3) + 2*Power(t2,4) + 
                        Power(s2,2)*(-1 + t1)*t2*
                        (t1*(-4 + t2) + 3*t2) - 
                        2*Power(t1,3)*
                        (1 + Power(s,2)*t2 - 4*Power(t2,2) + 
                        s*(-1 - 6*t2 + 3*Power(t2,2))) + 
                        Power(t1,2)*
                        (2 - 4*t2 + 7*Power(t2,2) - 8*Power(t2,3) + 
                        Power(s,2)*(1 + Power(t2,2)) + 
                        2*s*(1 - t2 - 5*Power(t2,2) + Power(t2,3))) - 
                        2*t1*
                        (1 - 4*t2 + Power(s,2)*t2 + 4*Power(t2,2) + 
                        Power(t2,3) - 2*Power(t2,4) - 
                        s*(2 - 4*t2 + Power(t2,2) + 2*Power(t2,3))) + 
                        2*s2*(-1 + t1)*
                        (-1 + s*Power(t1,3) + t2 + 2*s*t2 - 
                        s*Power(t2,2) + 
                        Power(t1,2)*(1 + s + t2 + Power(t2,2)) - 
                        t1*t2*(-2 + 4*t2 + Power(t2,2) + s*(2 + t2)))) \
+ Power(s1,2)*(Power(s2,2)*(-1 + t1)*(1 + Power(t1,2) - t2 - t1*t2) - 
                        2*s2*(-1 + t1)*
                        (-1 + t1 + Power(t1,3) + s*t2 - t1*t2 + 
                        Power(t2,2) - Power(t1,2)*(-1 + s + 2*t2)) + 
                        t1*(Power(t1,4) + Power(t1,3)*(1 - 3*t2) - 
                        4*(-1 + t2)*Power(t2,2) + 
                        Power(t1,2)*(3 - 6*t2 + 4*Power(t2,2)) + 
                        t1*(-5 + 3*t2 + 4*Power(t2,2) - 
                       2*Power(t2,3)) - 
                        2*s*(1 + Power(t1,3) + 
                        Power(t1,2)*(1 - 2*t2) - t2 + Power(t2,2) + 
                        t1*(-1 - t2 + Power(t2,2)))))))/
                  ((-1 + s1 + t1 - t2)*Power(-(s1*t1) + t2,2))))/
             ((-1 + s1)*(-s + s1 - t2))))/Power(-1 + t2,2))*
     B1(t1,1 - s1 - t1 + t2,t2))/(16.*Power(Pi,2)) + 
  (((2*s*(-2 + 2*t1 + s1*(1 + t1 + s1*(2 - s2 + t1)) + 
            Power(s,2)*(t1 - t2) + t2 - 
            (3*t1 + s1*(3 + 2*s1 + t1))*t2 + 
            (1 + 2*s1 + s2)*Power(t2,2) + 
            s*(1 - 3*t1 + t2 + (s2 + t1 - t2)*t2 + 
               s1*(-2 + s2 - 2*t1 + 3*t2)))*
          ((4 - Power(s1,4) + Power(s,3)*(3 + s1) + 15*s2 + 
               Power(s2,2) - 4*t1 - s2*t1 - 
               Power(s,2)*(2 + 3*Power(s1,2) + 4*s2 + 
                  s1*(3 + 2*s2 - t2) - 5*t2) - 4*t2 - s2*t2 + 
               2*Power(t2,2) - 3*s2*Power(t2,2) + 
               Power(s1,3)*(5 - 2*s2 + t2) + 
               Power(s1,2)*(9 - Power(s2,2) + t1 - 3*t2 + 
                  2*s2*(2 + t2)) + 
               s*(-19 + 3*Power(s1,3) + 3*s2 + t1 + 
                  Power(s1,2)*(-5 + 4*s2 - 2*t2) + 5*t2 - 7*s2*t2 + 
                  2*Power(t2,2) + 
                  s1*(-7 + 2*s2 + Power(s2,2) - t1 - 2*t2 - 2*s2*t2)) \
+ s1*(s2*(4 + t1 - 3*t2) + Power(s2,2)*(-1 + t2) + 
                  3*(5 + t1 - 3*t2 - Power(t2,2))))/
             ((-1 + s1)*(-1 + t2)) + 
            (-7 + Power(s,3)*s1 - Power(s1,4) + 11*s2 - 
               6*Power(s2,2) - t1 + 7*s2*t1 + 9*t2 - 5*s2*t2 + 
               Power(s2,2)*t2 - 2*s2*t1*t2 + 2*s2*Power(t2,2) + 
               Power(s1,3)*(5 - 2*s2 + t2) - 
               Power(s,2)*(5 + 3*Power(s1,2) + s2 - 2*t1 + 
                  s1*(-8 + 2*s2 - t2) + 2*t2) + 
               s1*(8 + s2*(8 + 3*t1 - 11*t2) + t1*(9 - 2*t2) + 
                  Power(s2,2)*(-3 + t2) - 7*t2 + 2*Power(t2,2)) + 
               Power(s1,2)*(17 - Power(s2,2) + 2*t1 - 11*t2 + 
                  s2*(3 + 2*t2)) + 
               s*(-2 + 3*Power(s1,3) + Power(s2,2) - 7*t1 + 
                  Power(s1,2)*(-13 + 4*s2 - 2*t2) + t2 + 2*t1*t2 - 
                  2*Power(t2,2) + s2*(7 - 2*t1 + t2) + 
                  s1*(-11 + Power(s2,2) - 5*t1 + 14*t2 - 
                     2*s2*(2 + t2))))/((s - s2 + t1)*(s - s1 + t2)) + 
            (-1 + Power(s,3)*(-1 + s1) - Power(s1,4) + 12*s2 + 
               2*Power(s2,2) + t1 - 2*s2*t1 + 3*t2 - 2*s2*t2 - 
               Power(s2,2)*t2 - 2*Power(t2,2) + 3*s2*Power(t2,2) + 
               Power(s1,3)*(7 - 2*s2 + t2) + 
               s1*(11 - 4*t1 + 2*s2*(4 + t1 - 5*t2) + 
                  Power(s2,2)*(-1 + t2) + 3*Power(t2,2)) + 
               Power(s,2)*(2 - 3*Power(s1,2) + 3*s2 - 3*t2 + 
                  s1*(6 - 2*s2 + t2)) + 
               Power(s1,2)*(5 - Power(s2,2) + 3*t1 - 10*t2 + 
                  s2*(5 + 2*t2)) + 
               s*(-9 + 3*Power(s1,3) - Power(s2,2) + 2*t1 + 
                  2*Power(s1,2)*(-6 + 2*s2 - t2) + 6*s2*(-1 + t2) - 
                  2*t2 - 2*Power(t2,2) + 
                  s1*(Power(s2,2) - 2*(4 + t1 - 6*t2) - 2*s2*(3 + t2))\
))/((-1 + s1)*(-s + s1 - t2)) + 
            (-9 - Power(s1,4) + Power(s,3)*(1 + s1) + 8*s2 - 
               2*Power(s2,2) - 2*Power(s2,3) + 3*t1 - s2*t1 + 
               6*Power(s2,2)*t1 - 4*s2*Power(t1,2) + 
               Power(s1,3)*(2 - 3*s2 + t1) + 
               Power(s,2)*(1 - 3*Power(s1,2) - 3*s2 + 4*t1 + 
                  s1*(2 - 3*s2 + t1) - t2) + 3*t2 - 3*s2*t2 - 
               2*Power(s2,2)*t2 + 3*s2*t1*t2 - s2*Power(t2,2) + 
               Power(s1,2)*(16 - 2*Power(s2,2) - 5*t1 + 
                  s2*(3 + t1 + t2)) + 
               s1*(12 - 4*Power(t1,2) + s2*(15 + t1 - 2*t2) + 
                  Power(s2,2)*(-1 + t2) - Power(t2,2) + 
                  t1*(-4 + 3*t2)) + 
               s*(-3 + 3*Power(s1,3) + 5*Power(s2,2) + 
                  Power(s1,2)*(-5 + 6*s2 - 2*t1) + t1 + 
                  4*Power(t1,2) - t2 - 3*t1*t2 + Power(t2,2) + 
                  s2*(-3 - 10*t1 + 4*t2) + 
                  s1*(-17 + 2*Power(s2,2) + t1 + t2 - 
                     s2*(2 + t1 + t2))))/((-1 + s2)*(-s + s2 - t1)) + 
            (-6 - 2*Power(s1,4) + Power(s,3)*(-3 + 2*s1) + 13*s2 + 
               7*Power(s2,2) - 2*t1 - 17*s2*t1 - 2*Power(s2,2)*t1 + 
               4*s2*Power(t1,2) + 2*t2 + 6*s2*t2 - Power(s2,2)*t2 - 
               3*s2*t1*t2 + s2*Power(t2,2) + 
               Power(s1,3)*(2 - 4*s2 + 2*t2) + 
               s1*(17 + 4*Power(t1,2) + s2*(17 - 4*t1 - 3*t2) + 
                  2*Power(s2,2)*t2 + Power(t2,2) - 3*t1*(5 + t2)) - 
               2*Power(s1,2)*
                (-8 + Power(s2,2) + t1 + t2 - s2*(1 + 2*t2)) + 
               Power(s,2)*(16 - 6*Power(s1,2) + 4*s2 - 5*t1 + 
                  s1*(7 - 4*s2 + 2*t2)) + 
               s*(-11 + 6*Power(s1,3) - Power(s2,2) + 17*t1 - 
                  4*Power(t1,2) + Power(s1,2)*(-6 + 8*s2 - 4*t2) - 
                  8*t2 + 3*t1*t2 - Power(t2,2) + 
                  s2*(-25 + 7*t1 + t2) + 
                  s1*(-32 + 2*Power(s2,2) + 7*t1 + 2*t2 - 
                     s2*(5 + 4*t2))))/((-1 + s2)*(-1 + t1)) + 
            (2 - Power(s,4) - Power(s1,4) + 19*s2 - 2*Power(s2,2) - 
               3*t1 - 5*s2*t1 + Power(t1,2) + 
               Power(s,3)*(4*s1 + 2*s2 - t1 - t2) - 8*t2 - 8*s2*t2 + 
               3*Power(s2,2)*t2 + t1*t2 + s2*t1*t2 - 
               2*s2*Power(t2,2) + Power(s1,3)*(1 - 2*s2 + t2) - 
               Power(s1,2)*(-17 + Power(s2,2) + t1 - 2*s2*(-1 + t2) + 
                  t2) + s1*(25 - Power(t1,2) + 
                  Power(s2,2)*(-3 + t2) - 12*t2 - 2*Power(t2,2) + 
                  t1*(-3 + 2*t2) + s2*(11 + 2*t2)) + 
               Power(s,2)*(16 - 6*Power(s1,2) - Power(s2,2) + t1 - 
                  2*t2 - t1*t2 + s2*(-3 + t1 + 2*t2) + 
                  s1*(2 - 6*s2 + 2*t1 + 3*t2)) + 
               s*(-27 + 4*Power(s1,3) + 6*t1 + 
                  s2*(-11 + t1*(-1 + t2) - t2) - 
                  Power(s2,2)*(-3 + t2) + 11*t2 - t1*t2 + 
                  2*Power(t2,2) + 
                  Power(s1,2)*(6*s2 - t1 - 3*(1 + t2)) + 
                  s1*(-33 + 2*Power(s2,2) + (3 + t1)*t2 - 
                     s2*(-4 + t1 + 4*t2))))/((-1 + t1)*(-1 + t2))))/
        (Power(s,3)*(-1 + s1)*(-1 + s2 - t1 + t2) - 
          Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
          Power(s,2)*(3 + s2*(-3 + t1) + 2*t1 - (3 + t1)*t2 + 
             Power(t2,2) + 2*Power(s1,2)*(-1 + s2 - t1 + t2) + 
             s1*(s2 - t2)*(-1 + s2 - t1 + t2)) + 
          s*(2*(-1 + s2)*(-1 + t1) - 
             (2*(2 + t1) + s2*(-5 + s2 + t1))*t2 + 
             (2 + s2)*Power(t2,2) + Power(s1,3)*(-1 + s2 - t1 + t2) + 
             Power(s1,2)*(1 + s2 - t2)*(-1 + s2 - t1 + t2) + 
             s1*(2 + s2*(-3 + s2 - t1) + 4*t1 - t2 - 
                (1 + s2)*(s2 - t1)*t2 - (1 + s2)*Power(t2,2)))) + 
       2*((16 - 16*t1 - 36*t2 - 4*s2*t2 + 16*t1*t2 + 27*Power(t2,2) + 
             8*s2*Power(t2,2) + 3*t1*Power(t2,2) - 11*Power(t2,3) - 
             7*s2*Power(t2,3) + 4*Power(t2,4) + Power(s,4)*(t1 + t2) + 
             Power(s1,4)*(-2 - s2 + t1 + 2*t2) + 
             Power(s,3)*(1 + s1*(s2 - 4*t1 - 5*t2) + 2*t1*(-1 + t2) + 
                10*t2 - 3*s2*t2 + 2*Power(t2,2)) - 
             Power(s1,3)*(-11 + s2*(3 - 5*t2) + (7 + 2*t1)*t2 + 
                4*Power(t2,2)) + 
             Power(s1,2)*(27 - 61*t2 + 32*Power(t2,2) + 
                2*Power(t2,3) + s2*(8 + 3*t2 - 7*Power(t2,2)) + 
                t1*(-5 + Power(t2,2))) + 
             s1*(12 + t1*(8 - 6*t2) - 46*t2 + 61*Power(t2,2) - 
                27*Power(t2,3) + 
                s2*(4 - 16*t2 + 7*Power(t2,2) + 3*Power(t2,3))) + 
             Power(s,2)*(-1 - 17*t2 + s2*t2 + 30*Power(t2,2) - 
                6*s2*Power(t2,2) + Power(t2,3) + 
                Power(s1,2)*(-2 - 3*s2 + 6*t1 + 9*t2) + 
                t1*(-9 - 6*t2 + Power(t2,2)) + 
                s1*(9 + 4*t1 - 31*t2 - 6*t1*t2 - 8*Power(t2,2) + 
                   s2*(-3 + 11*t2))) + 
             s*(-10 + 26*t1 + Power(s1,3)*(4 + 3*s2 - 4*t1 - 7*t2) + 
                30*t2 - 37*Power(t2,2) - 6*s2*Power(t2,2) - 
                4*t1*Power(t2,2) + 24*Power(t2,3) - 3*s2*Power(t2,3) + 
                Power(s1,2)*(-23 + s2*(6 - 13*t2) + 28*t2 + 
                   10*Power(t2,2) + t1*(-2 + 6*t2)) + 
                s1*(4 + 72*t2 - 64*Power(t2,2) - 3*Power(t2,3) - 
                   2*t1*(-5 - 3*t2 + Power(t2,2)) + 
                   s2*(-12 - 4*t2 + 13*Power(t2,2)))))/
           ((s - s2 + t1)*(s - s1 + t2)*(s1 - s1*t2 + t2*(-1 + s + t2))) \
+ (10 - 10*t1 - 43*t2 - 10*s2*t2 + 25*t1*t2 + 56*Power(t2,2) - 
             s2*Power(t2,2) - 9*t1*Power(t2,2) - 19*Power(t2,3) + 
             5*s2*Power(t2,3) - 4*Power(t2,4) + 
             2*Power(s,4)*(t1 + 2*t2) + 
             Power(s1,4)*(-2 - s2 + t1 + 2*t2) + 
             Power(s1,3)*(3 + t1*(3 - 2*t2) + t2 - 4*Power(t2,2) + 
                s2*(-6 + 5*t2)) + 
             Power(s,3)*(2 + s1*(2 + 2*s2 - 7*t1 - 13*t2) - t2 - 
                4*s2*t2 + 8*Power(t2,2) + t1*(-9 + 6*t2)) + 
             Power(s1,2)*(s2*(-3 + 17*t2 - 7*Power(t2,2)) + 
                t1*(-1 - 7*t2 + Power(t2,2)) + 
                2*(15 - 18*t2 + 2*Power(t2,2) + Power(t2,3))) + 
             s1*(31 - 84*t2 + 52*Power(t2,2) + Power(t2,3) + 
                t1*(-13 + 8*t2 + 4*Power(t2,2)) + 
                s2*(10 + 4*t2 - 16*Power(t2,2) + 3*Power(t2,3))) + 
             Power(s,2)*(-7 - 46*t2 + 13*s2*t2 + 8*Power(t2,2) - 
                8*s2*Power(t2,2) + 4*Power(t2,3) + 
                Power(s1,2)*(-6 - 5*s2 + 9*t1 + 16*t2) + 
                2*t1*(4 - 9*t2 + 3*Power(t2,2)) + 
                s1*(5 + 21*t1 - 14*t1*t2 - 18*Power(t2,2) + 
                   3*s2*(-3 + 5*t2))) + 
             s*(1 + Power(s1,3)*(6 + 4*s2 - 5*t1 - 9*t2) + 79*t2 - 
                5*s2*t2 - 77*Power(t2,2) + 18*s2*Power(t2,2) + 
                5*Power(t2,3) - 4*s2*Power(t2,3) + 
                t1*(9 - t2 - 9*Power(t2,2) + 2*Power(t2,3)) + 
                Power(s1,2)*(s2*(15 - 16*t2) + 5*t1*(-3 + 2*t2) + 
                   2*(-5 + 7*Power(t2,2))) + 
                s1*(-37 + 92*t2 - 15*Power(t2,2) - 5*Power(t2,3) + 
                   t1*(-5 + 25*t2 - 7*Power(t2,2)) + 
                   s2*(5 - 33*t2 + 16*Power(t2,2)))))/
           ((-1 + t1)*(-1 + t2)*(s1 - s1*t2 + t2*(-1 + s + t2))) - 
          (Power(s,4)*(t1 + 3*t2) + 
             Power(s,3)*(1 + s1*(2 + s2 - 3*t1 - 10*t2) - 3*t2 - 
                s2*t2 + 6*Power(t2,2) + t1*(-3 + 4*t2)) + 
             Power(s,2)*(-2 - 4*t2 + 4*s2*t2 + 4*Power(t2,2) - 
                2*s2*Power(t2,2) + 3*Power(t2,3) + 
                Power(s1,2)*(-6 - 2*s2 + 3*t1 + 13*t2) + 
                t1*(4 - 12*t2 + 5*Power(t2,2)) + 
                2*s1*(2*t1 + t2 - 4*t1*t2 - 7*Power(t2,2) + 
                   s2*(-1 + 3*t2))) + 
             s*(2 + Power(s1,3)*(6 + s2 - t1 - 8*t2) + 32*t2 - 
                8*s2*t2 - 17*Power(t2,2) + 6*s2*Power(t2,2) + 
                13*Power(t2,3) - s2*Power(t2,3) + 
                Power(s1,2)*
                 (-5 + t1 + s2*(2 - 7*t2) + 3*t2 + 4*t1*t2 + 
                   12*Power(t2,2)) + 
                t1*(-6 + 8*t2 - 13*Power(t2,2) + 2*Power(t2,3)) + 
                s1*(-8 + 20*t2 - 22*Power(t2,2) - 4*Power(t2,3) + 
                   t1*(-4 + 14*t2 - 5*Power(t2,2)) + 
                   s2*(4 - 12*t2 + 7*Power(t2,2)))) + 
             2*(-2 + Power(s1,4)*(-1 + t2) - 13*t2 + 4*s2*t2 + 
                18*Power(t2,2) - 4*s2*Power(t2,2) - 6*Power(t2,3) + 
                s2*Power(t2,3) + 3*Power(t2,4) - 
                Power(s1,3)*(-3 + t1 + t2 - s2*t2 + 2*Power(t2,2)) + 
                t1*(2 - 3*t2 + 2*Power(t2,2) - 2*Power(t2,3)) + 
                Power(s1,2)*(-(t1*(-2 + t2)) + 
                   s2*(-2 + 3*t2 - 2*Power(t2,2)) + 
                   t2*(-9 + 8*t2 + Power(t2,2))) + 
                s1*(17 - 20*t2 + 12*Power(t2,2) - 9*Power(t2,3) + 
                   t1*(-1 - 2*t2 + 4*Power(t2,2)) + 
                   s2*(-4 + 6*t2 - 4*Power(t2,2) + Power(t2,3)))))/
           ((-1 + s1)*(-s + s1 - t2)*(s1*(-1 + t2) - t2*(-1 + s + t2))) \
+ (Power(s,4)*(t1 - t2) + Power(s,3)*
              (1 + 7*t2 - s2*t2 - 2*Power(t2,2) + t1*(-5 + 4*t2) + 
                s1*(-2 + s2 - 3*t1 + 6*t2)) - 
             Power(s,2)*(4 + 4*t2 - 2*s2*t2 - 18*Power(t2,2) + 
                2*s2*Power(t2,2) + Power(t2,3) + 
                Power(s1,2)*(-6 + 2*s2 - 3*t1 + 11*t2) + 
                t1*(-6 + 14*t2 - 5*Power(t2,2)) - 
                2*s1*(4 + 6*t1 + s2*(-2 + t2) - 9*t2 - 4*t1*t2 + 
                   5*Power(t2,2))) + 
             s*(2 + 2*t1 - 36*t2 + 6*s2*t2 + 10*t1*t2 + Power(t2,2) - 
                13*t1*Power(t2,2) + 17*Power(t2,3) - s2*Power(t2,3) + 
                2*t1*Power(t2,3) + 
                Power(s1,3)*(-6 + s2 - t1 + 8*t2) + 
                Power(s1,2)*
                 (-7 + 13*t2 - 12*Power(t2,2) + s2*(4 + t2) + 
                   t1*(-9 + 4*t2)) - 
                s1*(-4 + 2*(-5 + s2)*t2 + (22 + s2)*Power(t2,2) - 
                   4*Power(t2,3) + 5*t1*(2 - 4*t2 + Power(t2,2)))) - 
             2*(-2 + Power(s1,4)*(-1 + t2) - 15*t2 + 5*s2*t2 + 
                23*Power(t2,2) - 5*s2*Power(t2,2) - 3*Power(t2,3) + 
                s2*Power(t2,3) - 3*Power(t2,4) + 
                Power(s1,3)*(1 - t1 + t2 + s2*t2 - 2*Power(t2,2)) - 
                Power(s1,2)*(-2 + s2 + 2*t1 + 2*t2 - 2*s2*t2 - 
                   3*t1*t2 + Power(t2,2) + 2*s2*Power(t2,2) - 
                   Power(t2,3)) + 
                t1*(2 - 3*t2 - 2*Power(t2,2) + 2*Power(t2,3)) + 
                s1*(19 - 27*t2 + 4*Power(t2,2) + 4*Power(t2,3) + 
                   t1*(-1 + 6*t2 - 4*Power(t2,2)) + 
                   s2*(-5 + 6*t2 - 3*Power(t2,2) + Power(t2,3)))))/
           ((-1 + s1)*(-1 + t2)*(s1*(-1 + t2) - t2*(-1 + s + t2))) + 
          (2*Power(s,6)*(Power(t1,2) + Power(t2,2)) - 
             4*(-1 + t2)*Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2)*
              (-s1 + Power(s1,2) + t2 - 2*s1*t2 + Power(t2,2)) - 
             Power(s,5)*(2*Power(t1,2)*(5 + 5*s1 - 3*t2) + 
                t1*(-4 - 7*t2 + s1*(4 - 4*s2 + 2*t2)) + 
                t2*(-4 + t2 + 12*s1*t2 - 6*Power(t2,2) + 
                   s2*(4 - 4*s1 + 4*t2))) + 
             Power(s,3)*(-6 + 12*t1 + 2*Power(t1,2) - 13*t2 + 
                21*t1*t2 + 22*Power(t1,2)*t2 + 77*Power(t2,2) + 
                81*s2*Power(t2,2) - 6*Power(s2,2)*Power(t2,2) - 
                111*t1*Power(t2,2) - 10*s2*t1*Power(t2,2) - 
                18*Power(t1,2)*Power(t2,2) - 41*Power(t2,3) - 
                10*s2*Power(t2,3) + 6*Power(s2,2)*Power(t2,3) + 
                37*t1*Power(t2,3) + 2*Power(t1,2)*Power(t2,3) - 
                11*Power(t2,4) - 12*s2*Power(t2,4) + 2*Power(t2,5) - 
                2*Power(s1,3)*
                 (2 + s2 + 3*Power(s2,2) + 3*t1 - 12*s2*t1 + 
                   10*Power(t1,2) - 14*t2 - 9*s2*t2 + 6*t1*t2 + 
                   16*Power(t2,2)) + 
                Power(s1,2)*
                 (36*Power(t1,2)*(-1 + t2) + 
                   6*Power(s2,2)*(-1 + 3*t2) + 
                   s2*(7 + t1*(32 - 48*t2) + 8*t2 - 
                      60*Power(t2,2)) + 
                   t1*(-26 + 47*t2 + 18*Power(t2,2)) + 
                   t2*(-5 - 31*t2 + 54*Power(t2,2))) + 
                s1*(1 + (-41 - 70*s2 + 12*Power(s2,2))*t2 + 
                   2*(43 + s2 - 9*Power(s2,2))*Power(t2,2) + 
                   2*(8 + 27*s2)*Power(t2,3) - 24*Power(t2,4) - 
                   2*Power(t1,2)*(10 - 31*t2 + 9*Power(t2,2)) + 
                   t1*(1 + (129 - 38*s2)*t2 + 
                      4*(-23 + 6*s2)*Power(t2,2) - 6*Power(t2,3)))) + 
             Power(s,4)*(2 - 9*t2 + 12*s2*t2 - 30*Power(t2,2) - 
                11*s2*Power(t2,2) + 2*Power(s2,2)*Power(t2,2) - 
                8*Power(t2,3) - 12*s2*Power(t2,3) + 6*Power(t2,4) + 
                2*Power(t1,2)*(7 - 13*t2 + 3*Power(t2,2)) + 
                2*Power(s1,2)*
                 (Power(s2,2) + 5*t1 - 8*s2*t1 + 10*Power(t1,2) - 
                   5*t2 - 7*s2*t2 + 4*t1*t2 + 14*Power(t2,2)) + 
                2*t1*(-8 - 11*t2 + 16*Power(t2,2)) + 
                s1*(-8*Power(t1,2)*(-4 + 3*t2) + 
                   t1*(7 - 35*t2 - 6*Power(t2,2) + 
                      4*s2*(-3 + 4*t2)) + 
                   t2*(-3 - 4*Power(s2,2) + 6*t2 - 30*Power(t2,2) + 
                      s2*(7 + 30*t2)))) + 
             Power(s,2)*(2*Power(t1,2)*(-8 + 11*t2 + Power(t2,3)) + 
                t1*(16 - 46*t2 + (107 + 22*s2)*Power(t2,2) - 
                   (111 + 20*s2)*Power(t2,3) + 6*Power(t2,4)) + 
                2*Power(s1,4)*
                 (4 + 3*Power(s2,2) + 5*Power(t1,2) + 
                   s2*(2 - 8*t1 - 5*t2) - 13*t2 + 9*Power(t2,2) + 
                   t1*(-1 + 4*t2)) - 
                Power(s1,3)*
                 (6 - 4*t2 - 44*Power(t2,2) + 42*Power(t2,3) + 
                   12*Power(s2,2)*(-1 + 2*t2) + 
                   8*Power(t1,2)*(-2 + 3*t2) + 
                   s2*(5 + t1*(28 - 48*t2) + 25*t2 - 
                      46*Power(t2,2)) + 
                   t1*(-19 + 17*t2 + 18*Power(t2,2))) - 
                t2*(-44 + 103*t2 - 56*Power(t2,2) + 
                   2*Power(s2,2)*(8 - 3*t2)*Power(t2,2) - 
                   5*Power(t2,3) + 2*Power(t2,4) + 
                   s2*(16 + 76*t2 - 127*Power(t2,2) + 
                      3*Power(t2,3) + 4*Power(t2,4))) + 
                Power(s1,2)*
                 (-17 + 102*t2 - 99*Power(t2,2) - 16*Power(t2,3) + 
                   30*Power(t2,4) + 4*Power(s2,2)*t2*(-10 + 9*t2) + 
                   2*Power(t1,2)*(3 - 25*t2 + 9*Power(t2,2)) + 
                   t1*(43 - 155*t2 + 74*Power(t2,2) + 
                      12*Power(t2,3)) + 
                   s2*(-28 + 81*t2 + 39*Power(t2,2) - 
                      66*Power(t2,3) + 
                      t1*(2 + 76*t2 - 48*Power(t2,2)))) - 
                s1*(12 - 12*(7 + 10*s2)*t2 + 
                   (148 + 203*s2 - 44*Power(s2,2))*Power(t2,2) + 
                   (-82 + 15*s2 + 24*Power(s2,2))*Power(t2,3) - 
                   34*s2*Power(t2,4) + 6*Power(t2,5) + 
                   2*Power(t1,2)*
                    (7 - 3*t2 - 14*Power(t2,2) + 2*Power(t2,3)) + 
                   t1*(-6 + 126*t2 - 247*Power(t2,2) + 
                      61*Power(t2,3) + 2*Power(t2,4) - 
                      4*s2*(4 - 10*t2 - 7*Power(t2,2) + 4*Power(t2,3))\
))) - 2*s*(-4 + 8*t1 - 4*Power(t1,2) + 18*t2 - 32*t1*t2 + 
                14*Power(t1,2)*t2 - 29*Power(t2,2) - 
                3*s2*Power(t2,2) - 4*Power(s2,2)*Power(t2,2) + 
                29*t1*Power(t2,2) + 6*s2*t1*Power(t2,2) - 
                10*Power(t1,2)*Power(t2,2) + 12*Power(t2,3) + 
                36*s2*Power(t2,3) - 21*t1*Power(t2,3) - 
                15*s2*t1*Power(t2,3) + 4*Power(t1,2)*Power(t2,3) + 
                10*Power(t2,4) - 33*s2*Power(t2,4) + 
                7*Power(s2,2)*Power(t2,4) + 13*t1*Power(t2,4) + 
                5*s2*t1*Power(t2,4) - 2*Power(t1,2)*Power(t2,4) - 
                6*Power(t2,5) - Power(s2,2)*Power(t2,5) + 
                3*t1*Power(t2,5) - Power(t2,6) + 
                Power(s1,5)*(Power(s2,2) + Power(t1,2) + 
                   t1*(-1 + t2) + 2*Power(-1 + t2,2) - 
                   s2*(-1 + 2*t1 + t2)) + 
                s1*(-8 + 2*(13 + 9*s2 + 4*Power(s2,2))*t2 + 
                   (-24 - 95*s2 + 2*Power(s2,2))*Power(t2,2) + 
                   (7 + 75*s2 - 24*Power(s2,2))*Power(t2,3) + 
                   (-6 + 5*s2 + 5*Power(s2,2))*Power(t2,4) + 
                   (5 - 3*s2)*Power(t2,5) + 
                   Power(t1,2)*
                    (-4 + 12*t2 - 18*Power(t2,2) + Power(t2,3)) + 
                   t1*(12 - 6*(3 + 4*s2)*t2 + 
                      (58 + 39*s2)*Power(t2,2) + 
                      (-52 + 5*s2)*Power(t2,3) - 2*s2*Power(t2,4))) + 
                Power(s1,4)*(Power(s2,2)*(3 - 5*t2) + 
                   Power(t1,2)*(1 - 3*t2) - 
                   3*Power(-1 + t2,2)*(1 + 2*t2) + 
                   t1*(2 + t2 - 3*Power(t2,2)) + 
                   s2*(1 - 7*t2 + 6*Power(t2,2) + t1*(-4 + 8*t2))) + 
                Power(s1,3)*(3*Power(-1 + t2,2)*
                    (-4 + 3*t2 + 2*Power(t2,2)) + 
                   Power(t1,2)*(2 - 9*t2 + 3*Power(t2,2)) + 
                   2*Power(s2,2)*(1 - 8*t2 + 5*Power(t2,2)) + 
                   t1*(16 - 26*t2 + 7*Power(t2,2) + 3*Power(t2,3)) + 
                   s2*(-11 + 7*t2 + 16*Power(t2,2) - 12*Power(t2,3) + 
                      t1*(-3 + 23*t2 - 12*Power(t2,2)))) - 
                Power(s1,2)*(Power(t1,2)*
                    (10 - 14*t2 - 9*Power(t2,2) + Power(t2,3)) + 
                   Power(-1 + t2,2)*
                    (5 - 16*t2 + 10*Power(t2,2) + 2*Power(t2,3)) + 
                   2*Power(s2,2)*
                    (2 + 2*t2 - 15*Power(t2,2) + 5*Power(t2,3)) + 
                   t1*(-5 + 57*t2 - 63*Power(t2,2) + 10*Power(t2,3) + 
                      Power(t2,4)) + 
                   s2*(t1*(-18 + 21*t2 + 29*Power(t2,2) - 
                        8*Power(t2,3)) + 
                      5*(3 - 14*t2 + 10*Power(t2,2) + 3*Power(t2,3) - 
                        2*Power(t2,4))))))/
           (s*(-1 + s2)*(-1 + t1)*
             Power(s1 - s1*t2 + t2*(-1 + s + t2),2)) - 
          (-2*(-3 + s1)*(s1 - t2)*(-1 + t2)*
              Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
             Power(s,6)*(Power(t1,2) + Power(t2,2)) + 
             Power(s,5)*(Power(t1,2)*(-3 - 5*s1 + 3*t2) + 
                t1*(2 + s1*(-2 + 2*s2 - t2) - s2*t2 + Power(t2,2)) + 
                t2*(2 + s2*(-2 + 2*s1 - 3*t2) + t2 - 6*s1*t2 + 
                   2*Power(t2,2))) + 
             Power(s,4)*(1 - 4*t2 + s2*t2 - 27*Power(t2,2) + 
                2*s2*Power(t2,2) + 2*Power(s2,2)*Power(t2,2) + 
                9*Power(t2,3) - 7*s2*Power(t2,3) + Power(t2,4) + 
                Power(t1,2)*(-5 - 7*t2 + 3*Power(t2,2)) + 
                Power(s1,2)*(Power(s2,2) + 5*t1 - 8*s2*t1 + 
                   10*Power(t1,2) - 5*t2 - 7*s2*t2 + 4*t1*t2 + 
                   14*Power(t2,2)) + 
                t1*(-4 + t2 + 6*s2*t2 + (2 - 4*s2)*Power(t2,2) + 
                   2*Power(t2,3)) + 
                s1*(Power(t1,2)*(9 - 12*t2) - 
                   t1*(4 + s2*(3 - 11*t2) + 5*t2 + 7*Power(t2,2)) + 
                   t2*(1 - 3*Power(s2,2) - 6*t2 - 10*Power(t2,2) + 
                      3*s2*(1 + 6*t2)))) - 
             Power(s,3)*(1 + 18*t2 - 21*s2*t2 - 74*Power(t2,2) - 
                14*s2*Power(t2,2) + 6*Power(s2,2)*Power(t2,2) + 
                70*Power(t2,3) - 5*Power(s2,2)*Power(t2,3) - 
                11*Power(t2,4) + 5*s2*Power(t2,4) + 
                Power(s1,3)*(2 + s2 + 3*Power(s2,2) + 3*t1 - 
                   12*s2*t1 + 10*Power(t1,2) - 14*t2 - 9*s2*t2 + 
                   6*t1*t2 + 16*Power(t2,2)) - 
                Power(t1,2)*
                 (27 - 21*t2 - 8*Power(t2,2) + Power(t2,3)) - 
                t1*(-14 + (16 - 13*s2)*t2 + 
                   (-15 + 17*s2)*Power(t2,2) + 
                   (11 - 5*s2)*Power(t2,3) + Power(t2,4)) - 
                Power(s1,2)*(3 - 16*t2 + 5*Power(t2,2) + 
                   18*Power(t2,3) + 9*Power(t1,2)*(-1 + 2*t2) + 
                   Power(s2,2)*(-2 + 11*t2) + 
                   s2*(t1*(8 - 27*t2) + 2*(4 - 17*t2)*t2) + 
                   t1*(1 + 7*t2 + 15*Power(t2,2))) + 
                s1*(3 + 29*t2 - 109*Power(t2,2) + 27*Power(t2,3) + 
                   4*Power(t2,4) + Power(s2,2)*t2*(-8 + 13*t2) + 
                   Power(t1,2)*(-15 - 16*t2 + 9*Power(t2,2)) + 
                   t1*(-19 - 21*t2 + 20*Power(t2,2) + 9*Power(t2,3)) + 
                   s2*(1 + 22*t2 + 9*Power(t2,2) - 28*Power(t2,3) + 
                      t1*(11 + 21*t2 - 20*Power(t2,2))))) + 
             Power(s,2)*(-8 + 58*t2 - 32*s2*t2 - 116*Power(t2,2) + 
                14*s2*Power(t2,2) + 103*Power(t2,3) + 
                48*s2*Power(t2,3) - 10*Power(s2,2)*Power(t2,3) - 
                38*Power(t2,4) - 10*s2*Power(t2,4) + 
                4*Power(s2,2)*Power(t2,4) + Power(t2,5) - 
                s2*Power(t2,5) - 
                Power(t1,2)*(32 - 65*t2 + 17*Power(t2,2) + 
                   8*Power(t2,3)) + 
                t1*(40 + 3*(-29 + 4*s2)*t2 + 
                   (47 - 27*s2)*Power(t2,2) + 
                   (-34 + 15*s2)*Power(t2,3) + (15 - 2*s2)*Power(t2,4)\
) + Power(s1,4)*(4 + 3*Power(s2,2) + 5*Power(t1,2) + 
                   s2*(2 - 8*t1 - 5*t2) - 13*t2 + 9*Power(t2,2) + 
                   t1*(-1 + 4*t2)) + 
                Power(s1,3)*(-8 + Power(s2,2)*(4 - 13*t2) + 
                   Power(t1,2)*(3 - 12*t2) + 19*t2 + 3*Power(t2,2) - 
                   14*Power(t2,3) + t1*(2 + t2 - 13*Power(t2,2)) + 
                   s2*(2 - 18*t2 + 26*Power(t2,2) + t1*(-7 + 25*t2))) \
+ Power(s1,2)*(-11 + 117*t2 - 135*Power(t2,2) + 24*Power(t2,3) + 
                   5*Power(t2,4) + 
                   Power(t1,2)*(-13 - 13*t2 + 9*Power(t2,2)) + 
                   Power(s2,2)*(-4 - 18*t2 + 21*Power(t2,2)) + 
                   t1*(8 - 65*t2 + 29*Power(t2,2) + 12*Power(t2,3)) + 
                   s2*(-6 + 34*t2 + 22*Power(t2,2) - 34*Power(t2,3) + 
                      t1*(19 + 27*t2 - 28*Power(t2,2)))) + 
                s1*(-4 + 65*t2 - 194*Power(t2,2) + 148*Power(t2,3) - 
                   15*Power(t2,4) + 
                   Power(s2,2)*t2*(4 + 24*t2 - 15*Power(t2,2)) + 
                   Power(t1,2)*
                    (-33 + 34*t2 + 14*Power(t2,2) - 2*Power(t2,3)) + 
                   t1*(1 + 3*t2 + 83*Power(t2,2) - 42*Power(t2,3) - 
                      3*Power(t2,4)) + 
                   s2*(4 + 28*t2 - 92*Power(t2,2) + 4*Power(t2,3) + 
                      14*Power(t2,4) + 
                      t1*(16 - 28*t2 - 27*Power(t2,2) + 13*Power(t2,3))\
))) + s*(12 - 24*t1 + 12*Power(t1,2) - 52*t2 + 4*s2*t2 + 98*t1*t2 - 
                4*s2*t1*t2 - 46*Power(t1,2)*t2 + 87*Power(t2,2) - 
                32*s2*Power(t2,2) + 8*Power(s2,2)*Power(t2,2) - 
                87*t1*Power(t2,2) + 18*s2*t1*Power(t2,2) + 
                36*Power(t1,2)*Power(t2,2) - 69*Power(t2,3) + 
                3*s2*Power(t2,3) - 8*Power(s2,2)*Power(t2,3) + 
                25*t1*Power(t2,3) - 12*s2*t1*Power(t2,3) - 
                Power(t1,2)*Power(t2,3) + 23*Power(t2,4) + 
                31*s2*Power(t2,4) - 4*Power(s2,2)*Power(t2,4) - 
                18*t1*Power(t2,4) + 4*s2*t1*Power(t2,4) - 
                4*Power(t1,2)*Power(t2,4) + Power(t2,5) - 
                6*s2*Power(t2,5) + Power(s2,2)*Power(t2,5) + 
                6*t1*Power(t2,5) - 2*Power(t2,6) - 
                Power(s1,5)*(Power(s2,2) + Power(t1,2) + 
                   t1*(-1 + t2) + 2*Power(-1 + t2,2) - 
                   s2*(-1 + 2*t1 + t2)) + 
                Power(s1,4)*(3*Power(t1,2)*t2 + 
                   Power(-1 + t2,2)*(5 + 4*t2) + 
                   Power(s2,2)*(-2 + 5*t2) + 
                   s2*(-2 + t1*(2 - 8*t2) + 9*t2 - 7*Power(t2,2)) + 
                   t1*(-1 - 3*t2 + 4*Power(t2,2))) + 
                Power(s1,3)*(Power(s2,2)*
                    (2 + 12*t2 - 10*Power(t2,2)) + 
                   Power(t1,2)*(1 + 6*t2 - 3*Power(t2,2)) - 
                   Power(-1 + t2,2)*(-31 + 10*t2 + 2*Power(t2,2)) - 
                   t1*(27 - 43*t2 + 11*Power(t2,2) + 5*Power(t2,3)) + 
                   s2*(17 - 15*t2 - 15*Power(t2,2) + 13*Power(t2,3) + 
                      4*t1*(-1 - 4*t2 + 3*Power(t2,2)))) + 
                Power(s1,2)*(Power(-1 + t2,2)*
                    (11 - 63*t2 + 3*Power(t2,2)) + 
                   Power(t1,2)*
                    (12 - 21*t2 - 8*Power(t2,2) + Power(t2,3)) + 
                   2*Power(s2,2)*
                    (4 - 6*t2 - 11*Power(t2,2) + 5*Power(t2,3)) + 
                   t1*(13 + 51*t2 - 97*Power(t2,2) + 31*Power(t2,3) + 
                      2*Power(t2,4)) + 
                   s2*(12 - 83*t2 + 79*Power(t2,2) + Power(t2,3) - 
                      9*Power(t2,4) + 
                      t1*(-26 + 48*t2 + 18*Power(t2,2) - 8*Power(t2,3))\
)) + s1*(Power(s2,2)*t2*(-16 + 18*t2 + 16*Power(t2,2) - 
                      5*Power(t2,3)) + 
                   Power(-1 + t2,2)*
                    (18 - 28*t2 + 35*Power(t2,2) + 4*Power(t2,3)) + 
                   Power(t1,2)*
                    (12 - 14*t2 + 7*Power(t2,2) + 8*Power(t2,3)) - 
                   3*t1*(10 - 2*t2 + 7*Power(t2,2) - 23*Power(t2,3) + 
                      8*Power(t2,4)) + 
                   s2*(-4 + 20*t2 + 63*Power(t2,2) - 93*Power(t2,3) + 
                      12*Power(t2,4) + 2*Power(t2,5) + 
                      2*t1*(2 + 4*t2 - 16*Power(t2,2) - 4*Power(t2,3) + 
                         Power(t2,4))))))/
           (s*(-1 + s2)*(s - s2 + t1)*
             Power(s1 - s1*t2 + t2*(-1 + s + t2),2))) + 
       ((2*(((-1 + s1)*(-3 + 2*Power(s,2) + Power(s1,2) - 
                    s1*(-3 + t2) + s*(-4 - 3*s1 + 2*t2))*
                  (2 - 2*t1 + s1*(1 + t1 - 2*t2)*(-1 + t2) - t2 + 
                    3*t1*t2 - Power(t2,2) - s2*Power(t2,2) + 
                    Power(s,2)*(-t1 + t2) + 
                    Power(s1,2)*(-2 + s2 - t1 + 2*t2) - 
                    s*(1 + t1*(-3 + t2) + t2 + s2*t2 - Power(t2,2) + 
                       s1*(-2 + s2 - 2*t1 + 3*t2))))/
                (s1*(-1 + t2) - t2*(-1 + s + t2)) + 
               2*(30 + 4*Power(s,3) - Power(s1,4) - 5*s2 + 
                  Power(s1,2)*(18 + 2*t1 + s2*(-2 + t2) - 7*t2) - 
                  8*t2 + 2*t1*t2 - 6*Power(t2,2) + 
                  Power(s1,3)*(2 - s2 + t2) - 
                  s1*(-3 + s2*(-8 + t2) + 2*t2 - 2*Power(t2,2) + 
                     2*t1*(1 + t2)) - 
                  Power(s,2)*
                   (13 + Power(s1,2) - 2*s2 + t1 - 7*t2 + 
                     s1*(2*s2 - t1 + t2)) + 
                  s*(-4 + 2*Power(s1,3) + 
                     Power(s1,2)*(-6 + 3*s2 - t1) + 3*t1 - 14*t2 - 
                     t1*t2 + 3*Power(t2,2) + s2*(-3 + 2*t2) + 
                     s1*(-10 + t1*(-2 + t2) + 6*t2 - 2*s2*t2 - 
                        Power(t2,2))))))/(-1 + t2) + 
          (2*(((-1 + s1)*(4 - Power(s1,2) + s*(1 + s1) + 
                    s1*(-2 + t2) - 2*t2)*
                  (-2 + 2*t1 + Power(s1,2)*(2 - s2 + t1 - 2*t2) + 
                    Power(s,2)*(t1 - t2) - 
                    s1*(1 + t1 - 2*t2)*(-1 + t2) + t2 - 3*t1*t2 + 
                    Power(t2,2) + s2*Power(t2,2) + 
                    s*(1 + t1*(-3 + t2) + t2 + s2*t2 - Power(t2,2) + 
                       s1*(-2 + s2 - 2*t1 + 3*t2))))/
                (s1*(-1 + t2) - t2*(-1 + s + t2)) - 
               2*(-27 + Power(s1,4) + 4*s2 + 3*t1 + 
                  Power(s1,3)*(-1 + s2 - t2) + 3*t2 - 3*s2*t2 + 
                  2*t1*t2 - 6*Power(t2,2) + 
                  Power(s,2)*
                   (-3 + Power(s1,2) - t1 + s1*(-4 + t1 - t2) + 3*t2) \
+ Power(s1,2)*(-20 + t1 - (-3 + s2)*t2) + 
                  s1*(3 + 11*t2 + 2*Power(t2,2) - 2*t1*(2 + t2) + 
                     s2*(-5 + 4*t2)) - 
                  s*(-15 + 2*Power(s1,3) + s2 + 
                     Power(s1,2)*(-7 + s2 + t1 - 2*t2) + 7*t2 + 
                     t1*t2 - 3*Power(t2,2) - 
                     s1*(14 + 2*s2 + t1 - 7*t2 + t1*t2 - Power(t2,2))))\
))/(-s + s1 - t2) + (2*(-1 + s1)*
             (((Power(s,3) + Power(s,2)*(-1 - 2*s1 + t2) + 
                    s*(-7 + s1 + Power(s1,2) + 2*t2 - s1*t2) + 
                    2*(8 + s1 - 6*t2 + Power(t2,2)))*
                  (2 - 2*t1 + s1*(1 + t1 - 2*t2)*(-1 + t2) - t2 + 
                    3*t1*t2 - Power(t2,2) - s2*Power(t2,2) + 
                    Power(s,2)*(-t1 + t2) + 
                    Power(s1,2)*(-2 + s2 - t1 + 2*t2) - 
                    s*(1 + t1*(-3 + t2) + t2 + s2*t2 - Power(t2,2) + 
                       s1*(-2 + s2 - 2*t1 + 3*t2))))/
                (s1 - s1*t2 + t2*(-1 + s + t2)) + 
               2*(Power(s,3)*(-3 + s1) + 
                  Power(s,2)*
                   (11 - 2*Power(s1,2) + 2*t1 - 8*t2 + s1*(4 + t2)) + 
                  s*(-14 + Power(s1,3) - 12*t1 - 
                     Power(s1,2)*(-1 + t2) + 26*t2 + 2*t1*t2 - 
                     5*Power(t2,2) + s2*(-1 + 2*t2) + 
                     s1*(-33 + 5*s2 - 2*t1 + 7*t2)) - 
                  2*(-8 + Power(s1,2)*(-5 + 2*s2 - t1) - 8*t1 + 4*t2 + 
                     t1*t2 - 3*Power(t2,2) + 
                     s2*(2 + 5*t2 - Power(t2,2)) - 
                     s1*(16 + t1 + 2*s2*(-1 + t2) - 13*t2 - t1*t2 + 
                        Power(t2,2))))))/((s - s2 + t1)*(s - s1 + t2)) + 
          (2*(-1 + s1)*(((Power(s,2)*(-2 + s1 - t2) + 
                    2*(5 + s1*(-2 + t2) + t2 - Power(t2,2)) - 
                    s*(-1 + Power(s1,2) + 2*t2 + Power(t2,2) - 
                       s1*(1 + 2*t2)))*
                  (2 - 2*t1 + s1*(1 + t1 - 2*t2)*(-1 + t2) - t2 + 
                    3*t1*t2 - Power(t2,2) - s2*Power(t2,2) + 
                    Power(s,2)*(-t1 + t2) + 
                    Power(s1,2)*(-2 + s2 - t1 + 2*t2) - 
                    s*(1 + t1*(-3 + t2) + t2 + s2*t2 - Power(t2,2) + 
                       s1*(-2 + s2 - 2*t1 + 3*t2))))/
                (s1 - s1*t2 + t2*(-1 + s + t2)) - 
               2*(-36 + 2*s2 + 8*t1 + Power(s1,2)*(s2 + 3*t1 - 2*t2) + 
                  Power(s,3)*(-2 + s1 + t1 - t2) + 19*t2 - 10*s2*t2 - 
                  t1*t2 + 7*Power(t2,2) + s2*Power(t2,2) + 
                  s1*(-29 + t1 - 2*s2*(-6 + t2) - 5*t2 - 3*t1*t2 + 
                     2*Power(t2,2)) - 
                  Power(s,2)*
                   (8 + 2*Power(s1,2) + 3*t1 + 
                     s1*(-4 + s2 + t1 - 3*t2) + 2*t2 - t1*t2 + 
                     Power(t2,2) - s2*(5 + t2)) + 
                  s*(50 + Power(s1,3) - 2*t1 - 4*t2 - t1*t2 - 
                     4*Power(t2,2) + s2*(-13 + 6*t2 + Power(t2,2)) + 
                     Power(s1,2)*(s2 - 2*(1 + t2)) + 
                     s1*(5 - t1 + 7*t2 + Power(t2,2) - s2*(5 + 2*t2))))\
))/((-1 + t1)*(-1 + t2)) + (2*(-1 + s1)*
             (-(((23 + Power(s,2)*(-3 + s1) + 2*Power(s1,3) + 
                      3*s2 - 10*t1 + 2*Power(s1,2)*(-4 + s2 - t2) + 
                      t2 + 7*s2*t2 - 2*t1*t2 + 
                      2*s1*(4 + t1 + 4*t2 - s2*(4 + t2)) + 
                      s*(-17 - 3*Power(s1,2) + 6*s2 + 2*t1 - 6*t2 + 
                       s1*(13 - 2*s2 + t2)))*
                    (-2 + 2*t1 + Power(s1,2)*(2 - s2 + t1 - 2*t2) + 
                      Power(s,2)*(t1 - t2) - 
                      s1*(1 + t1 - 2*t2)*(-1 + t2) + t2 - 3*t1*t2 + 
                      Power(t2,2) + s2*Power(t2,2) + 
                      s*(1 + t1*(-3 + t2) + t2 + s2*t2 - 
                        Power(t2,2) + s1*(-2 + s2 - 2*t1 + 3*t2))))/
                  (s1 - s1*t2 + t2*(-1 + s + t2))) + 
               2*(11 + Power(s1,4) + 13*s2 + Power(s2,2) - 23*t1 - 
                  4*s2*t1 + 6*Power(t1,2) + 
                  Power(s1,3)*(-4 + 2*s2 - t2) + 19*t2 - 4*s2*t2 + 
                  3*Power(s2,2)*t2 - 10*t1*t2 - 2*s2*t1*t2 + 
                  4*Power(t2,2) + 
                  Power(s,2)*
                   (14 + Power(s1,2) + s1*(-5 + s2) - 4*s2 - t1 + 
                     2*t2) + 
                  Power(s1,2)*
                   (6 + Power(s2,2) + t1 + 4*t2 - s2*(9 + 2*t2)) + 
                  s1*(-(Power(s2,2)*(4 + t2)) + 
                     2*s2*(2 + t1 + 4*t2) + 
                     2*(Power(t1,2) - 2*(-6 + t2) - t1*(5 + t2))) - 
                  s*(33 + 2*Power(s1,3) + 13*s2 - 3*Power(s2,2) - 
                     21*t1 + 4*Power(t1,2) + 
                     Power(s1,2)*(-9 + 3*s2 - t2) + 14*t2 + 4*s2*t2 - 
                     7*t1*t2 + 2*Power(t2,2) + 
                     s1*(16 + Power(s2,2) + 2*t1 + 3*t2 - 
                        s2*(15 + t2)))) - 
               ((s*(-2 + s1) - Power(s1,2) - 4*(1 + t2) + 
                    s1*(4 + t2))*
                  (Power(s,5)*Power(t1 - t2,2) + 
                    2*(s1 - t2)*(-1 + t2)*
                     Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                    2*Power(s,4)*
                     (Power(t1,2)*(3 + 2*s1 - t2) + 
                       t2*(-((-2 + t2)*t2) + s2*(1 + t2) + 
                       s1*(-1 + 2*t2)) - 
                       t1*(1 + (5 + s2)*t2 - 2*Power(t2,2) + 
                        s1*(-2 + s2 + 4*t2))) + 
                    Power(s,3)*
                     (1 - 2*t2 + 10*s2*t2 + 7*Power(t2,2) + 
                       Power(s2,2)*Power(t2,2) - 8*Power(t2,3) - 
                       4*s2*Power(t2,3) + Power(t2,4) + 
                       Power(t1,2)*(13 - 12*t2 + Power(t2,2)) - 
                       2*t1*
                       (5 + (6 + 4*s2)*t2 - 
                       2*(5 + s2)*Power(t2,2) + Power(t2,3)) + 
                       Power(s1,2)*
                       (2 + Power(s2,2) + 6*Power(t1,2) - 8*t2 + 
                       7*Power(t2,2) + s2*(-2 - 6*t1 + 2*t2) - 
                       2*t1*(-5 + 6*t2)) - 
                       2*s1*
                        (1 + t2 - 3*(2 + s2)*Power(t2,2) + 
                        3*Power(t2,3) + Power(t1,2)*(-7 + 3*t2) + 
                        t1*(-6 + 3*s2 + 17*t2 - 6*Power(t2,2)))) + 
                    s*(4 - 8*t1 + 4*Power(t1,2) + 
                       Power(s1,4)*
                       (Power(s2,2) + Power(t1,2) - 
                       2*s2*(1 + t1 - t2) - 2*t1*(-1 + t2) + 
                       2*Power(-1 + t2,2)) - 10*t2 + 4*s2*t2 + 
                       24*t1*t2 - 4*s2*t1*t2 - 14*Power(t1,2)*t2 + 
                       13*Power(t2,2) - 22*s2*Power(t2,2) + 
                       2*Power(s2,2)*Power(t2,2) - 
                       6*t1*Power(t2,2) + 14*s2*t1*Power(t2,2) + 
                       9*Power(t1,2)*Power(t2,2) - 12*Power(t2,3) + 
                       14*s2*Power(t2,3) - 
                       4*Power(s2,2)*Power(t2,3) - 
                       10*t1*Power(t2,3) - 8*s2*t1*Power(t2,3) + 
                       5*Power(t2,4) + 4*s2*Power(t2,4) + 
                       Power(s2,2)*Power(t2,4) - 
                       2*Power(s1,3)*(-1 + t2)*
                       (1 + Power(s2,2) + Power(t1,2) - 
                       2*t1*(-1 + t2) - 3*t2 + 2*Power(t2,2) + 
                       s2*(-1 - 2*t1 + t2)) + 
                       Power(s1,2)*
                       (-2*t1*(-7 + t2)*Power(-1 + t2,2) + 
                       Power(t1,2)*(5 - 10*t2 + Power(t2,2)) + 
                       2*Power(s2,2)*(1 - 4*t2 + Power(t2,2)) + 
                       Power(-1 + t2,2)*
                       (-3 - 4*t2 + 2*Power(t2,2)) - 
                       2*s2*
                       (Power(-1 + t2,2)*(1 + t2) + 
                       t1*(3 - 8*t2 + Power(t2,2)))) - 
                       2*s1*
                        (-(Power(-1 + t2,3)*t2) + 
                        Power(t1,2)*(-2 + 3*t2 - 3*Power(t2,2)) + 
                        Power(s2,2)*t2*(2 - 5*t2 + Power(t2,2)) + 
                        2*t1*
                       (1 + 6*t2 - 10*Power(t2,2) + 3*Power(t2,3)) + 
                        s2*(2 - 12*t2 + 7*Power(t2,2) + 
                        4*Power(t2,3) - Power(t2,4) + 
                        2*t1*(-1 + 2*t2 + Power(t2,2))))) - 
                    2*Power(s,2)*
                     (2 - 4*t2 + 7*s2*t2 + 6*Power(t2,2) - 
                       10*s2*Power(t2,2) + Power(s2,2)*Power(t2,2) - 
                       6*Power(t2,3) - 3*s2*Power(t2,3) - 
                       Power(s2,2)*Power(t2,3) + 2*Power(t2,4) + 
                       s2*Power(t2,4) + 
                       Power(t1,2)*(6 - 11*t2 + 3*Power(t2,2)) - 
                       t1*(8 + (-6 + 5*s2)*t2 - 
                        4*(3 + 2*s2)*Power(t2,2) + 
                        (5 + s2)*Power(t2,3)) + 
                       Power(s1,3)*
                        (2 + Power(s2,2) + 2*Power(t1,2) - 
                        4*t1*(-1 + t2) - 5*t2 + 3*Power(t2,2) + 
                        s2*(-2 - 3*t1 + 2*t2)) + 
                       Power(s1,2)*
                        (Power(t1,2)*(5 - 3*t2) - 
                        Power(s2,2)*(-1 + t2) + 
                        t2*(-5 + 9*t2 - 4*Power(t2,2)) + 
                        t1*(9 - 16*t2 + 6*Power(t2,2)) + 
                        s2*(-1 + t2 + Power(t2,2) + t1*(-5 + 3*t2))) \
+ s1*(-2 + 2*t2 + Power(s2,2)*(-2 + t2)*t2 + 5*Power(t2,2) - 
                        6*Power(t2,3) + Power(t2,4) + 
                        Power(t1,2)*(7 - 10*t2 + Power(t2,2)) + 
                        t1*(4 - 27*t2 + 21*Power(t2,2) - 
                        2*Power(t2,3)) + 
                        s2*(-1 + 5*t2 + 4*Power(t2,2) - 
                        4*Power(t2,3) + t1*(-1 + 3*t2 + Power(t2,2))))\
)))/(s*Power(s1 - s1*t2 + t2*(-1 + s + t2),2))))/((-1 + s2)*(-1 + t1)) + 
          ((-1 + s1)*((-2*(23 + Power(s,3) + 4*s2 - 5*t1 + 
                    s1*(16 + 2*t1 + s2*(-6 + t2) - 3*t2) - 
                    Power(s,2)*(1 + 3*s1 + s2 - t2) - 6*t2 - 
                    2*s2*t2 + 2*t1*t2 + 
                    Power(s1,2)*(-4 + s2 - t1 + 2*t2) + 
                    s*(-15 + 2*Power(s1,2) + 2*s2 + 
                       s1*(3 + t1 - 3*t2) + 3*t2 - s2*t2))*
                  (-2 + 2*t1 + Power(s1,2)*(2 - s2 + t1 - 2*t2) + 
                    Power(s,2)*(t1 - t2) - 
                    s1*(1 + t1 - 2*t2)*(-1 + t2) + t2 - 3*t1*t2 + 
                    Power(t2,2) + s2*Power(t2,2) + 
                    s*(1 + t1*(-3 + t2) + t2 + s2*t2 - Power(t2,2) + 
                       s1*(-2 + s2 - 2*t1 + 3*t2))))/
                (s1 - s1*t2 + t2*(-1 + s + t2)) + 
               4*(12 + Power(s,3)*(-1 + s1) + 4*s2 + Power(s2,2) - 
                  6*t1 + 4*s2*t1 - 6*Power(t1,2) + 
                  Power(s,2)*
                   (11 - 2*Power(s1,2) + t1 + s1*(-2 + t2) - 3*t2) - 
                  11*s2*t2 - Power(s2,2)*t2 + 10*t1*t2 + 2*s2*t1*t2 - 
                  4*Power(t2,2) + Power(s1,3)*(-2 + s2 - t1 + t2) + 
                  Power(s1,2)*
                   (7 + Power(s2,2) + 3*t1 - t2 + s2*(-7 - t1 + t2)) - 
                  s1*(-32 + 5*Power(s2,2) + 2*Power(t1,2) - 
                     2*t1*(-6 + t2) + t2 + s2*(-6 - 7*t1 + 2*t2)) + 
                  s*(-22 + Power(s1,3) + 2*Power(s2,2) + 4*t1 + 
                     4*Power(t1,2) + 
                     Power(s1,2)*(5 - s2 + t1 - 2*t2) + 7*t2 - 
                     7*t1*t2 + 2*Power(t2,2) + s2*(-7 - 5*t1 + 7*t2) - 
                     s1*(22 + Power(s2,2) + s2*(-6 - t1 + t2)))) - 
               (2*(-7 + Power(s,2) - s1*(-2 + t2) + 2*t2 + 
                    s*(-s1 + t2))*
                  (Power(s,5)*Power(t1 - t2,2) + 
                    2*(s1 - t2)*(-1 + t2)*
                     Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                    2*Power(s,4)*
                     (Power(t1,2)*(3 + 2*s1 - t2) + 
                       t2*(-((-2 + t2)*t2) + s2*(1 + t2) + 
                        s1*(-1 + 2*t2)) - 
                       t1*(1 + (5 + s2)*t2 - 2*Power(t2,2) + 
                        s1*(-2 + s2 + 4*t2))) + 
                    Power(s,3)*
                     (1 - 2*t2 + 10*s2*t2 + 7*Power(t2,2) + 
                       Power(s2,2)*Power(t2,2) - 8*Power(t2,3) - 
                       4*s2*Power(t2,3) + Power(t2,4) + 
                       Power(t1,2)*(13 - 12*t2 + Power(t2,2)) - 
                       2*t1*
                        (5 + (6 + 4*s2)*t2 - 
                        2*(5 + s2)*Power(t2,2) + Power(t2,3)) + 
                       Power(s1,2)*
                        (2 + Power(s2,2) + 6*Power(t1,2) - 8*t2 + 
                        7*Power(t2,2) + s2*(-2 - 6*t1 + 2*t2) - 
                        2*t1*(-5 + 6*t2)) - 
                       2*s1*(1 + t2 - 3*(2 + s2)*Power(t2,2) + 
                        3*Power(t2,3) + Power(t1,2)*(-7 + 3*t2) + 
                        t1*(-6 + 3*s2 + 17*t2 - 6*Power(t2,2)))) + 
                    s*(4 - 8*t1 + 4*Power(t1,2) + 
                       Power(s1,4)*
                        (Power(s2,2) + Power(t1,2) - 
                        2*s2*(1 + t1 - t2) - 2*t1*(-1 + t2) + 
                        2*Power(-1 + t2,2)) - 10*t2 + 4*s2*t2 + 
                       24*t1*t2 - 4*s2*t1*t2 - 14*Power(t1,2)*t2 + 
                       13*Power(t2,2) - 22*s2*Power(t2,2) + 
                       2*Power(s2,2)*Power(t2,2) - 6*t1*Power(t2,2) + 
                       14*s2*t1*Power(t2,2) + 
                       9*Power(t1,2)*Power(t2,2) - 12*Power(t2,3) + 
                       14*s2*Power(t2,3) - 
                       4*Power(s2,2)*Power(t2,3) - 
                       10*t1*Power(t2,3) - 8*s2*t1*Power(t2,3) + 
                       5*Power(t2,4) + 4*s2*Power(t2,4) + 
                       Power(s2,2)*Power(t2,4) - 
                       2*Power(s1,3)*(-1 + t2)*
                        (1 + Power(s2,2) + Power(t1,2) - 
                        2*t1*(-1 + t2) - 3*t2 + 2*Power(t2,2) + 
                        s2*(-1 - 2*t1 + t2)) + 
                       Power(s1,2)*
                        (-2*t1*(-7 + t2)*Power(-1 + t2,2) + 
                        Power(t1,2)*(5 - 10*t2 + Power(t2,2)) + 
                        2*Power(s2,2)*(1 - 4*t2 + Power(t2,2)) + 
                        Power(-1 + t2,2)*
                       (-3 - 4*t2 + 2*Power(t2,2)) - 
                        2*s2*
                        (Power(-1 + t2,2)*(1 + t2) + 
                        t1*(3 - 8*t2 + Power(t2,2)))) - 
                       2*s1*(-(Power(-1 + t2,3)*t2) + 
                        Power(t1,2)*(-2 + 3*t2 - 3*Power(t2,2)) + 
                        Power(s2,2)*t2*(2 - 5*t2 + Power(t2,2)) + 
                        2*t1*
                        (1 + 6*t2 - 10*Power(t2,2) + 3*Power(t2,3)) + 
                        s2*(2 - 12*t2 + 7*Power(t2,2) + 
                        4*Power(t2,3) - Power(t2,4) + 
                        2*t1*(-1 + 2*t2 + Power(t2,2))))) - 
                    2*Power(s,2)*
                     (2 - 4*t2 + 7*s2*t2 + 6*Power(t2,2) - 
                       10*s2*Power(t2,2) + Power(s2,2)*Power(t2,2) - 
                       6*Power(t2,3) - 3*s2*Power(t2,3) - 
                       Power(s2,2)*Power(t2,3) + 2*Power(t2,4) + 
                       s2*Power(t2,4) + 
                       Power(t1,2)*(6 - 11*t2 + 3*Power(t2,2)) - 
                       t1*(8 + (-6 + 5*s2)*t2 - 
                        4*(3 + 2*s2)*Power(t2,2) + 
                        (5 + s2)*Power(t2,3)) + 
                       Power(s1,3)*
                        (2 + Power(s2,2) + 2*Power(t1,2) - 
                        4*t1*(-1 + t2) - 5*t2 + 3*Power(t2,2) + 
                        s2*(-2 - 3*t1 + 2*t2)) + 
                       Power(s1,2)*
                        (Power(t1,2)*(5 - 3*t2) - 
                        Power(s2,2)*(-1 + t2) + 
                        t2*(-5 + 9*t2 - 4*Power(t2,2)) + 
                        t1*(9 - 16*t2 + 6*Power(t2,2)) + 
                        s2*(-1 + t2 + Power(t2,2) + t1*(-5 + 3*t2))) + 
                       s1*(-2 + 2*t2 + Power(s2,2)*(-2 + t2)*t2 + 
                        5*Power(t2,2) - 6*Power(t2,3) + Power(t2,4) + 
                        Power(t1,2)*(7 - 10*t2 + Power(t2,2)) + 
                        t1*(4 - 27*t2 + 21*Power(t2,2) - 
                        2*Power(t2,3)) + 
                        s2*(-1 + 5*t2 + 4*Power(t2,2) - 
                        4*Power(t2,3) + t1*(-1 + 3*t2 + Power(t2,2)))))\
))/(s*Power(s1 - s1*t2 + t2*(-1 + s + t2),2))))/
           ((-1 + s2)*(-s + s2 - t1)))/Power(-1 + s1,2))*
     B1(1 - s + s1 - t2,s,s1))/(4.*Power(Pi,2)) + 
  (((64*(1 - s1 - t1 + t2)*(2*s2 - 2*t1 - 2*s2*t1 + 2*Power(t1,2) + 
            Power(s1,2)*(-s2 + t1) + 
            s1*(-2*Power(s2,2) + 2*s2*t1 - (1 + t1)*(-1 + t2)) + 
            Power(s,2)*(2 - 2*s2 + t1 - t2) - t2 + 2*Power(s2,2)*t2 - 
            t1*t2 - 2*s2*t1*t2 + Power(t2,2) + s2*Power(t2,2) + 
            s*(-3 + 2*Power(s2,2) + t1 + 3*t2 + t1*t2 - Power(t2,2) + 
               s1*(-2 + 3*s2 - 2*t1 + t2) - s2*(2*t1 + t2)))*
          ((-Power(s,4) + 22*s2 - 5*Power(s2,2) - 4*Power(s2,3) + 
               Power(s2,4) - 2*Power(s1,3)*(s2 - t1) - t1 + 6*s2*t1 + 
               9*Power(s2,2)*t1 - 2*Power(s2,3)*t1 - 4*Power(t1,2) - 
               6*s2*Power(t1,2) + Power(s2,2)*Power(t1,2) + 
               Power(t1,3) + Power(s,3)*(-1 + 4*s1 + 2*s2 - 2*t2) - 
               17*s2*t2 - 6*Power(s2,2)*t2 + Power(s2,3)*t2 + 
               2*t1*t2 + 7*s2*t1*t2 - s2*Power(t1,2)*t2 + 
               2*Power(t2,2) + s2*Power(t2,2) - 
               Power(s2,2)*Power(t2,2) + t1*Power(t2,2) + 
               s2*t1*Power(t2,2) - 2*Power(t2,3) + 
               Power(s1,2)*(13 - 2*Power(s2,2) + Power(t1,2) - 
                  3*t2 - t1*(2 + 3*t2) + s2*(2 + t1 + 3*t2)) + 
               Power(s,2)*(16 - 5*Power(s1,2) + 4*t1 + Power(t1,2) - 
                  8*t2 - t1*t2 - Power(t2,2) + 
                  s2*(-2 - 2*t1 + 5*t2) + s1*(3 - 7*s2 + 2*t1 + 5*t2)\
) + s1*(21 + Power(s2,3) - 18*t2 + 5*Power(t2,2) - 
                  Power(t1,2)*(4 + t2) + 
                  Power(s2,2)*(-1 - 3*t1 + 4*t2) + 
                  s2*(14 + 2*Power(t1,2) + t1*(5 - 3*t2) - 6*t2 - 
                     Power(t2,2)) + t1*(-5 + 5*t2 + Power(t2,2))) + 
               s*(-22 + 2*Power(s1,3) - 2*Power(s2,3) - 2*t1 + 
                  6*Power(t1,2) + Power(s2,2)*(7 + 4*t1 - 4*t2) + 
                  Power(s1,2)*(-2 + 7*s2 - 4*t1 - 3*t2) + 20*t2 - 
                  5*t1*t2 + Power(t1,2)*t2 - 4*Power(t2,2) - 
                  t1*Power(t2,2) + 
                  s1*(-32 + 2*Power(s2,2) - t1 - 2*Power(t1,2) + 
                     s2*(-1 + t1 - 9*t2) + 11*t2 + 4*t1*t2 + 
                     Power(t2,2)) + 
                  s2*(-10 - 2*Power(t1,2) + t1*(-12 + t2) + 13*t2 + 
                     2*Power(t2,2))))/((-1 + s2)*(-s + s2 - t1)) - 
            (Power(s,4) + 2*Power(s1,4) - 23*s2 - 11*Power(s2,2) + 
               4*Power(s2,3) - Power(s2,4) + 11*t1 + 11*s2*t1 + 
               Power(s2,3)*t1 - Power(t1,2) - 3*s2*Power(t1,2) + 
               Power(s1,3)*(4*s2 + t1 - 5*t2) - 
               Power(s,3)*(-2 + 5*s1 + 2*s2 + t1 - 3*t2) + 6*t2 + 
               8*s2*t2 + 2*Power(s2,2)*t2 + t1*t2 + 2*s2*t1*t2 - 
               2*Power(s2,2)*t1*t2 - 6*Power(t2,2) - 
               3*s2*Power(t2,2) + 2*Power(s2,2)*Power(t2,2) - 
               2*t1*Power(t2,2) + s2*t1*Power(t2,2) + 2*Power(t2,3) - 
               s2*Power(t2,3) + 
               Power(s,2)*(-27 + 9*Power(s1,2) + 
                  s1*(-3 + 8*s2 + 3*t1 - 11*t2) + 
                  t1*(3 + 3*s2 - 2*t2) + 4*t2 - 6*s2*t2 + 
                  3*Power(t2,2)) + 
               Power(s1,2)*(-19 + Power(s2,2) + 
                  3*s2*(1 + t1 - 3*t2) + 5*t2 + 4*Power(t2,2) - 
                  2*t1*(1 + t2)) + 
               s1*(-18 - 2*Power(s2,3) - 4*Power(t1,2) + 
                  Power(s2,2)*(6 + 3*t1 - 4*t2) + 17*t2 - 
                  7*Power(t2,2) - Power(t2,3) + 
                  t1*(12 + 5*t2 + Power(t2,2)) + 
                  s2*(-30 - 4*(-1 + t1)*t2 + 6*Power(t2,2))) - 
               s*(-29 + 7*Power(s1,3) - 2*Power(s2,3) + 15*t1 - 
                  4*Power(t1,2) + 
                  Power(s1,2)*(-1 + 10*s2 + 3*t1 - 13*t2) + 
                  3*Power(s2,2)*(2 + t1 - t2) + 16*t2 + 4*t1*t2 - 
                  5*Power(t2,2) + t1*Power(t2,2) - Power(t2,3) + 
                  s2*(-36 + 3*t1 + 6*t2 - 4*t1*t2 + 5*Power(t2,2)) + 
                  s1*(-44 + Power(s2,2) + t1 + 
                     3*s2*(1 + 2*t1 - 5*t2) + 8*t2 - 4*t1*t2 + 
                     7*Power(t2,2))))/((-1 + s2)*(-1 + t1)) + 
            (3*Power(s1,3) + 36*s2 + 15*Power(s2,2) - 4*Power(s2,3) + 
               Power(s,3)*(-3 + t1) - 24*t1 - 19*s2*t1 + 
               5*Power(s2,2)*t1 + 4*Power(t1,2) - s2*Power(t1,2) + 
               Power(s1,2)*(-2 + s2 + 4*t1 - s2*t1 + Power(t1,2) - 
                  8*t2) + 12*t2 - 3*s2*t2 - 7*Power(s2,2)*t2 + 
               9*t1*t2 + 8*s2*t1*t2 + Power(s2,2)*t1*t2 - 
               3*Power(t1,2)*t2 - s2*Power(t1,2)*t2 - 
               12*Power(t2,2) - 2*s2*Power(t2,2) + 5*t1*Power(t2,2) - 
               2*Power(t2,3) + 
               Power(s,2)*(18 - 2*s1*(-3 + t1) - 2*s2*(-1 + t1) - 
                  t1 + Power(t1,2) - 7*t2 + t1*t2) + 
               s1*(-(Power(s2,2)*(1 + t1)) - Power(t1,2)*(-2 + t2) - 
                  5*t1*(4 + t2) + t2*(18 + 7*t2) + 
                  s2*(18 + t1 + Power(t1,2) - 3*t2 + t1*t2)) + 
               s*(-24 + Power(s1,2)*(-6 + t1) + 28*t1 - 
                  2*Power(t1,2) + Power(s2,2)*(5 + t1) - 9*t2 - 
                  t1*t2 + Power(t1,2)*t2 - 2*Power(t2,2) - 
                  s2*(33 + Power(t1,2) - 12*t2 + 2*t1*(1 + t2)) + 
                  s1*(3*s2*(-1 + t1) - 2*Power(t1,2) + 12*(-1 + t2) - 
                     t1*(3 + t2))))/((-1 + s1)*(-s + s1 - t2)) + 
            (Power(s1,3) + 20*s2 - 4*Power(s2,3) + 
               Power(s,3)*(-2 + t1) + 4*s2*t1 + 8*Power(s2,2)*t1 - 
               4*Power(t1,2) - 4*s2*Power(t1,2) + 
               Power(s1,2)*(2 - s2*(-7 + t1) - 5*t1 + Power(t1,2) - 
                  t2) - 10*s2*t2 - 6*Power(s2,2)*t2 + 3*t1*t2 + 
               2*s2*t1*t2 + Power(s2,2)*t1*t2 + 3*Power(t1,2)*t2 - 
               s2*Power(t1,2)*t2 - Power(t2,2) + 5*s2*Power(t2,2) - 
               4*t1*Power(t2,2) + Power(t2,3) + 
               Power(s,2)*(9 - 2*s1*(-4 + t1) + Power(t1,2) - 
                  10*t2 + t1*(5 - 2*s2 + t2)) + 
               s*(-20 + Power(s1,2)*(-7 + t1) - t1 + 7*Power(t1,2) + 
                  Power(s2,2)*(6 + t1) + 8*t2 - 7*t1*t2 + 
                  Power(t1,2)*t2 - 3*Power(t2,2) - 
                  s1*(15 + s2*(7 - 3*t1) + 2*Power(t1,2) - 14*t2 + 
                     t1*t2) - 
                  s2*(10 + Power(t1,2) - 17*t2 + 2*t1*(7 + t2))) - 
               s1*(-20 + Power(s2,2)*(2 + t1) + t1*(7 - 13*t2) + 
                  5*t2 + Power(t2,2) + Power(t1,2)*(7 + t2) - 
                  s2*(10 + Power(t1,2) - 16*t2 + t1*(10 + t2))))/
             ((s - s2 + t1)*(s - s1 + t2)) + 
            (-2*Power(s,4) + 31*s2 + 20*Power(s2,2) + Power(s2,3) - 
               Power(s2,4) - 18*t1 - 30*s2*t1 - 4*Power(s2,2)*t1 + 
               2*Power(s2,3)*t1 + 9*Power(t1,2) + 4*s2*Power(t1,2) - 
               Power(s2,2)*Power(t1,2) - Power(t1,3) + 
               Power(s1,3)*(4 - s2 + t1) + 
               Power(s,3)*(1 + 5*s1 + 7*s2 - 4*t1 - t2) + 
               Power(s1,2)*(1 - 4*Power(s2,2) + t1 - 2*Power(t1,2) + 
                  s2*(7 + 6*t1) - t2) + 10*t2 - 6*s2*t2 - 
               Power(s2,2)*t2 + 2*Power(s2,3)*t2 + 3*t1*t2 + 
               2*s2*t1*t2 - 4*Power(s2,2)*t1*t2 + 
               2*s2*Power(t1,2)*t2 - 5*Power(t2,2) + 
               Power(s2,2)*Power(t2,2) - 2*t1*Power(t2,2) - 
               s2*t1*Power(t2,2) + 3*Power(t2,3) + 
               Power(s,2)*(16 - 4*Power(s1,2) - 9*Power(s2,2) - 
                  3*t1 - 2*Power(t1,2) - t2 - 3*t1*t2 + 
                  Power(t2,2) + s1*(1 - 14*s2 + 9*t1 + t2) + 
                  2*s2*(-1 + 5*t1 + 2*t2)) + 
               s1*(3 - 4*Power(s2,3) + 3*t2 - 6*Power(t2,2) + 
                  2*Power(t1,2)*(1 + t2) + 
                  Power(s2,2)*(3 + 7*t1 + 2*t2) - 
                  s2*(-18 + t1 + 3*Power(t1,2) + 6*t2 + 4*t1*t2 - 
                     Power(t2,2)) - t1*(17 - 5*t2 + Power(t2,2))) + 
               s*(-21 + Power(s1,3) + 5*Power(s2,3) + 29*t1 - 
                  4*Power(t1,2) + Power(s1,2)*(8*s2 - 6*(1 + t1)) + 
                  2*t2 - 5*t1*t2 - 2*Power(t1,2)*t2 + 4*Power(t2,2) + 
                  t1*Power(t2,2) - Power(s2,2)*(8*t1 + 5*t2) + 
                  s2*(-35 + 6*t1 + 3*Power(t1,2) + 3*t2 + 7*t1*t2 - 
                     2*Power(t2,2)) + 
                  s1*(-14 + 13*Power(s2,2) + t1 + 4*Power(t1,2) + 
                     2*t2 + 3*t1*t2 - Power(t2,2) - 
                     s2*(4 + 16*t1 + 3*t2))))/((-1 + s1)*(-1 + t2)) + 
            (-2*Power(s,4) - Power(s1,4) + 18*s2 + 15*Power(s2,2) + 
               Power(s2,3) - Power(s2,4) - 6*t1 - 12*s2*t1 - 
               3*Power(s2,2)*t1 + Power(s2,3)*t1 + s2*Power(t1,2) + 
               Power(s,3)*(-1 + 7*s1 + 7*s2 - 2*t1 - 3*t2) - 5*t2 - 
               21*s2*t2 - 3*Power(s2,2)*t2 + 3*Power(s2,3)*t2 + 
               9*t1*t2 + 5*s2*t1*t2 - 3*Power(s2,2)*t1*t2 - 
               Power(s2,2)*Power(t2,2) - t1*Power(t2,2) + 
               2*s2*t1*Power(t2,2) + Power(t2,3) - s2*Power(t2,3) + 
               Power(s1,3)*(2 - 5*s2 + 2*t1 + t2) + 
               Power(s1,2)*(13 - 8*Power(s2,2) - 2*t2 + Power(t2,2) - 
                  t1*(3 + 4*t2) + s2*(6 + 5*t1 + 6*t2)) + 
               Power(s,2)*(21 - 9*Power(s1,2) - 9*Power(s2,2) - 
                  2*t1 - 4*t2 - 4*t1*t2 + 
                  s1*(5 - 19*s2 + 6*t1 + 7*t2) + s2*(4 + 5*t1 + 9*t2)) \
+ s*(-23 + 5*Power(s1,3) + 5*Power(s2,3) + 16*t1 + 
                  Power(s1,2)*(-6 + 17*s2 - 6*t1 - 5*t2) + 20*t2 - 
                  5*t1*t2 - 2*t1*Power(t2,2) + Power(t2,3) - 
                  Power(s2,2)*(4 + 4*t1 + 9*t2) + 
                  s2*(-37 + 6*t1 + 6*t2 + 7*t1*t2 + Power(t2,2)) + 
                  s1*(-36 + 17*Power(s2,2) + 5*t1 + 7*t2 + 8*t1*t2 - 
                     Power(t2,2) - 5*s2*(2 + 2*t1 + 3*t2))) - 
               s1*(-17 + 5*Power(s2,3) + 13*t2 + Power(t2,2) + 
                  Power(t2,3) - Power(s2,2)*(5 + 4*t1 + 8*t2) - 
                  2*t1*(-9 + 3*t2 + Power(t2,2)) + 
                  s2*(-33 + 5*t2 + t1*(5 + 7*t2))))/
             ((-1 + t1)*(-1 + t2))))/
        ((s - s2 + t1)*Power(-1 + s1 + t1 - t2,3)*(-1 + s2 - t1 + t2) - 
          Power(-1 + s2 + s*t1 - s2*t1 + Power(t1,2) + t2 - s*t2 - 
            t1*t2,2) - Power(-1 + s1 + t1 - t2,2)*
           (3 + Power(s2,2) + 2*t1 - s2*(3 + t1) + (-3 + t1)*t2 + 
             2*Power(-1 + s - s2 + t1,2)*(-1 + s2 - t1 + t2) + 
             (-1 + s - s2 + t1)*(s2 - t2)*(-1 + s2 - t1 + t2)) + 
          (1 - s1 - t1 + t2)*(2*(-1 + t1)*(-1 + t2) + 
             Power(s2,2)*(2 + t2) - 
             Power(-1 + s - s2 + t1,3)*(-1 + s2 - t1 + t2) + 
             Power(-1 + s - s2 + t1,2)*(1 - s2 + t2)*
              (-1 + s2 - t1 + t2) + 
             (1 - s + s2 - t1)*
              (2 - s2 + 4*t1 - Power(s2,2)*(1 + t2) + 
                t2*(-3 - t1 + t2) - s2*(1 + t2)*(-t1 + t2)) - 
             s2*(2*(2 + t1) + t2*(-5 + t1 + t2)))) + 
       8*((256*(9 + Power(s,3) + 2*s2 - 5*Power(s2,2) + 
                2*Power(s2,3) - 2*Power(s1,2)*(1 + s2 - t1) - 6*t1 + 
                7*s2*t1 - 3*Power(s2,2)*t1 - Power(t1,2) + 
                s2*Power(t1,2) - t2 - 4*s2*t2 + Power(s2,2)*t2 + 
                5*t1*t2 - Power(t1,2)*t2 - 3*Power(t2,2) - 
                s2*Power(t2,2) + t1*Power(t2,2) + 
                Power(s,2)*(-3 - 3*s1 + 2*t2) + 
                s1*(2 + Power(s2,2) - 2*s2*t1 + Power(t1,2) + 5*t2 + 
                   3*s2*t2 - 3*t1*t2) + 
                s*(-7 + 2*Power(s1,2) - 3*Power(s2,2) - 2*t1 - 
                   Power(t1,2) + s1*(4 + 2*s2 - 2*t1 - 3*t2) + 
                   s2*(8 + 3*t1 - 3*t2) + t1*t2 + Power(t2,2))) - 
             (128*(3 + 2*s1 + Power(s2,2) + 2*t1 - s2*(1 + s + t1) - 
                  2*t2)*(-2*s2 + Power(s1,2)*(s2 - t1) + 2*t1 + 
                  2*s2*t1 - 2*Power(t1,2) + 
                  s1*(2*Power(s2,2) - 2*s2*t1 + (1 + t1)*(-1 + t2)) + 
                  t2 - 2*Power(s2,2)*t2 + t1*t2 + 2*s2*t1*t2 - 
                  Power(t2,2) - s2*Power(t2,2) + 
                  Power(s,2)*(-2 + 2*s2 - t1 + t2) + 
                  s*(3 - 2*Power(s2,2) - t1 + 2*s2*t1 - 3*t2 + 
                     s2*t2 - t1*t2 + Power(t2,2) - 
                     s1*(-2 + 3*s2 - 2*t1 + t2))))/
              (-1 + s - s*s2 + s1*s2 + t1 - s2*t2))/
           (16.*(-1 + s2)*(-s + s2 - t1)) - 
          (8*(-20 + 4*s2 + 8*Power(s2,2) + 20*t1 - 16*s2*t1 - 
               8*Power(s2,2)*t1 + 4*Power(t1,2) + 12*s2*Power(t1,2) - 
               4*Power(t1,3) + 
               Power(s1,3)*(5*s2 + Power(s2,2) + t1 - Power(t1,2)) + 
               Power(s,4)*(t1 - t2) - 6*t2 - 16*s2*t2 + 
               4*Power(s2,2)*t2 + 8*Power(s2,3)*t2 - 2*t1*t2 - 
               4*s2*t1*t2 - 12*Power(s2,2)*t1*t2 + 
               4*s2*Power(t1,2)*t2 + Power(t2,2) + 3*s2*Power(t2,2) + 
               4*Power(s2,2)*Power(t2,2) + 4*t1*Power(t2,2) + 
               s2*t1*Power(t2,2) - Power(t1,2)*Power(t2,2) - 
               Power(t2,3) - 6*s2*Power(t2,3) - 
               Power(s2,2)*Power(t2,3) + t1*Power(t2,3) + 
               s2*t1*Power(t2,3) + 
               Power(s,3)*(-5 + Power(t1,2) - 
                  s2*(-6 + s1 + t1 - 2*t2) + 4*t2 + 2*s1*t2 - 
                  2*Power(t2,2) + t1*(-4 - 3*s1 + t2)) + 
               Power(s,2)*(3 + Power(s2,2)*(-14 + s1 - t2) - t2 - 
                  6*s1*t2 - Power(s1,2)*t2 + 6*Power(t2,2) + 
                  2*s1*Power(t2,2) - Power(t2,3) + 
                  Power(t1,2)*(-3 - 3*s1 + 2*t2) + 
                  t1*(-10 + 9*s1 + 3*Power(s1,2) - 3*t2 - 2*s1*t2 - 
                     Power(t2,2)) + 
                  s2*(13 + 2*Power(s1,2) + 17*t1 + 
                     s1*(-5 + 2*t1 - 6*t2) - 6*t2 - t1*t2 + 
                     4*Power(t2,2))) - 
               Power(s1,2)*(5 + Power(t1,2)*(3 - 2*t2) + 
                  t1*(-12 + t2) + t2 + Power(s2,2)*(4 + 3*t2) + 
                  s2*(7 + 16*t2 - t1*(11 + t2))) - 
               s1*(8*Power(s2,3) + 
                  Power(s2,2)*(4 - 12*t1 - 3*Power(t2,2)) + 
                  Power(t1,2)*(-8 - 4*t2 + Power(t2,2)) - 
                  2*(3 + 2*t2 + Power(t2,2)) + 
                  t1*(6 + 16*t2 + Power(t2,2)) + 
                  s2*(-24 + 4*Power(t1,2) - 4*t2 - 17*Power(t2,2) + 
                     2*t1*(2 + 6*t2 + Power(t2,2)))) + 
               s*(18 + 8*Power(s2,3) - 2*t1 - 16*Power(t1,2) - 
                  Power(s1,3)*(s2 + t1) + 8*t2 + 14*t1*t2 - 
                  4*Power(t1,2)*t2 - 5*Power(t2,2) + 
                  2*t1*Power(t2,2) + Power(t1,2)*Power(t2,2) + 
                  2*Power(t2,3) - t1*Power(t2,3) + 
                  Power(s1,2)*
                   (5 - 2*Power(s2,2) + 3*Power(t1,2) - 
                     s2*(6 + t1 - 4*t2) + t1*(-6 + t2) + 2*t2) - 
                  2*Power(s2,2)*(2 + 6*t1 + 5*t2 + Power(t2,2)) + 
                  s1*(Power(t1,2)*(6 - 4*t2) + 
                     2*Power(s2,2)*(9 + 2*t2) + 
                     s2*(2 - 28*t1 + 24*t2 - 5*Power(t2,2)) + 
                     t1*(-6 + 4*t2 + Power(t2,2)) - 
                     2*(5 + 2*Power(t2,2))) + 
                  s2*(4*Power(t1,2) + t1*(24 + 18*t2 + Power(t2,2)) + 
                     2*(-14 - 9*Power(t2,2) + Power(t2,3))))))/
           ((s - s2 + t1)*(s - s1 + t2)*
             (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2)) + 
          (8*(-42 - 22*s2 + 8*Power(s2,2) + 128*t1 + 64*s2*t1 - 
               24*Power(s2,2)*t1 - 132*Power(t1,2) - 
               60*s2*Power(t1,2) + 24*Power(s2,2)*Power(t1,2) + 
               48*Power(t1,3) + 16*s2*Power(t1,3) - 
               8*Power(s2,2)*Power(t1,3) - 2*Power(t1,4) + 
               2*s2*Power(t1,4) + Power(s1,6)*t1*(-s2 + t1) - 
               Power(s,5)*Power(t1 - t2,2)*(-1 + s1 + t1 - t2) - 
               47*t2 - 94*s2*t2 - 36*Power(s2,2)*t2 + 
               16*Power(s2,3)*t2 + 76*t1*t2 + 184*s2*t1*t2 + 
               68*Power(s2,2)*t1*t2 - 32*Power(s2,3)*t1*t2 - 
               10*Power(t1,2)*t2 - 84*s2*Power(t1,2)*t2 - 
               28*Power(s2,2)*Power(t1,2)*t2 + 
               16*Power(s2,3)*Power(t1,2)*t2 - 20*Power(t1,3)*t2 - 
               8*s2*Power(t1,3)*t2 - 4*Power(s2,2)*Power(t1,3)*t2 + 
               Power(t1,4)*t2 + 2*s2*Power(t1,4)*t2 + 10*Power(t2,2) - 
               62*s2*Power(t2,2) - 70*Power(s2,2)*Power(t2,2) - 
               6*Power(s2,3)*Power(t2,2) + 8*Power(s2,4)*Power(t2,2) - 
               42*t1*Power(t2,2) + 24*s2*t1*Power(t2,2) + 
               56*Power(s2,2)*t1*Power(t2,2) + 
               4*Power(s2,3)*t1*Power(t2,2) - 
               8*Power(s2,4)*t1*Power(t2,2) + 
               30*Power(t1,2)*Power(t2,2) + 
               42*s2*Power(t1,2)*Power(t2,2) + 
               18*Power(s2,2)*Power(t1,2)*Power(t2,2) + 
               2*Power(s2,3)*Power(t1,2)*Power(t2,2) + 
               2*Power(t1,3)*Power(t2,2) - 
               4*s2*Power(t1,3)*Power(t2,2) - 
               4*Power(s2,2)*Power(t1,3)*Power(t2,2) + 
               12*Power(t2,3) + 29*s2*Power(t2,3) - 
               13*Power(s2,2)*Power(t2,3) - 
               18*Power(s2,3)*Power(t2,3) + 
               8*Power(s2,4)*Power(t2,3) - 11*t1*Power(t2,3) - 
               48*s2*t1*Power(t2,3) - 24*Power(s2,2)*t1*Power(t2,3) - 
               8*Power(s2,3)*t1*Power(t2,3) - 
               5*s2*Power(t1,2)*Power(t2,3) + 
               5*Power(s2,2)*Power(t1,2)*Power(t2,3) + 
               2*Power(s2,3)*Power(t1,2)*Power(t2,3) - 
               Power(t1,3)*Power(t2,3) + Power(t2,4) + 
               5*s2*Power(t2,4) + 19*Power(s2,2)*Power(t2,4) + 
               2*Power(s2,3)*Power(t2,4) + 2*t1*Power(t2,4) + 
               3*s2*t1*Power(t2,4) + 3*Power(s2,2)*t1*Power(t2,4) - 
               2*Power(s2,3)*t1*Power(t2,4) + 
               Power(t1,2)*Power(t2,4) + 
               2*s2*Power(t1,2)*Power(t2,4) - 3*Power(t2,5) + 
               2*s2*Power(t2,5) - 3*Power(s2,2)*Power(t2,5) - 
               t1*Power(t2,5) - 2*s2*t1*Power(t2,5) - 
               Power(s2,2)*t1*Power(t2,5) + Power(t2,6) + 
               Power(s2,2)*Power(t2,6) + 
               Power(s1,5)*(t1*
                   (1 + Power(t1,2) + t1*(2 - 4*t2) - 2*t2) - 
                  Power(s2,2)*(4 + 2*t1 + t2) + 
                  s2*(-Power(t1,2) + t2 + t1*(-5 + 4*t2))) + 
               Power(s1,4)*(-2*Power(s2,3)*(3 + t1) - 
                  3*Power(t1,3)*(-1 + t2) + (-1 + t2)*t2 + 
                  Power(t1,2)*(-4 - 8*t2 + 6*Power(t2,2)) + 
                  t1*(5 - 8*t2 + 8*Power(t2,2)) + 
                  Power(s2,2)*
                   (21 - 2*Power(t1,2) + 13*t2 + 5*Power(t2,2) + 
                     t1*(-5 + 7*t2)) + 
                  s2*(8 + 7*t2 - 4*Power(t2,2) + 
                     Power(t1,2)*(-13 + 3*t2) + 
                     t1*(7 + 14*t2 - 6*Power(t2,2)))) - 
               Power(s,4)*(6 - 9*t1 + 7*Power(t1,2) - 6*Power(t1,3) + 
                  2*Power(s2,2)*(-3 + t1)*(-1 + t1 - t2) + 5*t2 - 
                  10*t1*t2 + 15*Power(t1,2)*t2 + 3*Power(t1,3)*t2 + 
                  3*Power(t2,2) - 12*t1*Power(t2,2) - 
                  9*Power(t1,2)*Power(t2,2) + 3*Power(t2,3) + 
                  9*t1*Power(t2,3) - 3*Power(t2,4) + 
                  Power(s1,2)*
                   (-5*Power(t1,2) + (1 + s2 - 3*t2)*t2 + 
                     t1*(-3 + s2 + 8*t2)) + 
                  s2*(-12 - 9*t2 + 7*Power(t2,2) + 2*Power(t2,3) + 
                     Power(t1,2)*(-1 + 2*t2) + 
                     t1*(15 - 6*t2 - 4*Power(t2,2))) + 
                  s1*(-6 + 2*Power(s2,2)*(-3 + t1) - 5*Power(t1,3) - 
                     4*Power(t2,2) + 6*Power(t2,3) + 
                     4*Power(t1,2)*(-1 + 4*t2) + 
                     t1*(6 + 8*t2 - 17*Power(t2,2)) + 
                     s2*(12 + Power(t1,2) + 2*t1*(-2 + t2) - 4*t2 - 
                        3*Power(t2,2)))) - 
               Power(s1,3)*(4 + 8*Power(s2,4) + 6*t2 - 
                  6*Power(t2,2) + 4*Power(t2,3) + 
                  2*Power(s2,3)*
                   (-13 + Power(t1,2) - 8*t2 - 4*t1*t2) + 
                  Power(t1,3)*(7 + 6*t2 - 3*Power(t2,2)) + 
                  4*Power(t1,2)*(-4 - 3*Power(t2,2) + Power(t2,3)) + 
                  t1*(5 + 4*t2 - 19*Power(t2,2) + 12*Power(t2,3)) + 
                  Power(s2,2)*
                   (-39 + Power(t1,2)*(11 - 6*t2) + 82*t2 + 
                     12*Power(t2,2) + 10*Power(t2,3) + 
                     4*t1*(-1 - 3*t2 + 2*Power(t2,2))) + 
                  s2*(42 + 5*Power(t1,3) + 29*t2 + 23*Power(t2,2) - 
                     6*Power(t2,3) + 
                     Power(t1,2)*(-2 - 29*t2 + 3*Power(t2,2)) + 
                     t1*(-69 + 16*t2 + 10*Power(t2,2) - 4*Power(t2,3))\
)) + s1*(27 - 23*t2 - 24*Power(t2,2) - 8*Power(t2,3) + 
                  10*Power(t2,4) - 4*Power(t2,5) + 
                  Power(t1,4)*(-5 + 3*t2) - 
                  8*Power(s2,4)*t2*(2 - 2*t1 + 3*t2) - 
                  Power(t1,3)*(-52 + 18*t2 + Power(t2,2)) + 
                  2*Power(t1,2)*
                   (-31 - 22*t2 + 6*Power(t2,2) - 5*Power(t2,3) + 
                     Power(t2,4)) + 
                  t1*(-12 + 82*t2 + 13*Power(t2,2) + 2*Power(t2,3) + 
                     8*Power(t2,4) - 2*Power(t2,5)) - 
                  2*Power(s2,3)*
                   (8 - 6*t2 - 31*Power(t2,2) + 
                     Power(t1,2)*(8 + 2*t2 + 3*Power(t2,2)) - 
                     4*t1*(4 - t2 + 2*Power(t2,2) + Power(t2,3))) - 
                  s2*(-110 + 2*Power(t1,4) - 122*t2 + 
                     92*Power(t2,2) + 23*Power(t2,3) + 
                     13*Power(t2,4) - Power(t2,5) + 
                     Power(t1,3)*(8 - 18*t2 + 5*Power(t2,2)) + 
                     Power(t1,2)*
                      (-132 + 106*t2 - 20*Power(t2,2) + Power(t2,3)) \
+ t1*(232 + 34*t2 - 149*Power(t2,2) + 8*Power(t2,3) - 7*Power(t2,4))) \
+ Power(s2,2)*(36 + 164*t2 + 65*Power(t2,2) - 78*Power(t2,3) + 
                     8*Power(t2,4) - 5*Power(t2,5) + 
                     Power(t1,3)*(4 + 8*t2) + 
                     Power(t1,2)*
                      (28 - 12*t2 - 21*Power(t2,2) + 2*Power(t2,3)) + 
                     2*t1*(-34 - 80*t2 + 26*Power(t2,2) - 
                        2*Power(t2,3) + Power(t2,4)))) + 
               Power(s1,2)*(21 - 3*Power(t1,4) - 
                  8*Power(s2,4)*(-1 + t1 - 3*t2) + 16*t2 + 
                  13*Power(t2,2) - 12*Power(t2,3) + 6*Power(t2,4) + 
                  Power(t1,3)*
                   (8 + 9*t2 + 3*Power(t2,2) - Power(t2,3)) + 
                  Power(t1,2)*
                   (38 - 28*t2 + 13*Power(t2,2) - 8*Power(t2,3) + 
                     Power(t2,4)) + 
                  t1*(-64 + 3*t2 - 5*Power(t2,2) - 19*Power(t2,3) + 
                     8*Power(t2,4)) + 
                  2*Power(s2,3)*
                   (-3 - 35*t2 - 6*Power(t2,2) + 
                     Power(t1,2)*(1 + 3*t2) + 
                     t1*(2 - 4*t2 - 6*Power(t2,2))) - 
                  Power(s2,2)*
                   (94 + 4*Power(t1,3) + 91*t2 - 120*Power(t2,2) + 
                     2*Power(t2,3) - 10*Power(t2,4) + 
                     3*Power(t1,2)*(2 - 9*t2 + 2*Power(t2,2)) - 
                     2*t1*(52 - 16*t2 - 3*Power(t2,2) + Power(t2,3))) \
+ s2*(-60 + 105*t2 + 39*Power(t2,2) + 27*Power(t2,3) - 
                     4*Power(t2,4) + 2*Power(t1,3)*(-7 + 5*t2) + 
                     Power(t1,2)*
                      (64 - 17*t2 - 17*Power(t2,2) + Power(t2,3)) - 
                     t1*(-10 + 170*t2 - 14*Power(t2,2) + 
                        4*Power(t2,3) + Power(t2,4)))) + 
               Power(s,3)*(44 - 73*t1 + 44*Power(t1,2) - 
                  15*Power(t1,3) + 
                  2*Power(s2,3)*(1 + t1)*(-1 + t1 - t2) + 49*t2 - 
                  39*t1*t2 + 18*Power(t1,2)*t2 + 14*Power(t1,3)*t2 + 
                  9*Power(t2,2) - 5*t1*Power(t2,2) - 
                  39*Power(t1,2)*Power(t2,2) - 
                  3*Power(t1,3)*Power(t2,2) + 2*Power(t2,3) + 
                  36*t1*Power(t2,3) + 9*Power(t1,2)*Power(t2,3) - 
                  11*Power(t2,4) - 9*t1*Power(t2,4) + 3*Power(t2,5) + 
                  Power(s1,3)*
                   (-10*Power(t1,2) + (2 - 3*t2)*t2 + 
                     s2*(-1 + 4*t1 + 3*t2) + 3*t1*(-3 + 4*t2)) + 
                  Power(s2,2)*
                   (36 + 29*t2 - 10*Power(t2,2) + Power(t2,3) - 
                     2*Power(t1,2)*(1 + 3*t2) + 
                     t1*(-34 + 19*t2 + 5*Power(t2,2))) + 
                  s2*(-77 + 5*Power(t1,3) - 66*t2 + 23*Power(t2,2) - 
                     10*Power(t2,3) - 6*Power(t2,4) + 
                     Power(t1,2)*(-25 + 6*t2 - 6*Power(t2,2)) + 
                     t1*(97 - 18*t2 - Power(t2,2) + 12*Power(t2,3))) \
+ Power(s1,2)*(-3 - 10*Power(t1,3) + 7*t2 - 15*Power(t2,2) + 
                     9*Power(t2,3) + Power(s2,2)*(-15 + 8*t1 + t2) + 
                     Power(t1,2)*(-20 + 34*t2) + 
                     t1*(1 + 44*t2 - 33*Power(t2,2)) + 
                     s2*(21 + 4*Power(t1,2) - 8*t2 - 
                        12*Power(t2,2) + t1*(-9 + 2*t2))) + 
                  s1*(-41 + 2*Power(s2,3)*(1 + t1) - 7*t2 - 
                     9*Power(t2,2) + 24*Power(t2,3) - 9*Power(t2,4) + 
                     3*Power(t1,3)*(-7 + 4*t2) + 
                     Power(t1,2)*(-10 + 68*t2 - 33*Power(t2,2)) + 
                     t1*(30 - 3*t2 - 71*Power(t2,2) + 
                       30*Power(t2,3)) + 
                     Power(s2,2)*
                      (-21 + 8*Power(t1,2) + 25*t2 - 2*Power(t2,2) - 
                        t1*(29 + 13*t2)) + 
                     s2*(57 - 33*t2 + 19*Power(t2,2) + 
                        15*Power(t2,3) + 3*Power(t1,2)*(3 + t2) - 
                        2*t1*(-6 + t2 + 9*Power(t2,2))))) + 
               Power(s,2)*(-98 + 229*t1 - 159*Power(t1,2) + 
                  31*Power(t1,3) - 3*Power(t1,4) - 131*t2 + 
                  141*t1*t2 + 19*Power(t1,2)*t2 - 21*Power(t1,3)*t2 - 
                  14*Power(t2,2) - 48*t1*Power(t2,2) + 
                  46*Power(t1,2)*Power(t2,2) + 
                  10*Power(t1,3)*Power(t2,2) + 16*Power(t2,3) - 
                  35*t1*Power(t2,3) - 29*Power(t1,2)*Power(t2,3) - 
                  Power(t1,3)*Power(t2,3) + 13*Power(t2,4) + 
                  28*t1*Power(t2,4) + 3*Power(t1,2)*Power(t2,4) - 
                  9*Power(t2,5) - 3*t1*Power(t2,5) + Power(t2,6) + 
                  Power(s2,4)*(8 - 8*t1 + 8*t2) + 
                  Power(s1,4)*
                   (10*Power(t1,2) + t1*(9 - 8*t2) + 
                     s2*(2 - 6*t1 - 3*t2) + (-1 + t2)*t2) - 
                  2*Power(s2,3)*
                   (23 + 23*t2 + Power(t2,2) - 
                     Power(t1,2)*(1 + 3*t2) + 
                     t1*(-22 + 4*t2 + 3*Power(t2,2))) + 
                  Power(s1,3)*
                   (-4 + 10*Power(t1,3) + Power(t1,2)*(26 - 36*t2) + 
                     Power(s2,2)*(8 - 12*t1 - 3*t2) - 9*t2 + 
                     12*Power(t2,2) - 4*Power(t2,3) + 
                     s2*(-1 + t1 - 6*Power(t1,2) + 3*t2 + 6*t1*t2 + 
                       15*Power(t2,2)) + 
                     t1*(17 - 58*t2 + 27*Power(t2,2))) + 
                  Power(s2,2)*
                   (44 - 4*Power(t1,3) + 101*t2 + 67*Power(t2,2) - 
                     5*Power(t2,3) + 3*Power(t2,4) + 
                     Power(t1,2)*(-12 + t2 - 6*Power(t2,2)) + 
                     t1*(-28 - 38*t2 + 17*Power(t2,2) + 
                       3*Power(t2,3))) + 
                  s2*(83 + 34*t2 - 82*Power(t2,2) + 16*Power(t2,3) + 
                     Power(t2,4) - 6*Power(t2,5) + 
                     Power(t1,3)*(1 + 10*t2) + 
                     Power(t1,2)*
                      (113 - 64*t2 + 11*Power(t2,2) - 6*Power(t2,3)) \
+ t1*(-197 - 12*t2 + 23*Power(t2,2) - 22*Power(t2,3) + 
                       12*Power(t2,4))) + 
                  Power(s1,2)*
                   (15 - 2*Power(s2,3)*(5 + 3*t1) + 8*t2 + 
                     31*Power(t2,2) - 30*Power(t2,3) + 
                     6*Power(t2,4) - 9*Power(t1,3)*(-3 + 2*t2) + 
                     Power(t1,2)*(37 - 99*t2 + 45*Power(t2,2)) - 
                     t1*(53 + 46*t2 - 117*Power(t2,2) + 
                       33*Power(t2,3)) + 
                     Power(s2,2)*
                      (85 - 12*Power(t1,2) - 21*t2 + 
                       9*Power(t2,2) + t1*(29 + 27*t2)) + 
                     s2*(-111 + 9*t2 - 11*Power(t2,2) - 
                        27*Power(t2,3) + Power(t1,2)*(-34 + 3*t2) + 
                        t1*(49 - 4*t2 + 18*Power(t2,2)))) + 
                  s1*(89 - 8*Power(s2,4) + 14*t2 - 20*Power(t2,2) - 
                     35*Power(t2,3) + 28*Power(t2,4) - 
                     4*Power(t2,5) - 
                     6*Power(s2,3)*
                      (-9 + Power(t1,2) - 2*t2 - 2*t1*t2) + 
                     Power(t1,3)*(23 - 34*t2 + 9*Power(t2,2)) - 
                     Power(t1,2)*
                      (33 + 82*t2 - 102*Power(t2,2) + 22*Power(t2,3)) \
+ t1*(-87 + 82*t2 + 64*Power(t2,2) - 96*Power(t2,3) + 
                       17*Power(t2,4)) + 
                     Power(s2,2)*
                      (-131 - 152*t2 + 18*Power(t2,2) - 
                        9*Power(t2,3) + Power(t1,2)*(-7 + 18*t2) - 
                        2*t1*(-37 + 23*t2 + 9*Power(t2,2))) + 
                     s2*(21 - 15*Power(t1,3) + 177*t2 - 
                        24*Power(t2,2) + 5*Power(t2,3) + 
                        21*Power(t2,4) + 
                        3*Power(t1,2)*(27 + 5*t2 + 3*Power(t2,2)) - 
                        t1*(55 + 48*t2 - 25*Power(t2,2) + 
                        30*Power(t2,3))))) + 
               s*(103 - 278*t1 + 248*Power(t1,2) - 74*Power(t1,3) + 
                  Power(t1,4) + 134*t2 - 181*t1*t2 + 
                  13*Power(t1,2)*t2 + 37*Power(t1,3)*t2 - 
                  3*Power(t1,4)*t2 - 5*Power(t2,2) + 
                  73*t1*Power(t2,2) - 37*Power(t1,2)*Power(t2,2) - 
                  7*Power(t1,3)*Power(t2,2) - 28*Power(t2,3) + 
                  14*t1*Power(t2,3) + 22*Power(t1,2)*Power(t2,3) + 
                  2*Power(t1,3)*Power(t2,3) - 5*Power(t2,4) - 
                  21*t1*Power(t2,4) - 6*Power(t1,2)*Power(t2,4) + 
                  9*Power(t2,5) + 6*t1*Power(t2,5) - 2*Power(t2,6) + 
                  16*Power(s2,4)*t2*(1 - t1 + t2) + 
                  Power(s1,5)*
                   (s2*(-1 + 4*t1 + t2) + t1*(-3 - 5*t1 + 2*t2)) + 
                  Power(s1,4)*
                   (1 - 5*Power(t1,3) + 2*t2 - 2*Power(t2,2) + 
                     Power(s2,2)*(5 + 8*t1 + 3*t2) + 
                     Power(t1,2)*(-13 + 19*t2) + 
                     t1*(-13 + 26*t2 - 8*Power(t2,2)) + 
                     s2*(-8 + 4*Power(t1,2) + t1*(9 - 10*t2) - 
                        6*Power(t2,2))) + 
                  2*Power(s2,3)*
                   (8 - 26*t2 - 31*Power(t2,2) + Power(t2,3) + 
                     Power(t1,2)*(8 + 2*t2 + 3*Power(t2,2)) - 
                     t1*(16 - 24*t2 + 8*Power(t2,2) + 3*Power(t2,3))) \
- Power(s2,2)*(68 + 66*t2 - 36*Power(t2,2) - 51*Power(t2,3) + 
                     4*Power(t2,4) - 3*Power(t2,5) + 
                     Power(t1,3)*(4 + 8*t2) + 
                     2*Power(t1,2)*
                      (30 - 3*t2 - 4*Power(t2,2) + Power(t2,3)) + 
                     t1*(-132 - 68*t2 + 28*Power(t2,2) - 
                        9*Power(t2,3) + Power(t2,4))) + 
                  s2*(-12 + 2*Power(t1,4) + 99*t2 + 138*Power(t2,2) - 
                     19*Power(t2,3) + 4*Power(t2,4) + 4*Power(t2,5) - 
                     2*Power(t2,6) + 
                     Power(t1,3)*(14 - 3*t2 + 5*Power(t2,2)) + 
                     Power(t1,2)*
                      (-46 + 97*t2 - 44*Power(t2,2) + 
                        8*Power(t2,3) - 2*Power(t2,4)) + 
                     t1*(42 - 193*t2 - 99*Power(t2,2) + 
                        29*Power(t2,3) - 17*Power(t2,4) + 
                        4*Power(t2,5))) + 
                  Power(s1,3)*
                   (8 + 2*Power(s2,3)*(7 + 3*t1) + 9*t2 - 
                     15*Power(t2,2) + 8*Power(t2,3) + 
                     3*Power(t1,3)*(-5 + 4*t2) + 
                     Power(t1,2)*(-16 + 54*t2 - 27*Power(t2,2)) + 
                     t1*(5 + 47*t2 - 66*Power(t2,2) + 
                       12*Power(t2,3)) + 
                     Power(s2,2)*
                      (-79 + 8*Power(t1,2) - 11*t2 - 
                        12*Power(t2,2) - t1*(3 + 23*t2)) + 
                     s2*(30 + Power(t1,2)*(37 - 7*t2) + 13*t2 + 
                        2*Power(t2,2) + 14*Power(t2,3) + 
                        t1*(-53 - 14*t2 + 2*Power(t2,2)))) + 
                  Power(s1,2)*
                   (-35 + 16*Power(s2,4) - 30*t2 - 26*Power(t2,2) + 
                     33*Power(t2,3) - 12*Power(t2,4) + 
                     Power(t1,3)*(-1 + 26*t2 - 9*Power(t2,2)) + 
                     Power(t1,2)*
                      (-39 + 64*t2 - 75*Power(t2,2) + 17*Power(t2,3)) \
+ t1*(99 - 14*t2 - 76*Power(t2,2) + 72*Power(t2,3) - 8*Power(t2,4)) + 
                     2*Power(s2,3)*
                      (3*Power(t1,2) - 9*t1*t2 - 13*(3 + t2)) + 
                     Power(s2,2)*
                      (40 + Power(t1,2)*(20 - 18*t2) + 209*t2 + 
                        3*Power(t2,2) + 18*Power(t2,3) + 
                        t1*(-44 + 15*t2 + 21*Power(t2,2))) + 
                     s2*(160 + 15*Power(t1,3) - 63*t2 + 
                        2*Power(t2,2) + 4*Power(t2,3) - 
                        16*Power(t2,4) - 2*Power(t1,2)*(29 + 25*t2) + 
                        t1*(-117 + 103*t2 - 16*Power(t2,2) + 
                        14*Power(t2,3)))) + 
                  s1*(-81 + 6*Power(t1,4) + 
                     16*Power(s2,4)*(-1 + t1 - 2*t2) + 20*t2 + 
                     50*Power(t2,2) + 21*Power(t2,3) - 
                     29*Power(t2,4) + 8*Power(t2,5) + 
                     Power(t1,3)*
                      (-51 + 12*t2 - 13*Power(t2,2) + 2*Power(t2,3)) + 
                     Power(t1,2)*
                      (59 + 48*t2 - 70*Power(t2,2) + 40*Power(t2,3) - 
                        4*Power(t2,4)) + 
                     t1*(67 - 128*t2 - 5*Power(t2,2) + 
                        63*Power(t2,3) - 35*Power(t2,4) + 
                        2*Power(t2,5)) + 
                     2*Power(s2,3)*
                      (26 + 70*t2 + 5*Power(t2,2) - 
                        2*Power(t1,2)*(1 + 3*t2) + 
                        t1*(-24 + 8*t2 + 9*Power(t2,2))) + 
                     Power(s2,2)*
                      (90 + 8*Power(t1,3) - 76*t2 - 181*Power(t2,2) + 
                        7*Power(t2,3) - 12*Power(t2,4) + 
                        2*Power(t1,2)*(9 - 14*t2 + 6*Power(t2,2)) - 
                        t1*(116 - 72*t2 + 21*Power(t2,2) + 
                        5*Power(t2,3))) + 
                     s2*(-157 + Power(t1,3)*(13 - 20*t2) - 290*t2 + 
                        52*Power(t2,2) - 11*Power(t2,3) - 
                        9*Power(t2,4) + 9*Power(t2,5) + 
                        5*Power(t1,2)*
                        (-35 + 22*t2 + Power(t2,2) + Power(t2,3)) + 
                        t1*(319 + 200*t2 - 79*Power(t2,2) + 
                        38*Power(t2,3) - 14*Power(t2,4)))))))/
           ((-1 + s1)*(-s + s1 - t2)*(-1 + s1 + t1 - t2)*
             Power(-1 + s - s*s2 + s1*s2 + t1 - s2*t2,2)) + 
          (4*(64*(-1 + s2)*(-1 + s1 + t1 - t2)*(s - s1 + t2) - 
               (2*(-11 + Power(s,2) + 2*Power(s1,2) + 6*s2 - 
                    3*Power(s2,2) + t1 + s2*t1 + 
                    s1*(3 - 2*s2 + t1 - 3*t2) - 
                    s*(3 + 3*s1 - 2*s2 + t1 - 2*t2) + t2 - t1*t2 + 
                    Power(t2,2))*
                  (-2*s2 + Power(s1,2)*(s2 - t1) + 2*t1 + 2*s2*t1 - 
                    2*Power(t1,2) + 
                    s1*(2*Power(s2,2) - 2*s2*t1 + 
                       (1 + t1)*(-1 + t2)) + t2 - 2*Power(s2,2)*t2 + 
                    t1*t2 + 2*s2*t1*t2 - Power(t2,2) - 
                    s2*Power(t2,2) + 
                    Power(s,2)*(-2 + 2*s2 - t1 + t2) + 
                    s*(3 - 2*Power(s2,2) - t1 + 2*s2*t1 - 3*t2 + 
                       s2*t2 - t1*t2 + Power(t2,2) - 
                       s1*(-2 + 3*s2 - 2*t1 + t2))))/
                (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2) - 
               4*(-1 + 2*Power(s,3) - 20*Power(s1,3) + 6*s2 + 
                  7*Power(s2,2) - Power(s2,3) - 38*t1 - 7*s2*t1 + 
                  14*Power(s2,2)*t1 + 20*Power(t1,2) - 
                  18*s2*Power(t1,2) + 11*t2 + 63*s2*t2 - 
                  42*Power(s2,2)*t2 - t1*t2 - 43*s2*t1*t2 + 
                  38*Power(t1,2)*t2 - 40*Power(t2,2) + 
                  75*s2*Power(t2,2) - 57*t1*Power(t2,2) + 
                  19*Power(t2,3) + 
                  Power(s,2)*(-20 - 4*s1 + 27*s2 - 2*t1 + 5*t2) + 
                  Power(s1,2)*(-50 + 103*s2 - 70*t1 + 58*t2) + 
                  s1*(-14 + 42*Power(s2,2) - 48*Power(t1,2) + 
                     3*s2*(-25 + 18*t1 - 59*t2) + 89*t2 - 
                     57*Power(t2,2) + t1*(11 + 125*t2)) + 
                  s*(-19 + 22*Power(s1,2) - 28*Power(s2,2) + 39*t1 + 
                     48*Power(t1,2) + 
                     s1*(100 - 142*s2 + 72*t1 - 55*t2) - 78*t2 - 
                     79*t1*t2 + 32*Power(t2,2) + 
                     s2*(51 - 84*t1 + 125*t2))) + 
               16*(-9 + 9*Power(s1,3) + 4*s2 - Power(s2,2) + 21*t1 - 
                  10*s2*t1 - 2*Power(s2,2)*t1 - 11*Power(t1,2) + 
                  8*s2*Power(t1,2) + 
                  Power(s,2)*(s1 - 16*s2 + 3*(4 + t1 - t2)) + 
                  Power(s1,2)*(38 - 68*s2 + 33*t1 - 26*t2) + 18*t2 - 
                  50*s2*t2 + 18*Power(s2,2)*t2 - 12*t1*t2 + 
                  41*s2*t1*t2 - 16*Power(t1,2)*t2 + 37*Power(t2,2) - 
                  57*s2*Power(t2,2) + 24*t1*Power(t2,2) - 
                  8*Power(t2,3) + 
                  s*(34 - 10*Power(s1,2) + 16*Power(s2,2) - 28*t1 - 
                     22*Power(t1,2) + s2*(-55 + 57*t1 - 81*t2) + 
                     54*t2 + 39*t1*t2 - 17*Power(t2,2) + 
                     s1*(-66 + 91*s2 - 36*t1 + 28*t2)) + 
                  s1*(-13 - 18*Power(s2,2) + t1 + 22*Power(t1,2) - 
                     76*t2 - 56*t1*t2 + 25*Power(t2,2) + 
                     s2*(51 - 42*t1 + 125*t2))) + 
               32*(2 - Power(s1,3) + 2*Power(s,2)*(-1 + s2) - s2 - 
                  4*t1 + 2*s2*t1 + 2*Power(t1,2) - s2*Power(t1,2) - 
                  7*t2 + 12*s2*t2 - 2*Power(s2,2)*t2 + 6*t1*t2 - 
                  11*s2*t1*t2 + 2*Power(t1,2)*t2 - 10*Power(t2,2) + 
                  13*s2*Power(t2,2) - 3*t1*Power(t2,2) + 
                  Power(t2,3) + 
                  Power(s1,2)*(-10 + 14*s2 - 4*t1 + 3*t2) + 
                  s*(-10 + Power(s1,2) - 2*Power(s2,2) + 8*t1 + 
                     3*Power(t1,2) + s1*(14 - 17*s2 + 4*t1 - 3*t2) - 
                     12*t2 - 5*t1*t2 + 2*Power(t2,2) + 
                     s2*(14 - 13*t1 + 16*t2)) + 
                  s1*(6 + 2*Power(s2,2) - 3*Power(t1,2) + 20*t2 - 
                     3*Power(t2,2) + t1*(-4 + 7*t2) + 
                     s2*(11*t1 - 3*(4 + 9*t2)))) + 
               8*(17 + 2*Power(s,3) - 26*Power(s1,3) + 10*s2 + 
                  Power(s2,2) + Power(s2,3) - 48*t1 + 8*s2*t1 + 
                  10*Power(s2,2)*t1 + 24*Power(t1,2) - 
                  21*s2*Power(t1,2) - 17*t2 + 91*s2*t2 - 
                  49*Power(s2,2)*t2 + 8*t1*t2 - 66*s2*t1*t2 + 
                  43*Power(t1,2)*t2 - 61*Power(t2,2) + 
                  106*s2*Power(t2,2) - 65*t1*Power(t2,2) + 
                  22*Power(t2,3) + 
                  Power(s,2)*(-26 - 8*s1 + 37*s2 - 8*t1 + 11*t2) + 
                  Power(s1,2)*(-71 + 138*s2 - 86*t1 + 73*t2) + 
                  s*(-61 + 32*Power(s1,2) - 40*Power(s2,2) + 53*t1 + 
                     56*Power(t1,2) + 
                     s1*(136 - 191*s2 + 94*t1 - 78*t2) - 106*t2 - 
                     99*t1*t2 + 44*Power(t2,2) + 
                     s2*(92 - 110*t1 + 166*t2)) + 
                  s1*(17 + 51*Power(s2,2) - 56*Power(t1,2) + 132*t2 - 
                     69*Power(t2,2) + 2*t1*(5 + 74*t2) + 
                     s2*(73*t1 - 3*(35 + 81*t2)))) + 
               (2*(-1 + s2)*(-2*Power(s2,2) + 
                    Power(s1,5)*Power(s2 - t1,2) + 4*s2*t1 + 
                    6*Power(s2,2)*t1 - 2*Power(t1,2) - 
                    12*s2*Power(t1,2) - 6*Power(s2,2)*Power(t1,2) + 
                    6*Power(t1,3) + 12*s2*Power(t1,3) + 
                    2*Power(s2,2)*Power(t1,3) - 6*Power(t1,4) - 
                    4*s2*Power(t1,4) + 2*Power(t1,5) + 4*s2*t2 - 
                    2*Power(s2,2)*t2 - 4*Power(s2,3)*t2 - 4*t1*t2 - 
                    2*s2*t1*t2 + 12*Power(s2,2)*t1*t2 + 
                    8*Power(s2,3)*t1*t2 + 4*Power(t1,2)*t2 - 
                    12*s2*Power(t1,2)*t2 - 
                    18*Power(s2,2)*Power(t1,2)*t2 - 
                    4*Power(s2,3)*Power(t1,2)*t2 + 4*Power(t1,3)*t2 + 
                    14*s2*Power(t1,3)*t2 + 
                    8*Power(s2,2)*Power(t1,3)*t2 - 4*Power(t1,4)*t2 - 
                    4*s2*Power(t1,4)*t2 + Power(t2,2) - 
                    4*s2*Power(t2,2) + 4*Power(s2,2)*Power(t2,2) - 
                    4*Power(s2,3)*Power(t2,2) - 
                    2*Power(s2,4)*Power(t2,2) + t1*Power(t2,2) + 
                    4*s2*t1*Power(t2,2) + 
                    8*Power(s2,2)*t1*Power(t2,2) + 
                    8*Power(s2,3)*t1*Power(t2,2) + 
                    2*Power(s2,4)*t1*Power(t2,2) - 
                    7*Power(t1,2)*Power(t2,2) - 
                    8*s2*Power(t1,2)*Power(t2,2) - 
                    14*Power(s2,2)*Power(t1,2)*Power(t2,2) - 
                    4*Power(s2,3)*Power(t1,2)*Power(t2,2) + 
                    5*Power(t1,3)*Power(t2,2) + 
                    8*s2*Power(t1,3)*Power(t2,2) + 
                    2*Power(s2,2)*Power(t1,3)*Power(t2,2) - 
                    3*Power(t2,3) + 4*s2*Power(t2,3) - 
                    6*Power(s2,2)*Power(t2,3) - 
                    2*Power(s2,4)*Power(t2,3) + 2*t1*Power(t2,3) + 
                    6*s2*t1*Power(t2,3) + 
                    6*Power(s2,2)*t1*Power(t2,3) + 
                    6*Power(s2,3)*t1*Power(t2,3) - 
                    3*Power(t1,2)*Power(t2,3) - 
                    8*s2*Power(t1,2)*Power(t2,3) - 
                    4*Power(s2,2)*Power(t1,2)*Power(t2,3) + 
                    3*Power(t2,4) - 4*s2*Power(t2,4) - 
                    Power(s2,2)*Power(t2,4) - 
                    2*Power(s2,3)*Power(t2,4) + t1*Power(t2,4) + 
                    4*s2*t1*Power(t2,4) + 
                    3*Power(s2,2)*t1*Power(t2,4) - Power(t2,5) - 
                    Power(s2,2)*Power(t2,5) + 
                    Power(s,4)*(-1 + s1 + t1 - t2)*
                     (2 + 2*Power(s2,2) + Power(t1,2) - 
                       2*s2*(2 + t1 - t2) - 2*t1*(-1 + t2) - 2*t2 + 
                       Power(t2,2)) + 
                    Power(s1,4)*(s2 - t1)*
                     (-2 + 2*Power(s2,2) - t1 - Power(t1,2) + 2*t2 + 
                       3*t1*t2 - s2*(-1 + t1 + 3*t2)) + 
                    Power(s1,3)*
                     (2*Power(s2,4) - 2*Power(t1,3)*(-2 + t2) + 
                       Power(-1 + t2,2) - 
                       2*Power(s2,3)*(2 + t1 + 2*t2) + 
                       Power(t1,2)*(-3 - 4*t2 + 3*Power(t2,2)) + 
                       t1*(2 - 8*t2 + 6*Power(t2,2)) - 
                       2*Power(s2,2)*
                       (2 + Power(t1,2) + 2*t2 - 2*Power(t2,2) - 
                       t1*(5 + 4*t2)) + 
                       2*s2*
                        (-1 + Power(t1,3) + 4*t2 - 3*Power(t2,2) - 
                        Power(t1,2)*(5 + t2) + 
                        t1*(4 + 3*t2 - 3*Power(t2,2)))) + 
                    Power(s1,2)*
                     (2*Power(t1,4) + 
                       2*Power(s2,4)*(-1 + t1 - 3*t2) - 
                       Power(-1 + t2,2)*(-1 + 3*t2) + 
                       Power(t1,3)*(-5 - 6*t2 + Power(t2,2)) + 
                       t1*(-5 + 11*Power(t2,2) - 6*Power(t2,3)) - 
                       Power(t1,2)*
                       (-7 + t2 - 5*Power(t2,2) + Power(t2,3)) + 
                       2*Power(s2,2)*
                       (4 - 5*Power(t1,2) + Power(t1,3) + t2 + 
                       2*Power(t2,2) - 2*Power(t2,3) - 
                       t1*t2*(7 + 2*t2)) + 
                       Power(s2,3)*
                       (-4 - 4*Power(t1,2) + 8*t2 + 2*t1*(4 + 5*t2)) \
- 2*s2*(-1 - 3*t2 + 7*Power(t2,2) - 3*Power(t2,3) + 
                        Power(t1,3)*(-1 + 2*t2) - 
                        Power(t1,2)*(5 + 5*t2 + 2*Power(t2,2)) + 
                        t1*(7 + 3*t2 + Power(t2,2) - Power(t2,3)))) + 
                    s1*(-2*Power(t1,4)*(-2 + t2) + 
                       Power(-1 + t2,2)*t2*(-2 + 3*t2) + 
                       2*Power(s2,4)*t2*(2 - 2*t1 + 3*t2) + 
                       2*Power(t1,3)*(-2 + Power(t2,2)) + 
                       Power(t1,2)*
                       (-4 + 7*Power(t2,2) - 2*Power(t2,3)) + 
                       2*t1*
                       (2 + 2*t2 - 2*Power(t2,2) - 3*Power(t2,3) + 
                       Power(t2,4)) + 
                       2*s2*
                       (-2 + 2*Power(t1,4) + t2 - 4*Power(t2,2) + 
                       6*Power(t2,3) - Power(t2,4) + 
                       Power(t1,3)*(-7 - 5*t2 + Power(t2,2)) + 
                       t1*(1 + 5*t2 - 4*Power(t2,2) - 
                       3*Power(t2,3)) - 
                       Power(t1,2)*
                       (-6 + t2 - 4*Power(t2,2) + Power(t2,3))) + 
                       2*Power(s2,3)*
                       (Power(t1,2)*(2 + 4*t2) - 
                       t1*(4 + 8*t2 + 7*Power(t2,2)) + 
                       2*(1 + 2*t2 - Power(t2,2) + Power(t2,3))) + 
                       Power(s2,2)*
                        (2 - 12*t2 + 8*Power(t2,2) + 3*Power(t2,4) - 
                        4*Power(t1,3)*(2 + t2) + 
                        6*Power(t1,2)*(3 + 4*t2 + Power(t2,2)) - 
                        2*t1*(6 + 4*t2 + Power(t2,2) + 2*Power(t2,3))\
)) - 2*Power(s,3)*(-3 + 3*t1 + Power(t1,3) + 
                       2*Power(s2,3)*(-1 + t1 - t2) - 
                       Power(s2,2)*(2 + 3*t1)*(-1 + t1 - t2) + 
                       2*t2 - 2*t1*t2 - 5*Power(t1,2)*t2 - 
                       Power(t1,3)*t2 + 2*Power(t2,2) + 
                       7*t1*Power(t2,2) + 
                       3*Power(t1,2)*Power(t2,2) - 3*Power(t2,3) - 
                       3*t1*Power(t2,3) + Power(t2,4) + 
                       Power(s1,2)*
                       (2 + 3*Power(s2,2) + 2*Power(t1,2) + 
                       t1*(4 - 3*t2) - 2*t2 + Power(t2,2) + 
                       s2*(-5 - 4*t1 + 2*t2)) + 
                       s2*(3 + Power(t1,3) - 
                       Power(t1,2)*(-1 + t2) - 2*t2 - 
                       5*Power(t2,2) + Power(t2,3) - 
                       t1*(6 - 4*t2 + Power(t2,2))) + 
                       s1*(1 + 2*Power(s2,3) + 2*Power(t1,3) + 
                        Power(t1,2)*(3 - 6*t2) - 5*t2 + 
                        5*Power(t2,2) - 2*Power(t2,3) - 
                        Power(s2,2)*(5 + 3*t2) + 
                        t1*(-2 - 8*t2 + 6*Power(t2,2)) + 
                        s2*(2 - 3*Power(t1,2) + 8*t2 - 
                        3*Power(t2,2) + t1*(2 + 6*t2)))) + 
                    Power(s,2)*
                     (-5 + 13*t1 - 15*Power(t1,2) + 5*Power(t1,3) + 
                       2*Power(t1,4) + 
                       2*Power(s2,4)*(-1 + t1 - t2) + 7*t2 - 
                       6*t1*t2 - 9*Power(t1,2)*t2 - 
                       8*Power(t1,3)*t2 + 3*Power(t2,2) + 
                       15*t1*Power(t2,2) + 
                       17*Power(t1,2)*Power(t2,2) + 
                       Power(t1,3)*Power(t2,2) - 11*Power(t2,3) - 
                       18*t1*Power(t2,3) - 
                       3*Power(t1,2)*Power(t2,3) + 7*Power(t2,4) + 
                       3*t1*Power(t2,4) - Power(t2,5) + 
                       Power(s1,3)*
                       (2 + 7*Power(s2,2) + 6*Power(t1,2) + 
                       t1*(10 - 6*t2) - 2*s2*(4 + 6*t1 - t2) - 
                       2*t2 + Power(t2,2)) + 
                       Power(s2,3)*
                       (-4*Power(t1,2) - 2*t1*(-2 + t2) + 
                       6*t2*(1 + t2)) - 
                       2*s2*
                       (2 + 2*t2 + 4*Power(t2,2) - 4*Power(t2,3) - 
                       Power(t2,4) + 2*Power(t1,3)*(2 + t2) + 
                       t1*t2*(-10 + t2 + 4*Power(t2,2)) - 
                       Power(t1,2)*(6 + t2 + 5*Power(t2,2))) + 
                       Power(s2,2)*
                       (10 + 2*Power(t1,3) - 8*t2 - 9*Power(t2,2) + 
                       5*Power(t2,3) + Power(t1,2)*(4 + 8*t2) - 
                       t1*(16 + 2*t2 + 15*Power(t2,2))) + 
                       Power(s1,2)*
                       (6 + 10*Power(s2,3) + 6*Power(t1,3) - 
                       14*t2 + 11*Power(t2,2) - 3*Power(t2,3) - 
                       Power(s2,2)*(15 + 9*t1 + 7*t2) - 
                       2*Power(t1,2)*(-5 + 9*t2) + 
                       t1*(6 - 34*t2 + 15*Power(t2,2)) - 
                       2*s2*
                       (2 + 3*Power(t1,2) - 10*t2 + Power(t2,2) - 
                       t1*(5 + 9*t2))) + 
                       s1*(-3 + 2*Power(s2,4) + 
                        2*Power(s2,3)*(-5 + 3*t1 - 8*t2) + 
                        Power(t1,3)*(8 - 6*t2) - 10*t2 + 
                        23*Power(t2,2) - 16*Power(t2,3) + 
                        3*Power(t2,4) + 
                        Power(t1,2)*(21 - 34*t2 + 15*Power(t2,2)) - 
                        2*t1*
                       (5 + 7*t2 - 21*Power(t2,2) + 6*Power(t2,3)) + 
                        Power(s2,2)*
                       (-2 - 14*Power(t1,2) + 22*t2 - 
                       5*Power(t2,2) + 2*t1*(9 + 13*t2)) + 
                        2*s2*
                        (8 + 3*Power(t1,3) + 3*t2 - 10*Power(t2,2) - 
                        Power(t2,3) - Power(t1,2)*(7 + 3*t2) + 
                        t1*(-11 + Power(t2,2))))) - 
                    2*s*(4*t1 - 9*Power(t1,2) + 6*Power(t1,3) - 
                       Power(t1,4) + 
                       Power(s1,4)*
                        (2*Power(s2,2) - s2*(1 + 4*t1) + 
                        t1*(2 + 2*t1 - t2)) + t2 - 2*t1*t2 + 
                       4*Power(t1,2)*t2 - 2*Power(t1,3)*t2 - 
                       Power(t1,4)*t2 + t1*Power(t2,2) + 
                       5*Power(t1,2)*Power(t2,2) + 
                       3*Power(t1,3)*Power(t2,2) - 4*Power(t2,3) - 
                       6*t1*Power(t2,3) - 4*Power(t1,2)*Power(t2,3) + 
                       4*Power(t2,4) + 3*t1*Power(t2,4) - 
                       Power(t2,5) + 2*Power(s2,4)*t2*(1 - t1 + t2) + 
                       Power(s2,3)*
                        (2 + 2*t2 - Power(t2,2) + 
                        Power(t1,2)*(2 + 4*t2) - 
                        2*t1*(2 + 3*t2 + 2*Power(t2,2))) + 
                       Power(s1,3)*
                        (4*Power(s2,3) + 2*Power(t1,3) + 
                        Power(t1,2)*(3 - 6*t2) + Power(-1 + t2,2) - 
                        Power(s2,2)*(3 + 5*t1 + 3*t2) + 
                        s2*(-3 + t1 - Power(t1,2) + 4*t2 + 
                       8*t1*t2 + Power(t2,2)) + 
                        t1*(6 - 11*t2 + 3*Power(t2,2))) + 
                       Power(s2,2)*
                        (-2*Power(t1,3)*(2 + t2) + 
                        Power(t1,2)*(8 + 5*t2 + Power(t2,2)) - 
                        t2*
                       (7 - 5*t2 - 4*Power(t2,2) + Power(t2,3)) + 
                        t1*(-4 + 4*t2 - 3*Power(t2,2) + 
                       2*Power(t2,3))) + 
                       s2*(-4 + 2*Power(t1,4) + 2*t2 - 
                        3*Power(t2,2) + 6*Power(t2,3) - 
                        Power(t2,5) + 
                        Power(t1,3)*(-3 + Power(t2,2)) - 
                        Power(t1,2)*
                       (4 + 4*t2 - 2*Power(t2,2) + 3*Power(t2,3)) + 
                        t1*(9 + 2*t2 - 7*Power(t2,2) - 
                        4*Power(t2,3) + 3*Power(t2,4))) + 
                       Power(s1,2)*
                        (1 + 2*Power(s2,4) + 
                        Power(t1,3)*(5 - 3*t2) - 6*t2 + 
                        8*Power(t2,2) - 3*Power(t2,3) - 
                        Power(s2,3)*(5 + 8*t2) + 
                        2*Power(t1,2)*(4 - 7*t2 + 3*Power(t2,2)) - 
                        t1*(5 + 15*t2 - 19*Power(t2,2) + 
                       3*Power(t2,3)) - 
                        Power(s2,2)*
                       (5 + 5*Power(t1,2) - 8*t2 + Power(t2,2) - 
                       t1*(13 + 14*t2)) + 
                        s2*(4 + 3*Power(t1,3) + 8*t2 - 
                        5*Power(t2,2) - 3*Power(t2,3) - 
                        Power(t1,2)*(13 + 3*t2) - 
                        t1*(1 + Power(t2,2)))) + 
                       s1*(-1 + 2*Power(t1,4) + 
                        2*Power(s2,4)*(-1 + t1 - 2*t2) - t2 + 
                        9*Power(t2,2) - 10*Power(t2,3) + 
                        3*Power(t2,4) + 
                        Power(t1,3)*(-3 - 7*t2 + Power(t2,2)) + 
                        Power(t1,2)*
                        (3 - 15*t2 + 15*Power(t2,2) - 2*Power(t2,3)) \
+ t1*(-1 + 5*t2 + 15*Power(t2,2) - 13*Power(t2,3) + Power(t2,4)) + 
                        Power(s2,3)*
                        (-2 - 4*Power(t1,2) + 6*t2 + 4*Power(t2,2) + 
                        t1*(6 + 4*t2)) + 
                        s2*(1 - 2*t2 - 11*Power(t2,2) + 
                        2*Power(t2,3) + 3*Power(t2,4) - 
                        Power(t1,3)*(3 + 4*t2) + 
                        Power(t1,2)*(13 + 10*t2 + 7*Power(t2,2)) + 
                        t1*(-11 + 10*t2 + 3*Power(t2,2) - 
                        6*Power(t2,3))) + 
                        Power(s2,2)*
                        (2*Power(t1,3) + Power(t1,2)*(-3 + 4*t2) - 
                        t1*(8 + 10*t2 + 11*Power(t2,2)) + 
                        3*(3 - 3*Power(t2,2) + Power(t2,3)))))))/
                ((-1 + s1 + t1 - t2)*
                  Power(-1 + s - s*s2 + s1*s2 + t1 - s2*t2,2))))/
           ((-1 + s2)*(-1 + t1)) + 
          (8*(32*(-1 + s1 + t1 - t2)*(s - s1 + t2)*
                (-1 + s2 - t1 + t2) - 
               ((-2*s2 + Power(s1,2)*(s2 - t1) + 2*t1 + 2*s2*t1 - 
                    2*Power(t1,2) + 
                    s1*(2*Power(s2,2) - 2*s2*t1 + 
                       (1 + t1)*(-1 + t2)) + t2 - 
                    2*Power(s2,2)*t2 + t1*t2 + 2*s2*t1*t2 - 
                    Power(t2,2) - s2*Power(t2,2) + 
                    Power(s,2)*(-2 + 2*s2 - t1 + t2) + 
                    s*(3 - 2*Power(s2,2) - t1 + 2*s2*t1 - 3*t2 + 
                       s2*t2 - t1*t2 + Power(t2,2) - 
                       s1*(-2 + 3*s2 - 2*t1 + t2)))*
                  (-15 + 6*Power(s,2) + 3*Power(s1,2) - 6*s2 + 
                    3*Power(s2,2) + 6*t1 - 2*s2*t1 - 6*s2*t2 + 
                    4*t1*t2 - Power(t2,2) + 
                    s*(5 - 9*s1 - 9*s2 + 4*t1 + 5*t2) + 
                    s1*(8*s2 - 2*(2 + 2*t1 + t2))))/
                (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2) + 
               16*(2 + 2*Power(s,3) - s2 - 2*t1 + 3*s2*t1 - 
                  2*Power(t1,2) + 
                  Power(s1,2)*(12*s2 - 12*t1 + 11*(-1 + t2)) - 
                  11*t2 + 12*s2*t2 + 5*t1*t2 - 14*s2*t1*t2 + 
                  13*Power(t1,2)*t2 - 5*Power(t2,2) + 
                  12*s2*Power(t2,2) - 24*t1*Power(t2,2) + 
                  12*Power(t2,3) + 
                  Power(s,2)*(-2 - 3*s1 - 2*s2 + 2*t1 + 3*t2) + 
                  s1*(8 - 13*Power(t1,2) + 
                     s2*(-13 + 13*t1 - 24*t2) + 17*t2 + 35*t1*t2 - 
                     23*Power(t2,2)) + 
                  s*(-12 + Power(s1,2) + 13*Power(t1,2) + 
                     s1*(15 - 11*s2 + 10*t1 - 14*t2) - 4*t2 - 
                     22*t1*t2 + 13*Power(t2,2) + 
                     s2*(15 - 14*t1 + 10*t2))) + 
               4*(30 + 48*Power(s,3) - Power(s1,3) - 22*s2 - 
                  11*Power(s2,2) - 2*Power(s2,3) - 21*t1 + 76*s2*t1 + 
                  8*Power(s2,2)*t1 - 52*Power(t1,2) - 
                  6*s2*Power(t1,2) - 60*t2 + 119*s2*t2 - 
                  Power(s2,2)*t2 + 64*t1*t2 - 138*s2*t1*t2 + 
                  125*Power(t1,2)*t2 - 77*Power(t2,2) + 
                  104*s2*Power(t2,2) - 205*t1*Power(t2,2) + 
                  97*Power(t2,3) + 
                  Power(s,2)*(-88 - 69*s1 - 54*s2 + 54*t1 + 55*t2) + 
                  Power(s1,2)*(-89 + 112*s2 - 105*t1 + 86*t2) + 
                  s1*(36 + 7*Power(s2,2) - 127*Power(t1,2) + 
                     s2*(-171 + 128*t1 - 208*t2) + 171*t2 - 
                     182*Power(t2,2) + 9*t1*(2 + 33*t2)) + 
                  s*(-67 + 22*Power(s1,2) + 8*Power(s2,2) - 46*t1 + 
                     127*Power(t1,2) + 
                     s1*(199 - 73*s2 + 51*t1 - 128*t2) - 103*t2 - 
                     155*t1*t2 + 106*Power(t2,2) + 
                     s2*(197 - 153*t1 + 59*t2))) - 
               8*(13 + 17*Power(s,3) - 6*s2 - 3*Power(s2,2) - 13*t1 + 
                  26*s2*t1 + 2*Power(s2,2)*t1 - 17*Power(t1,2) - 
                  2*s2*Power(t1,2) - 41*t2 + 56*s2*t2 - 
                  Power(s2,2)*t2 + 32*t1*t2 - 66*s2*t1*t2 + 
                  61*Power(t1,2)*t2 - 35*Power(t2,2) + 
                  53*s2*Power(t2,2) - 106*t1*Power(t2,2) + 
                  52*Power(t2,3) + 
                  Power(s,2)*(-23 - 25*s1 - 18*s2 + 18*t1 + 22*t2) + 
                  Power(s1,2)*(-47 + 56*s2 - 54*t1 + 46*t2) + 
                  s*(-49 + 8*Power(s1,2) + Power(s2,2) - 7*t1 + 
                     62*Power(t1,2) + 
                     s1*(83 - 45*s2 + 36*t1 - 66*t2) - 36*t2 - 
                     90*t1*t2 + 58*Power(t2,2) + 
                     s2*(82 - 71*t1 + 38*t2)) + 
                  s1*(26 + 3*Power(s2,2) - 62*Power(t1,2) + 
                     s2*(-73 + 61*t1 - 107*t2) + 86*t2 - 
                     98*Power(t2,2) + 2*t1*(2 + 77*t2))) - 
               2*(8 + 38*Power(s,3) - Power(s1,3) - 62*s2 - 
                  14*Power(s2,2) - Power(s2,3) + 17*t1 + 81*s2*t1 + 
                  4*Power(s2,2)*t1 - 56*Power(t1,2) - 
                  3*s2*Power(t1,2) - 27*t2 + 96*s2*t2 - 
                  4*Power(s2,2)*t2 + 39*t1*t2 - 100*s2*t1*t2 + 
                  92*Power(t1,2)*t2 - 54*Power(t2,2) + 
                  74*s2*Power(t2,2) - 145*t1*Power(t2,2) + 
                  66*Power(t2,3) + 
                  Power(s,2)*(-101 - 50*s1 - 37*s2 + 44*t1 + 42*t2) + 
                  Power(s1,2)*(-67 + 87*s2 - 79*t1 + 60*t2) + 
                  s*(13 + 13*Power(s1,2) - 73*t1 + 92*Power(t1,2) - 
                     95*t2 - 101*t1*t2 + 70*Power(t2,2) - 
                     5*s1*(-35 + 12*s2 - 7*t1 + 17*t2) + 
                     s2*(174 - 106*t1 + 40*t2)) + 
                  s1*(2 + 12*Power(s2,2) - 92*Power(t1,2) + 
                     s2*(-145 + 89*t1 - 152*t2) + 122*t2 - 
                     125*Power(t2,2) + t1*(24 + 214*t2))) + 
               ((3 + 2*s - 2*s1 - s2 + 2*t2)*
                  (-2*Power(s2,2) + Power(s1,5)*Power(s2 - t1,2) + 
                    4*s2*t1 + 6*Power(s2,2)*t1 - 2*Power(t1,2) - 
                    12*s2*Power(t1,2) - 6*Power(s2,2)*Power(t1,2) + 
                    6*Power(t1,3) + 12*s2*Power(t1,3) + 
                    2*Power(s2,2)*Power(t1,3) - 6*Power(t1,4) - 
                    4*s2*Power(t1,4) + 2*Power(t1,5) + 4*s2*t2 - 
                    2*Power(s2,2)*t2 - 4*Power(s2,3)*t2 - 4*t1*t2 - 
                    2*s2*t1*t2 + 12*Power(s2,2)*t1*t2 + 
                    8*Power(s2,3)*t1*t2 + 4*Power(t1,2)*t2 - 
                    12*s2*Power(t1,2)*t2 - 
                    18*Power(s2,2)*Power(t1,2)*t2 - 
                    4*Power(s2,3)*Power(t1,2)*t2 + 4*Power(t1,3)*t2 + 
                    14*s2*Power(t1,3)*t2 + 
                    8*Power(s2,2)*Power(t1,3)*t2 - 4*Power(t1,4)*t2 - 
                    4*s2*Power(t1,4)*t2 + Power(t2,2) - 
                    4*s2*Power(t2,2) + 4*Power(s2,2)*Power(t2,2) - 
                    4*Power(s2,3)*Power(t2,2) - 
                    2*Power(s2,4)*Power(t2,2) + t1*Power(t2,2) + 
                    4*s2*t1*Power(t2,2) + 
                    8*Power(s2,2)*t1*Power(t2,2) + 
                    8*Power(s2,3)*t1*Power(t2,2) + 
                    2*Power(s2,4)*t1*Power(t2,2) - 
                    7*Power(t1,2)*Power(t2,2) - 
                    8*s2*Power(t1,2)*Power(t2,2) - 
                    14*Power(s2,2)*Power(t1,2)*Power(t2,2) - 
                    4*Power(s2,3)*Power(t1,2)*Power(t2,2) + 
                    5*Power(t1,3)*Power(t2,2) + 
                    8*s2*Power(t1,3)*Power(t2,2) + 
                    2*Power(s2,2)*Power(t1,3)*Power(t2,2) - 
                    3*Power(t2,3) + 4*s2*Power(t2,3) - 
                    6*Power(s2,2)*Power(t2,3) - 
                    2*Power(s2,4)*Power(t2,3) + 2*t1*Power(t2,3) + 
                    6*s2*t1*Power(t2,3) + 
                    6*Power(s2,2)*t1*Power(t2,3) + 
                    6*Power(s2,3)*t1*Power(t2,3) - 
                    3*Power(t1,2)*Power(t2,3) - 
                    8*s2*Power(t1,2)*Power(t2,3) - 
                    4*Power(s2,2)*Power(t1,2)*Power(t2,3) + 
                    3*Power(t2,4) - 4*s2*Power(t2,4) - 
                    Power(s2,2)*Power(t2,4) - 
                    2*Power(s2,3)*Power(t2,4) + t1*Power(t2,4) + 
                    4*s2*t1*Power(t2,4) + 
                    3*Power(s2,2)*t1*Power(t2,4) - Power(t2,5) - 
                    Power(s2,2)*Power(t2,5) + 
                    Power(s,4)*(-1 + s1 + t1 - t2)*
                     (2 + 2*Power(s2,2) + Power(t1,2) - 
                       2*s2*(2 + t1 - t2) - 2*t1*(-1 + t2) - 2*t2 + 
                       Power(t2,2)) + 
                    Power(s1,4)*(s2 - t1)*
                     (-2 + 2*Power(s2,2) - t1 - Power(t1,2) + 2*t2 + 
                       3*t1*t2 - s2*(-1 + t1 + 3*t2)) + 
                    Power(s1,3)*
                     (2*Power(s2,4) - 2*Power(t1,3)*(-2 + t2) + 
                       Power(-1 + t2,2) - 
                       2*Power(s2,3)*(2 + t1 + 2*t2) + 
                       Power(t1,2)*(-3 - 4*t2 + 3*Power(t2,2)) + 
                       t1*(2 - 8*t2 + 6*Power(t2,2)) - 
                       2*Power(s2,2)*
                       (2 + Power(t1,2) + 2*t2 - 2*Power(t2,2) - 
                       t1*(5 + 4*t2)) + 
                       2*s2*
                        (-1 + Power(t1,3) + 4*t2 - 3*Power(t2,2) - 
                        Power(t1,2)*(5 + t2) + 
                        t1*(4 + 3*t2 - 3*Power(t2,2)))) + 
                    Power(s1,2)*
                     (2*Power(t1,4) + 
                       2*Power(s2,4)*(-1 + t1 - 3*t2) - 
                       Power(-1 + t2,2)*(-1 + 3*t2) + 
                       Power(t1,3)*(-5 - 6*t2 + Power(t2,2)) + 
                       t1*(-5 + 11*Power(t2,2) - 6*Power(t2,3)) - 
                       Power(t1,2)*
                       (-7 + t2 - 5*Power(t2,2) + Power(t2,3)) + 
                       2*Power(s2,2)*
                       (4 - 5*Power(t1,2) + Power(t1,3) + t2 + 
                       2*Power(t2,2) - 2*Power(t2,3) - 
                       t1*t2*(7 + 2*t2)) + 
                       Power(s2,3)*
                       (-4 - 4*Power(t1,2) + 8*t2 + 2*t1*(4 + 5*t2)) \
- 2*s2*(-1 - 3*t2 + 7*Power(t2,2) - 3*Power(t2,3) + 
                        Power(t1,3)*(-1 + 2*t2) - 
                        Power(t1,2)*(5 + 5*t2 + 2*Power(t2,2)) + 
                        t1*(7 + 3*t2 + Power(t2,2) - Power(t2,3)))) + 
                    s1*(-2*Power(t1,4)*(-2 + t2) + 
                       Power(-1 + t2,2)*t2*(-2 + 3*t2) + 
                       2*Power(s2,4)*t2*(2 - 2*t1 + 3*t2) + 
                       2*Power(t1,3)*(-2 + Power(t2,2)) + 
                       Power(t1,2)*
                       (-4 + 7*Power(t2,2) - 2*Power(t2,3)) + 
                       2*t1*
                       (2 + 2*t2 - 2*Power(t2,2) - 3*Power(t2,3) + 
                       Power(t2,4)) + 
                       2*s2*
                       (-2 + 2*Power(t1,4) + t2 - 4*Power(t2,2) + 
                       6*Power(t2,3) - Power(t2,4) + 
                       Power(t1,3)*(-7 - 5*t2 + Power(t2,2)) + 
                       t1*(1 + 5*t2 - 4*Power(t2,2) - 
                       3*Power(t2,3)) - 
                       Power(t1,2)*
                       (-6 + t2 - 4*Power(t2,2) + Power(t2,3))) + 
                       2*Power(s2,3)*
                       (Power(t1,2)*(2 + 4*t2) - 
                       t1*(4 + 8*t2 + 7*Power(t2,2)) + 
                       2*(1 + 2*t2 - Power(t2,2) + Power(t2,3))) + 
                       Power(s2,2)*
                        (2 - 12*t2 + 8*Power(t2,2) + 3*Power(t2,4) - 
                        4*Power(t1,3)*(2 + t2) + 
                        6*Power(t1,2)*(3 + 4*t2 + Power(t2,2)) - 
                        2*t1*(6 + 4*t2 + Power(t2,2) + 2*Power(t2,3))\
)) - 2*Power(s,3)*(-3 + 3*t1 + Power(t1,3) + 
                       2*Power(s2,3)*(-1 + t1 - t2) - 
                       Power(s2,2)*(2 + 3*t1)*(-1 + t1 - t2) + 
                       2*t2 - 2*t1*t2 - 5*Power(t1,2)*t2 - 
                       Power(t1,3)*t2 + 2*Power(t2,2) + 
                       7*t1*Power(t2,2) + 
                       3*Power(t1,2)*Power(t2,2) - 3*Power(t2,3) - 
                       3*t1*Power(t2,3) + Power(t2,4) + 
                       Power(s1,2)*
                       (2 + 3*Power(s2,2) + 2*Power(t1,2) + 
                       t1*(4 - 3*t2) - 2*t2 + Power(t2,2) + 
                       s2*(-5 - 4*t1 + 2*t2)) + 
                       s2*(3 + Power(t1,3) - 
                       Power(t1,2)*(-1 + t2) - 2*t2 - 
                       5*Power(t2,2) + Power(t2,3) - 
                       t1*(6 - 4*t2 + Power(t2,2))) + 
                       s1*(1 + 2*Power(s2,3) + 2*Power(t1,3) + 
                        Power(t1,2)*(3 - 6*t2) - 5*t2 + 
                        5*Power(t2,2) - 2*Power(t2,3) - 
                        Power(s2,2)*(5 + 3*t2) + 
                        t1*(-2 - 8*t2 + 6*Power(t2,2)) + 
                        s2*(2 - 3*Power(t1,2) + 8*t2 - 
                        3*Power(t2,2) + t1*(2 + 6*t2)))) + 
                    Power(s,2)*
                     (-5 + 13*t1 - 15*Power(t1,2) + 5*Power(t1,3) + 
                       2*Power(t1,4) + 
                       2*Power(s2,4)*(-1 + t1 - t2) + 7*t2 - 
                       6*t1*t2 - 9*Power(t1,2)*t2 - 
                       8*Power(t1,3)*t2 + 3*Power(t2,2) + 
                       15*t1*Power(t2,2) + 
                       17*Power(t1,2)*Power(t2,2) + 
                       Power(t1,3)*Power(t2,2) - 11*Power(t2,3) - 
                       18*t1*Power(t2,3) - 
                       3*Power(t1,2)*Power(t2,3) + 7*Power(t2,4) + 
                       3*t1*Power(t2,4) - Power(t2,5) + 
                       Power(s1,3)*
                       (2 + 7*Power(s2,2) + 6*Power(t1,2) + 
                       t1*(10 - 6*t2) - 2*s2*(4 + 6*t1 - t2) - 
                       2*t2 + Power(t2,2)) + 
                       Power(s2,3)*
                       (-4*Power(t1,2) - 2*t1*(-2 + t2) + 
                       6*t2*(1 + t2)) - 
                       2*s2*
                       (2 + 2*t2 + 4*Power(t2,2) - 4*Power(t2,3) - 
                       Power(t2,4) + 2*Power(t1,3)*(2 + t2) + 
                       t1*t2*(-10 + t2 + 4*Power(t2,2)) - 
                       Power(t1,2)*(6 + t2 + 5*Power(t2,2))) + 
                       Power(s2,2)*
                       (10 + 2*Power(t1,3) - 8*t2 - 9*Power(t2,2) + 
                       5*Power(t2,3) + Power(t1,2)*(4 + 8*t2) - 
                       t1*(16 + 2*t2 + 15*Power(t2,2))) + 
                       Power(s1,2)*
                       (6 + 10*Power(s2,3) + 6*Power(t1,3) - 
                       14*t2 + 11*Power(t2,2) - 3*Power(t2,3) - 
                       Power(s2,2)*(15 + 9*t1 + 7*t2) - 
                       2*Power(t1,2)*(-5 + 9*t2) + 
                       t1*(6 - 34*t2 + 15*Power(t2,2)) - 
                       2*s2*
                       (2 + 3*Power(t1,2) - 10*t2 + Power(t2,2) - 
                       t1*(5 + 9*t2))) + 
                       s1*(-3 + 2*Power(s2,4) + 
                        2*Power(s2,3)*(-5 + 3*t1 - 8*t2) + 
                        Power(t1,3)*(8 - 6*t2) - 10*t2 + 
                        23*Power(t2,2) - 16*Power(t2,3) + 
                        3*Power(t2,4) + 
                        Power(t1,2)*(21 - 34*t2 + 15*Power(t2,2)) - 
                        2*t1*
                       (5 + 7*t2 - 21*Power(t2,2) + 6*Power(t2,3)) + 
                        Power(s2,2)*
                       (-2 - 14*Power(t1,2) + 22*t2 - 
                       5*Power(t2,2) + 2*t1*(9 + 13*t2)) + 
                        2*s2*
                        (8 + 3*Power(t1,3) + 3*t2 - 10*Power(t2,2) - 
                        Power(t2,3) - Power(t1,2)*(7 + 3*t2) + 
                        t1*(-11 + Power(t2,2))))) - 
                    2*s*(4*t1 - 9*Power(t1,2) + 6*Power(t1,3) - 
                       Power(t1,4) + 
                       Power(s1,4)*
                        (2*Power(s2,2) - s2*(1 + 4*t1) + 
                        t1*(2 + 2*t1 - t2)) + t2 - 2*t1*t2 + 
                       4*Power(t1,2)*t2 - 2*Power(t1,3)*t2 - 
                       Power(t1,4)*t2 + t1*Power(t2,2) + 
                       5*Power(t1,2)*Power(t2,2) + 
                       3*Power(t1,3)*Power(t2,2) - 4*Power(t2,3) - 
                       6*t1*Power(t2,3) - 4*Power(t1,2)*Power(t2,3) + 
                       4*Power(t2,4) + 3*t1*Power(t2,4) - 
                       Power(t2,5) + 2*Power(s2,4)*t2*(1 - t1 + t2) + 
                       Power(s2,3)*
                        (2 + 2*t2 - Power(t2,2) + 
                        Power(t1,2)*(2 + 4*t2) - 
                        2*t1*(2 + 3*t2 + 2*Power(t2,2))) + 
                       Power(s1,3)*
                        (4*Power(s2,3) + 2*Power(t1,3) + 
                        Power(t1,2)*(3 - 6*t2) + Power(-1 + t2,2) - 
                        Power(s2,2)*(3 + 5*t1 + 3*t2) + 
                        s2*(-3 + t1 - Power(t1,2) + 4*t2 + 
                       8*t1*t2 + Power(t2,2)) + 
                        t1*(6 - 11*t2 + 3*Power(t2,2))) + 
                       Power(s2,2)*
                        (-2*Power(t1,3)*(2 + t2) + 
                        Power(t1,2)*(8 + 5*t2 + Power(t2,2)) - 
                        t2*
                       (7 - 5*t2 - 4*Power(t2,2) + Power(t2,3)) + 
                        t1*(-4 + 4*t2 - 3*Power(t2,2) + 
                       2*Power(t2,3))) + 
                       s2*(-4 + 2*Power(t1,4) + 2*t2 - 
                        3*Power(t2,2) + 6*Power(t2,3) - 
                        Power(t2,5) + 
                        Power(t1,3)*(-3 + Power(t2,2)) - 
                        Power(t1,2)*
                       (4 + 4*t2 - 2*Power(t2,2) + 3*Power(t2,3)) + 
                        t1*(9 + 2*t2 - 7*Power(t2,2) - 
                        4*Power(t2,3) + 3*Power(t2,4))) + 
                       Power(s1,2)*
                        (1 + 2*Power(s2,4) + 
                        Power(t1,3)*(5 - 3*t2) - 6*t2 + 
                        8*Power(t2,2) - 3*Power(t2,3) - 
                        Power(s2,3)*(5 + 8*t2) + 
                        2*Power(t1,2)*(4 - 7*t2 + 3*Power(t2,2)) - 
                        t1*(5 + 15*t2 - 19*Power(t2,2) + 
                       3*Power(t2,3)) - 
                        Power(s2,2)*
                       (5 + 5*Power(t1,2) - 8*t2 + Power(t2,2) - 
                       t1*(13 + 14*t2)) + 
                        s2*(4 + 3*Power(t1,3) + 8*t2 - 
                        5*Power(t2,2) - 3*Power(t2,3) - 
                        Power(t1,2)*(13 + 3*t2) - 
                        t1*(1 + Power(t2,2)))) + 
                       s1*(-1 + 2*Power(t1,4) + 
                        2*Power(s2,4)*(-1 + t1 - 2*t2) - t2 + 
                        9*Power(t2,2) - 10*Power(t2,3) + 
                        3*Power(t2,4) + 
                        Power(t1,3)*(-3 - 7*t2 + Power(t2,2)) + 
                        Power(t1,2)*
                        (3 - 15*t2 + 15*Power(t2,2) - 2*Power(t2,3)) \
+ t1*(-1 + 5*t2 + 15*Power(t2,2) - 13*Power(t2,3) + Power(t2,4)) + 
                        Power(s2,3)*
                        (-2 - 4*Power(t1,2) + 6*t2 + 4*Power(t2,2) + 
                        t1*(6 + 4*t2)) + 
                        s2*(1 - 2*t2 - 11*Power(t2,2) + 
                        2*Power(t2,3) + 3*Power(t2,4) - 
                        Power(t1,3)*(3 + 4*t2) + 
                        Power(t1,2)*(13 + 10*t2 + 7*Power(t2,2)) + 
                        t1*(-11 + 10*t2 + 3*Power(t2,2) - 
                        6*Power(t2,3))) + 
                        Power(s2,2)*
                        (2*Power(t1,3) + Power(t1,2)*(-3 + 4*t2) - 
                        t1*(8 + 10*t2 + 11*Power(t2,2)) + 
                        3*(3 - 3*Power(t2,2) + Power(t2,3)))))))/
                ((-1 + s1 + t1 - t2)*
                  Power(-1 + s - s*s2 + s1*s2 + t1 - s2*t2,2))))/
           ((-1 + s1)*(-1 + t2)) + 
          (8*(-32*(s2 - t1)*(-1 + s1 + t1 - t2)*(s - s1 + t2) - 
               ((-2*s2 + Power(s1,2)*(s2 - t1) + 2*t1 + 2*s2*t1 - 
                    2*Power(t1,2) + 
                    s1*(2*Power(s2,2) - 2*s2*t1 + 
                       (1 + t1)*(-1 + t2)) + t2 - 2*Power(s2,2)*t2 + 
                    t1*t2 + 2*s2*t1*t2 - Power(t2,2) - 
                    s2*Power(t2,2) + 
                    Power(s,2)*(-2 + 2*s2 - t1 + t2) + 
                    s*(3 - 2*Power(s2,2) - t1 + 2*s2*t1 - 3*t2 + 
                       s2*t2 - t1*t2 + Power(t2,2) - 
                       s1*(-2 + 3*s2 - 2*t1 + t2)))*
                  (-15 + 6*Power(s,2) + 5*Power(s1,2) - 7*s2 + 
                    3*Power(s2,2) + 6*t1 - s2*t1 + 2*t2 - 7*s2*t2 + 
                    2*t1*t2 + Power(t2,2) + 
                    s*(7 - 11*s1 - 9*s2 + 2*t1 + 7*t2) + 
                    s1*(9*s2 - 2*(3 + t1 + 3*t2))))/
                (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2) - 
               16*(Power(s1,3) - 2*s2 + 2*t1 + 2*s2*t1 - 
                  2*Power(t1,2) + 
                  Power(s1,2)*(-1 + 15*s2 - 13*t1 - 2*t2) + 
                  Power(s,2)*(-3 + 2*s1 + 2*s2 + t1 - 2*t2) - 3*t2 + 
                  15*s2*t2 - 2*Power(s2,2)*t2 - 11*t1*t2 - 
                  12*s2*t1*t2 + 13*Power(t1,2)*t2 + Power(t2,2) + 
                  14*s2*Power(t2,2) - 14*t1*Power(t2,2) - 
                  s*(3 + 3*Power(s1,2) + 2*Power(s2,2) + 11*t1 - 
                     13*Power(t1,2) + s2*(-16 + 13*t1 - 17*t2) + 
                     2*s1*(-3 + 9*s2 - 6*t1 - 2*t2) + 3*t2 + 
                     14*t1*t2 + Power(t2,2)) + 
                  s1*(3 + 2*Power(s2,2) - 13*Power(t1,2) + 
                     s2*(-15 + 12*t1 - 29*t2) + Power(t2,2) + 
                     t1*(11 + 27*t2))) - 
               8*(-3 + Power(s,3) - 8*Power(s1,3) + 17*s2 - 
                  Power(s2,2) - 11*t1 - 18*s2*t1 + 15*Power(t1,2) + 
                  s2*Power(t1,2) + 20*t2 - 75*s2*t2 + 
                  19*Power(s2,2)*t2 + 43*t1*t2 + 50*s2*t1*t2 - 
                  60*Power(t1,2)*t2 - 7*Power(t2,2) - 
                  67*s2*Power(t2,2) + 68*t1*Power(t2,2) - 
                  Power(t2,3) + 
                  Power(s,2)*(23 - 18*s1 - 18*s2 - 7*t1 + 15*t2) + 
                  Power(s1,2)*(4 - 79*s2 + 62*t1 + 15*t2) + 
                  s*(17 + 25*Power(s1,2) + 17*Power(s2,2) + 46*t1 - 
                     60*Power(t1,2) + s2*(-78 + 58*t1 - 93*t2) + 
                     s1*(-43 + 104*s2 - 55*t1 - 31*t2) + 23*t2 + 
                     67*t1*t2 + 7*Power(t2,2)) - 
                  s1*(17 + 19*Power(s2,2) - 60*Power(t1,2) + 
                     s2*(-76 + 49*t1 - 146*t2) - 2*t2 + 
                     6*Power(t2,2) + 3*t1*(16 + 43*t2))) + 
               4*(-5 + 5*Power(s,3) - 22*Power(s1,3) + 42*s2 - 
                  7*Power(s2,2) - Power(s2,3) - 17*t1 - 44*s2*t1 - 
                  2*Power(s2,2)*t1 + 35*Power(t1,2) + 
                  5*s2*Power(t1,2) + 36*t2 - 148*s2*t2 + 
                  60*Power(s2,2)*t2 + 67*t1*t2 + 85*s2*t1*t2 - 
                  118*Power(t1,2)*t2 - 14*Power(t2,2) - 
                  135*s2*Power(t2,2) + 142*t1*Power(t2,2) - 
                  7*Power(t2,3) + 
                  Power(s1,2)*(5 - 180*s2 + 129*t1 + 36*t2) + 
                  Power(s,2)*(53 - 54*s1 - 55*s2 - 16*t1 + 40*t2) + 
                  s*(26 + 71*Power(s1,2) + 51*Power(s2,2) + 80*t1 - 
                     119*Power(t1,2) + s2*(-144 + 109*t1 - 213*t2) + 
                     s1*(-97 + 251*s2 - 113*t1 - 82*t2) + 57*t2 + 
                     140*t1*t2 + 15*Power(t2,2)) - 
                  s1*(25 + 62*Power(s2,2) - 119*Power(t1,2) + 
                     2*s2*(-76 + 42*t1 - 157*t2) - 5*t2 + 
                     7*Power(t2,2) + t1*(86 + 268*t2))) + 
               2*(11 + 16*Power(s1,3) - 2*s2 + 15*Power(s2,2) - 
                  Power(s2,3) - 5*t1 + 23*s2*t1 + 6*Power(s2,2)*t1 - 
                  25*Power(t1,2) - 6*s2*Power(t1,2) + 
                  Power(s,2)*(-33 + 36*s1 + 39*s2 + 16*t1 - 28*t2) + 
                  Power(s1,2)*(4 + 135*s2 - 94*t1 - 23*t2) - 32*t2 + 
                  95*s2*t2 - 52*Power(s2,2)*t2 - 31*t1*t2 - 
                  55*s2*t1*t2 + 84*Power(t1,2)*t2 + 9*Power(t2,2) + 
                  96*s2*Power(t2,2) - 106*t1*Power(t2,2) + 
                  9*Power(t2,3) + 
                  s1*(35 + 52*Power(s2,2) - 86*Power(t1,2) + 
                     s2*(-93 + 58*t1 - 230*t2) - 8*t2 - 
                     2*Power(t2,2) + 9*t1*(5 + 22*t2)) - 
                  s*(44 + 52*Power(s1,2) + 38*Power(s2,2) + 40*t1 - 
                     86*Power(t1,2) + s2*(-77 + 80*t1 - 158*t2) + 
                     40*t2 + 102*t1*t2 + 9*Power(t2,2) + 
                     s1*(186*s2 - 78*t1 - 59*(1 + t2)))) + 
               ((3 + 2*s - 2*s1 - s2 + 2*t2)*
                  (-2*Power(s2,2) + Power(s1,5)*Power(s2 - t1,2) + 
                    4*s2*t1 + 6*Power(s2,2)*t1 - 2*Power(t1,2) - 
                    12*s2*Power(t1,2) - 6*Power(s2,2)*Power(t1,2) + 
                    6*Power(t1,3) + 12*s2*Power(t1,3) + 
                    2*Power(s2,2)*Power(t1,3) - 6*Power(t1,4) - 
                    4*s2*Power(t1,4) + 2*Power(t1,5) + 4*s2*t2 - 
                    2*Power(s2,2)*t2 - 4*Power(s2,3)*t2 - 4*t1*t2 - 
                    2*s2*t1*t2 + 12*Power(s2,2)*t1*t2 + 
                    8*Power(s2,3)*t1*t2 + 4*Power(t1,2)*t2 - 
                    12*s2*Power(t1,2)*t2 - 
                    18*Power(s2,2)*Power(t1,2)*t2 - 
                    4*Power(s2,3)*Power(t1,2)*t2 + 4*Power(t1,3)*t2 + 
                    14*s2*Power(t1,3)*t2 + 
                    8*Power(s2,2)*Power(t1,3)*t2 - 4*Power(t1,4)*t2 - 
                    4*s2*Power(t1,4)*t2 + Power(t2,2) - 
                    4*s2*Power(t2,2) + 4*Power(s2,2)*Power(t2,2) - 
                    4*Power(s2,3)*Power(t2,2) - 
                    2*Power(s2,4)*Power(t2,2) + t1*Power(t2,2) + 
                    4*s2*t1*Power(t2,2) + 
                    8*Power(s2,2)*t1*Power(t2,2) + 
                    8*Power(s2,3)*t1*Power(t2,2) + 
                    2*Power(s2,4)*t1*Power(t2,2) - 
                    7*Power(t1,2)*Power(t2,2) - 
                    8*s2*Power(t1,2)*Power(t2,2) - 
                    14*Power(s2,2)*Power(t1,2)*Power(t2,2) - 
                    4*Power(s2,3)*Power(t1,2)*Power(t2,2) + 
                    5*Power(t1,3)*Power(t2,2) + 
                    8*s2*Power(t1,3)*Power(t2,2) + 
                    2*Power(s2,2)*Power(t1,3)*Power(t2,2) - 
                    3*Power(t2,3) + 4*s2*Power(t2,3) - 
                    6*Power(s2,2)*Power(t2,3) - 
                    2*Power(s2,4)*Power(t2,3) + 2*t1*Power(t2,3) + 
                    6*s2*t1*Power(t2,3) + 
                    6*Power(s2,2)*t1*Power(t2,3) + 
                    6*Power(s2,3)*t1*Power(t2,3) - 
                    3*Power(t1,2)*Power(t2,3) - 
                    8*s2*Power(t1,2)*Power(t2,3) - 
                    4*Power(s2,2)*Power(t1,2)*Power(t2,3) + 
                    3*Power(t2,4) - 4*s2*Power(t2,4) - 
                    Power(s2,2)*Power(t2,4) - 
                    2*Power(s2,3)*Power(t2,4) + t1*Power(t2,4) + 
                    4*s2*t1*Power(t2,4) + 
                    3*Power(s2,2)*t1*Power(t2,4) - Power(t2,5) - 
                    Power(s2,2)*Power(t2,5) + 
                    Power(s,4)*(-1 + s1 + t1 - t2)*
                     (2 + 2*Power(s2,2) + Power(t1,2) - 
                       2*s2*(2 + t1 - t2) - 2*t1*(-1 + t2) - 2*t2 + 
                       Power(t2,2)) + 
                    Power(s1,4)*(s2 - t1)*
                     (-2 + 2*Power(s2,2) - t1 - Power(t1,2) + 2*t2 + 
                       3*t1*t2 - s2*(-1 + t1 + 3*t2)) + 
                    Power(s1,3)*
                     (2*Power(s2,4) - 2*Power(t1,3)*(-2 + t2) + 
                       Power(-1 + t2,2) - 
                       2*Power(s2,3)*(2 + t1 + 2*t2) + 
                       Power(t1,2)*(-3 - 4*t2 + 3*Power(t2,2)) + 
                       t1*(2 - 8*t2 + 6*Power(t2,2)) - 
                       2*Power(s2,2)*
                        (2 + Power(t1,2) + 2*t2 - 2*Power(t2,2) - 
                        t1*(5 + 4*t2)) + 
                       2*s2*(-1 + Power(t1,3) + 4*t2 - 
                        3*Power(t2,2) - Power(t1,2)*(5 + t2) + 
                        t1*(4 + 3*t2 - 3*Power(t2,2)))) + 
                    Power(s1,2)*
                     (2*Power(t1,4) + 
                       2*Power(s2,4)*(-1 + t1 - 3*t2) - 
                       Power(-1 + t2,2)*(-1 + 3*t2) + 
                       Power(t1,3)*(-5 - 6*t2 + Power(t2,2)) + 
                       t1*(-5 + 11*Power(t2,2) - 6*Power(t2,3)) - 
                       Power(t1,2)*
                        (-7 + t2 - 5*Power(t2,2) + Power(t2,3)) + 
                       2*Power(s2,2)*
                        (4 - 5*Power(t1,2) + Power(t1,3) + t2 + 
                        2*Power(t2,2) - 2*Power(t2,3) - 
                        t1*t2*(7 + 2*t2)) + 
                       Power(s2,3)*
                        (-4 - 4*Power(t1,2) + 8*t2 + 2*t1*(4 + 5*t2)) \
- 2*s2*(-1 - 3*t2 + 7*Power(t2,2) - 3*Power(t2,3) + 
                        Power(t1,3)*(-1 + 2*t2) - 
                        Power(t1,2)*(5 + 5*t2 + 2*Power(t2,2)) + 
                        t1*(7 + 3*t2 + Power(t2,2) - Power(t2,3)))) + 
                    s1*(-2*Power(t1,4)*(-2 + t2) + 
                       Power(-1 + t2,2)*t2*(-2 + 3*t2) + 
                       2*Power(s2,4)*t2*(2 - 2*t1 + 3*t2) + 
                       2*Power(t1,3)*(-2 + Power(t2,2)) + 
                       Power(t1,2)*
                        (-4 + 7*Power(t2,2) - 2*Power(t2,3)) + 
                       2*t1*
                        (2 + 2*t2 - 2*Power(t2,2) - 3*Power(t2,3) + 
                        Power(t2,4)) + 
                       2*s2*
                        (-2 + 2*Power(t1,4) + t2 - 4*Power(t2,2) + 
                        6*Power(t2,3) - Power(t2,4) + 
                        Power(t1,3)*(-7 - 5*t2 + Power(t2,2)) + 
                        t1*(1 + 5*t2 - 4*Power(t2,2) - 
                       3*Power(t2,3)) - 
                        Power(t1,2)*
                        (-6 + t2 - 4*Power(t2,2) + Power(t2,3))) + 
                       2*Power(s2,3)*
                        (Power(t1,2)*(2 + 4*t2) - 
                        t1*(4 + 8*t2 + 7*Power(t2,2)) + 
                        2*(1 + 2*t2 - Power(t2,2) + Power(t2,3))) + 
                       Power(s2,2)*
                        (2 - 12*t2 + 8*Power(t2,2) + 3*Power(t2,4) - 
                        4*Power(t1,3)*(2 + t2) + 
                        6*Power(t1,2)*(3 + 4*t2 + Power(t2,2)) - 
                        2*t1*(6 + 4*t2 + Power(t2,2) + 2*Power(t2,3)))\
) - 2*Power(s,3)*(-3 + 3*t1 + Power(t1,3) + 
                       2*Power(s2,3)*(-1 + t1 - t2) - 
                       Power(s2,2)*(2 + 3*t1)*(-1 + t1 - t2) + 2*t2 - 
                       2*t1*t2 - 5*Power(t1,2)*t2 - Power(t1,3)*t2 + 
                       2*Power(t2,2) + 7*t1*Power(t2,2) + 
                       3*Power(t1,2)*Power(t2,2) - 3*Power(t2,3) - 
                       3*t1*Power(t2,3) + Power(t2,4) + 
                       Power(s1,2)*
                        (2 + 3*Power(s2,2) + 2*Power(t1,2) + 
                        t1*(4 - 3*t2) - 2*t2 + Power(t2,2) + 
                        s2*(-5 - 4*t1 + 2*t2)) + 
                       s2*(3 + Power(t1,3) - Power(t1,2)*(-1 + t2) - 
                        2*t2 - 5*Power(t2,2) + Power(t2,3) - 
                        t1*(6 - 4*t2 + Power(t2,2))) + 
                       s1*(1 + 2*Power(s2,3) + 2*Power(t1,3) + 
                        Power(t1,2)*(3 - 6*t2) - 5*t2 + 
                        5*Power(t2,2) - 2*Power(t2,3) - 
                        Power(s2,2)*(5 + 3*t2) + 
                        t1*(-2 - 8*t2 + 6*Power(t2,2)) + 
                        s2*(2 - 3*Power(t1,2) + 8*t2 - 
                        3*Power(t2,2) + t1*(2 + 6*t2)))) + 
                    Power(s,2)*
                     (-5 + 13*t1 - 15*Power(t1,2) + 5*Power(t1,3) + 
                       2*Power(t1,4) + 2*Power(s2,4)*(-1 + t1 - t2) + 
                       7*t2 - 6*t1*t2 - 9*Power(t1,2)*t2 - 
                       8*Power(t1,3)*t2 + 3*Power(t2,2) + 
                       15*t1*Power(t2,2) + 
                       17*Power(t1,2)*Power(t2,2) + 
                       Power(t1,3)*Power(t2,2) - 11*Power(t2,3) - 
                       18*t1*Power(t2,3) - 
                       3*Power(t1,2)*Power(t2,3) + 7*Power(t2,4) + 
                       3*t1*Power(t2,4) - Power(t2,5) + 
                       Power(s1,3)*
                        (2 + 7*Power(s2,2) + 6*Power(t1,2) + 
                        t1*(10 - 6*t2) - 2*s2*(4 + 6*t1 - t2) - 
                        2*t2 + Power(t2,2)) + 
                       Power(s2,3)*
                        (-4*Power(t1,2) - 2*t1*(-2 + t2) + 
                        6*t2*(1 + t2)) - 
                       2*s2*
                        (2 + 2*t2 + 4*Power(t2,2) - 4*Power(t2,3) - 
                        Power(t2,4) + 2*Power(t1,3)*(2 + t2) + 
                        t1*t2*(-10 + t2 + 4*Power(t2,2)) - 
                        Power(t1,2)*(6 + t2 + 5*Power(t2,2))) + 
                       Power(s2,2)*
                        (10 + 2*Power(t1,3) - 8*t2 - 9*Power(t2,2) + 
                        5*Power(t2,3) + Power(t1,2)*(4 + 8*t2) - 
                        t1*(16 + 2*t2 + 15*Power(t2,2))) + 
                       Power(s1,2)*
                        (6 + 10*Power(s2,3) + 6*Power(t1,3) - 
                        14*t2 + 11*Power(t2,2) - 3*Power(t2,3) - 
                        Power(s2,2)*(15 + 9*t1 + 7*t2) - 
                        2*Power(t1,2)*(-5 + 9*t2) + 
                        t1*(6 - 34*t2 + 15*Power(t2,2)) - 
                        2*s2*
                        (2 + 3*Power(t1,2) - 10*t2 + Power(t2,2) - 
                        t1*(5 + 9*t2))) + 
                       s1*(-3 + 2*Power(s2,4) + 
                        2*Power(s2,3)*(-5 + 3*t1 - 8*t2) + 
                        Power(t1,3)*(8 - 6*t2) - 10*t2 + 
                        23*Power(t2,2) - 16*Power(t2,3) + 
                        3*Power(t2,4) + 
                        Power(t1,2)*(21 - 34*t2 + 15*Power(t2,2)) - 
                        2*t1*
                        (5 + 7*t2 - 21*Power(t2,2) + 6*Power(t2,3)) + 
                        Power(s2,2)*
                        (-2 - 14*Power(t1,2) + 22*t2 - 
                        5*Power(t2,2) + 2*t1*(9 + 13*t2)) + 
                        2*s2*
                        (8 + 3*Power(t1,3) + 3*t2 - 10*Power(t2,2) - 
                        Power(t2,3) - Power(t1,2)*(7 + 3*t2) + 
                        t1*(-11 + Power(t2,2))))) - 
                    2*s*(4*t1 - 9*Power(t1,2) + 6*Power(t1,3) - 
                       Power(t1,4) + 
                       Power(s1,4)*
                        (2*Power(s2,2) - s2*(1 + 4*t1) + 
                        t1*(2 + 2*t1 - t2)) + t2 - 2*t1*t2 + 
                       4*Power(t1,2)*t2 - 2*Power(t1,3)*t2 - 
                       Power(t1,4)*t2 + t1*Power(t2,2) + 
                       5*Power(t1,2)*Power(t2,2) + 
                       3*Power(t1,3)*Power(t2,2) - 4*Power(t2,3) - 
                       6*t1*Power(t2,3) - 4*Power(t1,2)*Power(t2,3) + 
                       4*Power(t2,4) + 3*t1*Power(t2,4) - 
                       Power(t2,5) + 2*Power(s2,4)*t2*(1 - t1 + t2) + 
                       Power(s2,3)*
                        (2 + 2*t2 - Power(t2,2) + 
                        Power(t1,2)*(2 + 4*t2) - 
                        2*t1*(2 + 3*t2 + 2*Power(t2,2))) + 
                       Power(s1,3)*
                        (4*Power(s2,3) + 2*Power(t1,3) + 
                        Power(t1,2)*(3 - 6*t2) + Power(-1 + t2,2) - 
                        Power(s2,2)*(3 + 5*t1 + 3*t2) + 
                        s2*(-3 + t1 - Power(t1,2) + 4*t2 + 8*t1*t2 + 
                        Power(t2,2)) + t1*(6 - 11*t2 + 3*Power(t2,2))) \
+ Power(s2,2)*(-2*Power(t1,3)*(2 + t2) + 
                        Power(t1,2)*(8 + 5*t2 + Power(t2,2)) - 
                        t2*(7 - 5*t2 - 4*Power(t2,2) + Power(t2,3)) + 
                        t1*(-4 + 4*t2 - 3*Power(t2,2) + 2*Power(t2,3))\
) + s2*(-4 + 2*Power(t1,4) + 2*t2 - 3*Power(t2,2) + 6*Power(t2,3) - 
                        Power(t2,5) + 
                        Power(t1,3)*(-3 + Power(t2,2)) - 
                        Power(t1,2)*
                        (4 + 4*t2 - 2*Power(t2,2) + 3*Power(t2,3)) + 
                        t1*(9 + 2*t2 - 7*Power(t2,2) - 
                        4*Power(t2,3) + 3*Power(t2,4))) + 
                       Power(s1,2)*
                        (1 + 2*Power(s2,4) + Power(t1,3)*(5 - 3*t2) - 
                        6*t2 + 8*Power(t2,2) - 3*Power(t2,3) - 
                        Power(s2,3)*(5 + 8*t2) + 
                        2*Power(t1,2)*(4 - 7*t2 + 3*Power(t2,2)) - 
                        t1*(5 + 15*t2 - 19*Power(t2,2) + 
                        3*Power(t2,3)) - 
                        Power(s2,2)*
                        (5 + 5*Power(t1,2) - 8*t2 + Power(t2,2) - 
                        t1*(13 + 14*t2)) + 
                        s2*(4 + 3*Power(t1,3) + 8*t2 - 
                        5*Power(t2,2) - 3*Power(t2,3) - 
                        Power(t1,2)*(13 + 3*t2) - t1*(1 + Power(t2,2))\
)) + s1*(-1 + 2*Power(t1,4) + 2*Power(s2,4)*(-1 + t1 - 2*t2) - t2 + 
                        9*Power(t2,2) - 10*Power(t2,3) + 
                        3*Power(t2,4) + 
                        Power(t1,3)*(-3 - 7*t2 + Power(t2,2)) + 
                        Power(t1,2)*
                        (3 - 15*t2 + 15*Power(t2,2) - 2*Power(t2,3)) + 
                        t1*(-1 + 5*t2 + 15*Power(t2,2) - 
                        13*Power(t2,3) + Power(t2,4)) + 
                        Power(s2,3)*
                        (-2 - 4*Power(t1,2) + 6*t2 + 4*Power(t2,2) + 
                        t1*(6 + 4*t2)) + 
                        s2*(1 - 2*t2 - 11*Power(t2,2) + 
                        2*Power(t2,3) + 3*Power(t2,4) - 
                        Power(t1,3)*(3 + 4*t2) + 
                        Power(t1,2)*(13 + 10*t2 + 7*Power(t2,2)) + 
                        t1*(-11 + 10*t2 + 3*Power(t2,2) - 
                        6*Power(t2,3))) + 
                        Power(s2,2)*
                        (2*Power(t1,3) + Power(t1,2)*(-3 + 4*t2) - 
                        t1*(8 + 10*t2 + 11*Power(t2,2)) + 
                        3*(3 - 3*Power(t2,2) + Power(t2,3)))))))/
                ((-1 + s1 + t1 - t2)*
                  Power(-1 + s - s*s2 + s1*s2 + t1 - s2*t2,2))))/
           ((-1 + t1)*(-1 + t2))) + 
       (((-64*(6 + Power(s,3) + Power(s1,2) - 8*s2 - 
                  2*Power(s2,2) + 10*t1 + 2*s2*t1 - 5*t2 + 
                  5*s2*t2 + Power(s2,2)*t2 - 6*t1*t2 - 2*s2*t1*t2 + 
                  Power(t1,2)*t2 + Power(t2,2) + 
                  Power(s,2)*(-5 - s1 - 2*s2 + 2*t1 + t2) - 
                  s1*(-5 + Power(s2,2) + s2*(5 - 2*t1) - 6*t1 + 
                     Power(t1,2) + 2*t2) + 
                  s*(5 + Power(s2,2) + 2*s1*(2 + s2 - t1) - 6*t1 + 
                     Power(t1,2) + s2*(7 - 2*t1 - 2*t2) - 4*t2 + 
                     2*t1*t2))*
                (-2*s2 + Power(s1,2)*(s2 - t1) + 2*t1 + 2*s2*t1 - 
                  2*Power(t1,2) + 
                  s1*(2*Power(s2,2) - 2*s2*t1 + 
                     (1 + t1)*(-1 + t2)) + t2 - 2*Power(s2,2)*t2 + 
                  t1*t2 + 2*s2*t1*t2 - Power(t2,2) - 
                  s2*Power(t2,2) + 
                  Power(s,2)*(-2 + 2*s2 - t1 + t2) + 
                  s*(3 - 2*Power(s2,2) - t1 + 2*s2*t1 - 3*t2 + 
                     s2*t2 - t1*t2 + Power(t2,2) - 
                     s1*(-2 + 3*s2 - 2*t1 + t2))))/
              (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2) + 
             128*(-20 + Power(s,4) - 10*s2 + 2*Power(s2,3) - 10*t1 - 
                8*s2*t1 - 2*Power(s2,2)*t1 + 8*Power(t1,2) - 
                2*s2*Power(t1,2) + 2*Power(t1,3) + 
                Power(s1,2)*(-7 + Power(s2,2) - t1 - s2*t1) + 29*t2 + 
                21*s2*t2 - 3*Power(s2,2)*t2 - Power(s2,3)*t2 - 
                9*t1*t2 + 9*s2*t1*t2 + 2*Power(s2,2)*t1*t2 - 
                5*Power(t1,2)*t2 - 2*s2*Power(t1,2)*t2 + 
                Power(t1,3)*t2 - 9*Power(t2,2) - 7*s2*Power(t2,2) + 
                6*t1*Power(t2,2) + s2*t1*Power(t2,2) - 
                Power(t1,2)*Power(t2,2) + 
                Power(s,3)*(-1 - 2*s1 - 3*s2 + 2*t1 + t2) + 
                Power(s,2)*(-11 + Power(s1,2) + 3*Power(s2,2) - t1 + 
                   2*Power(t1,2) + s2*(4 - 4*t1 - 3*t2) + 
                   s1*(3 + 5*s2 - 3*t1 - t2) + 2*t2 + t1*t2) + 
                s1*(-25 + Power(s2,3) - Power(t1,3) + t1*(7 - 5*t2) + 
                   16*t2 + Power(t1,2)*(3 + t2) - 
                   Power(s2,2)*(-5 + 2*t1 + t2) + 
                   s2*(-15 - 9*t1 + 2*Power(t1,2) + 7*t2)) + 
                s*(19 - Power(s2,3) + 7*t1 + Power(t1,2) + 
                   Power(t1,3) + Power(s1,2)*(-2 - 2*s2 + t1) - 
                   28*t2 - t1*t2 + Power(t1,2)*t2 + 5*Power(t2,2) - 
                   t1*Power(t2,2) + Power(s2,2)*(-5 + 2*t1 + 3*t2) - 
                   s2*(-13 - 5*t1 + 2*Power(t1,2) + t2 + 3*t1*t2) + 
                   s1*(24 - 4*Power(s2,2) + 4*t1 - 2*Power(t1,2) - 
                      3*t2 + s2*(-6 + 5*t1 + 2*t2)))))/
           ((s - s2 + t1)*(s - s1 + t2)) - 
          (64*(42 + 14*s2 + 12*Power(s2,2) + 
               Power(s1,4)*(1 + 2*s2)*(s2 - t1) - 20*t1 - 34*s2*t1 - 
               14*Power(s2,2)*t1 - 20*Power(t1,2) + 
               22*s2*Power(t1,2) + 2*Power(s2,2)*Power(t1,2) - 
               2*Power(t1,3) - 2*s2*Power(t1,3) - 60*t2 + 16*s2*t2 + 
               2*Power(s2,2)*t2 + 12*Power(s2,3)*t2 + 86*t1*t2 + 
               55*s2*t1*t2 - 13*Power(s2,2)*t1*t2 - 
               2*Power(s2,3)*t1*t2 - Power(t1,2)*t2 - 
               14*s2*Power(t1,2)*t2 + Power(s2,2)*Power(t1,2)*t2 + 
               5*Power(t1,3)*t2 + s2*Power(t1,3)*t2 + 17*Power(t2,2) - 
               55*s2*Power(t2,2) - 10*Power(s2,2)*Power(t2,2) - 
               4*Power(s2,3)*Power(t2,2) - 53*t1*Power(t2,2) - 
               9*s2*t1*Power(t2,2) + 12*Power(s2,2)*t1*Power(t2,2) + 
               Power(s2,3)*t1*Power(t2,2) + 
               2*Power(t1,2)*Power(t2,2) - 
               3*s2*Power(t1,2)*Power(t2,2) - 
               Power(s2,2)*Power(t1,2)*Power(t2,2) - 
               Power(t1,3)*Power(t2,2) + 2*Power(t2,3) + 
               28*s2*Power(t2,3) + Power(s2,2)*Power(t2,3) + 
               11*t1*Power(t2,3) + s2*t1*Power(t2,3) - 
               Power(s2,2)*t1*Power(t2,3) - Power(t1,2)*Power(t2,3) + 
               s2*Power(t1,2)*Power(t2,3) - Power(t2,4) - 
               5*s2*Power(t2,4) + Power(s,5)*(-t1 + t2) + 
               Power(s,4)*(9 + (5 + 3*s1)*t1 - 2*Power(t1,2) + 
                  s2*(-10 + s1 + 4*t1 - 3*t2) - 7*t2 - 2*s1*t2 + 
                  2*Power(t2,2)) + 
               Power(s1,3)*(-1 - 8*Power(t1,2) + Power(t1,3) + 
                  Power(s2,2)*(-1 + 3*t1 - 6*t2) + t2 + 
                  t1*(-4 + 3*t2) + 
                  s2*(-15 - 4*Power(t1,2) + 2*t2 + t1*(11 + 6*t2))) + 
               Power(s,3)*(-46 + 9*t1 + 6*Power(t1,2) - Power(t1,3) + 
                  35*t2 - 2*t1*t2 - 3*Power(t1,2)*t2 - 
                  12*Power(t2,2) + 3*t1*Power(t2,2) + Power(t2,3) + 
                  Power(s1,2)*(2 - 4*s2 - 3*t1 + t2) + 
                  Power(s2,2)*(20 - 5*t1 + 3*t2) + 
                  s2*(24 - 26*t1 + 5*Power(t1,2) - 10*t2 + 4*t1*t2 - 
                     5*Power(t2,2)) + 
                  s1*(-2*Power(s2,2) + s2*(22 - 11*t1 + 9*t2) + 
                     2*(-8 - 6*t1 + 3*Power(t1,2) + 5*t2 - 
                        Power(t2,2)))) + 
               Power(s1,2)*(13 + Power(s2,3)*(-8 + t1) + 4*t2 - 
                  3*Power(t2,2) - 2*Power(t1,3)*(1 + t2) + 
                  Power(t1,2)*(-13 + 15*t2) + 
                  t1*(-33 + 19*t2 - 3*Power(t2,2)) + 
                  Power(s2,2)*
                   (-20 - 2*Power(t1,2) + t1*(19 - 7*t2) + 3*t2 + 
                     6*Power(t2,2)) + 
                  s2*(-53 + Power(t1,3) + 58*t2 - 12*Power(t2,2) + 
                     Power(t1,2)*(-7 + 9*t2) - 
                     3*t1*(-8 + 7*t2 + 2*Power(t2,2)))) + 
               Power(s,2)*(77 - 63*t1 - Power(t1,2) + 
                  Power(s1,3)*(-2 + 5*s2 + t1) + 
                  Power(s1,2)*
                   (3 + 5*Power(s2,2) + 14*t1 - 6*Power(t1,2) + 
                     4*s2*(-4 + 2*t1 - 3*t2)) + 
                  Power(s2,3)*(-10 + 2*t1 - t2) - 83*t2 + 50*t1*t2 + 
                  8*Power(t1,2)*t2 - 2*Power(t1,3)*t2 + 
                  29*Power(t2,2) - 16*t1*Power(t2,2) - 
                  4*Power(t2,3) + 2*t1*Power(t2,3) + 
                  Power(s2,2)*
                   (-49 - 3*Power(t1,2) + t1*(29 - 7*t2) + 25*t2 + 
                     4*Power(t2,2)) + 
                  s2*(20 + Power(t1,3) + 38*t2 - 6*Power(t2,2) - 
                     2*Power(t2,3) + 4*Power(t1,2)*(-3 + 2*t2) - 
                     t1*(-13 + 30*t2 + Power(t2,2))) + 
                  s1*(64 + Power(s2,3) + 3*Power(t1,3) + 
                     Power(s2,2)*(-34 + 11*t1 - 9*t2) - 32*t2 + 
                     6*Power(t2,2) + 2*Power(t1,2)*(-7 + 3*t2) + 
                     t1*(-19 + 2*t2 - 3*Power(t2,2)) + 
                     s2*(-64 - 14*Power(t1,2) - 7*t1*(-8 + t2) + 
                        22*t2 + 9*Power(t2,2)))) + 
               s1*(54 - 2*Power(s2,3)*(-6 + t1)*(-1 + t2) - 30*t2 - 
                  5*Power(t2,2) + 3*Power(t2,3) + 
                  Power(t1,2)*(1 + 11*t2 - 6*Power(t2,2)) + 
                  Power(t1,3)*(-7 + 3*t2 + Power(t2,2)) + 
                  t1*(-78 + 86*t2 - 26*Power(t2,2) + Power(t2,3)) + 
                  Power(s2,2)*
                   (-6 + 3*Power(t1,2)*(-1 + t2) + 30*t2 - 
                     3*Power(t2,2) - 2*Power(t2,3) + 
                     t1*(19 - 31*t2 + 5*Power(t2,2))) + 
                  s2*(-22 - Power(t1,3)*(-1 + t2) + 108*t2 - 
                     71*Power(t2,2) + 14*Power(t2,3) + 
                     Power(t1,2)*(8 + 10*t2 - 6*Power(t2,2)) + 
                     t1*(-45 - 15*t2 + 9*Power(t2,2) + 2*Power(t2,3)))\
) + s*(-76 - 2*Power(s1,4)*s2 + 60*t1 + 3*Power(t1,2) + 
                  7*Power(t1,3) + 100*t2 - 130*t1*t2 - 
                  5*Power(t1,2)*t2 - Power(t1,3)*t2 - 27*Power(t2,2) + 
                  58*t1*Power(t2,2) + Power(t1,2)*Power(t2,2) - 
                  Power(t1,3)*Power(t2,2) - 9*t1*Power(t2,3) + 
                  Power(t1,2)*Power(t2,3) + Power(t2,4) - 
                  Power(s2,3)*
                   (-20 + t1*(4 - 3*t2) + 12*t2 + Power(t2,2)) + 
                  Power(s1,3)*
                   (4 - 5*Power(s2,2) - 6*t1 + 2*Power(t1,2) - t2 + 
                     s2*(3 + t1 + 6*t2)) + 
                  Power(s2,2)*
                   (2 + Power(t1,2)*(5 - 4*t2) - 45*t2 + 
                     10*Power(t2,2) + Power(t2,3) + 
                     t1*(-25 + 33*t2 - 3*Power(t2,2))) + 
                  s2*(-32 + Power(t1,3)*(-1 + t2) - 29*t2 + 
                     46*Power(t2,2) - 9*Power(t2,3) + 
                     Power(t1,2)*(-14 - 9*t2 + 4*Power(t2,2)) - 
                     t1*(-71 + 2*t2 + 9*Power(t2,2) + Power(t2,3))) - 
                  Power(s1,2)*
                   (9 + Power(s2,3) + 3*Power(t1,3) + 
                     Power(s2,2)*(-17 + 9*t1 - 11*t2) + 8*t2 - 
                     3*Power(t2,2) + Power(t1,2)*(-16 + 3*t2) - 
                     t1*(14 + 3*t2) + 
                     s2*(-53 - 13*Power(t1,2) + 15*t2 + 
                        6*Power(t2,2) + t1*(41 + 3*t2))) + 
                  s1*(-88 + Power(t1,2)*(14 - 17*t2) + 36*t2 + 
                     4*Power(t2,2) - 3*Power(t2,3) + 
                     Power(t1,3)*(2 + 4*t2) + 
                     2*t1*(55 - 36*t2 + 6*Power(t2,2)) + 
                     Power(s2,3)*(-3*t1 + 2*(8 + t2)) + 
                     Power(s2,2)*
                      (61 + 5*Power(t1,2) - 27*t2 - 7*Power(t2,2) + 
                        4*t1*(-11 + 3*t2)) + 
                     s2*(37 - 2*Power(t1,3) + 
                        Power(t1,2)*(19 - 17*t2) - 99*t2 + 
                        21*Power(t2,2) + 2*Power(t2,3) + 
                        t1*(-43 + 50*t2 + 3*Power(t2,2)))))))/
           ((-1 + s2)*(-s + s2 - t1)*
             (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2)) + 
          (32*((-2*(-2*s2 + Power(s1,2)*(s2 - t1) + 2*t1 + 
                    2*s2*t1 - 2*Power(t1,2) + 
                    s1*(2*Power(s2,2) - 2*s2*t1 + 
                       (1 + t1)*(-1 + t2)) + t2 - 
                    2*Power(s2,2)*t2 + t1*t2 + 2*s2*t1*t2 - 
                    Power(t2,2) - s2*Power(t2,2) + 
                    Power(s,2)*(-2 + 2*s2 - t1 + t2) + 
                    s*(3 - 2*Power(s2,2) - t1 + 2*s2*t1 - 3*t2 + 
                       s2*t2 - t1*t2 + Power(t2,2) - 
                       s1*(-2 + 3*s2 - 2*t1 + t2)))*
                  (23 + 2*Power(s,3) + 5*s2 + 2*Power(s2,2) - 
                    2*Power(s2,3) + Power(s1,2)*(-7 + t1) - 8*t1 - 
                    8*s2*t1 + 4*Power(s2,2)*t1 + Power(t1,2) - 
                    s2*Power(t1,2) - 7*t2 + 8*s2*t2 + 
                    Power(s2,2)*t2 + 4*t1*t2 - 4*s2*t1*t2 + 
                    Power(t1,2)*t2 - 4*Power(t2,2) + 
                    s2*Power(t2,2) + 
                    Power(s,2)*(-3*s1 - 6*s2 + 3*t1 + 2*t2) - 
                    s1*(-7 + 3*Power(s2,2) + Power(t1,2) - 11*t2 + 
                       t1*(8 + t2) + s2*(4 - 6*t1 + t2)) + 
                    s*(-12 + Power(s1,2) + 6*Power(s2,2) + 9*t1 + 
                       Power(t1,2) + s1*(7 + 6*s2 - 4*t1 - t2) - 
                       8*t2 + 3*t1*t2 - s2*(2 + 7*t1 + 3*t2))))/
                (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2) + 
               4*(-12 + Power(s,4) - 26*s2 - 6*Power(s2,2) - 
                  Power(s2,3) + Power(s2,4) + 13*t1 + 11*s2*t1 + 
                  7*Power(s2,2)*t1 - 3*Power(s2,3)*t1 - 
                  5*Power(t1,2) - 5*s2*Power(t1,2) + 
                  2*Power(s2,2)*Power(t1,2) + Power(s1,3)*(3 + t1) + 
                  17*t2 + 17*s2*t2 - 6*Power(s2,2)*t2 - 12*t1*t2 + 
                  2*Power(s2,2)*t1*t2 + 4*Power(t1,2)*t2 - 
                  3*s2*Power(t1,2)*t2 + s2*Power(t2,2) - 
                  Power(s2,2)*Power(t2,2) - 2*t1*Power(t2,2) + 
                  2*s2*t1*Power(t2,2) + Power(t1,2)*Power(t2,2) + 
                  Power(t2,3) - t1*Power(t2,3) + 
                  Power(s,3)*(1 - 2*s1 - 4*s2 + 2*t1 + t2) - 
                  Power(s1,2)*
                   (7 + s2*(-6 + t1) - 2*Power(t1,2) + 
                     3*t1*(-1 + t2) + 4*t2) + 
                  Power(s,2)*
                   (-14 + Power(s1,2) + 6*Power(s2,2) + 8*t1 + 
                     2*Power(t1,2) + s1*(1 + 6*s2 - 3*t1 - t2) - 
                     3*t2 + t1*t2 - s2*(3 + 7*t1 + 2*t2)) + 
                  s*(32 - 4*Power(s2,3) - Power(s1,2)*(5 + 2*s2) - 
                     14*t1 + 5*Power(t1,2) - 17*t2 + 3*t1*t2 + 
                     3*Power(t1,2)*t2 - Power(t2,2) - 
                     2*t1*Power(t2,2) + 
                     Power(s2,2)*(3 + 8*t1 + t2) + 
                     s1*(19 - 6*Power(s2,2) - 11*t1 - 
                       4*Power(t1,2) + s2*(-2 + 8*t1) + 8*t2 + 
                       2*t1*t2) + 
                     s2*(20 - 4*Power(t1,2) + 9*t2 + Power(t2,2) - 
                        3*t1*(5 + t2))) + 
                  s1*(-25 + 2*Power(s2,3) + 7*t2 + 
                     Power(s2,2)*(-6*t1 + 2*t2) - 
                     Power(t1,2)*(5 + 3*t2) + 
                     t1*(7 + 2*t2 + 3*Power(t2,2)) + 
                     s2*(4*Power(t1,2) + t1*(11 - 2*t2) - 
                        5*(3 + 2*t2)))) + 
               (2*(Power(s,2) - s2 + Power(s2,2) + 
                    s1*(2 + s2 - t1) + t1 - s2*t1 - 2*t2 - s2*t2 + 
                    t1*t2 + s*(-1 - s1 - 2*s2 + t1 + t2))*
                  (-2*Power(s2,2) + Power(s1,5)*Power(s2 - t1,2) + 
                    4*s2*t1 + 6*Power(s2,2)*t1 - 2*Power(t1,2) - 
                    12*s2*Power(t1,2) - 6*Power(s2,2)*Power(t1,2) + 
                    6*Power(t1,3) + 12*s2*Power(t1,3) + 
                    2*Power(s2,2)*Power(t1,3) - 6*Power(t1,4) - 
                    4*s2*Power(t1,4) + 2*Power(t1,5) + 4*s2*t2 - 
                    2*Power(s2,2)*t2 - 4*Power(s2,3)*t2 - 4*t1*t2 - 
                    2*s2*t1*t2 + 12*Power(s2,2)*t1*t2 + 
                    8*Power(s2,3)*t1*t2 + 4*Power(t1,2)*t2 - 
                    12*s2*Power(t1,2)*t2 - 
                    18*Power(s2,2)*Power(t1,2)*t2 - 
                    4*Power(s2,3)*Power(t1,2)*t2 + 4*Power(t1,3)*t2 + 
                    14*s2*Power(t1,3)*t2 + 
                    8*Power(s2,2)*Power(t1,3)*t2 - 4*Power(t1,4)*t2 - 
                    4*s2*Power(t1,4)*t2 + Power(t2,2) - 
                    4*s2*Power(t2,2) + 4*Power(s2,2)*Power(t2,2) - 
                    4*Power(s2,3)*Power(t2,2) - 
                    2*Power(s2,4)*Power(t2,2) + t1*Power(t2,2) + 
                    4*s2*t1*Power(t2,2) + 
                    8*Power(s2,2)*t1*Power(t2,2) + 
                    8*Power(s2,3)*t1*Power(t2,2) + 
                    2*Power(s2,4)*t1*Power(t2,2) - 
                    7*Power(t1,2)*Power(t2,2) - 
                    8*s2*Power(t1,2)*Power(t2,2) - 
                    14*Power(s2,2)*Power(t1,2)*Power(t2,2) - 
                    4*Power(s2,3)*Power(t1,2)*Power(t2,2) + 
                    5*Power(t1,3)*Power(t2,2) + 
                    8*s2*Power(t1,3)*Power(t2,2) + 
                    2*Power(s2,2)*Power(t1,3)*Power(t2,2) - 
                    3*Power(t2,3) + 4*s2*Power(t2,3) - 
                    6*Power(s2,2)*Power(t2,3) - 
                    2*Power(s2,4)*Power(t2,3) + 2*t1*Power(t2,3) + 
                    6*s2*t1*Power(t2,3) + 
                    6*Power(s2,2)*t1*Power(t2,3) + 
                    6*Power(s2,3)*t1*Power(t2,3) - 
                    3*Power(t1,2)*Power(t2,3) - 
                    8*s2*Power(t1,2)*Power(t2,3) - 
                    4*Power(s2,2)*Power(t1,2)*Power(t2,3) + 
                    3*Power(t2,4) - 4*s2*Power(t2,4) - 
                    Power(s2,2)*Power(t2,4) - 
                    2*Power(s2,3)*Power(t2,4) + t1*Power(t2,4) + 
                    4*s2*t1*Power(t2,4) + 
                    3*Power(s2,2)*t1*Power(t2,4) - Power(t2,5) - 
                    Power(s2,2)*Power(t2,5) + 
                    Power(s,4)*(-1 + s1 + t1 - t2)*
                     (2 + 2*Power(s2,2) + Power(t1,2) - 
                       2*s2*(2 + t1 - t2) - 2*t1*(-1 + t2) - 2*t2 + 
                       Power(t2,2)) + 
                    Power(s1,4)*(s2 - t1)*
                     (-2 + 2*Power(s2,2) - t1 - Power(t1,2) + 2*t2 + 
                       3*t1*t2 - s2*(-1 + t1 + 3*t2)) + 
                    Power(s1,3)*
                     (2*Power(s2,4) - 2*Power(t1,3)*(-2 + t2) + 
                       Power(-1 + t2,2) - 
                       2*Power(s2,3)*(2 + t1 + 2*t2) + 
                       Power(t1,2)*(-3 - 4*t2 + 3*Power(t2,2)) + 
                       t1*(2 - 8*t2 + 6*Power(t2,2)) - 
                       2*Power(s2,2)*
                       (2 + Power(t1,2) + 2*t2 - 2*Power(t2,2) - 
                       t1*(5 + 4*t2)) + 
                       2*s2*
                        (-1 + Power(t1,3) + 4*t2 - 3*Power(t2,2) - 
                        Power(t1,2)*(5 + t2) + 
                        t1*(4 + 3*t2 - 3*Power(t2,2)))) + 
                    Power(s1,2)*
                     (2*Power(t1,4) + 
                       2*Power(s2,4)*(-1 + t1 - 3*t2) - 
                       Power(-1 + t2,2)*(-1 + 3*t2) + 
                       Power(t1,3)*(-5 - 6*t2 + Power(t2,2)) + 
                       t1*(-5 + 11*Power(t2,2) - 6*Power(t2,3)) - 
                       Power(t1,2)*
                       (-7 + t2 - 5*Power(t2,2) + Power(t2,3)) + 
                       2*Power(s2,2)*
                       (4 - 5*Power(t1,2) + Power(t1,3) + t2 + 
                       2*Power(t2,2) - 2*Power(t2,3) - 
                       t1*t2*(7 + 2*t2)) + 
                       Power(s2,3)*
                       (-4 - 4*Power(t1,2) + 8*t2 + 2*t1*(4 + 5*t2)) \
- 2*s2*(-1 - 3*t2 + 7*Power(t2,2) - 3*Power(t2,3) + 
                        Power(t1,3)*(-1 + 2*t2) - 
                        Power(t1,2)*(5 + 5*t2 + 2*Power(t2,2)) + 
                        t1*(7 + 3*t2 + Power(t2,2) - Power(t2,3)))) + 
                    s1*(-2*Power(t1,4)*(-2 + t2) + 
                       Power(-1 + t2,2)*t2*(-2 + 3*t2) + 
                       2*Power(s2,4)*t2*(2 - 2*t1 + 3*t2) + 
                       2*Power(t1,3)*(-2 + Power(t2,2)) + 
                       Power(t1,2)*
                       (-4 + 7*Power(t2,2) - 2*Power(t2,3)) + 
                       2*t1*
                       (2 + 2*t2 - 2*Power(t2,2) - 3*Power(t2,3) + 
                       Power(t2,4)) + 
                       2*s2*
                       (-2 + 2*Power(t1,4) + t2 - 4*Power(t2,2) + 
                       6*Power(t2,3) - Power(t2,4) + 
                       Power(t1,3)*(-7 - 5*t2 + Power(t2,2)) + 
                       t1*(1 + 5*t2 - 4*Power(t2,2) - 
                       3*Power(t2,3)) - 
                       Power(t1,2)*
                       (-6 + t2 - 4*Power(t2,2) + Power(t2,3))) + 
                       2*Power(s2,3)*
                       (Power(t1,2)*(2 + 4*t2) - 
                       t1*(4 + 8*t2 + 7*Power(t2,2)) + 
                       2*(1 + 2*t2 - Power(t2,2) + Power(t2,3))) + 
                       Power(s2,2)*
                        (2 - 12*t2 + 8*Power(t2,2) + 3*Power(t2,4) - 
                        4*Power(t1,3)*(2 + t2) + 
                        6*Power(t1,2)*(3 + 4*t2 + Power(t2,2)) - 
                        2*t1*(6 + 4*t2 + Power(t2,2) + 2*Power(t2,3))\
)) - 2*Power(s,3)*(-3 + 3*t1 + Power(t1,3) + 
                       2*Power(s2,3)*(-1 + t1 - t2) - 
                       Power(s2,2)*(2 + 3*t1)*(-1 + t1 - t2) + 
                       2*t2 - 2*t1*t2 - 5*Power(t1,2)*t2 - 
                       Power(t1,3)*t2 + 2*Power(t2,2) + 
                       7*t1*Power(t2,2) + 
                       3*Power(t1,2)*Power(t2,2) - 3*Power(t2,3) - 
                       3*t1*Power(t2,3) + Power(t2,4) + 
                       Power(s1,2)*
                       (2 + 3*Power(s2,2) + 2*Power(t1,2) + 
                       t1*(4 - 3*t2) - 2*t2 + Power(t2,2) + 
                       s2*(-5 - 4*t1 + 2*t2)) + 
                       s2*(3 + Power(t1,3) - 
                       Power(t1,2)*(-1 + t2) - 2*t2 - 
                       5*Power(t2,2) + Power(t2,3) - 
                       t1*(6 - 4*t2 + Power(t2,2))) + 
                       s1*(1 + 2*Power(s2,3) + 2*Power(t1,3) + 
                        Power(t1,2)*(3 - 6*t2) - 5*t2 + 
                        5*Power(t2,2) - 2*Power(t2,3) - 
                        Power(s2,2)*(5 + 3*t2) + 
                        t1*(-2 - 8*t2 + 6*Power(t2,2)) + 
                        s2*(2 - 3*Power(t1,2) + 8*t2 - 
                        3*Power(t2,2) + t1*(2 + 6*t2)))) + 
                    Power(s,2)*
                     (-5 + 13*t1 - 15*Power(t1,2) + 5*Power(t1,3) + 
                       2*Power(t1,4) + 
                       2*Power(s2,4)*(-1 + t1 - t2) + 7*t2 - 
                       6*t1*t2 - 9*Power(t1,2)*t2 - 
                       8*Power(t1,3)*t2 + 3*Power(t2,2) + 
                       15*t1*Power(t2,2) + 
                       17*Power(t1,2)*Power(t2,2) + 
                       Power(t1,3)*Power(t2,2) - 11*Power(t2,3) - 
                       18*t1*Power(t2,3) - 
                       3*Power(t1,2)*Power(t2,3) + 7*Power(t2,4) + 
                       3*t1*Power(t2,4) - Power(t2,5) + 
                       Power(s1,3)*
                       (2 + 7*Power(s2,2) + 6*Power(t1,2) + 
                       t1*(10 - 6*t2) - 2*s2*(4 + 6*t1 - t2) - 
                       2*t2 + Power(t2,2)) + 
                       Power(s2,3)*
                       (-4*Power(t1,2) - 2*t1*(-2 + t2) + 
                       6*t2*(1 + t2)) - 
                       2*s2*
                       (2 + 2*t2 + 4*Power(t2,2) - 4*Power(t2,3) - 
                       Power(t2,4) + 2*Power(t1,3)*(2 + t2) + 
                       t1*t2*(-10 + t2 + 4*Power(t2,2)) - 
                       Power(t1,2)*(6 + t2 + 5*Power(t2,2))) + 
                       Power(s2,2)*
                       (10 + 2*Power(t1,3) - 8*t2 - 9*Power(t2,2) + 
                       5*Power(t2,3) + Power(t1,2)*(4 + 8*t2) - 
                       t1*(16 + 2*t2 + 15*Power(t2,2))) + 
                       Power(s1,2)*
                       (6 + 10*Power(s2,3) + 6*Power(t1,3) - 
                       14*t2 + 11*Power(t2,2) - 3*Power(t2,3) - 
                       Power(s2,2)*(15 + 9*t1 + 7*t2) - 
                       2*Power(t1,2)*(-5 + 9*t2) + 
                       t1*(6 - 34*t2 + 15*Power(t2,2)) - 
                       2*s2*
                       (2 + 3*Power(t1,2) - 10*t2 + Power(t2,2) - 
                       t1*(5 + 9*t2))) + 
                       s1*(-3 + 2*Power(s2,4) + 
                        2*Power(s2,3)*(-5 + 3*t1 - 8*t2) + 
                        Power(t1,3)*(8 - 6*t2) - 10*t2 + 
                        23*Power(t2,2) - 16*Power(t2,3) + 
                        3*Power(t2,4) + 
                        Power(t1,2)*(21 - 34*t2 + 15*Power(t2,2)) - 
                        2*t1*
                       (5 + 7*t2 - 21*Power(t2,2) + 6*Power(t2,3)) + 
                        Power(s2,2)*
                       (-2 - 14*Power(t1,2) + 22*t2 - 
                       5*Power(t2,2) + 2*t1*(9 + 13*t2)) + 
                        2*s2*
                        (8 + 3*Power(t1,3) + 3*t2 - 10*Power(t2,2) - 
                        Power(t2,3) - Power(t1,2)*(7 + 3*t2) + 
                        t1*(-11 + Power(t2,2))))) - 
                    2*s*(4*t1 - 9*Power(t1,2) + 6*Power(t1,3) - 
                       Power(t1,4) + 
                       Power(s1,4)*
                        (2*Power(s2,2) - s2*(1 + 4*t1) + 
                        t1*(2 + 2*t1 - t2)) + t2 - 2*t1*t2 + 
                       4*Power(t1,2)*t2 - 2*Power(t1,3)*t2 - 
                       Power(t1,4)*t2 + t1*Power(t2,2) + 
                       5*Power(t1,2)*Power(t2,2) + 
                       3*Power(t1,3)*Power(t2,2) - 4*Power(t2,3) - 
                       6*t1*Power(t2,3) - 4*Power(t1,2)*Power(t2,3) + 
                       4*Power(t2,4) + 3*t1*Power(t2,4) - 
                       Power(t2,5) + 2*Power(s2,4)*t2*(1 - t1 + t2) + 
                       Power(s2,3)*
                        (2 + 2*t2 - Power(t2,2) + 
                        Power(t1,2)*(2 + 4*t2) - 
                        2*t1*(2 + 3*t2 + 2*Power(t2,2))) + 
                       Power(s1,3)*
                        (4*Power(s2,3) + 2*Power(t1,3) + 
                        Power(t1,2)*(3 - 6*t2) + Power(-1 + t2,2) - 
                        Power(s2,2)*(3 + 5*t1 + 3*t2) + 
                        s2*(-3 + t1 - Power(t1,2) + 4*t2 + 
                       8*t1*t2 + Power(t2,2)) + 
                        t1*(6 - 11*t2 + 3*Power(t2,2))) + 
                       Power(s2,2)*
                        (-2*Power(t1,3)*(2 + t2) + 
                        Power(t1,2)*(8 + 5*t2 + Power(t2,2)) - 
                        t2*
                       (7 - 5*t2 - 4*Power(t2,2) + Power(t2,3)) + 
                        t1*(-4 + 4*t2 - 3*Power(t2,2) + 
                       2*Power(t2,3))) + 
                       s2*(-4 + 2*Power(t1,4) + 2*t2 - 
                        3*Power(t2,2) + 6*Power(t2,3) - 
                        Power(t2,5) + 
                        Power(t1,3)*(-3 + Power(t2,2)) - 
                        Power(t1,2)*
                       (4 + 4*t2 - 2*Power(t2,2) + 3*Power(t2,3)) + 
                        t1*(9 + 2*t2 - 7*Power(t2,2) - 
                        4*Power(t2,3) + 3*Power(t2,4))) + 
                       Power(s1,2)*
                        (1 + 2*Power(s2,4) + 
                        Power(t1,3)*(5 - 3*t2) - 6*t2 + 
                        8*Power(t2,2) - 3*Power(t2,3) - 
                        Power(s2,3)*(5 + 8*t2) + 
                        2*Power(t1,2)*(4 - 7*t2 + 3*Power(t2,2)) - 
                        t1*(5 + 15*t2 - 19*Power(t2,2) + 
                       3*Power(t2,3)) - 
                        Power(s2,2)*
                       (5 + 5*Power(t1,2) - 8*t2 + Power(t2,2) - 
                       t1*(13 + 14*t2)) + 
                        s2*(4 + 3*Power(t1,3) + 8*t2 - 
                        5*Power(t2,2) - 3*Power(t2,3) - 
                        Power(t1,2)*(13 + 3*t2) - 
                        t1*(1 + Power(t2,2)))) + 
                       s1*(-1 + 2*Power(t1,4) + 
                        2*Power(s2,4)*(-1 + t1 - 2*t2) - t2 + 
                        9*Power(t2,2) - 10*Power(t2,3) + 
                        3*Power(t2,4) + 
                        Power(t1,3)*(-3 - 7*t2 + Power(t2,2)) + 
                        Power(t1,2)*
                        (3 - 15*t2 + 15*Power(t2,2) - 2*Power(t2,3)) \
+ t1*(-1 + 5*t2 + 15*Power(t2,2) - 13*Power(t2,3) + Power(t2,4)) + 
                        Power(s2,3)*
                        (-2 - 4*Power(t1,2) + 6*t2 + 4*Power(t2,2) + 
                        t1*(6 + 4*t2)) + 
                        s2*(1 - 2*t2 - 11*Power(t2,2) + 
                        2*Power(t2,3) + 3*Power(t2,4) - 
                        Power(t1,3)*(3 + 4*t2) + 
                        Power(t1,2)*(13 + 10*t2 + 7*Power(t2,2)) + 
                        t1*(-11 + 10*t2 + 3*Power(t2,2) - 
                        6*Power(t2,3))) + 
                        Power(s2,2)*
                        (2*Power(t1,3) + Power(t1,2)*(-3 + 4*t2) - 
                        t1*(8 + 10*t2 + 11*Power(t2,2)) + 
                        3*(3 - 3*Power(t2,2) + Power(t2,3)))))))/
                ((-1 + s1 + t1 - t2)*
                  Power(-1 + s - s*s2 + s1*s2 + t1 - s2*t2,2))))/
           ((-1 + s2)*(-1 + t1)) + 
          (32*((2*(-2*s2 + Power(s1,2)*(s2 - t1) + 2*t1 + 2*s2*t1 - 
                    2*Power(t1,2) + 
                    s1*(2*Power(s2,2) - 2*s2*t1 + 
                       (1 + t1)*(-1 + t2)) + t2 - 
                    2*Power(s2,2)*t2 + t1*t2 + 2*s2*t1*t2 - 
                    Power(t2,2) - s2*Power(t2,2) + 
                    Power(s,2)*(-2 + 2*s2 - t1 + t2) + 
                    s*(3 - 2*Power(s2,2) - t1 + 2*s2*t1 - 3*t2 + 
                       s2*t2 - t1*t2 + Power(t2,2) - 
                       s1*(-2 + 3*s2 - 2*t1 + t2)))*
                  (-19 + 2*Power(s,3) + 2*Power(s1,3) - 5*s2 + 
                    4*Power(s2,2) - 2*Power(s2,3) - 5*t1 - 7*s2*t1 + 
                    4*Power(s2,2)*t1 + 5*Power(t1,2) - 
                    s2*Power(t1,2) + 
                    Power(s1,2)*(1 + 3*s2 + 3*t1 - 5*t2) + 19*t2 + 
                    2*s2*t2 + Power(s2,2)*t2 - 3*t1*t2 - 
                    8*s2*t1*t2 + 2*Power(t1,2)*t2 - 3*Power(t2,2) + 
                    5*s2*Power(t2,2) - Power(t2,3) - 
                    Power(s,2)*(-4 + s1 + 6*s2 - 6*t1 + t2) + 
                    s*(8 - 3*Power(s1,2) + 6*Power(s2,2) + 7*t1 + 
                       2*Power(t1,2) - 4*t2 + 6*t1*t2 - 
                       4*Power(t2,2) + s2*(-7 - 11*t1 + t2) + 
                       s1*(-4 + 3*s2 - 9*t1 + 7*t2)) - 
                    s1*(3*Power(s2,2) + t1 + 2*Power(t1,2) + 
                       3*t1*t2 + s2*(-3 - 10*t1 + 8*t2) - 
                       2*(-9 + t2 + 2*Power(t2,2)))))/
                (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2) - 
               4*(12 + Power(s,4) - Power(s1,4) + 15*s2 + 
                  Power(s2,2) - Power(s2,3) + Power(s2,4) + 
                  12*s2*t1 + 3*Power(s2,2)*t1 - 3*Power(s2,3)*t1 - 
                  11*Power(t1,2) - 4*s2*Power(t1,2) + 
                  2*Power(s2,2)*Power(t1,2) + Power(t1,3) - 
                  3*Power(s1,3)*(s2 - t2) - 27*t2 - 24*s2*t2 - 
                  4*Power(s2,2)*t2 + 12*t1*t2 + 9*s2*t1*t2 + 
                  5*Power(s2,2)*t1*t2 - 5*s2*Power(t1,2)*t2 + 
                  12*Power(t2,2) + 3*s2*Power(t2,2) - 
                  4*Power(s2,2)*Power(t2,2) - 5*t1*Power(t2,2) + 
                  3*s2*t1*Power(t2,2) + Power(t1,2)*Power(t2,2) - 
                  Power(t2,3) + s2*Power(t2,3) - t1*Power(t2,3) - 
                  Power(s,3)*(-1 + s1 + 4*s2 - 4*t1 + t2) + 
                  Power(s,2)*
                   (4 - 2*Power(s1,2) + 6*Power(s2,2) + 2*t1 + 
                     3*Power(t1,2) - 5*t2 + 3*t1*t2 - 
                     3*Power(t2,2) + s2*(-2 - 12*t1 + 3*t2) + 
                     s1*(-2 + 3*s2 - 8*t1 + 6*t2)) - 
                  Power(s1,2)*
                   (2*Power(s2,2) - 3*Power(t1,2) + 
                     4*s2*(t1 - 2*t2) + 2*t1*(1 + t2) + 
                     3*(-5 + t2 + Power(t2,2))) + 
                  s1*(27 + 2*Power(s2,3) - 27*t2 + 4*Power(t2,2) + 
                     Power(t2,3) - Power(t1,2)*(3 + 4*t2) + 
                     Power(s2,2)*(-1 - 9*t1 + 7*t2) + 
                     s2*(19 + t1 + 6*Power(t1,2) - 5*t2 - 
                       6*Power(t2,2)) + 
                     t1*(-4 + 9*t2 + 3*Power(t2,2))) + 
                  s*(-21 + 3*Power(s1,3) - 4*Power(s2,3) - 10*t1 + 
                     5*Power(t1,2) + 
                     Power(s1,2)*(1 + 4*s2 + 4*t1 - 8*t2) + 
                     Power(s2,2)*(2 + 11*t1 - 2*t2) + 27*t2 - 
                     8*t1*t2 + 4*Power(t1,2)*t2 - 6*Power(t2,2) - 
                     2*t1*Power(t2,2) - Power(t2,3) - 
                     s2*(6 + 4*t1 + 6*Power(t1,2) - 8*t2 + 7*t1*t2 - 
                        7*Power(t2,2)) + 
                     s1*(-22 - 4*Power(s2,2) + t1 - 6*Power(t1,2) + 
                        s2*(3 + 16*t1 - 13*t2) + 9*t2 - t1*t2 + 
                        6*Power(t2,2)))) - 
               (2*(5 + Power(s,2) - Power(s1,2) - 3*s2 + 
                    Power(s2,2) + 2*t1 - s2*t1 + 
                    s*(3 - 2*s2 + 2*t1) + t2 - s2*t2 + 2*t1*t2 - 
                    Power(t2,2) + s1*(-1 + s2 - 2*t1 + 2*t2))*
                  (-2*Power(s2,2) + Power(s1,5)*Power(s2 - t1,2) + 
                    4*s2*t1 + 6*Power(s2,2)*t1 - 2*Power(t1,2) - 
                    12*s2*Power(t1,2) - 6*Power(s2,2)*Power(t1,2) + 
                    6*Power(t1,3) + 12*s2*Power(t1,3) + 
                    2*Power(s2,2)*Power(t1,3) - 6*Power(t1,4) - 
                    4*s2*Power(t1,4) + 2*Power(t1,5) + 4*s2*t2 - 
                    2*Power(s2,2)*t2 - 4*Power(s2,3)*t2 - 4*t1*t2 - 
                    2*s2*t1*t2 + 12*Power(s2,2)*t1*t2 + 
                    8*Power(s2,3)*t1*t2 + 4*Power(t1,2)*t2 - 
                    12*s2*Power(t1,2)*t2 - 
                    18*Power(s2,2)*Power(t1,2)*t2 - 
                    4*Power(s2,3)*Power(t1,2)*t2 + 4*Power(t1,3)*t2 + 
                    14*s2*Power(t1,3)*t2 + 
                    8*Power(s2,2)*Power(t1,3)*t2 - 4*Power(t1,4)*t2 - 
                    4*s2*Power(t1,4)*t2 + Power(t2,2) - 
                    4*s2*Power(t2,2) + 4*Power(s2,2)*Power(t2,2) - 
                    4*Power(s2,3)*Power(t2,2) - 
                    2*Power(s2,4)*Power(t2,2) + t1*Power(t2,2) + 
                    4*s2*t1*Power(t2,2) + 
                    8*Power(s2,2)*t1*Power(t2,2) + 
                    8*Power(s2,3)*t1*Power(t2,2) + 
                    2*Power(s2,4)*t1*Power(t2,2) - 
                    7*Power(t1,2)*Power(t2,2) - 
                    8*s2*Power(t1,2)*Power(t2,2) - 
                    14*Power(s2,2)*Power(t1,2)*Power(t2,2) - 
                    4*Power(s2,3)*Power(t1,2)*Power(t2,2) + 
                    5*Power(t1,3)*Power(t2,2) + 
                    8*s2*Power(t1,3)*Power(t2,2) + 
                    2*Power(s2,2)*Power(t1,3)*Power(t2,2) - 
                    3*Power(t2,3) + 4*s2*Power(t2,3) - 
                    6*Power(s2,2)*Power(t2,3) - 
                    2*Power(s2,4)*Power(t2,3) + 2*t1*Power(t2,3) + 
                    6*s2*t1*Power(t2,3) + 
                    6*Power(s2,2)*t1*Power(t2,3) + 
                    6*Power(s2,3)*t1*Power(t2,3) - 
                    3*Power(t1,2)*Power(t2,3) - 
                    8*s2*Power(t1,2)*Power(t2,3) - 
                    4*Power(s2,2)*Power(t1,2)*Power(t2,3) + 
                    3*Power(t2,4) - 4*s2*Power(t2,4) - 
                    Power(s2,2)*Power(t2,4) - 
                    2*Power(s2,3)*Power(t2,4) + t1*Power(t2,4) + 
                    4*s2*t1*Power(t2,4) + 
                    3*Power(s2,2)*t1*Power(t2,4) - Power(t2,5) - 
                    Power(s2,2)*Power(t2,5) + 
                    Power(s,4)*(-1 + s1 + t1 - t2)*
                     (2 + 2*Power(s2,2) + Power(t1,2) - 
                       2*s2*(2 + t1 - t2) - 2*t1*(-1 + t2) - 2*t2 + 
                       Power(t2,2)) + 
                    Power(s1,4)*(s2 - t1)*
                     (-2 + 2*Power(s2,2) - t1 - Power(t1,2) + 2*t2 + 
                       3*t1*t2 - s2*(-1 + t1 + 3*t2)) + 
                    Power(s1,3)*
                     (2*Power(s2,4) - 2*Power(t1,3)*(-2 + t2) + 
                       Power(-1 + t2,2) - 
                       2*Power(s2,3)*(2 + t1 + 2*t2) + 
                       Power(t1,2)*(-3 - 4*t2 + 3*Power(t2,2)) + 
                       t1*(2 - 8*t2 + 6*Power(t2,2)) - 
                       2*Power(s2,2)*
                       (2 + Power(t1,2) + 2*t2 - 2*Power(t2,2) - 
                       t1*(5 + 4*t2)) + 
                       2*s2*
                        (-1 + Power(t1,3) + 4*t2 - 3*Power(t2,2) - 
                        Power(t1,2)*(5 + t2) + 
                        t1*(4 + 3*t2 - 3*Power(t2,2)))) + 
                    Power(s1,2)*
                     (2*Power(t1,4) + 
                       2*Power(s2,4)*(-1 + t1 - 3*t2) - 
                       Power(-1 + t2,2)*(-1 + 3*t2) + 
                       Power(t1,3)*(-5 - 6*t2 + Power(t2,2)) + 
                       t1*(-5 + 11*Power(t2,2) - 6*Power(t2,3)) - 
                       Power(t1,2)*
                       (-7 + t2 - 5*Power(t2,2) + Power(t2,3)) + 
                       2*Power(s2,2)*
                       (4 - 5*Power(t1,2) + Power(t1,3) + t2 + 
                       2*Power(t2,2) - 2*Power(t2,3) - 
                       t1*t2*(7 + 2*t2)) + 
                       Power(s2,3)*
                       (-4 - 4*Power(t1,2) + 8*t2 + 2*t1*(4 + 5*t2)) \
- 2*s2*(-1 - 3*t2 + 7*Power(t2,2) - 3*Power(t2,3) + 
                        Power(t1,3)*(-1 + 2*t2) - 
                        Power(t1,2)*(5 + 5*t2 + 2*Power(t2,2)) + 
                        t1*(7 + 3*t2 + Power(t2,2) - Power(t2,3)))) + 
                    s1*(-2*Power(t1,4)*(-2 + t2) + 
                       Power(-1 + t2,2)*t2*(-2 + 3*t2) + 
                       2*Power(s2,4)*t2*(2 - 2*t1 + 3*t2) + 
                       2*Power(t1,3)*(-2 + Power(t2,2)) + 
                       Power(t1,2)*
                       (-4 + 7*Power(t2,2) - 2*Power(t2,3)) + 
                       2*t1*
                       (2 + 2*t2 - 2*Power(t2,2) - 3*Power(t2,3) + 
                       Power(t2,4)) + 
                       2*s2*
                       (-2 + 2*Power(t1,4) + t2 - 4*Power(t2,2) + 
                       6*Power(t2,3) - Power(t2,4) + 
                       Power(t1,3)*(-7 - 5*t2 + Power(t2,2)) + 
                       t1*(1 + 5*t2 - 4*Power(t2,2) - 
                       3*Power(t2,3)) - 
                       Power(t1,2)*
                       (-6 + t2 - 4*Power(t2,2) + Power(t2,3))) + 
                       2*Power(s2,3)*
                       (Power(t1,2)*(2 + 4*t2) - 
                       t1*(4 + 8*t2 + 7*Power(t2,2)) + 
                       2*(1 + 2*t2 - Power(t2,2) + Power(t2,3))) + 
                       Power(s2,2)*
                        (2 - 12*t2 + 8*Power(t2,2) + 3*Power(t2,4) - 
                        4*Power(t1,3)*(2 + t2) + 
                        6*Power(t1,2)*(3 + 4*t2 + Power(t2,2)) - 
                        2*t1*(6 + 4*t2 + Power(t2,2) + 2*Power(t2,3))\
)) - 2*Power(s,3)*(-3 + 3*t1 + Power(t1,3) + 
                       2*Power(s2,3)*(-1 + t1 - t2) - 
                       Power(s2,2)*(2 + 3*t1)*(-1 + t1 - t2) + 
                       2*t2 - 2*t1*t2 - 5*Power(t1,2)*t2 - 
                       Power(t1,3)*t2 + 2*Power(t2,2) + 
                       7*t1*Power(t2,2) + 
                       3*Power(t1,2)*Power(t2,2) - 3*Power(t2,3) - 
                       3*t1*Power(t2,3) + Power(t2,4) + 
                       Power(s1,2)*
                       (2 + 3*Power(s2,2) + 2*Power(t1,2) + 
                       t1*(4 - 3*t2) - 2*t2 + Power(t2,2) + 
                       s2*(-5 - 4*t1 + 2*t2)) + 
                       s2*(3 + Power(t1,3) - 
                       Power(t1,2)*(-1 + t2) - 2*t2 - 
                       5*Power(t2,2) + Power(t2,3) - 
                       t1*(6 - 4*t2 + Power(t2,2))) + 
                       s1*(1 + 2*Power(s2,3) + 2*Power(t1,3) + 
                        Power(t1,2)*(3 - 6*t2) - 5*t2 + 
                        5*Power(t2,2) - 2*Power(t2,3) - 
                        Power(s2,2)*(5 + 3*t2) + 
                        t1*(-2 - 8*t2 + 6*Power(t2,2)) + 
                        s2*(2 - 3*Power(t1,2) + 8*t2 - 
                        3*Power(t2,2) + t1*(2 + 6*t2)))) + 
                    Power(s,2)*
                     (-5 + 13*t1 - 15*Power(t1,2) + 5*Power(t1,3) + 
                       2*Power(t1,4) + 
                       2*Power(s2,4)*(-1 + t1 - t2) + 7*t2 - 
                       6*t1*t2 - 9*Power(t1,2)*t2 - 
                       8*Power(t1,3)*t2 + 3*Power(t2,2) + 
                       15*t1*Power(t2,2) + 
                       17*Power(t1,2)*Power(t2,2) + 
                       Power(t1,3)*Power(t2,2) - 11*Power(t2,3) - 
                       18*t1*Power(t2,3) - 
                       3*Power(t1,2)*Power(t2,3) + 7*Power(t2,4) + 
                       3*t1*Power(t2,4) - Power(t2,5) + 
                       Power(s1,3)*
                       (2 + 7*Power(s2,2) + 6*Power(t1,2) + 
                       t1*(10 - 6*t2) - 2*s2*(4 + 6*t1 - t2) - 
                       2*t2 + Power(t2,2)) + 
                       Power(s2,3)*
                       (-4*Power(t1,2) - 2*t1*(-2 + t2) + 
                       6*t2*(1 + t2)) - 
                       2*s2*
                       (2 + 2*t2 + 4*Power(t2,2) - 4*Power(t2,3) - 
                       Power(t2,4) + 2*Power(t1,3)*(2 + t2) + 
                       t1*t2*(-10 + t2 + 4*Power(t2,2)) - 
                       Power(t1,2)*(6 + t2 + 5*Power(t2,2))) + 
                       Power(s2,2)*
                       (10 + 2*Power(t1,3) - 8*t2 - 9*Power(t2,2) + 
                       5*Power(t2,3) + Power(t1,2)*(4 + 8*t2) - 
                       t1*(16 + 2*t2 + 15*Power(t2,2))) + 
                       Power(s1,2)*
                       (6 + 10*Power(s2,3) + 6*Power(t1,3) - 
                       14*t2 + 11*Power(t2,2) - 3*Power(t2,3) - 
                       Power(s2,2)*(15 + 9*t1 + 7*t2) - 
                       2*Power(t1,2)*(-5 + 9*t2) + 
                       t1*(6 - 34*t2 + 15*Power(t2,2)) - 
                       2*s2*
                       (2 + 3*Power(t1,2) - 10*t2 + Power(t2,2) - 
                       t1*(5 + 9*t2))) + 
                       s1*(-3 + 2*Power(s2,4) + 
                        2*Power(s2,3)*(-5 + 3*t1 - 8*t2) + 
                        Power(t1,3)*(8 - 6*t2) - 10*t2 + 
                        23*Power(t2,2) - 16*Power(t2,3) + 
                        3*Power(t2,4) + 
                        Power(t1,2)*(21 - 34*t2 + 15*Power(t2,2)) - 
                        2*t1*
                       (5 + 7*t2 - 21*Power(t2,2) + 6*Power(t2,3)) + 
                        Power(s2,2)*
                       (-2 - 14*Power(t1,2) + 22*t2 - 
                       5*Power(t2,2) + 2*t1*(9 + 13*t2)) + 
                        2*s2*
                        (8 + 3*Power(t1,3) + 3*t2 - 10*Power(t2,2) - 
                        Power(t2,3) - Power(t1,2)*(7 + 3*t2) + 
                        t1*(-11 + Power(t2,2))))) - 
                    2*s*(4*t1 - 9*Power(t1,2) + 6*Power(t1,3) - 
                       Power(t1,4) + 
                       Power(s1,4)*
                        (2*Power(s2,2) - s2*(1 + 4*t1) + 
                        t1*(2 + 2*t1 - t2)) + t2 - 2*t1*t2 + 
                       4*Power(t1,2)*t2 - 2*Power(t1,3)*t2 - 
                       Power(t1,4)*t2 + t1*Power(t2,2) + 
                       5*Power(t1,2)*Power(t2,2) + 
                       3*Power(t1,3)*Power(t2,2) - 4*Power(t2,3) - 
                       6*t1*Power(t2,3) - 4*Power(t1,2)*Power(t2,3) + 
                       4*Power(t2,4) + 3*t1*Power(t2,4) - 
                       Power(t2,5) + 2*Power(s2,4)*t2*(1 - t1 + t2) + 
                       Power(s2,3)*
                        (2 + 2*t2 - Power(t2,2) + 
                        Power(t1,2)*(2 + 4*t2) - 
                        2*t1*(2 + 3*t2 + 2*Power(t2,2))) + 
                       Power(s1,3)*
                        (4*Power(s2,3) + 2*Power(t1,3) + 
                        Power(t1,2)*(3 - 6*t2) + Power(-1 + t2,2) - 
                        Power(s2,2)*(3 + 5*t1 + 3*t2) + 
                        s2*(-3 + t1 - Power(t1,2) + 4*t2 + 
                       8*t1*t2 + Power(t2,2)) + 
                        t1*(6 - 11*t2 + 3*Power(t2,2))) + 
                       Power(s2,2)*
                        (-2*Power(t1,3)*(2 + t2) + 
                        Power(t1,2)*(8 + 5*t2 + Power(t2,2)) - 
                        t2*
                       (7 - 5*t2 - 4*Power(t2,2) + Power(t2,3)) + 
                        t1*(-4 + 4*t2 - 3*Power(t2,2) + 
                       2*Power(t2,3))) + 
                       s2*(-4 + 2*Power(t1,4) + 2*t2 - 
                        3*Power(t2,2) + 6*Power(t2,3) - 
                        Power(t2,5) + 
                        Power(t1,3)*(-3 + Power(t2,2)) - 
                        Power(t1,2)*
                       (4 + 4*t2 - 2*Power(t2,2) + 3*Power(t2,3)) + 
                        t1*(9 + 2*t2 - 7*Power(t2,2) - 
                        4*Power(t2,3) + 3*Power(t2,4))) + 
                       Power(s1,2)*
                        (1 + 2*Power(s2,4) + 
                        Power(t1,3)*(5 - 3*t2) - 6*t2 + 
                        8*Power(t2,2) - 3*Power(t2,3) - 
                        Power(s2,3)*(5 + 8*t2) + 
                        2*Power(t1,2)*(4 - 7*t2 + 3*Power(t2,2)) - 
                        t1*(5 + 15*t2 - 19*Power(t2,2) + 
                       3*Power(t2,3)) - 
                        Power(s2,2)*
                       (5 + 5*Power(t1,2) - 8*t2 + Power(t2,2) - 
                       t1*(13 + 14*t2)) + 
                        s2*(4 + 3*Power(t1,3) + 8*t2 - 
                        5*Power(t2,2) - 3*Power(t2,3) - 
                        Power(t1,2)*(13 + 3*t2) - 
                        t1*(1 + Power(t2,2)))) + 
                       s1*(-1 + 2*Power(t1,4) + 
                        2*Power(s2,4)*(-1 + t1 - 2*t2) - t2 + 
                        9*Power(t2,2) - 10*Power(t2,3) + 
                        3*Power(t2,4) + 
                        Power(t1,3)*(-3 - 7*t2 + Power(t2,2)) + 
                        Power(t1,2)*
                        (3 - 15*t2 + 15*Power(t2,2) - 2*Power(t2,3)) \
+ t1*(-1 + 5*t2 + 15*Power(t2,2) - 13*Power(t2,3) + Power(t2,4)) + 
                        Power(s2,3)*
                        (-2 - 4*Power(t1,2) + 6*t2 + 4*Power(t2,2) + 
                        t1*(6 + 4*t2)) + 
                        s2*(1 - 2*t2 - 11*Power(t2,2) + 
                        2*Power(t2,3) + 3*Power(t2,4) - 
                        Power(t1,3)*(3 + 4*t2) + 
                        Power(t1,2)*(13 + 10*t2 + 7*Power(t2,2)) + 
                        t1*(-11 + 10*t2 + 3*Power(t2,2) - 
                        6*Power(t2,3))) + 
                        Power(s2,2)*
                        (2*Power(t1,3) + Power(t1,2)*(-3 + 4*t2) - 
                        t1*(8 + 10*t2 + 11*Power(t2,2)) + 
                        3*(3 - 3*Power(t2,2) + Power(t2,3)))))))/
                ((-1 + s1 + t1 - t2)*
                  Power(-1 + s - s*s2 + s1*s2 + t1 - s2*t2,2))))/
           ((-1 + t1)*(-1 + t2)) + 
          (32*((2*(-2*s2 + Power(s1,2)*(s2 - t1) + 2*t1 + 2*s2*t1 - 
                    2*Power(t1,2) + 
                    s1*(2*Power(s2,2) - 2*s2*t1 + 
                       (1 + t1)*(-1 + t2)) + t2 - 
                    2*Power(s2,2)*t2 + t1*t2 + 2*s2*t1*t2 - 
                    Power(t2,2) - s2*Power(t2,2) + 
                    Power(s,2)*(-2 + 2*s2 - t1 + t2) + 
                    s*(3 - 2*Power(s2,2) - t1 + 2*s2*t1 - 3*t2 + 
                       s2*t2 - t1*t2 + Power(t2,2) - 
                       s1*(-2 + 3*s2 - 2*t1 + t2)))*
                  (-22 + 2*Power(s,3) + 2*Power(s1,3) - 6*s2 + 
                    4*Power(s2,2) - 2*Power(s2,3) + t1 - 5*s2*t1 + 
                    4*Power(s2,2)*t1 + 2*Power(t1,2) - 
                    2*s2*Power(t1,2) + 
                    Power(s1,2)*(2*s2 + t1 - 5*t2) + 
                    Power(s,2)*(2 + s1 - 6*s2 + 8*t1 - 3*t2) + 
                    15*t2 + Power(s2,2)*t2 - t1*t2 - 6*s2*t1*t2 + 
                    4*Power(t1,2)*t2 - 4*Power(t2,2) + 
                    4*s2*Power(t2,2) - 2*t1*Power(t2,2) - 
                    Power(t2,3) + 
                    s1*(-14 - 3*Power(s2,2) - 4*Power(t1,2) + 
                       s2*(5 + 8*t1 - 6*t2) + t1*(-3 + t2) + 4*t2 + 
                       4*Power(t2,2)) + 
                    s*(10 - 5*Power(s1,2) + 6*Power(s2,2) + 3*t1 + 
                       4*Power(t1,2) - 4*t2 + 6*t1*t2 - 
                       6*Power(t2,2) + s2*(-7 - 11*t1 + t2) + 
                       s1*(-4 + 3*s2 - 9*t1 + 11*t2))))/
                (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2) - 
               4*(13 + Power(s,4) + 24*s2 + 2*Power(s2,2) - 
                  Power(s2,3) + Power(s2,4) - 12*t1 + 
                  Power(s2,2)*t1 - 3*Power(s2,3)*t1 - 3*Power(t1,2) + 
                  3*Power(s2,2)*Power(t1,2) - s2*Power(t1,3) + 
                  Power(s1,3)*(3 - 2*s2 + 2*t1) + 
                  Power(s,3)*(-1 + s1 - 4*s2 + 6*t1 - 3*t2) - 16*t2 - 
                  21*s2*t2 - Power(s2,2)*t2 + 10*t1*t2 + 4*s2*t1*t2 + 
                  3*Power(s2,2)*t1*t2 - 2*Power(t1,2)*t2 - 
                  5*s2*Power(t1,2)*t2 + 2*Power(t1,3)*t2 + 
                  11*Power(t2,2) + 5*s2*Power(t2,2) - 
                  3*Power(s2,2)*Power(t2,2) - 3*t1*Power(t2,2) + 
                  4*s2*t1*Power(t2,2) - Power(t1,2)*Power(t2,2) + 
                  s2*Power(t2,3) - t1*Power(t2,3) - 
                  Power(s1,2)*
                   (-14 + Power(s2,2) - 2*Power(t1,2) + 
                     s2*(t1 - 6*t2) + 8*t2 + t1*(-7 + 6*t2)) + 
                  Power(s,2)*
                   (8 - 4*Power(s1,2) + 6*Power(s2,2) - 3*t1 - 
                     14*s2*t1 + 7*Power(t1,2) - 4*t2 + 5*s2*t2 + 
                     t1*t2 - 5*Power(t2,2) + 
                     s1*(-4 + s2 - 6*t1 + 10*t2)) + 
                  s1*(2*Power(s2,3) - 2*Power(t1,3) - 
                     Power(t1,2)*(-1 + t2) + 
                     Power(s2,2)*(-4 - 7*t1 + 5*t2) + 
                     s2*(18 + 7*Power(t1,2) + t1*(7 - 4*t2) - 
                       7*t2 - 5*Power(t2,2)) + 
                     5*(3 - 5*t2 + Power(t2,2)) + 
                     t1*(-3 - 2*t2 + 5*Power(t2,2))) + 
                  s*(-20 + 2*Power(s1,3) - 4*Power(s2,3) - 
                     Power(t1,2) + 2*Power(t1,3) + 
                     Power(s1,2)*(2 + 5*s2 - 2*t1 - 6*t2) + 
                     Power(s2,2)*(2 + 11*t1 - 2*t2) + 26*t2 - 
                     6*t1*t2 + 6*Power(t1,2)*t2 - 4*Power(t2,2) - 
                     6*t1*Power(t2,2) - Power(t2,3) + 
                     s1*(-22 - 4*Power(s2,2) - 4*t1 - 
                        9*Power(t1,2) + s2*(8 + 14*t1 - 15*t2) + 
                        6*t2 + 9*t1*t2 + 5*Power(t2,2)) + 
                     s2*(-9 + t1 - 9*Power(t1,2) + 6*t2 - 5*t1*t2 + 
                        8*Power(t2,2)))) - 
               (2*(5 + Power(s,2) - Power(s1,2) - 3*s2 + 
                    Power(s2,2) + 2*t1 - s2*t1 + 
                    s*(3 - 2*s2 + 2*t1) + t2 - s2*t2 + 2*t1*t2 - 
                    Power(t2,2) + s1*(-1 + s2 - 2*t1 + 2*t2))*
                  (-2*Power(s2,2) + Power(s1,5)*Power(s2 - t1,2) + 
                    4*s2*t1 + 6*Power(s2,2)*t1 - 2*Power(t1,2) - 
                    12*s2*Power(t1,2) - 6*Power(s2,2)*Power(t1,2) + 
                    6*Power(t1,3) + 12*s2*Power(t1,3) + 
                    2*Power(s2,2)*Power(t1,3) - 6*Power(t1,4) - 
                    4*s2*Power(t1,4) + 2*Power(t1,5) + 4*s2*t2 - 
                    2*Power(s2,2)*t2 - 4*Power(s2,3)*t2 - 4*t1*t2 - 
                    2*s2*t1*t2 + 12*Power(s2,2)*t1*t2 + 
                    8*Power(s2,3)*t1*t2 + 4*Power(t1,2)*t2 - 
                    12*s2*Power(t1,2)*t2 - 
                    18*Power(s2,2)*Power(t1,2)*t2 - 
                    4*Power(s2,3)*Power(t1,2)*t2 + 4*Power(t1,3)*t2 + 
                    14*s2*Power(t1,3)*t2 + 
                    8*Power(s2,2)*Power(t1,3)*t2 - 4*Power(t1,4)*t2 - 
                    4*s2*Power(t1,4)*t2 + Power(t2,2) - 
                    4*s2*Power(t2,2) + 4*Power(s2,2)*Power(t2,2) - 
                    4*Power(s2,3)*Power(t2,2) - 
                    2*Power(s2,4)*Power(t2,2) + t1*Power(t2,2) + 
                    4*s2*t1*Power(t2,2) + 
                    8*Power(s2,2)*t1*Power(t2,2) + 
                    8*Power(s2,3)*t1*Power(t2,2) + 
                    2*Power(s2,4)*t1*Power(t2,2) - 
                    7*Power(t1,2)*Power(t2,2) - 
                    8*s2*Power(t1,2)*Power(t2,2) - 
                    14*Power(s2,2)*Power(t1,2)*Power(t2,2) - 
                    4*Power(s2,3)*Power(t1,2)*Power(t2,2) + 
                    5*Power(t1,3)*Power(t2,2) + 
                    8*s2*Power(t1,3)*Power(t2,2) + 
                    2*Power(s2,2)*Power(t1,3)*Power(t2,2) - 
                    3*Power(t2,3) + 4*s2*Power(t2,3) - 
                    6*Power(s2,2)*Power(t2,3) - 
                    2*Power(s2,4)*Power(t2,3) + 2*t1*Power(t2,3) + 
                    6*s2*t1*Power(t2,3) + 
                    6*Power(s2,2)*t1*Power(t2,3) + 
                    6*Power(s2,3)*t1*Power(t2,3) - 
                    3*Power(t1,2)*Power(t2,3) - 
                    8*s2*Power(t1,2)*Power(t2,3) - 
                    4*Power(s2,2)*Power(t1,2)*Power(t2,3) + 
                    3*Power(t2,4) - 4*s2*Power(t2,4) - 
                    Power(s2,2)*Power(t2,4) - 
                    2*Power(s2,3)*Power(t2,4) + t1*Power(t2,4) + 
                    4*s2*t1*Power(t2,4) + 
                    3*Power(s2,2)*t1*Power(t2,4) - Power(t2,5) - 
                    Power(s2,2)*Power(t2,5) + 
                    Power(s,4)*(-1 + s1 + t1 - t2)*
                     (2 + 2*Power(s2,2) + Power(t1,2) - 
                       2*s2*(2 + t1 - t2) - 2*t1*(-1 + t2) - 2*t2 + 
                       Power(t2,2)) + 
                    Power(s1,4)*(s2 - t1)*
                     (-2 + 2*Power(s2,2) - t1 - Power(t1,2) + 2*t2 + 
                       3*t1*t2 - s2*(-1 + t1 + 3*t2)) + 
                    Power(s1,3)*
                     (2*Power(s2,4) - 2*Power(t1,3)*(-2 + t2) + 
                       Power(-1 + t2,2) - 
                       2*Power(s2,3)*(2 + t1 + 2*t2) + 
                       Power(t1,2)*(-3 - 4*t2 + 3*Power(t2,2)) + 
                       t1*(2 - 8*t2 + 6*Power(t2,2)) - 
                       2*Power(s2,2)*
                       (2 + Power(t1,2) + 2*t2 - 2*Power(t2,2) - 
                       t1*(5 + 4*t2)) + 
                       2*s2*
                        (-1 + Power(t1,3) + 4*t2 - 3*Power(t2,2) - 
                        Power(t1,2)*(5 + t2) + 
                        t1*(4 + 3*t2 - 3*Power(t2,2)))) + 
                    Power(s1,2)*
                     (2*Power(t1,4) + 
                       2*Power(s2,4)*(-1 + t1 - 3*t2) - 
                       Power(-1 + t2,2)*(-1 + 3*t2) + 
                       Power(t1,3)*(-5 - 6*t2 + Power(t2,2)) + 
                       t1*(-5 + 11*Power(t2,2) - 6*Power(t2,3)) - 
                       Power(t1,2)*
                       (-7 + t2 - 5*Power(t2,2) + Power(t2,3)) + 
                       2*Power(s2,2)*
                       (4 - 5*Power(t1,2) + Power(t1,3) + t2 + 
                       2*Power(t2,2) - 2*Power(t2,3) - 
                       t1*t2*(7 + 2*t2)) + 
                       Power(s2,3)*
                       (-4 - 4*Power(t1,2) + 8*t2 + 2*t1*(4 + 5*t2)) \
- 2*s2*(-1 - 3*t2 + 7*Power(t2,2) - 3*Power(t2,3) + 
                        Power(t1,3)*(-1 + 2*t2) - 
                        Power(t1,2)*(5 + 5*t2 + 2*Power(t2,2)) + 
                        t1*(7 + 3*t2 + Power(t2,2) - Power(t2,3)))) + 
                    s1*(-2*Power(t1,4)*(-2 + t2) + 
                       Power(-1 + t2,2)*t2*(-2 + 3*t2) + 
                       2*Power(s2,4)*t2*(2 - 2*t1 + 3*t2) + 
                       2*Power(t1,3)*(-2 + Power(t2,2)) + 
                       Power(t1,2)*
                       (-4 + 7*Power(t2,2) - 2*Power(t2,3)) + 
                       2*t1*
                       (2 + 2*t2 - 2*Power(t2,2) - 3*Power(t2,3) + 
                       Power(t2,4)) + 
                       2*s2*
                       (-2 + 2*Power(t1,4) + t2 - 4*Power(t2,2) + 
                       6*Power(t2,3) - Power(t2,4) + 
                       Power(t1,3)*(-7 - 5*t2 + Power(t2,2)) + 
                       t1*(1 + 5*t2 - 4*Power(t2,2) - 
                       3*Power(t2,3)) - 
                       Power(t1,2)*
                       (-6 + t2 - 4*Power(t2,2) + Power(t2,3))) + 
                       2*Power(s2,3)*
                       (Power(t1,2)*(2 + 4*t2) - 
                       t1*(4 + 8*t2 + 7*Power(t2,2)) + 
                       2*(1 + 2*t2 - Power(t2,2) + Power(t2,3))) + 
                       Power(s2,2)*
                        (2 - 12*t2 + 8*Power(t2,2) + 3*Power(t2,4) - 
                        4*Power(t1,3)*(2 + t2) + 
                        6*Power(t1,2)*(3 + 4*t2 + Power(t2,2)) - 
                        2*t1*(6 + 4*t2 + Power(t2,2) + 2*Power(t2,3))\
)) - 2*Power(s,3)*(-3 + 3*t1 + Power(t1,3) + 
                       2*Power(s2,3)*(-1 + t1 - t2) - 
                       Power(s2,2)*(2 + 3*t1)*(-1 + t1 - t2) + 
                       2*t2 - 2*t1*t2 - 5*Power(t1,2)*t2 - 
                       Power(t1,3)*t2 + 2*Power(t2,2) + 
                       7*t1*Power(t2,2) + 
                       3*Power(t1,2)*Power(t2,2) - 3*Power(t2,3) - 
                       3*t1*Power(t2,3) + Power(t2,4) + 
                       Power(s1,2)*
                       (2 + 3*Power(s2,2) + 2*Power(t1,2) + 
                       t1*(4 - 3*t2) - 2*t2 + Power(t2,2) + 
                       s2*(-5 - 4*t1 + 2*t2)) + 
                       s2*(3 + Power(t1,3) - 
                       Power(t1,2)*(-1 + t2) - 2*t2 - 
                       5*Power(t2,2) + Power(t2,3) - 
                       t1*(6 - 4*t2 + Power(t2,2))) + 
                       s1*(1 + 2*Power(s2,3) + 2*Power(t1,3) + 
                        Power(t1,2)*(3 - 6*t2) - 5*t2 + 
                        5*Power(t2,2) - 2*Power(t2,3) - 
                        Power(s2,2)*(5 + 3*t2) + 
                        t1*(-2 - 8*t2 + 6*Power(t2,2)) + 
                        s2*(2 - 3*Power(t1,2) + 8*t2 - 
                        3*Power(t2,2) + t1*(2 + 6*t2)))) + 
                    Power(s,2)*
                     (-5 + 13*t1 - 15*Power(t1,2) + 5*Power(t1,3) + 
                       2*Power(t1,4) + 
                       2*Power(s2,4)*(-1 + t1 - t2) + 7*t2 - 
                       6*t1*t2 - 9*Power(t1,2)*t2 - 
                       8*Power(t1,3)*t2 + 3*Power(t2,2) + 
                       15*t1*Power(t2,2) + 
                       17*Power(t1,2)*Power(t2,2) + 
                       Power(t1,3)*Power(t2,2) - 11*Power(t2,3) - 
                       18*t1*Power(t2,3) - 
                       3*Power(t1,2)*Power(t2,3) + 7*Power(t2,4) + 
                       3*t1*Power(t2,4) - Power(t2,5) + 
                       Power(s1,3)*
                       (2 + 7*Power(s2,2) + 6*Power(t1,2) + 
                       t1*(10 - 6*t2) - 2*s2*(4 + 6*t1 - t2) - 
                       2*t2 + Power(t2,2)) + 
                       Power(s2,3)*
                       (-4*Power(t1,2) - 2*t1*(-2 + t2) + 
                       6*t2*(1 + t2)) - 
                       2*s2*
                       (2 + 2*t2 + 4*Power(t2,2) - 4*Power(t2,3) - 
                       Power(t2,4) + 2*Power(t1,3)*(2 + t2) + 
                       t1*t2*(-10 + t2 + 4*Power(t2,2)) - 
                       Power(t1,2)*(6 + t2 + 5*Power(t2,2))) + 
                       Power(s2,2)*
                       (10 + 2*Power(t1,3) - 8*t2 - 9*Power(t2,2) + 
                       5*Power(t2,3) + Power(t1,2)*(4 + 8*t2) - 
                       t1*(16 + 2*t2 + 15*Power(t2,2))) + 
                       Power(s1,2)*
                       (6 + 10*Power(s2,3) + 6*Power(t1,3) - 
                       14*t2 + 11*Power(t2,2) - 3*Power(t2,3) - 
                       Power(s2,2)*(15 + 9*t1 + 7*t2) - 
                       2*Power(t1,2)*(-5 + 9*t2) + 
                       t1*(6 - 34*t2 + 15*Power(t2,2)) - 
                       2*s2*
                       (2 + 3*Power(t1,2) - 10*t2 + Power(t2,2) - 
                       t1*(5 + 9*t2))) + 
                       s1*(-3 + 2*Power(s2,4) + 
                        2*Power(s2,3)*(-5 + 3*t1 - 8*t2) + 
                        Power(t1,3)*(8 - 6*t2) - 10*t2 + 
                        23*Power(t2,2) - 16*Power(t2,3) + 
                        3*Power(t2,4) + 
                        Power(t1,2)*(21 - 34*t2 + 15*Power(t2,2)) - 
                        2*t1*
                       (5 + 7*t2 - 21*Power(t2,2) + 6*Power(t2,3)) + 
                        Power(s2,2)*
                       (-2 - 14*Power(t1,2) + 22*t2 - 
                       5*Power(t2,2) + 2*t1*(9 + 13*t2)) + 
                        2*s2*
                        (8 + 3*Power(t1,3) + 3*t2 - 10*Power(t2,2) - 
                        Power(t2,3) - Power(t1,2)*(7 + 3*t2) + 
                        t1*(-11 + Power(t2,2))))) - 
                    2*s*(4*t1 - 9*Power(t1,2) + 6*Power(t1,3) - 
                       Power(t1,4) + 
                       Power(s1,4)*
                        (2*Power(s2,2) - s2*(1 + 4*t1) + 
                        t1*(2 + 2*t1 - t2)) + t2 - 2*t1*t2 + 
                       4*Power(t1,2)*t2 - 2*Power(t1,3)*t2 - 
                       Power(t1,4)*t2 + t1*Power(t2,2) + 
                       5*Power(t1,2)*Power(t2,2) + 
                       3*Power(t1,3)*Power(t2,2) - 4*Power(t2,3) - 
                       6*t1*Power(t2,3) - 4*Power(t1,2)*Power(t2,3) + 
                       4*Power(t2,4) + 3*t1*Power(t2,4) - 
                       Power(t2,5) + 2*Power(s2,4)*t2*(1 - t1 + t2) + 
                       Power(s2,3)*
                        (2 + 2*t2 - Power(t2,2) + 
                        Power(t1,2)*(2 + 4*t2) - 
                        2*t1*(2 + 3*t2 + 2*Power(t2,2))) + 
                       Power(s1,3)*
                        (4*Power(s2,3) + 2*Power(t1,3) + 
                        Power(t1,2)*(3 - 6*t2) + Power(-1 + t2,2) - 
                        Power(s2,2)*(3 + 5*t1 + 3*t2) + 
                        s2*(-3 + t1 - Power(t1,2) + 4*t2 + 
                       8*t1*t2 + Power(t2,2)) + 
                        t1*(6 - 11*t2 + 3*Power(t2,2))) + 
                       Power(s2,2)*
                        (-2*Power(t1,3)*(2 + t2) + 
                        Power(t1,2)*(8 + 5*t2 + Power(t2,2)) - 
                        t2*
                       (7 - 5*t2 - 4*Power(t2,2) + Power(t2,3)) + 
                        t1*(-4 + 4*t2 - 3*Power(t2,2) + 
                       2*Power(t2,3))) + 
                       s2*(-4 + 2*Power(t1,4) + 2*t2 - 
                        3*Power(t2,2) + 6*Power(t2,3) - 
                        Power(t2,5) + 
                        Power(t1,3)*(-3 + Power(t2,2)) - 
                        Power(t1,2)*
                       (4 + 4*t2 - 2*Power(t2,2) + 3*Power(t2,3)) + 
                        t1*(9 + 2*t2 - 7*Power(t2,2) - 
                        4*Power(t2,3) + 3*Power(t2,4))) + 
                       Power(s1,2)*
                        (1 + 2*Power(s2,4) + 
                        Power(t1,3)*(5 - 3*t2) - 6*t2 + 
                        8*Power(t2,2) - 3*Power(t2,3) - 
                        Power(s2,3)*(5 + 8*t2) + 
                        2*Power(t1,2)*(4 - 7*t2 + 3*Power(t2,2)) - 
                        t1*(5 + 15*t2 - 19*Power(t2,2) + 
                       3*Power(t2,3)) - 
                        Power(s2,2)*
                       (5 + 5*Power(t1,2) - 8*t2 + Power(t2,2) - 
                       t1*(13 + 14*t2)) + 
                        s2*(4 + 3*Power(t1,3) + 8*t2 - 
                        5*Power(t2,2) - 3*Power(t2,3) - 
                        Power(t1,2)*(13 + 3*t2) - 
                        t1*(1 + Power(t2,2)))) + 
                       s1*(-1 + 2*Power(t1,4) + 
                        2*Power(s2,4)*(-1 + t1 - 2*t2) - t2 + 
                        9*Power(t2,2) - 10*Power(t2,3) + 
                        3*Power(t2,4) + 
                        Power(t1,3)*(-3 - 7*t2 + Power(t2,2)) + 
                        Power(t1,2)*
                        (3 - 15*t2 + 15*Power(t2,2) - 2*Power(t2,3)) \
+ t1*(-1 + 5*t2 + 15*Power(t2,2) - 13*Power(t2,3) + Power(t2,4)) + 
                        Power(s2,3)*
                        (-2 - 4*Power(t1,2) + 6*t2 + 4*Power(t2,2) + 
                        t1*(6 + 4*t2)) + 
                        s2*(1 - 2*t2 - 11*Power(t2,2) + 
                        2*Power(t2,3) + 3*Power(t2,4) - 
                        Power(t1,3)*(3 + 4*t2) + 
                        Power(t1,2)*(13 + 10*t2 + 7*Power(t2,2)) + 
                        t1*(-11 + 10*t2 + 3*Power(t2,2) - 
                        6*Power(t2,3))) + 
                        Power(s2,2)*
                        (2*Power(t1,3) + Power(t1,2)*(-3 + 4*t2) - 
                        t1*(8 + 10*t2 + 11*Power(t2,2)) + 
                        3*(3 - 3*Power(t2,2) + Power(t2,3)))))))/
                ((-1 + s1 + t1 - t2)*
                  Power(-1 + s - s*s2 + s1*s2 + t1 - s2*t2,2))))/
           ((-1 + s1)*(-1 + t2)) + 
          (64*(-(((-2*s2 + Power(s1,2)*(s2 - t1) + 2*t1 + 2*s2*t1 - 
                      2*Power(t1,2) + 
                      s1*(2*Power(s2,2) - 2*s2*t1 + 
                       (1 + t1)*(-1 + t2)) + t2 - 2*Power(s2,2)*t2 + 
                      t1*t2 + 2*s2*t1*t2 - Power(t2,2) - 
                      s2*Power(t2,2) + 
                      Power(s,2)*(-2 + 2*s2 - t1 + t2) + 
                      s*(3 - 2*Power(s2,2) - t1 + 2*s2*t1 - 3*t2 + 
                        s2*t2 - t1*t2 + Power(t2,2) - 
                        s1*(-2 + 3*s2 - 2*t1 + t2)))*
                    (25 + 4*Power(s,3) + 5*s2 - 2*Power(s2,2) - 
                      7*t1 - s2*t1 - 7*t2 + 4*Power(s2,2)*t2 + 
                      7*t1*t2 - 5*s2*t1*t2 + 2*Power(t1,2)*t2 - 
                      7*Power(t2,2) - s2*Power(t2,2) + 
                      t1*Power(t2,2) - Power(t2,3) - 
                      Power(s1,2)*(7 + 4*s2 - 3*t1 + t2) + 
                      Power(s,2)*(1 - 7*s1 - 8*s2 + 6*t1 + 4*t2) + 
                      s*(-14 + 3*Power(s1,2) + 4*Power(s2,2) + 
                        6*t1 - 5*s2*t1 + 2*Power(t1,2) + 
                        s1*(3 + 12*s2 - 9*t1 - 2*t2) - 5*t2 - 
                        9*s2*t2 + 7*t1*t2 - Power(t2,2)) + 
                      s1*(5 - 4*Power(s2,2) - 9*t1 - 2*Power(t1,2) + 
                        14*t2 - 4*t1*t2 + 2*Power(t2,2) + 
                        5*s2*(t1 + t2))))/
                  (1 + s*(-1 + s2) - s1*s2 - t1 + s2*t2)) - 
               2*(12 - 2*Power(s,4) + Power(s1,4) + 30*s2 + 
                  4*Power(s2,2) - 18*t1 - 9*s2*t1 - 2*Power(s2,2)*t1 + 
                  5*Power(t1,2) + s2*Power(t1,2) + Power(t1,3) + 
                  Power(s1,3)*(-1 + t1 - 3*t2) + 
                  Power(s,3)*(4*s1 + 6*s2 - 5*t1 - t2) - 4*t2 - 
                  23*s2*t2 + 2*Power(s2,3)*t2 + 17*t1*t2 + 
                  5*s2*t1*t2 - 4*Power(s2,2)*t1*t2 - 
                  7*Power(t1,2)*t2 + 3*s2*Power(t1,2)*t2 - 
                  Power(t1,3)*t2 - 9*Power(t2,2) - 3*s2*Power(t2,2) + 
                  5*t1*Power(t2,2) - Power(t2,3) - s2*Power(t2,3) + 
                  t1*Power(t2,3) + 
                  Power(s1,2)*
                   (5 - 3*Power(s2,2) + 2*t1 - 2*Power(t1,2) + t2 + 
                     3*Power(t2,2) + s2*(5*t1 - 2*(1 + t2))) + 
                  Power(s,2)*
                   (13 - Power(s1,2) + s2 - 6*Power(s2,2) - 6*t1 + 
                     9*s2*t1 - 4*Power(t1,2) + 4*t2 + 5*s2*t2 - 
                     4*t1*t2 + 2*Power(t2,2) - 
                     s1*(11*s2 + 2*(1 - 5*t1 + t2))) - 
                  s*(22 + 2*Power(s1,3) - 2*Power(s2,3) - 14*t1 + 
                     6*Power(t1,2) + Power(t1,3) - 12*t2 + 
                     4*Power(t1,2)*t2 - 2*Power(t2,2) - 
                     2*t1*Power(t2,2) - Power(t2,3) - 
                     Power(s1,2)*(3 + 5*s2 - 6*t1 + 6*t2) + 
                     Power(s2,2)*(1 + 4*t1 + 6*t2) + 
                     s2*(19 - 8*t1 - 3*Power(t1,2) + 4*t2 - 
                        9*t1*t2 + 2*Power(t2,2)) + 
                     s1*(17 - 9*Power(s2,2) - 5*t1 - 6*Power(t1,2) + 
                        7*t2 - 3*t1*t2 + 5*Power(t2,2) + 
                        2*s2*(7*t1 + t2))) + 
                  s1*(8 - 2*Power(s2,3) + Power(t1,3) + 4*t2 + 
                     Power(t2,2) - Power(t2,3) + 
                     2*Power(t1,2)*(4 + t2) + 
                     Power(s2,2)*(2 + 4*t1 + 3*t2) - 
                     t1*(8 + 9*t2 + 2*Power(t2,2)) + 
                     s2*(18 - 3*Power(t1,2) + 7*t2 + 3*Power(t2,2) - 
                        t1*(8 + 5*t2)))) + 
               ((-3 + 2*Power(s,2) + Power(s1,2) + 2*s2 + t1 + 
                    s1*(2*s2 - t1 - 2*t2) - 2*s2*t2 + t1*t2 + 
                    Power(t2,2) + s*(1 - 3*s1 - 2*s2 + t1 + 3*t2))*
                  (-2*Power(s2,2) + Power(s1,5)*Power(s2 - t1,2) + 
                    4*s2*t1 + 6*Power(s2,2)*t1 - 2*Power(t1,2) - 
                    12*s2*Power(t1,2) - 6*Power(s2,2)*Power(t1,2) + 
                    6*Power(t1,3) + 12*s2*Power(t1,3) + 
                    2*Power(s2,2)*Power(t1,3) - 6*Power(t1,4) - 
                    4*s2*Power(t1,4) + 2*Power(t1,5) + 4*s2*t2 - 
                    2*Power(s2,2)*t2 - 4*Power(s2,3)*t2 - 4*t1*t2 - 
                    2*s2*t1*t2 + 12*Power(s2,2)*t1*t2 + 
                    8*Power(s2,3)*t1*t2 + 4*Power(t1,2)*t2 - 
                    12*s2*Power(t1,2)*t2 - 
                    18*Power(s2,2)*Power(t1,2)*t2 - 
                    4*Power(s2,3)*Power(t1,2)*t2 + 4*Power(t1,3)*t2 + 
                    14*s2*Power(t1,3)*t2 + 
                    8*Power(s2,2)*Power(t1,3)*t2 - 4*Power(t1,4)*t2 - 
                    4*s2*Power(t1,4)*t2 + Power(t2,2) - 
                    4*s2*Power(t2,2) + 4*Power(s2,2)*Power(t2,2) - 
                    4*Power(s2,3)*Power(t2,2) - 
                    2*Power(s2,4)*Power(t2,2) + t1*Power(t2,2) + 
                    4*s2*t1*Power(t2,2) + 
                    8*Power(s2,2)*t1*Power(t2,2) + 
                    8*Power(s2,3)*t1*Power(t2,2) + 
                    2*Power(s2,4)*t1*Power(t2,2) - 
                    7*Power(t1,2)*Power(t2,2) - 
                    8*s2*Power(t1,2)*Power(t2,2) - 
                    14*Power(s2,2)*Power(t1,2)*Power(t2,2) - 
                    4*Power(s2,3)*Power(t1,2)*Power(t2,2) + 
                    5*Power(t1,3)*Power(t2,2) + 
                    8*s2*Power(t1,3)*Power(t2,2) + 
                    2*Power(s2,2)*Power(t1,3)*Power(t2,2) - 
                    3*Power(t2,3) + 4*s2*Power(t2,3) - 
                    6*Power(s2,2)*Power(t2,3) - 
                    2*Power(s2,4)*Power(t2,3) + 2*t1*Power(t2,3) + 
                    6*s2*t1*Power(t2,3) + 
                    6*Power(s2,2)*t1*Power(t2,3) + 
                    6*Power(s2,3)*t1*Power(t2,3) - 
                    3*Power(t1,2)*Power(t2,3) - 
                    8*s2*Power(t1,2)*Power(t2,3) - 
                    4*Power(s2,2)*Power(t1,2)*Power(t2,3) + 
                    3*Power(t2,4) - 4*s2*Power(t2,4) - 
                    Power(s2,2)*Power(t2,4) - 
                    2*Power(s2,3)*Power(t2,4) + t1*Power(t2,4) + 
                    4*s2*t1*Power(t2,4) + 
                    3*Power(s2,2)*t1*Power(t2,4) - Power(t2,5) - 
                    Power(s2,2)*Power(t2,5) + 
                    Power(s,4)*(-1 + s1 + t1 - t2)*
                     (2 + 2*Power(s2,2) + Power(t1,2) - 
                       2*s2*(2 + t1 - t2) - 2*t1*(-1 + t2) - 2*t2 + 
                       Power(t2,2)) + 
                    Power(s1,4)*(s2 - t1)*
                     (-2 + 2*Power(s2,2) - t1 - Power(t1,2) + 2*t2 + 
                       3*t1*t2 - s2*(-1 + t1 + 3*t2)) + 
                    Power(s1,3)*
                     (2*Power(s2,4) - 2*Power(t1,3)*(-2 + t2) + 
                       Power(-1 + t2,2) - 
                       2*Power(s2,3)*(2 + t1 + 2*t2) + 
                       Power(t1,2)*(-3 - 4*t2 + 3*Power(t2,2)) + 
                       t1*(2 - 8*t2 + 6*Power(t2,2)) - 
                       2*Power(s2,2)*
                        (2 + Power(t1,2) + 2*t2 - 2*Power(t2,2) - 
                        t1*(5 + 4*t2)) + 
                       2*s2*(-1 + Power(t1,3) + 4*t2 - 
                        3*Power(t2,2) - Power(t1,2)*(5 + t2) + 
                        t1*(4 + 3*t2 - 3*Power(t2,2)))) + 
                    Power(s1,2)*
                     (2*Power(t1,4) + 
                       2*Power(s2,4)*(-1 + t1 - 3*t2) - 
                       Power(-1 + t2,2)*(-1 + 3*t2) + 
                       Power(t1,3)*(-5 - 6*t2 + Power(t2,2)) + 
                       t1*(-5 + 11*Power(t2,2) - 6*Power(t2,3)) - 
                       Power(t1,2)*
                        (-7 + t2 - 5*Power(t2,2) + Power(t2,3)) + 
                       2*Power(s2,2)*
                        (4 - 5*Power(t1,2) + Power(t1,3) + t2 + 
                        2*Power(t2,2) - 2*Power(t2,3) - 
                        t1*t2*(7 + 2*t2)) + 
                       Power(s2,3)*
                        (-4 - 4*Power(t1,2) + 8*t2 + 2*t1*(4 + 5*t2)) \
- 2*s2*(-1 - 3*t2 + 7*Power(t2,2) - 3*Power(t2,3) + 
                        Power(t1,3)*(-1 + 2*t2) - 
                        Power(t1,2)*(5 + 5*t2 + 2*Power(t2,2)) + 
                        t1*(7 + 3*t2 + Power(t2,2) - Power(t2,3)))) + 
                    s1*(-2*Power(t1,4)*(-2 + t2) + 
                       Power(-1 + t2,2)*t2*(-2 + 3*t2) + 
                       2*Power(s2,4)*t2*(2 - 2*t1 + 3*t2) + 
                       2*Power(t1,3)*(-2 + Power(t2,2)) + 
                       Power(t1,2)*
                        (-4 + 7*Power(t2,2) - 2*Power(t2,3)) + 
                       2*t1*
                        (2 + 2*t2 - 2*Power(t2,2) - 3*Power(t2,3) + 
                        Power(t2,4)) + 
                       2*s2*
                        (-2 + 2*Power(t1,4) + t2 - 4*Power(t2,2) + 
                        6*Power(t2,3) - Power(t2,4) + 
                        Power(t1,3)*(-7 - 5*t2 + Power(t2,2)) + 
                        t1*(1 + 5*t2 - 4*Power(t2,2) - 
                       3*Power(t2,3)) - 
                        Power(t1,2)*
                        (-6 + t2 - 4*Power(t2,2) + Power(t2,3))) + 
                       2*Power(s2,3)*
                        (Power(t1,2)*(2 + 4*t2) - 
                        t1*(4 + 8*t2 + 7*Power(t2,2)) + 
                        2*(1 + 2*t2 - Power(t2,2) + Power(t2,3))) + 
                       Power(s2,2)*
                        (2 - 12*t2 + 8*Power(t2,2) + 3*Power(t2,4) - 
                        4*Power(t1,3)*(2 + t2) + 
                        6*Power(t1,2)*(3 + 4*t2 + Power(t2,2)) - 
                        2*t1*(6 + 4*t2 + Power(t2,2) + 2*Power(t2,3)))\
) - 2*Power(s,3)*(-3 + 3*t1 + Power(t1,3) + 
                       2*Power(s2,3)*(-1 + t1 - t2) - 
                       Power(s2,2)*(2 + 3*t1)*(-1 + t1 - t2) + 2*t2 - 
                       2*t1*t2 - 5*Power(t1,2)*t2 - Power(t1,3)*t2 + 
                       2*Power(t2,2) + 7*t1*Power(t2,2) + 
                       3*Power(t1,2)*Power(t2,2) - 3*Power(t2,3) - 
                       3*t1*Power(t2,3) + Power(t2,4) + 
                       Power(s1,2)*
                        (2 + 3*Power(s2,2) + 2*Power(t1,2) + 
                        t1*(4 - 3*t2) - 2*t2 + Power(t2,2) + 
                        s2*(-5 - 4*t1 + 2*t2)) + 
                       s2*(3 + Power(t1,3) - Power(t1,2)*(-1 + t2) - 
                        2*t2 - 5*Power(t2,2) + Power(t2,3) - 
                        t1*(6 - 4*t2 + Power(t2,2))) + 
                       s1*(1 + 2*Power(s2,3) + 2*Power(t1,3) + 
                        Power(t1,2)*(3 - 6*t2) - 5*t2 + 
                        5*Power(t2,2) - 2*Power(t2,3) - 
                        Power(s2,2)*(5 + 3*t2) + 
                        t1*(-2 - 8*t2 + 6*Power(t2,2)) + 
                        s2*(2 - 3*Power(t1,2) + 8*t2 - 
                        3*Power(t2,2) + t1*(2 + 6*t2)))) + 
                    Power(s,2)*
                     (-5 + 13*t1 - 15*Power(t1,2) + 5*Power(t1,3) + 
                       2*Power(t1,4) + 2*Power(s2,4)*(-1 + t1 - t2) + 
                       7*t2 - 6*t1*t2 - 9*Power(t1,2)*t2 - 
                       8*Power(t1,3)*t2 + 3*Power(t2,2) + 
                       15*t1*Power(t2,2) + 
                       17*Power(t1,2)*Power(t2,2) + 
                       Power(t1,3)*Power(t2,2) - 11*Power(t2,3) - 
                       18*t1*Power(t2,3) - 
                       3*Power(t1,2)*Power(t2,3) + 7*Power(t2,4) + 
                       3*t1*Power(t2,4) - Power(t2,5) + 
                       Power(s1,3)*
                        (2 + 7*Power(s2,2) + 6*Power(t1,2) + 
                        t1*(10 - 6*t2) - 2*s2*(4 + 6*t1 - t2) - 
                        2*t2 + Power(t2,2)) + 
                       Power(s2,3)*
                        (-4*Power(t1,2) - 2*t1*(-2 + t2) + 
                        6*t2*(1 + t2)) - 
                       2*s2*
                        (2 + 2*t2 + 4*Power(t2,2) - 4*Power(t2,3) - 
                        Power(t2,4) + 2*Power(t1,3)*(2 + t2) + 
                        t1*t2*(-10 + t2 + 4*Power(t2,2)) - 
                        Power(t1,2)*(6 + t2 + 5*Power(t2,2))) + 
                       Power(s2,2)*
                        (10 + 2*Power(t1,3) - 8*t2 - 9*Power(t2,2) + 
                        5*Power(t2,3) + Power(t1,2)*(4 + 8*t2) - 
                        t1*(16 + 2*t2 + 15*Power(t2,2))) + 
                       Power(s1,2)*
                        (6 + 10*Power(s2,3) + 6*Power(t1,3) - 
                        14*t2 + 11*Power(t2,2) - 3*Power(t2,3) - 
                        Power(s2,2)*(15 + 9*t1 + 7*t2) - 
                        2*Power(t1,2)*(-5 + 9*t2) + 
                        t1*(6 - 34*t2 + 15*Power(t2,2)) - 
                        2*s2*
                        (2 + 3*Power(t1,2) - 10*t2 + Power(t2,2) - 
                        t1*(5 + 9*t2))) + 
                       s1*(-3 + 2*Power(s2,4) + 
                        2*Power(s2,3)*(-5 + 3*t1 - 8*t2) + 
                        Power(t1,3)*(8 - 6*t2) - 10*t2 + 
                        23*Power(t2,2) - 16*Power(t2,3) + 
                        3*Power(t2,4) + 
                        Power(t1,2)*(21 - 34*t2 + 15*Power(t2,2)) - 
                        2*t1*
                        (5 + 7*t2 - 21*Power(t2,2) + 6*Power(t2,3)) + 
                        Power(s2,2)*
                        (-2 - 14*Power(t1,2) + 22*t2 - 
                        5*Power(t2,2) + 2*t1*(9 + 13*t2)) + 
                        2*s2*
                        (8 + 3*Power(t1,3) + 3*t2 - 10*Power(t2,2) - 
                        Power(t2,3) - Power(t1,2)*(7 + 3*t2) + 
                        t1*(-11 + Power(t2,2))))) - 
                    2*s*(4*t1 - 9*Power(t1,2) + 6*Power(t1,3) - 
                       Power(t1,4) + 
                       Power(s1,4)*
                        (2*Power(s2,2) - s2*(1 + 4*t1) + 
                        t1*(2 + 2*t1 - t2)) + t2 - 2*t1*t2 + 
                       4*Power(t1,2)*t2 - 2*Power(t1,3)*t2 - 
                       Power(t1,4)*t2 + t1*Power(t2,2) + 
                       5*Power(t1,2)*Power(t2,2) + 
                       3*Power(t1,3)*Power(t2,2) - 4*Power(t2,3) - 
                       6*t1*Power(t2,3) - 4*Power(t1,2)*Power(t2,3) + 
                       4*Power(t2,4) + 3*t1*Power(t2,4) - 
                       Power(t2,5) + 2*Power(s2,4)*t2*(1 - t1 + t2) + 
                       Power(s2,3)*
                        (2 + 2*t2 - Power(t2,2) + 
                        Power(t1,2)*(2 + 4*t2) - 
                        2*t1*(2 + 3*t2 + 2*Power(t2,2))) + 
                       Power(s1,3)*
                        (4*Power(s2,3) + 2*Power(t1,3) + 
                        Power(t1,2)*(3 - 6*t2) + Power(-1 + t2,2) - 
                        Power(s2,2)*(3 + 5*t1 + 3*t2) + 
                        s2*(-3 + t1 - Power(t1,2) + 4*t2 + 8*t1*t2 + 
                        Power(t2,2)) + t1*(6 - 11*t2 + 3*Power(t2,2))) \
+ Power(s2,2)*(-2*Power(t1,3)*(2 + t2) + 
                        Power(t1,2)*(8 + 5*t2 + Power(t2,2)) - 
                        t2*(7 - 5*t2 - 4*Power(t2,2) + Power(t2,3)) + 
                        t1*(-4 + 4*t2 - 3*Power(t2,2) + 2*Power(t2,3))\
) + s2*(-4 + 2*Power(t1,4) + 2*t2 - 3*Power(t2,2) + 6*Power(t2,3) - 
                        Power(t2,5) + 
                        Power(t1,3)*(-3 + Power(t2,2)) - 
                        Power(t1,2)*
                        (4 + 4*t2 - 2*Power(t2,2) + 3*Power(t2,3)) + 
                        t1*(9 + 2*t2 - 7*Power(t2,2) - 
                        4*Power(t2,3) + 3*Power(t2,4))) + 
                       Power(s1,2)*
                        (1 + 2*Power(s2,4) + Power(t1,3)*(5 - 3*t2) - 
                        6*t2 + 8*Power(t2,2) - 3*Power(t2,3) - 
                        Power(s2,3)*(5 + 8*t2) + 
                        2*Power(t1,2)*(4 - 7*t2 + 3*Power(t2,2)) - 
                        t1*(5 + 15*t2 - 19*Power(t2,2) + 
                        3*Power(t2,3)) - 
                        Power(s2,2)*
                        (5 + 5*Power(t1,2) - 8*t2 + Power(t2,2) - 
                        t1*(13 + 14*t2)) + 
                        s2*(4 + 3*Power(t1,3) + 8*t2 - 
                        5*Power(t2,2) - 3*Power(t2,3) - 
                        Power(t1,2)*(13 + 3*t2) - t1*(1 + Power(t2,2))\
)) + s1*(-1 + 2*Power(t1,4) + 2*Power(s2,4)*(-1 + t1 - 2*t2) - t2 + 
                        9*Power(t2,2) - 10*Power(t2,3) + 
                        3*Power(t2,4) + 
                        Power(t1,3)*(-3 - 7*t2 + Power(t2,2)) + 
                        Power(t1,2)*
                        (3 - 15*t2 + 15*Power(t2,2) - 2*Power(t2,3)) + 
                        t1*(-1 + 5*t2 + 15*Power(t2,2) - 
                        13*Power(t2,3) + Power(t2,4)) + 
                        Power(s2,3)*
                        (-2 - 4*Power(t1,2) + 6*t2 + 4*Power(t2,2) + 
                        t1*(6 + 4*t2)) + 
                        s2*(1 - 2*t2 - 11*Power(t2,2) + 
                        2*Power(t2,3) + 3*Power(t2,4) - 
                        Power(t1,3)*(3 + 4*t2) + 
                        Power(t1,2)*(13 + 10*t2 + 7*Power(t2,2)) + 
                        t1*(-11 + 10*t2 + 3*Power(t2,2) - 
                        6*Power(t2,3))) + 
                        Power(s2,2)*
                        (2*Power(t1,3) + Power(t1,2)*(-3 + 4*t2) - 
                        t1*(8 + 10*t2 + 11*Power(t2,2)) + 
                        3*(3 - 3*Power(t2,2) + Power(t2,3)))))))/
                ((-1 + s1 + t1 - t2)*
                  Power(-1 + s - s*s2 + s1*s2 + t1 - s2*t2,2))))/
           ((-1 + s1)*(-s + s1 - t2)))/(s - s2 + t1))*
     B1(1 - s + s1 - t2,1 - s1 - t1 + t2,1 - s + s2 - t1))/
   (128.*Power(Pi,2)) + (((-64*s*
          (t1*(1 + s*(-1 + t2) - 3*t2) - 2*Power(s1,2)*(-1 + t2) - 
            (1 + t2)*(1 + (-1 + s - s2)*t2) - 
            s1*(t1*(-3 + t2) + s2*(1 + t2) - 2*t2*(-1 + s + t2)))*
          (-((5 + Power(s1,4) - 2*s2 + 3*Power(s2,2) + 7*t1 - 
                 8*s2*t1 - Power(s2,2)*t1 - 2*Power(t1,2) + 
                 2*s2*Power(t1,2) - 2*t2 + 2*Power(s2,2)*t2 - 
                 3*t1*t2 - s2*t1*t2 + Power(t2,2) + 
                 Power(s1,3)*(-4 + 2*s2 + t1 + t2) + 
                 Power(s,2)*
                  (1 + 2*Power(s1,2) - 6*t1 + Power(t1,2) + 
                    s1*(-5 + 3*t1) + 2*t2) + 
                 s1*(5 + Power(s2,2)*t1 + 4*Power(t1,2) + 
                    s2*(t1*(-6 + t2) + 2*(-2 + t2)) - 8*t2 + 
                    Power(t2,2) - 4*t1*(3 + t2)) + 
                 Power(s1,2)*
                  (Power(s2,2) + t1*(-2 + t2) - 2*(5 + t2) + 
                    s2*(-5 + 2*t1 + t2)) - 
                 s*(1 + 3*Power(s1,3) - 7*t1 + 2*Power(t1,2) - t2 - 
                    2*t1*t2 + Power(s1,2)*(-10 + 3*s2 + 4*t1 + t2) + 
                    s1*(-7 + Power(t1,2) + s2*(-5 + 4*t1) + 
                       t1*(-11 + t2) + t2) + 
                    s2*(4 - 7*t1 + Power(t1,2) + 4*t2)))/
               ((-1 + s2)*(-s + s2 - t1))) + 
            (-11 + 3*Power(s1,3) + 6*s2 + Power(s2,2) + 11*t1 - 
               s2*t1 - 4*Power(s,2)*(-2 + t2) + 8*t2 + s2*t2 - 
               2*Power(s2,2)*t2 + 3*t1*t2 + Power(t2,2) - 
               Power(s1,2)*(7 - 3*t1 + s2*(-3 + t2) + 3*t2) + 
               s1*(8 - 14*t1 + s2*(7 + t1) - Power(s2,2)*(-2 + t2) + 
                  14*t2 - 3*t1*t2 + 3*Power(t2,2)) - 
               s*(6 + Power(s1,2) + 9*s2 - t1 + t2 - 6*s2*t2 + 
                  s1*(14 + 2*s2 + t1 - 2*t2 - s2*t2)))/
             ((-1 + s1)*(-s + s1 - t2)) + 
            (-2 + Power(s1,3) + 11*s2 + 2*t1 - 2*s2*t1 + 3*t2 - 
               4*Power(s,2)*t2 + 4*s2*t2 + Power(s2,2)*t2 - 3*t1*t2 - 
               Power(t2,2) + 
               Power(s1,2)*(-3 + t1 - s2*(-2 + t2) + 4*t2) + 
               s1*(12 + 3*t1*(-1 + t2) + 5*t2 - Power(s2,2)*t2 - 
                  3*Power(t2,2) + s2*(-1 + 2*t1 + t2)) + 
               s*(-11 - 2*Power(s1,2) + 2*t1 - 4*t2 + 3*s2*t2 + 
                  s1*(1 - 2*t1 + 4*t2 + s2*t2)))/((-1 + s1)*(-1 + t2)) \
- (1 - 8*s2 + Power(s2,2) - 4*t1 + s2*t1 + Power(t1,2) + 
               2*Power(s1,3)*(-1 + t2) - 4*t2 + s2*t2 + 2*t1*t2 - 
               s2*t1*t2 + 3*Power(t2,2) + s2*Power(t2,2) + 
               Power(s,2)*(-1 + (3 + s1 - t1)*t2) - 
               Power(s1,2)*
                (6 + 3*t1 - 3*(-1 + s2)*t2 + Power(t2,2)) - 
               s1*(6 + Power(t1,2) + t1*(-6 + t2) + 4*t2 - 
                  Power(s2,2)*t2 - 2*Power(t2,2) + 
                  s2*(8 + 4*t2 + Power(t2,2))) + 
               s*(8 - 3*Power(s1,2)*t2 - 3*s2*t2 - 2*Power(t2,2) + 
                  t1*(-1 + t2 + s2*t2) + 
                  s1*(10 + (1 - 2*s2 + t1)*t2 + Power(t2,2))))/
             ((-1 + t1)*(-1 + t2)) + 
            (-4 + 2*s2 - 3*Power(s2,2) - 2*t1 + 7*s2*t1 + 
               2*Power(t1,2) - Power(s1,3)*(-5 + t2) + 9*t2 + 
               4*s2*t2 - 4*Power(s2,2)*t2 + 3*t1*t2 - 2*s2*t1*t2 + 
               3*Power(t2,2) + s2*Power(t2,2) - 
               Power(s,2)*(3 + t1*(-6 + t2) + 2*s1*(-3 + t2) + 
                  2*t2) - Power(s1,2)*
                (-5 - 2*t1 + 2*s2*(-4 + t2) + 4*t2 + Power(t2,2)) - 
               s1*(11 + 2*Power(t1,2) + t1*(-10 + t2) - 9*t2 + 
                  Power(s2,2)*t2 - 2*Power(t2,2) + 
                  s2*(-2 - 7*t1 + 10*t2 + Power(t2,2))) + 
               s*(1 - 7*t1 - 8*t2 + 2*t1*t2 - 2*Power(t2,2) + 
                  Power(s1,2)*(-14 + 3*t2) + 
                  s2*(t1*(-6 + t2) + 6*(1 + t2)) + 
                  s1*(t1*(-13 + t2) + 3*s2*(-2 + t2) + t2*(8 + t2))))/
             ((s - s2 + t1)*(s - s1 + t2)) + 
            (-11 - 2*Power(s1,4) + 3*s2 + 4*Power(s2,2) + 3*t1 - 
               4*s2*t1 - Power(s2,2)*t1 - 2*Power(t1,2) + 
               2*s2*Power(t1,2) + 8*t2 + 6*s2*t2 - 3*Power(s2,2)*t2 - 
               3*t1*t2 - s2*t1*t2 + Power(t2,2) + 
               Power(s1,3)*(4 - 3*s2 - 2*t1 + t2) + 
               Power(s,2)*(9 - s1 - Power(s1,2) - 5*t1 + 
                  Power(t1,2) + t2) + 
               Power(s1,2)*(11 - Power(s2,2) - t2 + t1*(4 + t2) + 
                  s2*(4 - 3*t1 + t2)) + 
               s1*(1 - Power(s2,2)*(-2 + t1) + 4*Power(t1,2) + 
                  s2*(12 + t1*(-1 + t2) - t2) + 7*t2 + Power(t2,2) - 
                  t1*(5 + 4*t2)) + 
               s*(-3 + 3*Power(s1,3) + 3*t1 - 2*Power(t1,2) - 
                  s2*(13 - 6*t1 + Power(t1,2) - 2*t2) + 
                  Power(s1,2)*(-1 + 2*s2 + 2*t1 - t2) - 6*t2 + 
                  2*t1*t2 + s1*
                   (s2*(-1 + t1) - Power(t1,2) - t1*(-5 + t2) - 
                     2*(9 + t2))))/((-1 + s2)*(-1 + t1))))/
        (Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
          Power(s,2)*(-1 + s1)*(t1*(-1 + t2) + (s1 - t2)*t2) - 
          s*(-1 + s1)*(Power(s1,2)*(-1 + t2) + 2*(1 + t1 - t2)*t2 + 
             s2*(1 + t1*(-1 + t2) - 2*t2 - Power(t2,2)) + 
             s1*(t1*(-3 + t2) - Power(-1 + t2,2) + s2*(1 + t2)))) + 
       64*((-12 + 12*Power(s,2) + 10*s1 - 36*s*s1 + 16*Power(s1,2) + 
             2*s2 - 12*s*s2 + 14*s1*s2 + 30*t1 - 12*s*t1 + 2*s1*t1 + 
             20*t2 + 14*s*t2 - 4*Power(s,2)*t2 - 10*s1*t2 + 
             8*s*s1*t2 - 4*Power(s1,2)*t2 - 18*s2*t2 + 4*s*s2*t2 - 
             4*s1*s2*t2 + 4*t1*t2 + 4*Power(t2,2) + 2*s*Power(t2,2) - 
             4*s1*Power(t2,2) - 2*s2*Power(t2,2) + 
             ((-4 - s*(-4 + t2) + s1*(-3 + t2) + t2 + Power(t2,2))*
                (2*Power(s1,2)*(-1 + t2) + 
                  t1*(-1 + s + 3*t2 - s*t2) + 
                  (1 + t2)*(1 + (-1 + s - s2)*t2) + 
                  s1*(t1*(-3 + t2) + s2*(1 + t2) - 2*t2*(-1 + s + t2))\
))/(s1*(-1 + t2) - t2*(-1 + s + t2)))/((s - s2 + t1)*(s - s1 + t2)) + 
          (1 - t1 + 6*Power(s1,3)*(-1 + t2) + 19*t2 - 5*s2*t2 - 
             2*t1*t2 + 3*Power(t2,2) - 2*s2*Power(t2,2) + 
             3*t1*Power(t2,2) - 17*Power(t2,3) + 3*s2*Power(t2,3) + 
             4*t1*Power(t2,3) - 6*Power(t2,4) - 
             Power(s1,2)*(-14 - 6*t2 + 20*Power(t2,2) + 
                s2*Power(1 + t2,2) + t1*(5 - 10*t2 + Power(t2,2))) + 
             Power(s,2)*(-(t1*Power(-1 + t2,2)) + 
                t2*(13 + 2*s1 + 2*s2*(-2 + t2) - 10*t2 + Power(t2,2))) \
+ s1*(-19 - 19*t2 + 17*Power(t2,2) + 21*Power(t2,3) + 
                t1*(2 + 4*t2 - 14*Power(t2,2)) + 
                s2*(5 + 3*t2 - Power(t2,2) + Power(t2,3))) + 
             s*(-1 + 2*t1 + Power(s1,2)*(2 - 8*t2) - 30*t2 + 9*s2*t2 + 
                5*Power(t2,2) - 2*s2*Power(t2,2) + 6*t1*Power(t2,2) - 
                16*Power(t2,3) + s2*Power(t2,3) - 
                s1*(-14 + 13*t2 - 24*Power(t2,2) + Power(t2,3) - 
                   2*t1*(1 - 6*t2 + Power(t2,2)) + 
                   s2*(5 - 6*t2 + Power(t2,2)))))/
           ((-1 + s1)*(-s + s1 - t2)*(s1*(-1 + t2) - t2*(-1 + s + t2))) \
+ (1 - t1 + 2*Power(s1,3)*(-1 + t2) + 23*t2 - 7*s2*t2 - 2*t1*t2 - 
             7*Power(t2,2) + 11*t1*Power(t2,2) - 23*Power(t2,3) + 
             3*s2*Power(t2,3) - 4*t1*Power(t2,3) + 6*Power(t2,4) - 
             Power(s1,2)*(s2*(3 + Power(t2,2)) + 
                t1*(-3 - 2*t2 + Power(t2,2)) + 
                2*(-5 + 4*t2 + Power(t2,2))) + 
             Power(s,2)*(t1 - t1*Power(t2,2) + 
                t2*(-7 + 4*s1 + 2*(-2 + s2)*t2 + Power(t2,2))) + 
             s1*(-23 - 5*t2 + 33*Power(t2,2) - 5*Power(t2,3) + 
                2*t1*(1 - 6*t2 + Power(t2,2)) + 
                s2*(7 + 3*t2 - 3*Power(t2,2) + Power(t2,3))) + 
             s*(1 + Power(s1,2)*(4 - 6*t2) + (5*s2 + 6*(-2 + t1))*t2 - 
                (17 + 2*t1)*Power(t2,2) + s2*Power(t2,3) + 
                s1*(-8 + 7*t2 + 8*Power(t2,2) - Power(t2,3) + 
                   s2*(1 + 2*t2 - Power(t2,2)) + 
                   2*t1*(-2 - t2 + Power(t2,2)))))/
           ((-1 + s1)*(-1 + t2)*(s1*(-1 + t2) - t2*(-1 + s + t2))) + 
          (2 - 2*t1 + 4*Power(s1,3)*Power(-1 + t2,2) - 12*t2 - 
             8*s2*t2 + 4*t1*t2 + 7*Power(t2,2) - 
             2*Power(s,3)*Power(t2,2) - 4*s2*Power(t2,2) + 
             9*t1*Power(t2,2) + 8*Power(t2,3) + s2*Power(t2,3) - 
             t1*Power(t2,3) - 5*Power(t2,4) + s2*Power(t2,4) - 
             Power(s,2)*(t1*(-1 + t2 + 2*Power(t2,2)) + 
                t2*(7 + s1*(2 - 8*t2) + t2 - 4*s2*t2 + 2*Power(t2,2))) \
+ Power(s1,2)*(4 - 2*t2 + 4*Power(t2,2) - 6*Power(t2,3) - 
                t1*(-5 + t2 + 2*Power(t2,2)) + 
                s2*(1 - 7*t2 + 4*Power(t2,2))) + 
             s1*(15 - 14*t2 - 13*Power(t2,2) + 10*Power(t2,3) + 
                2*Power(t2,4) + 
                s2*(8 + 3*t2 + 6*Power(t2,2) - 5*Power(t2,3)) + 
                t1*(-7 - 11*t2 + 5*Power(t2,2) + Power(t2,3))) + 
             s*(1 + 26*t2 + 5*s2*t2 - 10*Power(s1,2)*(-1 + t2)*t2 + 
                7*Power(t2,2) - s2*Power(t2,2) - 7*Power(t2,3) + 
                4*s2*Power(t2,3) - Power(t2,4) - 
                t1*(-1 - 3*t2 + Power(t2,2) + Power(t2,3)) + 
                s1*(-8 + 7*t2 + 3*Power(t2,2) + 8*Power(t2,3) + 
                   s2*(1 + 7*t2 - 8*Power(t2,2)) + 
                   2*t1*(-2 + t2 + 2*Power(t2,2)))))/
           ((-1 + t1)*(-1 + t2)*(s1 - s1*t2 + t2*(-1 + s + t2))) - 
          (-2*(-1 + t2)*Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2)*
              (-4*s1 + Power(s1,2) - (-4 + t2)*t2) - 
             2*Power(s,5)*t2*(t1 + 2*t1*t2 + t2*(-4 + 2*s1 + t2)) + 
             Power(s,4)*(-(Power(t1,2)*(-1 + t2)) + 
                2*Power(s1,2)*t2*(-4 + 7*t2) + 
                2*t1*t2*(s2 + 2*s2*t2 - 5*Power(t2,2)) + 
                t2*(-2 - (31 + 8*s2)*t2 + (19 + 4*s2)*Power(t2,2)) + 
                s1*(-(t2*(-18 + 2*s2 + 29*t2 - 4*s2*t2 + 
                       Power(t2,2))) + t1*(-2 + t2 + 13*Power(t2,2)))) \
- Power(s,3)*(2*Power(s1,3)*(2 - 10*t2 + 9*Power(t2,2)) + 
                Power(t1,2)*
                 (-2 - 3*t2 + 2*Power(t2,2) + Power(t2,3)) + 
                t1*(-2 + (-10 + 3*s2)*t2 + 
                   (-21 + 5*s2)*Power(t2,2) + 
                   (2 - 10*s2)*Power(t2,3) + Power(t2,4)) + 
                Power(t2,2)*
                 (-28 + 18*t2 + 2*Power(s2,2)*t2 - 8*Power(t2,2) + 
                   s2*(-31 + 27*t2)) + 
                Power(s1,2)*
                 (-10 + 49*t2 - 33*Power(t2,2) - 10*Power(t2,3) + 
                   s2*(2 - 15*t2 + 11*Power(t2,2)) + 
                   t1*(-3 - 4*t2 + 15*Power(t2,2))) + 
                s1*(2 + Power(t1,2)*(6 - 4*t2) + 
                   (39 + 20*s2 - 2*Power(s2,2))*t2 - 
                   (121 + 12*s2)*Power(t2,2) - 
                   2*(-24 + s2)*Power(t2,3) + 
                   t1*(4 - 4*t2 + 12*Power(t2,2) - 22*Power(t2,3) + 
                      s2*(-4 + t2 + 9*Power(t2,2))))) + 
             Power(s,2)*(1 + 15*t2 - 9*s2*t2 + 3*Power(t2,2) - 
                31*s2*Power(t2,2) - 2*Power(s2,2)*Power(t2,2) - 
                27*Power(t2,3) + 29*s2*Power(t2,3) + 
                4*Power(s2,2)*Power(t2,3) + 6*Power(t2,4) - 
                17*s2*Power(t2,4) + 2*Power(t2,5) + 
                2*Power(s1,4)*(3 - 8*t2 + 5*Power(t2,2)) + 
                Power(t1,2)*
                 (-7 + 15*t2 + 5*Power(t2,2) - 2*Power(t2,3) + 
                   Power(t2,4)) + 
                t1*(6 - 9*t2 - 8*(4 + s2)*Power(t2,2) + 
                   (68 - 7*s2)*Power(t2,3) + (-6 + s2)*Power(t2,4) + 
                   Power(t2,5)) + 
                Power(s1,3)*
                 (-20 + 45*t2 - 16*Power(t2,2) - 9*Power(t2,3) + 
                   t1*(-2 - 3*t2 + 7*Power(t2,2)) + 
                   s2*(9 - 21*t2 + 10*Power(t2,2))) + 
                Power(s1,2)*
                 (-7 + Power(t1,2)*(11 - 7*t2) + 108*t2 - 
                   145*Power(t2,2) + 46*Power(t2,3) - 
                   2*Power(t2,4) + 
                   Power(s2,2)*(3 - 6*t2 + Power(t2,2)) + 
                   t1*(13 - 18*t2 + 15*Power(t2,2) - 
                      14*Power(t2,3)) + 
                   s2*(-16 + 19*t2 + 12*Power(t2,2) - 
                      11*Power(t2,3) + 
                      t1*(-11 + 3*t2 + 6*Power(t2,2)))) + 
                s1*(-4 + t2 - Power(s2,2)*Power(-1 + t2,2)*t2 - 
                   58*Power(t2,2) + 70*Power(t2,3) - 10*Power(t2,4) + 
                   Power(t2,5) + 
                   Power(t1,2)*
                    (-11 - 8*t2 + Power(t2,2) + 2*Power(t2,3)) + 
                   t1*(-6 + 14*t2 - 43*Power(t2,2) + 9*Power(t2,3)) + 
                   s2*(2 + 64*t2 - 74*Power(t2,2) + 33*Power(t2,3) + 
                      Power(t2,4) + 
                      t1*(7 + 2*t2 + 21*Power(t2,2) - 14*Power(t2,3)))\
)) + s*(4 - 8*t1 + 4*Power(t1,2) - 2*Power(s1,5)*Power(-1 + t2,2) - 
                12*t2 - s2*t2 + 37*t1*t2 + s2*t1*t2 - 
                25*Power(t1,2)*t2 + 10*Power(t2,2) - 
                33*s2*Power(t2,2) + 7*Power(s2,2)*Power(t2,2) + 
                12*t1*Power(t2,2) + 15*s2*t1*Power(t2,2) + 
                17*Power(t1,2)*Power(t2,2) + 2*Power(t2,3) + 
                29*s2*Power(t2,3) - 4*Power(s2,2)*Power(t2,3) - 
                64*t1*Power(t2,3) - 32*s2*t1*Power(t2,3) + 
                13*Power(t1,2)*Power(t2,3) - 8*Power(t2,4) + 
                13*s2*Power(t2,4) + 3*Power(s2,2)*Power(t2,4) + 
                16*t1*Power(t2,4) + 5*s2*t1*Power(t2,4) - 
                3*Power(t1,2)*Power(t2,4) + 6*Power(t2,5) - 
                8*s2*Power(t2,5) + 7*t1*Power(t2,5) - 
                s2*t1*Power(t2,5) - 2*Power(t2,6) - 
                Power(s1,4)*(-1 + t2)*
                 (8 + t1 - 6*t2 + t1*t2 - 2*Power(t2,2) + 
                   s2*(-5 + 3*t2)) - 
                Power(s1,3)*(Power(t1,2)*(8 - 6*t2) + 
                   Power(s2,2)*(7 - 6*t2 + Power(t2,2)) - 
                   Power(-1 + t2,2)*(23 - 13*t2 + 2*Power(t2,2)) + 
                   t1*(13 - 14*t2 + 3*Power(t2,2) - 2*Power(t2,3)) + 
                   s2*(-12 + 5*t2 + 12*Power(t2,2) - 5*Power(t2,3) + 
                      t1*(-13 + 8*t2 + Power(t2,2)))) + 
                Power(s1,2)*(Power(t1,2)*
                    (17 - 7*t2 + Power(t2,2) - Power(t2,3)) - 
                   2*Power(-1 + t2,2)*
                    (4 + 15*t2 + Power(t2,2) + Power(t2,3)) + 
                   Power(s2,2)*
                    (7 + 10*t2 - 9*Power(t2,2) + 2*Power(t2,3)) + 
                   t1*(30 - 32*t2 + 5*Power(t2,2) - 4*Power(t2,3) + 
                      Power(t2,4)) + 
                   s2*(11 - 65*t2 + 66*Power(t2,2) - 11*Power(t2,3) - 
                      Power(t2,4) + 
                      t1*(-29 + 12*t2 - 7*Power(t2,2) + 4*Power(t2,3))\
)) + s1*(Power(s2,2)*t2*(-14 + t2 - Power(t2,3)) + 
                   Power(-1 + t2,2)*
                    (2 + 14*t2 + 3*Power(t2,2) + 9*Power(t2,3)) - 
                   Power(t1,2)*
                    (-15 + 22*t2 + 10*Power(t2,2) - 4*Power(t2,3) + 
                      Power(t2,4)) - 
                   t1*(17 + 66*t2 - 133*Power(t2,2) + 
                      52*Power(t2,3) - 4*Power(t2,4) + 2*Power(t2,5)) \
+ s2*(1 + 22*t2 + 24*Power(t2,2) - 69*Power(t2,3) + 23*Power(t2,4) - 
                      Power(t2,5) + 
                      t1*(-1 + 14*t2 + 7*Power(t2,2) + 
                        10*Power(t2,3) - 2*Power(t2,4))))))/
           (s*(-1 + s2)*(s - s2 + t1)*
             Power(s1 - s1*t2 + t2*(-1 + s + t2),2)) - 
          (2*(-1 + t2)*Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2)*
              (2*Power(s1,2) - 3*s1*(1 + t2) + t2*(3 + t2)) + 
             Power(s,5)*t2*(t1 - t1*t2 + t2*(3 + 2*s1 + t2)) + 
             Power(s,4)*(2*Power(t1,2)*(-1 + t2) + 
                2*Power(s1,2)*(2 - 5*t2)*t2 - 
                2*t1*t2*(s2 + t2 - Power(t2,2)) + 
                t2*(1 - (-15 + s2)*t2 + (10 - 3*s2)*Power(t2,2)) + 
                s1*(t1*(1 - 8*t2 + Power(t2,2)) + 
                   t2*(5 + s2 - 9*t2 - 3*s2*t2 + 2*Power(t2,2)))) + 
             Power(s,3)*(2*Power(s1,3)*(1 - 8*t2 + 9*Power(t2,2)) + 
                Power(t1,2)*(1 - 8*t2 + Power(t2,2) + 2*Power(t2,3)) + 
                t1*(-4 + (-9 + 2*s2)*t2 + (10 + s2)*Power(t2,2) - 
                   (7 + 5*s2)*Power(t2,3) + 2*Power(t2,4)) + 
                Power(s1,2)*(2 - 7*t2 + 16*Power(t2,2) - 
                   17*Power(t2,3) + t1*(-6 + 13*t2 + 3*Power(t2,2)) + 
                   s2*(1 - 11*t2 + 10*Power(t2,2))) - 
                t2*(3 + 35*t2 + 10*Power(t2,2) - 
                   2*Power(s2,2)*Power(t2,2) - 12*Power(t2,3) + 
                   s2*(-2 + 21*t2 - 8*Power(t2,2) + Power(t2,3))) + 
                s1*(1 + Power(t1,2)*(12 - 8*t2) - 
                   2*(-17 + s2 + Power(s2,2))*t2 + 
                   (-43 + 5*s2)*Power(t2,2) - 3*(7 + s2)*Power(t2,3) + 
                   3*Power(t2,4) + 
                   t1*(5 - 14*t2 + 10*Power(t2,2) - 11*Power(t2,3) + 
                      s2*(-6 + 11*t2 + Power(t2,2))))) + 
             Power(s,2)*(-2 - 12*t2 + 6*s2*t2 + 55*s2*Power(t2,2) + 
                2*Power(s2,2)*Power(t2,2) + 42*Power(t2,3) - 
                37*s2*Power(t2,3) - Power(s2,2)*Power(t2,3) - 
                30*Power(t2,4) + 10*s2*Power(t2,4) + 
                Power(s2,2)*Power(t2,4) + 2*Power(t2,5) - 
                2*Power(s1,4)*(3 - 10*t2 + 7*Power(t2,2)) + 
                Power(t1,2)*(4 - 6*t2 - 3*Power(t2,2) - 
                   4*Power(t2,3) + Power(t2,4)) + 
                t1*(-2 + (3 + 2*s2)*t2 - 5*(3 + s2)*Power(t2,2) + 
                   (-25 + 13*s2)*Power(t2,3) - 
                   4*(-1 + s2)*Power(t2,4) + Power(t2,5)) - 
                Power(s1,3)*(-2 + 4*t2 + 22*Power(t2,2) - 
                   24*Power(t2,3) + t1*(-7 + 6*t2 + 5*Power(t2,2)) + 
                   s2*(7 - 22*t2 + 11*Power(t2,2))) + 
                Power(s1,2)*(17 - 83*t2 + 70*Power(t2,2) + 
                   7*Power(t2,3) - 11*Power(t2,4) + 
                   2*Power(t1,2)*(-11 + 7*t2) - 
                   Power(s2,2)*(4 - 7*t2 + Power(t2,2)) + 
                   t1*(-27 + 31*t2 - 14*Power(t2,2) + 
                      16*Power(t2,3)) + 
                   s2*(4 + 6*t2 - 31*Power(t2,2) + 15*Power(t2,3) + 
                      t1*(25 - 17*t2 - 2*Power(t2,2)))) + 
                s1*(2 - 18*t2 + 2*Power(s2,2)*(1 - 3*t2)*t2 + 
                   44*Power(t2,2) - 19*Power(t2,3) - 10*Power(t2,4) + 
                   Power(t2,5) + 
                   Power(t1,2)*
                    (3 + 16*t2 + Power(t2,2) - 4*Power(t2,3)) + 
                   t1*(10 + 52*t2 - 29*Power(t2,2) + 8*Power(t2,3) - 
                      9*Power(t2,4)) - 
                   s2*(2 + 71*t2 - 50*Power(t2,2) + 5*Power(t2,3) + 
                      4*Power(t2,4) + 
                      t1*(6 + 8*t2 + 8*Power(t2,2) - 10*Power(t2,3))))) \
+ s*(-3 + 6*t1 - 3*Power(t1,2) + 4*Power(s1,5)*Power(-1 + t2,2) + 
                9*t2 + 2*s2*t2 - 27*t1*t2 - 2*s2*t1*t2 + 
                18*Power(t1,2)*t2 - Power(t2,2) + 7*s2*Power(t2,2) - 
                7*Power(s2,2)*Power(t2,2) - 9*t1*Power(t2,2) + 
                2*s2*t1*Power(t2,2) - 12*Power(t1,2)*Power(t2,2) - 
                35*Power(t2,3) + 28*s2*Power(t2,3) - 
                2*Power(s2,2)*Power(t2,3) + 36*t1*Power(t2,3) + 
                7*s2*t1*Power(t2,3) - 8*Power(t1,2)*Power(t2,3) + 
                50*Power(t2,4) - 39*s2*Power(t2,4) + 
                Power(s2,2)*Power(t2,4) - 13*t1*Power(t2,4) + 
                10*s2*t1*Power(t2,4) - 3*Power(t1,2)*Power(t2,4) - 
                18*Power(t2,5) + 2*s2*Power(t2,5) + 7*t1*Power(t2,5) - 
                s2*t1*Power(t2,5) - 2*Power(t2,6) + 
                2*Power(s1,4)*(-1 + t2)*
                 (4 + t1 + 2*s2*(-2 + t2) + t2 + t1*t2 - 5*Power(t2,2)) \
+ Power(s1,3)*(-4*Power(t1,2)*(-4 + 3*t2) + 
                   Power(s2,2)*(11 - 8*t2 + Power(t2,2)) + 
                   4*Power(-1 + t2,2)*(-5 + 5*t2 + 2*Power(t2,2)) + 
                   t1*(18 - 17*t2 + 6*Power(t2,2) - 7*Power(t2,3)) + 
                   s2*(-12 - 11*t2 + 32*Power(t2,2) - 9*Power(t2,3) + 
                      t1*(-25 + 16*t2 + Power(t2,2)))) - 
                Power(s1,2)*(-2*Power(t1,2)*
                    (-11 + 2*Power(t2,2) + Power(t2,3)) + 
                   Power(s2,2)*
                    (7 + 24*t2 - 17*Power(t2,2) + 2*Power(t2,3)) + 
                   Power(-1 + t2,2)*
                    (-6 - 11*t2 + 17*Power(t2,2) + 2*Power(t2,3)) + 
                   t1*(6 + 34*t2 - 31*Power(t2,2) - 2*Power(t2,3) - 
                      7*Power(t2,4)) + 
                   s2*(25 - 100*t2 + 57*Power(t2,2) + 24*Power(t2,3) - 
                      6*Power(t2,4) + 
                      t1*(-34 - 9*t2 + 6*Power(t2,2) + 5*Power(t2,3)))) \
+ s1*(Power(s2,2)*t2*(14 + 15*t2 - 10*Power(t2,2) + Power(t2,3)) + 
                   Power(-1 + t2,2)*
                    (-2 - 17*t2 + 32*Power(t2,2) + 7*Power(t2,3)) - 
                   Power(t1,2)*
                    (11 - 26*t2 - 6*Power(t2,3) + Power(t2,4)) + 
                   t1*(13 + 31*t2 - 36*Power(t2,2) + 11*Power(t2,3) - 
                      17*Power(t2,4) - 2*Power(t2,5)) + 
                   s2*(-2 + 18*t2 - 116*Power(t2,2) + 99*Power(t2,3) + 
                      2*Power(t2,4) - Power(t2,5) + 
                      t1*(2 - 36*t2 + 9*Power(t2,2) - 20*Power(t2,3) + 
                         5*Power(t2,4))))))/
           (s*(-1 + s2)*(-1 + t1)*Power(s1 - s1*t2 + t2*(-1 + s + t2),2))\
) + ((64*(((-1 + s1)*(-2 - 2*s*(-1 + t2) + t2 + s1*(-1 + 2*t2))*
                  (2*Power(s1,2)*(-1 + t2) + 
                    t1*(-1 + s + 3*t2 - s*t2) + 
                    (1 + t2)*(1 + (-1 + s - s2)*t2) + 
                    s1*(t1*(-3 + t2) + s2*(1 + t2) - 
                       2*t2*(-1 + s + t2))))/
                (s1*(-1 + t2) - t2*(-1 + s + t2)) - 
               2*(-20 + 2*s2 + 8*t1 + 4*Power(s,2)*(-2 + t2) - 20*t2 - 
                  3*s2*t2 + 2*t1*t2 - 6*Power(t2,2) + 
                  Power(s1,3)*(1 + 2*t2) + 
                  s1*(18 + 11*t2 + 2*Power(t2,2) + s2*(1 + t2) - 
                     2*t1*(5 + t2)) + 
                  Power(s1,2)*(-11 + 2*t1 - 9*t2 + s2*(-3 + 2*t2)) + 
                  s*(24 - t1 + 2*s2*(-2 + t2) - 
                     2*Power(s1,2)*(-1 + t2) - 11*t2 + 
                     s1*(6 + t1 - 2*s2*(-2 + t2) + 5*t2)))))/
           (-s + s1 - t2) + ((-64*(-1 + s1)*
                (1 + s + t2 + 2*s*t2 - 2*s1*t2)*
                (2*Power(s1,2)*(-1 + t2) + 
                  t1*(-1 + s + 3*t2 - s*t2) + 
                  (1 + t2)*(1 + (-1 + s - s2)*t2) + 
                  s1*(t1*(-3 + t2) + s2*(1 + t2) - 
                     2*t2*(-1 + s + t2))))/
              (s1*(-1 + t2) - t2*(-1 + s + t2)) - 
             128*(-21 + 3*s2 + 3*t1 - 17*t2 + 4*Power(s,2)*t2 + 
                2*Power(s1,3)*t2 + 2*s2*t2 - 2*t1*t2 + 6*Power(t2,2) + 
                Power(s1,2)*(-9 + s2 - 3*t1 - 5*t2 + 2*s2*t2) - 
                2*s1*(-5 - (2 + t1)*t2 + Power(t2,2) + 
                   2*s2*(1 + t2)) - 
                s*(-5 + 2*t1 + 2*t2 - 2*s2*t2 + 
                   Power(s1,2)*(1 + 2*t2) - 2*s1*(4 + t1 - s2*t2))))/
           (-1 + t2) + (64*(-1 + s1)*
             (-(((8 + Power(s1,2) + 
                      s*(s1*(-5 + t2) - Power(-3 + t2,2)) - 
                      Power(s,2)*(-4 + t2) - Power(t2,2) + 
                      s1*(3 - 4*t2 + Power(t2,2)))*
                    (2*Power(s1,2)*(-1 + t2) + 
                      t1*(-1 + s + 3*t2 - s*t2) + 
                      (1 + t2)*(1 + (-1 + s - s2)*t2) + 
                      s1*(t1*(-3 + t2) + s2*(1 + t2) - 
                        2*t2*(-1 + s + t2))))/
                  (s1*(-1 + t2) - t2*(-1 + s + t2))) - 
               2*(6 - Power(s1,3) + 2*s2 - 22*t1 - 23*t2 + 18*s2*t2 - 
                  t1*t2 - 5*Power(t2,2) + 
                  Power(s1,2)*
                   (-9 + 2*s2 - 3*t1 + 4*t2 - Power(t2,2)) + 
                  Power(s,2)*
                   (-2 + t1*(-6 + t2) + s1*(-4 + t2) + 7*t2 - 
                     Power(t2,2)) + 
                  s1*(-5 + t1 + 9*t2 + t1*t2 + Power(t2,2) - 
                     s2*(9 - 4*t2 + Power(t2,2))) + 
                  s*(-4 + 19*t1 - Power(s1,2)*(-3 + t2) - 6*t2 - 
                     3*t1*t2 + 2*Power(t2,2) + 
                     s2*(2 - 9*t2 + Power(t2,2)) + 
                     s1*(26 - 2*s2 + 7*t1 - 12*t2 - t1*t2 + 
                        2*Power(t2,2))))))/((s - s2 + t1)*(s - s1 + t2)) \
+ ((-1 + s1)*((-64*(5 + (5 + s)*t2 + Power(s1,2)*t2 + 
                    (1 + s)*Power(t2,2) - s1*(1 + t2)*(1 + s + t2))*
                  (2*Power(s1,2)*(-1 + t2) + 
                    t1*(-1 + s + 3*t2 - s*t2) + 
                    (1 + t2)*(1 + (-1 + s - s2)*t2) + 
                    s1*(t1*(-3 + t2) + s2*(1 + t2) - 
                       2*t2*(-1 + s + t2))))/
                (s1*(-1 + t2) - t2*(-1 + s + t2)) + 
               128*(16 - s2 + 8*t2 + Power(s1,3)*t2 + 3*s2*t2 - 
                  2*t1*t2 - 6*Power(t2,2) + s2*Power(t2,2) + 
                  Power(s,2)*(1 + (-3 + t1)*t2 - Power(t2,2)) - 
                  Power(s1,2)*
                   (-1 + 2*t1 - (-2 + s2)*t2 + Power(t2,2)) + 
                  s1*(5 + 13*t2 + Power(t2,2) + 2*t1*(1 + t2) - 
                     s2*(3 + 4*t2 + Power(t2,2))) + 
                  s*(-9 - t2 - 2*t1*t2 + Power(t2,2) - 
                     Power(s1,2)*(1 + t2) + 
                     s2*(1 + t2 + Power(t2,2)) + 
                     s1*(-4 + t1 + t2 - s2*t2 - t1*t2 + 2*Power(t2,2)))\
)))/((-1 + t1)*(-1 + t2)) + (64*(-1 + s1)*
             (-(((11 + 2*Power(s1,3) + s2 - t1 + 
                      Power(s,2)*(-2 + t1 - t2) + 8*t2 + 2*s2*t2 + 
                      t1*t2 + 
                      s*(-4 - 2*Power(s1,2) - 
                       s1*(-7 + s2 + 2*t1 - 3*t2) - 3*t2 + t1*t2 + 
                       s2*(2 + t2)) + 
                      Power(s1,2)*(s2 + t1 - 2*(4 + t2)) - 
                      s1*(-10 - 7*t2 + t1*(2 + t2) + s2*(5 + t2)))*
                    (t1*(1 + s*(-1 + t2) - 3*t2) - 
                      2*Power(s1,2)*(-1 + t2) - 
                      (1 + t2)*(1 + (-1 + s - s2)*t2) - 
                      s1*(t1*(-3 + t2) + s2*(1 + t2) - 
                        2*t2*(-1 + s + t2))))/
                  (s1 - s1*t2 + t2*(-1 + s + t2))) + 
               2*(-5 + Power(s1,4) + 6*s2 + Power(s2,2) - 4*t1 - 
                  2*s2*t1 + 3*Power(t1,2) + 
                  Power(s1,3)*(-4 + s2 + t1 - t2) + 23*t2 - 4*s2*t2 + 
                  Power(s2,2)*t2 - 9*t1*t2 + 4*Power(t2,2) + 
                  Power(s1,2)*
                   (8 + s2*(-6 + t1 - t2) + 3*t2 - t1*(3 + t2)) + 
                  Power(s,2)*
                   (6 + Power(t1,2) + s1*(-2 + t1 - t2) + 5*t2 - 
                     t1*(5 + t2)) - 
                  s*(15 + Power(s1,3) - 9*t1 + 3*Power(t1,2) + 
                     Power(s1,2)*(-5 + s2 + 2*t1 - 2*t2) + 17*t2 - 
                     2*t1*t2 + 
                     s1*(5 - 5*t1 + Power(t1,2) + 
                       s2*(-3 + t1 - t2) + 5*t2 - 2*t1*t2) - 
                     s2*(-5 + t1*(2 + t2))) + 
                  s1*(8 - 2*Power(s2,2) + Power(t1,2) + 3*t2 + 
                     t1*(3 + t2) + s2*(-(t1*(2 + t2)) + 2*(3 + t2)))) \
- ((-Power(s1,2) + s*(-2 + s1 - t2) - 4*(1 + t2) + s1*(4 + t2))*
                  (2*(s1 - t2)*(-1 + t2)*
                     Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                    Power(s,3)*
                     (Power(t1,2)*Power(-1 + t2,2) + 
                       2*t1*(s1 - t2)*(-1 + t2)*t2 + 
                       Power(t2,2)*
                        (1 + 2*Power(s1,2) + Power(t2,2) - 
                        2*s1*(1 + t2))) - 
                    2*Power(s,2)*
                     (Power(t1,2)*(-1 + t2)*
                       (-1 + s1*(-3 + t2) + 3*t2) + 
                       t1*(-1 + t2 + (3 + s2)*Power(t2,2) - 
                       (5 + s2)*Power(t2,3) + 
                       Power(s1,2)*(1 - 5*t2 + 2*Power(t2,2)) + 
                       s1*(1 + s2*(-1 + t2) - 5*t2 + 
                       10*Power(t2,2) - 2*Power(t2,3))) + 
                       t2*(-1 + 2*Power(s1,3)*(-1 + t2) + t2 - 
                        2*Power(t2,2) + 2*Power(t2,3) + 
                        Power(s1,2)*
                       (1 + s2 + 2*t2 + s2*t2 - 3*Power(t2,2)) + 
                        s2*(1 - t2 + Power(t2,2) + Power(t2,3)) - 
                        s1*(-1 + 2*Power(t2,2) - Power(t2,3) + 
                        s2*(1 + t2 + 2*Power(t2,2))))) + 
                    s*(1 - 2*t1 + Power(t1,2) + 
                       2*Power(s1,4)*Power(-1 + t2,2) - 2*t2 + 
                       10*t1*t2 - 8*Power(t1,2)*t2 + 6*Power(t2,2) - 
                       12*s2*Power(t2,2) + Power(s2,2)*Power(t2,2) + 
                       2*t1*Power(t2,2) + 4*s2*t1*Power(t2,2) + 
                       9*Power(t1,2)*Power(t2,2) - 10*Power(t2,3) + 
                       8*s2*Power(t2,3) - 10*t1*Power(t2,3) - 
                       8*s2*t1*Power(t2,3) + 5*Power(t2,4) + 
                       4*s2*Power(t2,4) + Power(s2,2)*Power(t2,4) + 
                       2*Power(s1,3)*(-1 + t2)*
                        (t1*(-3 + t2) - 2*(-1 + t2)*t2 + s2*(1 + t2)) \
+ Power(s1,2)*(Power(s2,2)*(1 + Power(t2,2)) + 
                        Power(t1,2)*(9 - 8*t2 + Power(t2,2)) + 
                        2*Power(-1 + t2,2)*(-1 - t2 + Power(t2,2)) + 
                        4*s2*
                       (t1*(-2 + t2) - (-1 + t2)*Power(t2,2)) - 
                        2*t1*
                        (-5 + 13*t2 - 9*Power(t2,2) + Power(t2,3))) + 
                       2*s1*(Power(t1,2)*(3 - 8*t2 + 3*Power(t2,2)) + 
                        t1*(-3 + 2*(-4 + s2)*t2 + 
                        (17 + 2*s2)*Power(t2,2) - 6*Power(t2,3)) + 
                        t2*(Power(-1 + t2,3) - 
                        Power(s2,2)*(1 + Power(t2,2)) + 
                        s2*(6 - 3*t2 - 4*Power(t2,2) + Power(t2,3)))))\
))/(s*Power(s1 - s1*t2 + t2*(-1 + s + t2),2))))/((-1 + s2)*(-1 + t1)) + 
          (32*(-1 + s1)*((-2*
                  (8 - s2 + 10*t1 + 
                    s1*(19 + t1*(-4 + t2) + s2*(-2 + t2) - 5*t2) + 
                    2*Power(s1,2)*(-2 + t2) + 4*t2 - 4*s2*t2 - 
                    t1*t2 + Power(s,2)*(2 - 2*s1 - 2*t1 + t2) + 
                    s*(-12 + 2*Power(s1,2) + s2 + 3*t1 + 
                       s1*(2 + 2*t1 - 3*t2) + 3*t2 - s2*t2 - t1*t2))*
                  (2*Power(s1,2)*(-1 + t2) + 
                    t1*(-1 + s + 3*t2 - s*t2) + 
                    (1 + t2)*(1 + (-1 + s - s2)*t2) + 
                    s1*(t1*(-3 + t2) + s2*(1 + t2) - 
                       2*t2*(-1 + s + t2))))/
                (s1*(-1 + t2) - t2*(-1 + s + t2)) - 
               4*(4 + 4*s2 - 2*Power(s2,2) - 11*t1 - 2*s2*t1 + 
                  3*Power(t1,2) - 
                  Power(s1,2)*
                   (10 + s2*(-3 + t2) + t1*(-3 + t2) - 3*t2) - 
                  Power(s1,3)*(-2 + t2) - 4*t2 + 6*s2*t2 + 
                  2*Power(s2,2)*t2 - 9*t1*t2 + 4*Power(t2,2) + 
                  Power(s,2)*
                   (-3 + Power(s1,2) + Power(t1,2) + 
                     s1*(-2 + 2*t1 - t2) + 5*t2 - t1*(4 + t2)) + 
                  s1*(-13 + Power(t1,2) + t1*(-7 + t2) - 9*t2 + 
                     s2*(-6 + 3*t1 + 6*t2 - t1*t2)) - 
                  s*(Power(s1,3) - 18*t1 + 3*Power(t1,2) + 
                     Power(s1,2)*(1 + 2*t1 - 2*t2) + 6*t2 - 2*t1*t2 + 
                     s1*(-23 + 2*t1 + Power(t1,2) - s2*(-3 + t2) + 
                        5*t2 - 2*t1*t2) + s2*(-3 + t1 + 7*t2 - t1*t2))) \
- (2*(-7 + Power(s,2) - s1*(-2 + t2) + 2*t2 + s*(-1 - s1 + t2))*
                  (2*(s1 - t2)*(-1 + t2)*
                     Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                    Power(s,3)*
                     (Power(t1,2)*Power(-1 + t2,2) + 
                       2*t1*(s1 - t2)*(-1 + t2)*t2 + 
                       Power(t2,2)*
                        (1 + 2*Power(s1,2) + Power(t2,2) - 
                        2*s1*(1 + t2))) - 
                    2*Power(s,2)*
                     (Power(t1,2)*(-1 + t2)*
                        (-1 + s1*(-3 + t2) + 3*t2) + 
                       t1*(-1 + t2 + (3 + s2)*Power(t2,2) - 
                        (5 + s2)*Power(t2,3) + 
                        Power(s1,2)*(1 - 5*t2 + 2*Power(t2,2)) + 
                        s1*(1 + s2*(-1 + t2) - 5*t2 + 
                        10*Power(t2,2) - 2*Power(t2,3))) + 
                       t2*(-1 + 2*Power(s1,3)*(-1 + t2) + t2 - 
                        2*Power(t2,2) + 2*Power(t2,3) + 
                        Power(s1,2)*
                        (1 + s2 + 2*t2 + s2*t2 - 3*Power(t2,2)) + 
                        s2*(1 - t2 + Power(t2,2) + Power(t2,3)) - 
                        s1*(-1 + 2*Power(t2,2) - Power(t2,3) + 
                        s2*(1 + t2 + 2*Power(t2,2))))) + 
                    s*(1 - 2*t1 + Power(t1,2) + 
                       2*Power(s1,4)*Power(-1 + t2,2) - 2*t2 + 
                       10*t1*t2 - 8*Power(t1,2)*t2 + 6*Power(t2,2) - 
                       12*s2*Power(t2,2) + Power(s2,2)*Power(t2,2) + 
                       2*t1*Power(t2,2) + 4*s2*t1*Power(t2,2) + 
                       9*Power(t1,2)*Power(t2,2) - 10*Power(t2,3) + 
                       8*s2*Power(t2,3) - 10*t1*Power(t2,3) - 
                       8*s2*t1*Power(t2,3) + 5*Power(t2,4) + 
                       4*s2*Power(t2,4) + Power(s2,2)*Power(t2,4) + 
                       2*Power(s1,3)*(-1 + t2)*
                        (t1*(-3 + t2) - 2*(-1 + t2)*t2 + s2*(1 + t2)) \
+ Power(s1,2)*(Power(s2,2)*(1 + Power(t2,2)) + 
                        Power(t1,2)*(9 - 8*t2 + Power(t2,2)) + 
                        2*Power(-1 + t2,2)*(-1 - t2 + Power(t2,2)) + 
                        4*s2*(t1*(-2 + t2) - (-1 + t2)*Power(t2,2)) - 
                        2*t1*
                        (-5 + 13*t2 - 9*Power(t2,2) + Power(t2,3))) + 
                       2*s1*(Power(t1,2)*(3 - 8*t2 + 3*Power(t2,2)) + 
                        t1*(-3 + 2*(-4 + s2)*t2 + 
                        (17 + 2*s2)*Power(t2,2) - 6*Power(t2,3)) + 
                        t2*(Power(-1 + t2,3) - 
                        Power(s2,2)*(1 + Power(t2,2)) + 
                        s2*(6 - 3*t2 - 4*Power(t2,2) + Power(t2,3))))))\
)/(s*Power(s1 - s1*t2 + t2*(-1 + s + t2),2))))/((-1 + s2)*(-s + s2 - t1))\
)/Power(-1 + s1,2))*B1(t2,s,s1))/(128.*Power(Pi,2)) + 
  (((4*(-1 + s2 - t1 + t2)*(-2*Power(t1,2) + s1*(s2 + t1)*(-1 + t2) + 
            (1 + t2)*(1 + (-1 + s - s2)*t2) + 
            t1*(1 + t2 + 2*s2*t2 - s*(1 + t2)))*
          ((-2 + 3*Power(s2,2) - Power(s1,3)*(-3 + t1) + t1 + 
               3*s2*t1 + 3*Power(s2,2)*t1 + 4*Power(t1,2) + 
               s2*Power(t1,2) - Power(s2,2)*Power(t1,2) + 
               Power(t1,3) - 2*t2 + 4*s2*t2 - 4*Power(s2,2)*t2 + 
               t1*t2 + 6*s2*t1*t2 - Power(t1,2)*t2 - 
               s2*Power(t1,2)*t2 - 5*s2*Power(t2,2) + 
               s2*t1*Power(t2,2) - 2*Power(t2,3) - 
               Power(s,2)*(1 + s1*(-1 + t1) + t1 + Power(t1,2) + 
                  t1*t2) + Power(s1,2)*
                (-3 - Power(t1,2) - 2*t1*(-2 + t2) + 2*t2 + 
                  Power(t2,2) + s2*(6 - 2*t1 + t2)) + 
               s*(2*Power(s1,2)*(-2 + t1) - t1 - Power(t1,2) - 3*t2 - 
                  6*t1*t2 + Power(t1,2)*t2 + 2*Power(t2,2) - 
                  t1*Power(t2,2) + 
                  s2*(-1 + 2*Power(t1,2) + t1*(-2 + t2) + 3*t2) + 
                  s1*(3 - 2*t1 + 2*Power(t1,2) + 
                     s2*(-4 + 2*t1 - t2) - 2*t2 + 3*t1*t2 - 
                     Power(t2,2))) + 
               s1*(-1 + s2*(-2*Power(t1,2) - 2*t1*(-3 + t2) + 
                     Power(-1 + t2,2)) - Power(t1,2)*(-3 + t2) + t2 - 
                  3*Power(t2,2) + Power(s2,2)*(3 - t1 + t2) + 
                  t1*(4 + 5*t2 + Power(t2,2))))/((-1 + s1)*(-1 + t2)) + 
            (-9 - 2*s2 + 14*Power(s2,2) - Power(s2,3) + 6*t1 - 
               7*s2*t1 - 3*Power(s2,2)*t1 - 6*Power(t1,2) + 
               s2*Power(t1,2) + Power(t1,3) - t2 + 11*s2*t2 + 
               5*Power(s2,2)*t2 - Power(s2,3)*t2 + 7*t1*t2 - 
               2*s2*t1*t2 - 5*Power(t1,2)*t2 + 3*Power(t2,2) + 
               6*s2*Power(t2,2) - Power(s2,2)*Power(t2,2) + 
               3*t1*Power(t2,2) + Power(t2,3) + 
               Power(s,2)*(4 + s2 - 2*t1 + 3*t2 - s2*t2 - 
                  Power(t2,2)) - 
               Power(s1,2)*(-1 - t1 + 2*t2 + Power(t2,2) + 
                  s2*(2 + t2)) + 
               s1*(2 + 3*Power(t1,2) + 12*t2 + 5*Power(t2,2) - 
                  5*t1*(2 + t2) - Power(s2,2)*(3 + 2*t2) + 
                  s2*(16 - 3*t1 + 3*t2 - 2*Power(t2,2))) + 
               s*(3 + 5*t1 - Power(t1,2) - 11*t2 + 2*Power(s2,2)*t2 + 
                  2*t1*t2 - 5*Power(t2,2) + 
                  s2*(-18 + 5*t1 - 8*t2 + 2*Power(t2,2)) + 
                  s1*(-6 + s2 + 2*t1 - t2 + 2*s2*t2 + 2*Power(t2,2))))/
             ((-1 + s2)*(-1 + t1)) - 
            (6 + 5*s2 - 6*Power(s2,2) - 2*Power(s2,3) + 
               2*Power(s1,3)*(-3 + t1) + 2*t1 - s2*t1 + 
               Power(s2,2)*t1 + Power(s2,3)*t1 - 4*Power(t1,2) - 
               5*s2*Power(t1,2) - 4*Power(t1,3) + s2*Power(t1,3) + 
               t2 - 4*s2*t2 - 4*Power(s2,2)*t2 + t1*t2 + 6*s2*t1*t2 + 
               Power(s2,2)*t1*t2 + 2*Power(t1,2)*t2 - 
               s2*Power(t1,2)*t2 - 4*Power(t2,2) - 4*s2*Power(t2,2) + 
               t1*Power(t2,2) - Power(t2,3) + 
               Power(s,2)*(s1 + t1)*(-5 + s2 + t1 + t2) + 
               Power(s1,2)*(Power(s2,2) + 3*Power(t1,2) + 
                  3*(2 + t2) - t1*(12 + t2) + s2*(-7 + 2*t1 + t2)) + 
               s1*(4 + Power(s2,3) + Power(t1,3) - 3*t2 - 
                  Power(t1,2)*(11 + t2) + 
                  Power(s2,2)*(-3 + t1 + t2) + 
                  s2*(-5 - 9*t1 + 3*Power(t1,2) + 2*t2) + 
                  t1*(2 + 5*t2)) - 
               s*(4 + 2*Power(s2,2)*(-1 + t1) - 5*Power(t1,2) + 
                  Power(t1,3) - 4*t2 + 3*t1*t2 - Power(t1,2)*t2 - 
                  Power(t2,2) + Power(s1,2)*(-11 + s2 + 3*t1 + t2) + 
                  s2*(Power(t1,2) + 2*t1*(-2 + t2) - 2*(4 + t2)) + 
                  s1*(4 + 2*Power(s2,2) - 15*t1 + 4*Power(t1,2) + 
                     3*t2 + 2*s2*(-3 + t1 + t2))))/
             ((-1 + s1)*(-s + s1 - t2)) + 
            (2 - 7*s2 + 2*Power(s2,2) - 2*Power(s2,3) - 10*t1 - 
               2*s2*t1 + 9*Power(s2,2)*t1 - Power(s2,3)*t1 + 
               8*Power(t1,2) + 4*s2*Power(t1,2) - 
               2*Power(s2,2)*Power(t1,2) + s2*Power(t1,3) - 4*t2 + 
               7*s2*t2 - 2*Power(s2,2)*t2 + 2*t1*t2 + 2*s2*t1*t2 - 
               Power(s2,2)*t1*t2 + 2*Power(t1,2)*t2 - 
               s2*Power(t1,2)*t2 - s2*Power(t2,2) - 
               2*Power(s1,3)*(-3 + s2 + t2) + 
               Power(s1,2)*(-1 - 3*Power(s2,2) + Power(t1,2) + 
                  t1*(7 - 3*t2) + s2*(15 - 4*t1 - 3*t2) + t2) - 
               Power(s,2)*(4 - 7*t1 + Power(t1,2) + s2*(2 + t1) + 
                  2*t2 + t1*t2 + s1*(-7 + s2 + t1 + t2)) + 
               s1*(-15 - Power(s2,3) + Power(t1,3) - 
                  Power(t1,2)*(-2 + t2) + 5*t2 + Power(t2,2) + 
                  2*t1*(4 + t2) - Power(s2,2)*(-7 + 5*t1 + t2) + 
                  s2*(-3 - Power(t1,2) + t1*(18 - 4*t2) + t2)) + 
               s*(3 + 2*t1 - 4*Power(t1,2) - Power(t1,3) + 
                  2*Power(s2,2)*(2 + t1) - 4*t2 + Power(t1,2)*t2 + 
                  Power(s1,2)*(-13 + 3*s2 + t1 + 3*t2) + 
                  s2*(3*Power(t1,2) + 2*t1*(-8 + t2) + 3*(1 + t2)) + 
                  s1*(7 + 2*Power(s2,2) + 4*t1*(-4 + t2) + t2 + 
                     2*s2*(-7 + 3*t1 + t2))))/
             ((s - s2 + t1)*(s - s1 + t2)) + 
            (-3 + 3*s2 + 5*Power(s2,2) + 9*t1 + s2*t1 - 
               4*Power(s2,2)*t1 - Power(t1,2) + 2*s2*Power(t1,2) - 
               Power(t1,3) - 5*t2 + 4*s2*t2 + 2*Power(s2,2)*t2 + 
               4*t1*t2 + s2*t1*t2 - Power(s2,2)*t1*t2 + 
               3*Power(t1,2)*t2 - 2*Power(t2,2) + s2*Power(t2,2) - 
               4*t1*Power(t2,2) + 
               Power(s1,2)*(3 + 2*t2 - t1*(2 + t2)) - 
               Power(s,2)*(-1 + t1*(2 + t2)) + 
               s*(s1*(-1 + t1)*(3 + 2*t2) - 
                  2*(1 + t1 + Power(t1,2) + 2*t2 + t1*t2) + 
                  s2*(-7 - t2 + 2*t1*(3 + t2))) + 
               s1*(2 + 3*t2 + t1*(2 + 5*t2) + 
                  s2*(8 + 3*t2 - t1*(5 + 2*t2))))/((-1 + t1)*(-1 + t2)) \
+ (-3 - 4*s2 - 2*Power(s2,3) - t1 + 10*Power(s2,2)*t1 - 
               Power(s2,3)*t1 + Power(t1,2) - s2*Power(t1,2) - 
               3*Power(t1,3) - 4*t2 + 3*s2*t2 + 3*t1*t2 + 7*s2*t1*t2 - 
               3*Power(s2,2)*t1*t2 + Power(t1,2)*t2 + 
               s2*Power(t1,2)*t2 + 4*Power(t2,2) + s2*Power(t2,2) + 
               t1*Power(t2,2) - s2*t1*Power(t2,2) + Power(t2,3) - 
               Power(s1,3)*(-2 + s2 + t2) - 
               Power(s,2)*(1 - 3*t1 + s2*(2 + t1) + t2 + 2*t1*t2 + 
                  s1*(-4 + s2 + t2)) - 
               Power(s1,2)*(-1 + 2*Power(s2,2) + (-4 + t1)*t2 + 
                  Power(t2,2) + s2*(-9 + t1 + 3*t2)) + 
               s*(2 - 2*t1 + Power(t1,2) + 2*Power(s2,2)*(2 + t1) - 
                  t2 - 7*t1*t2 - Power(t1,2)*t2 + Power(t2,2) + 
                  t1*Power(t2,2) + 2*Power(s1,2)*(-3 + s2 + t2) + 
                  s2*(-13*t1 + 2*t2 + 5*t1*t2) + 
                  s1*(1 + 2*Power(s2,2) - 4*t1 - 3*t2 + 3*t1*t2 + 
                     Power(t2,2) + s2*(-9 + 2*t1 + 3*t2))) - 
               s1*(4 + Power(s2,3) - Power(t1,2)*(-5 + t2) - 3*t2 + 
                  Power(s2,2)*(-5 + 2*t1 + 2*t2) + 
                  t1*(-4 - 8*t2 + Power(t2,2)) + 
                  s2*((-4 + t2)*t2 + t1*(-11 + 4*t2))))/
             ((-1 + s2)*(-s + s2 - t1))))/
        (-Power(s*(-1 + s2) - (s2 - t1)*(-1 + s1 + t1 - t2),2) + 
          (-1 + t1)*(-1 + s1 + t1 - t2)*Power(1 - s2 + t1 - t2,3) - 
          Power(-1 + s2 - t1 + t2,2)*
           (3 + Power(s2,2) - s2*(4 - s + s2 - t1) - 
             2*(-1 + s - s2 + t1) + 2*Power(t1,2)*(-1 + s1 + t1 - t2) + 
             (2 + s - s2 + t1)*(-1 + s - s1 + t2) - 
             t1*(-1 + s1 + t1 - t2)*(-1 + s - s1 + s2 + t2)) + 
          (1 - s2 + t1 - t2)*(-2*Power(s,2)*(-1 + s2) + 
             s*(t1 - Power(t1,3) - 
                s1*(2 + t1 + Power(t1,2) - s2*(3 + t1)) + 2*t2 + t1*t2 + 
                Power(t1,2)*t2 + s2*(-1 + Power(t1,2) - 3*t2 - t1*t2)) - 
             (s2 - t1)*(3 + Power(t1,3) + Power(s1,2)*(1 + t1) + 
                Power(t1,2)*(1 - 2*t2) + 4*t2 + Power(t2,2) + 
                2*s1*(-2 + t1 + Power(t1,2) - t2 - t1*t2) + 
                t1*(-5 - 2*t2 + Power(t2,2))))) + 
       4*((2 - 2*s + 4*s1 + 8*s2 + 2*s*s2 - 2*Power(s2,2) + 6*t1 - 
             4*s1*t1 - 4*Power(t1,2) + 4*t2 - 2*s*t2 + 6*s1*t2 + 
             8*t1*t2 + 2*s*t1*t2 - 2*s1*t1*t2 - 2*s2*t1*t2 - 
             8*Power(t2,2) + 
             ((2 + s2 - t1 + 2*t2)*
                (-2*Power(t1,2) + s1*(s2 + t1)*(-1 + t2) + 
                  (1 + t2)*(1 + (-1 + s - s2)*t2) + 
                  t1*(1 + t2 + 2*s2*t2 - s*(1 + t2))))/(t1 - s2*t2))/
           ((-1 + t1)*(-1 + t2)) - 
          (3 + 2*Power(t1,3) + t2 + 3*s*t2 + 5*s2*t2 + 4*s*s2*t2 - 
             16*Power(s2,2)*t2 + 4*Power(s2,3)*t2 - 3*Power(t2,2) + 
             4*s*Power(t2,2) - 22*s2*Power(t2,2) + 4*s*s2*Power(t2,2) - 
             2*Power(s2,2)*Power(t2,2) - 2*s*Power(s2,2)*Power(t2,2) + 
             2*Power(s2,3)*Power(t2,2) - Power(t2,3) + s*Power(t2,3) - 
             9*s2*Power(t2,3) - 2*s*s2*Power(t2,3) + 
             2*Power(s2,2)*Power(t2,3) - 
             2*Power(t1,2)*(4 + s2 + 6*t2 + s2*t2) + 
             t1*(-5 - 4*Power(s2,2) + 22*t2 + 9*Power(t2,2) + 
                s*(-7 + 2*(-4 + s2)*t2 + Power(t2,2)) + 
                2*s2*(8 + 5*t2 + 5*Power(t2,2))) + 
             s1*(2*Power(s2,2)*t2*(2 + t2) + 
                t1*(5 + 2*t1 - Power(t2,2)) - 
                s2*(3 + 6*t2 - 3*Power(t2,2) - 2*Power(t2,3) + 
                   4*t1*(1 + t2))))/((-1 + s2)*(-1 + t1)*(-t1 + s2*t2)) + 
          (-2 - 4*s1*t1 - 19*Power(t1,2) + 7*s1*Power(t1,2) - 
             8*Power(s1,2)*Power(t1,2) - 14*Power(t1,3) - 
             11*s1*Power(t1,3) - Power(s1,2)*Power(t1,3) - 
             11*Power(t1,4) - 6*s1*Power(t1,4) + 
             3*Power(s1,2)*Power(t1,4) - 2*Power(t1,5) + 2*t2 - 
             5*t1*t2 + s1*t1*t2 + 9*Power(t1,2)*t2 - 
             16*s1*Power(t1,2)*t2 + 5*Power(s1,2)*Power(t1,2)*t2 - 
             Power(t1,3)*t2 + 6*s1*Power(t1,3)*t2 - 
             5*Power(s1,2)*Power(t1,3)*t2 + Power(t1,4)*t2 + 
             2*s1*Power(t1,4)*t2 + Power(s1,2)*Power(t1,4)*t2 + 
             4*Power(t2,2) - 2*Power(s2,5)*(-1 + t1)*Power(t2,2) + 
             9*t1*Power(t2,2) + 7*s1*t1*Power(t2,2) + 
             17*Power(t1,2)*Power(t2,2) + 
             5*s1*Power(t1,2)*Power(t2,2) + 
             2*Power(s1,2)*Power(t1,2)*Power(t2,2) + 
             7*Power(t1,3)*Power(t2,2) - 3*s1*Power(t1,3)*Power(t2,2) - 
             2*Power(s1,2)*Power(t1,3)*Power(t2,2) - 4*Power(t2,3) - 
             3*t1*Power(t2,3) - s1*t1*Power(t2,3) - 
             7*Power(t1,2)*Power(t2,3) + 4*s1*Power(t1,2)*Power(t2,3) + 
             Power(s1,2)*Power(t1,2)*Power(t2,3) - 2*Power(t2,4) - 
             t1*Power(t2,4) - 3*s1*t1*Power(t2,4) + 2*Power(t2,5) + 
             Power(s2,4)*t2*(2*Power(t1,2)*(2 + t2) + 
                t1*(-4 + (5 + 4*s)*t2 - 5*Power(t2,2)) + 
                t2*(4 + 2*s*(-1 + t2) + t2 + Power(t2,2)) - 
                s1*(-2 + t1 - 9*t2 + 5*t1*t2 + Power(t2,2))) + 
             Power(s,2)*(2*(-1 + t2)*Power(t2,2)*Power(1 + t2,2) + 
                Power(t1,4)*(3 + t2) - 
                3*Power(t1,3)*(-1 + 2*t2 + Power(t2,2)) + 
                t1*t2*(6 + 5*t2 - 10*Power(t2,2) - 5*Power(t2,3)) + 
                Power(t1,2)*(-4 - 2*t2 + 11*Power(t2,2) + 
                   5*Power(t2,3))) + 
             s*(-4*t2*Power(-1 + Power(t2,2),2) + 
                Power(t1,2)*(21 + 7*(1 + s1)*t2 - 
                   7*(4 + s1)*Power(t2,2) - 2*(5 + 3*s1)*Power(t2,3)) \
- 2*Power(t1,4)*(-2 + t2 + s1*(3 + t2)) + 
                Power(t1,3)*(7 + 7*t2 + 6*Power(t2,2) + 
                   s1*(-2 + 11*t2 + 5*Power(t2,2))) + 
                t1*(2 - 9*t2 - 20*Power(t2,2) + 17*Power(t2,3) + 
                   10*Power(t2,4) + 
                   s1*(4 - 7*Power(t2,2) + 3*Power(t2,4)))) + 
             Power(s2,3)*(-2*Power(t1,3)*(1 + 2*t2) + 
                Power(t1,2)*(2 - (10 + 9*s)*t2 + 
                   (7 - 5*s)*Power(t2,2) + 3*Power(t2,3)) + 
                t2*(-2 + (3 - 6*s)*t2 - 
                   2*(-7 + 5*s + Power(s,2))*Power(t2,2) + 
                   2*(-4 + s)*Power(t2,3) + Power(t2,4)) - 
                Power(s1,2)*(-2 - 7*Power(t2,2) + Power(t2,3) + 
                   t1*t2*(1 + 3*t2)) - 
                t1*t2*(7 + 4*t2 + 2*Power(s,2)*t2 - 7*Power(t2,2) + 
                   4*Power(t2,3) + s*(-6 + 7*t2 - 5*Power(t2,2))) + 
                s1*(Power(t1,2)*(1 + 11*t2 + 6*Power(t2,2)) + 
                   t1*(-2 + (-19 + s)*t2 + (1 + 5*s)*Power(t2,2) - 
                      6*Power(t2,3)) + 
                   t2*(2*t2*(-1 + 5*t2) + 
                      s*(-2 - 7*t2 + 3*Power(t2,2))))) + 
             Power(s2,2)*(2*Power(t1,4) + 
                t1*(2 + (-5 + 13*s - 2*Power(s,2))*t2 + 
                   (-39 + 30*s + 4*Power(s,2))*Power(t2,2) + 
                   2*(7 + s)*Power(t2,3) + (4 + 3*s)*Power(t2,4)) + 
                Power(t1,3)*(5 + t2 - 8*Power(t2,2) + 
                   s*(5 + 10*t2 + Power(t2,2))) + 
                Power(s1,2)*(-2 - 2*t2 - 3*Power(t2,2) + 
                   8*Power(t2,3) - Power(t2,4) + 
                   Power(t1,2)*(1 + 7*t2 + 4*Power(t2,2)) - 
                   t1*(2 + 14*t2 + 7*Power(t2,2) + Power(t2,3))) + 
                Power(t1,2)*(3 + 4*t2 - 29*Power(t2,2) + 
                   10*Power(t2,3) + Power(s,2)*t2*(5 + 3*t2) - 
                   2*s*(2 - 5*t2 + 4*Power(t2,2) + 2*Power(t2,3))) + 
                t2*(Power(s,2)*t2*(2 + 3*t2 - 3*Power(t2,2)) + 
                   2*t2*(-9 + 4*t2 + 8*Power(t2,2) - 3*Power(t2,3)) - 
                   s*(2 - 9*t2 + 9*Power(t2,2) + 8*Power(t2,3))) + 
                s1*(-4 + 2*t2 + 2*Power(t2,2) + Power(t2,3) - 
                   2*Power(t2,4) + Power(t2,5) - 
                   Power(t1,3)*(6 + 11*t2 + Power(t2,2)) + 
                   t1*t2*(4 - 31*t2 + 6*Power(t2,2) - 3*Power(t2,3)) + 
                   Power(t1,2)*
                    (10 - 2*t2 + 7*Power(t2,2) + 3*Power(t2,3)) + 
                   s*(-(Power(t1,2)*(1 + 12*t2 + 7*Power(t2,2))) + 
                      t1*(6 + 16*t2 + 3*Power(t2,2) + Power(t2,3)) + 
                      t2*(2 - 3*t2 - 9*Power(t2,2) + 4*Power(t2,3))))) \
- s2*(-2 - 2*Power(t1,2) + Power(t1,3) + 3*Power(t1,4) + 2*t2 - 
                37*t1*t2 - 39*Power(t1,2)*t2 - 32*Power(t1,3)*t2 - 
                7*Power(t1,4)*t2 - Power(t2,2) + 17*t1*Power(t2,2) + 
                5*Power(t1,2)*Power(t2,2) + 7*Power(t1,3)*Power(t2,2) + 
                5*Power(t2,3) + 33*t1*Power(t2,3) + 
                12*Power(t1,2)*Power(t2,3) - 5*Power(t2,4) - 
                13*t1*Power(t2,4) + Power(t2,5) + 
                Power(s1,2)*t1*
                 (-4 - 15*t2 + 17*Power(t2,2) + 3*Power(t2,3) - 
                   Power(t2,4) + 
                   Power(t1,2)*(4 + 7*t2 + Power(t2,2)) - 
                   t1*(3 + 13*t2 + 8*Power(t2,2))) + 
                Power(s,2)*(-(t1*(-4 + t2)*t2*Power(1 + t2,2)) + 
                   Power(t2,3)*(-1 - 4*t2 + Power(t2,2)) + 
                   Power(t1,3)*(3 + 6*t2 + Power(t2,2)) - 
                   Power(t1,2)*(4 - 3*t2 + 4*Power(t2,2) + Power(t2,3))\
) + s*(Power(t1,4)*(5 + t2) + Power(t1,3)*(5 + t2 - 6*Power(t2,2)) + 
                   Power(t1,2)*
                    (7 + 27*t2 + 11*Power(t2,2) + 9*Power(t2,3)) + 
                   t2*(-6 - 9*t2 - 4*Power(t2,2) + 17*Power(t2,3) + 
                      2*Power(t2,4)) - 
                   2*t1*(-1 - 15*t2 + 5*Power(t2,2) + 18*Power(t2,3) + 
                      3*Power(t2,4))) + 
                s1*(-(Power(t1,4)*(5 + t2)) - 
                   (-1 + Power(t2,2))*
                    (-4 + 5*Power(t2,2) + (-1 + s)*Power(t2,3)) - 
                   Power(t1,3)*
                    (2 + 5*t2 - 5*Power(t2,2) + 
                      s*(7 + 13*t2 + 2*Power(t2,2))) + 
                   Power(t1,2)*
                    (2 - 32*t2 + 12*Power(t2,2) - 6*Power(t2,3) + 
                      s*(9 + 12*t2 + 12*Power(t2,2) + Power(t2,3))) + 
                   t1*(-2 + t2 - 7*Power(t2,2) + 7*Power(t2,3) + 
                      Power(t2,4) + 
                      s*(10 - 3*t2 - 18*Power(t2,2) - 3*Power(t2,3) + 
                         2*Power(t2,4))))))/
           ((-1 + s1)*(-1 + t2)*(-1 + s2 - t1 + t2)*Power(t1 - s2*t2,2)) \
+ (6*Power(t1,5) - 2*Power(s2,6)*Power(t2,2) + 
             2*(-1 + t2)*Power(1 + t2,2)*Power(1 + (-1 + s)*t2,2) + 
             Power(t1,4)*(9 - 9*t2 - s*(3 + t2) + s1*(5 + t2)) + 
             Power(s2,5)*t2*((13 + 4*s - 9*t2)*t2 + 2*t1*(2 + t2) - 
                s1*(1 + 5*t2)) - 
             Power(t1,3)*(-18 + 24*t2 + 8*Power(t2,2) + 
                Power(s1,2)*(3 + t2) + Power(s,2)*(7 + t2) + 
                s1*(4 + 13*t2 + Power(t2,2)) - 
                2*s*(-2 + 11*t2 + 4*Power(t2,2) + s1*(5 + t2))) + 
             Power(s2,4)*(-2*Power(t1,2)*(1 + 2*t2) - 
                Power(s1,2)*t2*(1 + 3*t2) + 
                t2*(1 - (5 + 20*s + 2*Power(s,2))*t2 + 
                   (29 + 16*s)*Power(t2,2) - 11*Power(t2,3)) + 
                t1*t2*(-26 + 9*t2 + 9*Power(t2,2) - s*(9 + 5*t2)) + 
                s1*(t2*(7 + s + 23*t2 + 5*s*t2 - 16*Power(t2,2)) + 
                   t1*(1 + 11*t2 + 6*Power(t2,2)))) + 
             t1*(Power(s,2)*t2*
                 (7 + 4*t2 - 11*Power(t2,2) - 4*Power(t2,3)) - 
                Power(-1 + t2,2)*
                 (-4 + 6*t2 + 6*Power(t2,2) + 
                   s1*(5 + 7*t2 + 2*Power(t2,2))) + 
                s*(-1 + t2)*(-3 + 4*t2 + 29*Power(t2,2) + 
                   14*Power(t2,3) + 
                   s1*(-4 - 3*t2 + 3*Power(t2,2) + 2*Power(t2,3)))) + 
             Power(t1,2)*(Power(s,2)*
                 (-11 + 9*t2 + 10*Power(t2,2) + 2*Power(t2,3)) + 
                (-1 + t2)*(Power(s1,2)*(7 + t2) + 
                   2*s1*(-1 + 8*t2 + Power(t2,2)) + 
                   5*(-1 + 4*t2 + 3*Power(t2,2))) - 
                s*(-18 - 21*t2 + 34*Power(t2,2) + 17*Power(t2,3) + 
                   s1*(-6 + 5*t2 + 5*Power(t2,2) + 2*Power(t2,3)))) + 
             Power(s2,3)*(2*Power(t1,3) + 
                Power(s1,2)*(2 + t1 + 2*t2 + 7*t1*t2 + 
                   7*Power(t2,2) + 4*t1*Power(t2,2) - 7*Power(t2,3)) + 
                t1*(-1 + (9 + 48*s + 5*Power(s,2))*t2 + 
                   (-70 - 11*s + 3*Power(s,2))*Power(t2,2) + 
                   (14 - 17*s)*Power(t2,3) + 6*Power(t2,4)) + 
                Power(t1,2)*(13 + 9*t2 - 26*Power(t2,2) - 
                   2*Power(t2,3) + s*(5 + 10*t2 + Power(t2,2))) - 
                t2*(7 + 11*t2 + 7*Power(s,2)*(-1 + t2)*t2 - 
                   15*Power(t2,2) - 21*Power(t2,3) + 4*Power(t2,4) + 
                   s*(1 - 5*t2 + 38*Power(t2,2) - 18*Power(t2,3))) - 
                s1*(Power(t1,2)*(6 + 11*t2 + Power(t2,2)) + 
                   t1*(7 + s + 58*t2 + 12*s*t2 - Power(t2,2) + 
                      7*s*Power(t2,2) - 16*Power(t2,3)) + 
                   t2*(7 - 40*Power(t2,2) + 15*Power(t2,3) + 
                      s*(4 + 14*t2 - 14*Power(t2,2))))) + 
             Power(s2,2)*(t2*
                 (8 - 2*(6 - 7*s + Power(s,2))*t2 + 
                   (-24 - 10*s + 11*Power(s,2))*Power(t2,2) + 
                   (22 - 22*s - 7*Power(s,2))*Power(t2,3) + 
                   6*(1 + s)*Power(t2,4)) - 
                Power(t1,3)*
                 (9 - 25*t2 - 10*Power(t2,2) + s*(5 + t2)) - 
                Power(s1,2)*(2 + 3*t2 - 9*Power(t2,3) + 
                   4*Power(t2,4) + 
                   Power(t1,2)*(4 + 7*t2 + Power(t2,2)) + 
                   t1*(4 + 18*t2 - 3*Power(t2,2) - 7*Power(t2,3))) - 
                Power(t1,2)*(4 - 52*t2 - 15*Power(t2,2) + 
                   21*Power(t2,3) + 
                   Power(s,2)*(3 + 6*t2 + Power(t2,2)) + 
                   s*(28 + 31*t2 - 23*Power(t2,2) - 4*Power(t2,3))) + 
                t1*(7 + 34*t2 - 10*Power(t2,2) - 78*Power(t2,3) + 
                   5*Power(t2,4) + 
                   Power(s,2)*t2*(-19 + 7*t2 + 8*Power(t2,2)) + 
                   s*(1 - 20*t2 + 91*Power(t2,2) - 4*Power(t2,3) - 
                      10*Power(t2,4))) + 
                s1*(-4 + (4 + 3*s)*t2 - 2*(4 + s)*Power(t2,2) - 
                   18*s*Power(t2,3) + (12 + 11*s)*Power(t2,4) - 
                   4*Power(t2,5) + Power(t1,3)*(5 + t2) + 
                   Power(t1,2)*
                    (35 + 44*t2 - 15*Power(t2,2) - 4*Power(t2,3) + 
                      s*(7 + 13*t2 + 2*Power(t2,2))) + 
                   t1*(7 + 16*t2 - 86*Power(t2,2) + Power(t2,3) + 
                      8*Power(t2,4) + 
                      s*(8 + 37*t2 - 10*Power(t2,2) - 15*Power(t2,3))))\
) + s2*(2 - 6*t2 + 5*s*t2 + 10*Power(t2,2) + 7*s*Power(t2,2) - 
                Power(s,2)*Power(t2,2) - 6*Power(t2,3) + 
                9*s*Power(t2,3) + 2*Power(s,2)*Power(t2,3) - 
                4*Power(t2,4) - 15*s*Power(t2,4) + 
                5*Power(s,2)*Power(t2,4) + 4*Power(t2,5) - 
                6*s*Power(t2,5) - 2*Power(s,2)*Power(t2,5) - 
                2*Power(t1,4)*(4 + 7*t2) + 
                Power(t1,3)*(-11 - 27*t2 + 24*Power(t2,2) + 
                   Power(s,2)*(3 + t2) + s*(26 - 3*t2 - 3*Power(t2,2))) \
+ Power(s1,2)*t1*(5 + 11*t2 - 19*Power(t2,2) + Power(t2,3) + 
                   2*Power(t2,4) + Power(t1,2)*(3 + t2) + 
                   t1*(7 + 11*t2 - 4*Power(t2,2) - 2*Power(t2,3))) + 
                Power(t1,2)*(-23 - 23*t2 + 81*Power(t2,2) + 
                   7*Power(t2,3) + 
                   Power(s,2)*
                    (14 + 9*t2 - 9*Power(t2,2) - 2*Power(t2,3)) + 
                   s*(15 - 49*t2 - 36*Power(t2,2) + 2*Power(t2,3))) + 
                t1*(-8 + 7*t2 + 49*Power(t2,2) - 27*Power(t2,3) - 
                   21*Power(t2,4) + 
                   Power(s,2)*t2*
                    (7 - 26*t2 + 3*Power(t2,2) + 4*Power(t2,3)) + 
                   s*(-4 - 32*t2 - 3*Power(t2,2) + 56*Power(t2,3) + 
                      7*Power(t2,4))) + 
                s1*(4 + t2 - 11*Power(t2,2) + Power(t2,3) + 
                   7*Power(t2,4) - 2*Power(t2,5) + 
                   Power(t1,3)*(-29 - 6*t2 + 3*Power(t2,2)) + 
                   Power(t1,2)*
                    (-16 + 50*t2 + 27*Power(t2,2) - 7*Power(t2,3)) + 
                   2*t1*t2*(7 + 5*t2 - 15*Power(t2,2) + 
                      3*Power(t2,3)) - 
                   s*(2*Power(t1,3)*(3 + t2) + 
                      Power(t1,2)*
                       (23 + 22*t2 - 13*Power(t2,2) - 4*Power(t2,3)) + 
                      Power(t2,2)*
                       (-1 + 2*t2 + Power(t2,2) - 2*Power(t2,3)) + 
                      t1*(11 + 4*t2 - 39*Power(t2,2) + 6*Power(t2,3) + 
                         6*Power(t2,4))))))/
           ((-1 + s2)*(-s + s2 - t1)*(-1 + s2 - t1 + t2)*
             Power(t1 - s2*t2,2)) + 
          (2*(-105 + 2*Power(s1,3) + 132*s2 - 15*Power(s2,2) - 
               Power(s2,3) + 22*t1 - 73*s2*t1 + 6*Power(s2,2)*t1 + 
               114*Power(t1,2) - 45*s2*Power(t1,2) + 
               Power(s,2)*(5 + 2*s1 + s2 + t1 - 3*t2) - 
               16*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)*(-1 + t2) + 
               173*t2 - 128*s2*t2 + 33*Power(s2,2)*t2 - 133*t1*t2 + 
               157*s2*t1*t2 - 126*Power(t1,2)*t2 - 11*Power(t2,2) - 
               92*s2*Power(t2,2) + 202*t1*Power(t2,2) - 
               90*Power(t2,3) - 
               2*Power(s1,2)*(-14 + 7*s2 - 7*t1 + 11*t2) - 
               s*(-26 + 4*Power(s1,2) + 8*t1 - 21*Power(t1,2) + 
                  s1*(33 - 13*s2 + 15*t1 - 25*t2) + 
                  s2*(22 + t1 - 6*t2) + 31*t2 + 35*t1*t2 - 
                  38*Power(t2,2)) - 
               s1*(-39 + 11*Power(s2,2) + 21*Power(t1,2) + 
                  s2*(63 + 10*t1 - 117*t2) + 132*t2 - 98*Power(t2,2) + 
                  t1*(-74 + 79*t2)) - 
               ((-2 + Power(s,2) + 2*Power(s1,2) - 3*s2 + 
                    Power(s2,2) - 7*t1 + 2*s2*t1 - 
                    s*(-3 + 3*s1 + 2*s2 + 2*t1) + 
                    s1*(-6 + 3*s2 + 2*t1))*
                  (2*Power(t1,2) - s1*(s2 + t1)*(-1 + t2) - 
                    (1 + t2)*(1 + (-1 + s - s2)*t2) + 
                    t1*(-1 + s - t2 + s*t2 - 2*s2*t2)))/(-t1 + s2*t2) - 
               8*(-13 + 15*s2 - Power(s2,2) - 2*t1 - 11*s2*t1 + 
                  15*Power(t1,2) - 2*s2*Power(t1,2) + 17*t2 - 
                  7*s2*t2 + 2*Power(s2,2)*t2 - 23*t1*t2 + 
                  15*s2*t1*t2 - 14*Power(t1,2)*t2 + 9*Power(t2,2) - 
                  13*s2*Power(t2,2) + 26*t1*Power(t2,2) - 
                  13*Power(t2,3) - Power(s1,2)*(-1 + s2 - t1 + t2) + 
                  s*(2 - s2 + t1 + Power(t1,2) - 3*t2 - 2*t1*t2 + 
                     2*Power(t2,2) + s1*(-1 + s2 - t1 + t2)) - 
                  s1*(-9 + Power(s2,2) + Power(t1,2) + 22*t2 - 
                     13*Power(t2,2) + 2*t1*(-5 + 6*t2) - 
                     2*s2*(-5 + 7*t2))) - 
               2*(-132 + 2*Power(s1,3) + 165*s2 - 15*Power(s2,2) - 
                  2*Power(s2,3) - 7*t1 - 105*s2*t1 + 
                  8*Power(s2,2)*t1 + 153*Power(t1,2) - 
                  53*s2*Power(t1,2) + 
                  Power(s1,2)*(22 - 15*s2 + 18*t1 - 23*t2) + 
                  Power(s,2)*(2 + 2*s1 + s2 + 2*t1 - 3*t2) + 215*t2 - 
                  140*s2*t2 + 37*Power(s2,2)*t2 - 189*t1*t2 + 
                  196*s2*t1*t2 - 159*Power(t1,2)*t2 + 14*Power(t2,2) - 
                  128*s2*Power(t2,2) + 263*t1*Power(t2,2) - 
                  123*Power(t2,3) + 
                  s*(38 - 4*Power(s1,2) + Power(s2,2) + 7*t1 + 
                     27*Power(t1,2) - 50*t2 - 47*t1*t2 + 
                     47*Power(t2,2) + s2*(-24 - 7*t1 + 8*t2) + 
                     2*s1*(-12 + 7*s2 - 10*t1 + 13*t2)) - 
                  s1*(-52 + 16*Power(s2,2) + 27*Power(t1,2) + 
                     s2*(88 + 3*t1 - 149*t2) + 175*t2 - 
                     127*Power(t2,2) + 3*t1*(-27 + 34*t2))) - 
               4*(64 - 77*s2 + 5*Power(s2,2) + Power(s2,3) + 11*t1 + 
                  50*s2*t1 - 2*Power(s2,2)*t1 - 77*Power(t1,2) + 
                  20*s2*Power(t1,2) + Power(s,2)*(-1 + t2) - 95*t2 + 
                  52*s2*t2 - 14*Power(s2,2)*t2 + 100*t1*t2 - 
                  85*s2*t1*t2 + 72*Power(t1,2)*t2 - 24*Power(t2,2) + 
                  63*s2*Power(t2,2) - 125*t1*Power(t2,2) + 
                  61*Power(t2,3) + 
                  Power(s1,2)*(-9 + 7*s2 - 6*t1 + 8*t2) - 
                  s*(16 + Power(s2,2) + 5*t1 + 11*Power(t1,2) - 
                     23*t2 - 19*t1*t2 + 18*Power(t2,2) + 
                     s2*(-10 - 2*t1 + 3*t2) + 
                     s1*(-10 + 7*s2 - 6*t1 + 9*t2)) + 
                  s1*(-30 + 8*Power(s2,2) + 11*Power(t1,2) + 
                     s2*(41 + t1 - 69*t2) + 91*t2 - 61*Power(t2,2) + 
                     t1*(-43 + 51*t2))) + 
               ((3 + s - s1 - s2)*
                  (-2*Power(t1,5) + 
                    (-1 + t2)*Power(1 + t2,2)*
                     Power(1 + (-1 + s)*t2,2) - 
                    2*Power(t1,4)*(s + s1 - 2*t2 + s*t2 - s1*t2) + 
                    Power(s2,3)*
                     (Power(s1,2)*Power(-1 + t2,2) + 
                       2*s1*(t1 - t2)*(-1 + t2)*t2 + 
                       Power(t2,2)*
                        (1 + 2*Power(t1,2) + Power(t2,2) - 
                        2*t1*(1 + t2))) + 
                    Power(t1,3)*
                     (5 - 2*t2 - 5*Power(t2,2) - 
                       Power(s1,2)*(-1 + Power(t2,2)) - 
                       2*s1*(2 - 3*t2 + Power(t2,2)) - 
                       Power(s,2)*(-1 + 2*t2 + Power(t2,2)) + 
                       2*s*(-2 + 3*t2 + 3*Power(t2,2) + 
                        s1*(-2 + t2 + Power(t2,2)))) + 
                    t1*(-(Power(-1 + t2,2)*
                        (-1 + 2*(1 + s1)*t2 + (1 + 2*s1)*Power(t2,2))\
) + Power(s,2)*t2*(2 + 3*t2 - 4*Power(t2,2) - 3*Power(t2,3)) + 
                       2*s*(-1 + t2)*
                        (t2*(1 + 6*t2 + 3*Power(t2,2)) + 
                        s1*(-1 - 2*t2 + Power(t2,3)))) + 
                    Power(t1,2)*
                     (Power(s,2)*
                        (-1 - 3*t2 + 5*Power(t2,2) + 3*Power(t2,3)) + 
                       (-1 + t2)*
                        (3 + 2*t2 + 3*Power(t2,2) + 
                        2*s1*(-3 + Power(t2,2)) + 
                        Power(s1,2)*(1 + Power(t2,2))) - 
                       2*s*(-3 - 4*t2 + 5*Power(t2,2) + 
                        4*Power(t2,3) + s1*(3 - 5*t2 + 2*Power(t2,3)))\
) + s2*(1 - 2*t2 + 4*s*t2 + 2*Power(t2,2) + 2*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) + 
                       2*s*Power(t2,3) + Power(t2,4) - 
                       6*s*Power(t2,4) + Power(s,2)*Power(t2,4) - 
                       2*s*Power(t2,5) + Power(t1,4)*(2 + 4*t2) + 
                       2*Power(t1,3)*
                        (-1 - t2 - 4*Power(t2,2) + s*Power(1 + t2,2)) \
+ Power(s1,2)*t1*(-1 + t2)*(-4*t2 + t1*(3 + t2)) + 
                       Power(t1,2)*
                        (-1 - 8*t2 + 7*Power(t2,2) + 8*Power(t2,3) + 
                        Power(s,2)*(1 + Power(t2,2)) - 
                        2*s*(1 - t2 + 5*Power(t2,2) + 3*Power(t2,3))) \
+ 2*t1*t2*(-(Power(s,2)*(1 + Power(t2,2))) - 
                        2*(-2 + t2 + Power(t2,3)) + 
                        s*(-4 - t2 + 6*Power(t2,2) + 3*Power(t2,3))) \
+ 2*s1*(-1 + t2)*(-(Power(t1,3)*t2) + 
                        (1 + t2)*(-1 + t2 + s*Power(t2,2)) - 
                        Power(t1,2)*(-((-4 + t2)*t2) + s*(1 + t2)) + 
                        t1*(-2*s*(-1 + t2) + Power(1 + t2,2)))) + 
                    Power(s2,2)*
                     (Power(s1,2)*(-1 + t2)*
                        (1 + Power(t2,2) - t1*(1 + t2)) - 
                       2*s1*(-1 + t2)*
                        (-1 + Power(t1,2) + t2 + Power(t2,2) + 
                        Power(t2,3) - t1*t2*(1 + 2*t2) + 
                        s*(t1 - Power(t2,2))) + 
                       t2*(-2*Power(t1,3)*(2 + t2) + 
                        4*Power(t1,2)*(1 + t2 + Power(t2,2)) - 
                        3*t1*t2*(-1 + 2*t2 + Power(t2,2)) + 
                        t2*(-5 + 3*t2 + Power(t2,2) + Power(t2,3)) - 
                        2*s*(1 - t2 + Power(t2,2) + Power(t2,3) + 
                        Power(t1,2)*(1 + t2) - 
                        t1*(1 + t2 + 2*Power(t2,2)))))))/
                ((-1 + s2 - t1 + t2)*Power(t1 - s2*t2,2))))/
           ((-1 + s1)*(-s + s1 - t2)) + 
          (2*(98 + 2*Power(s1,3) - 151*s2 + 19*Power(s2,2) + 
               3*Power(s2,3) - 31*t1 + 146*s2*t1 - 14*Power(s2,2)*t1 - 
               82*Power(t1,2) + 11*s2*Power(t1,2) + 
               Power(s,2)*(11 + 2*s1 - 7*s2 + 9*t1 - 11*t2) + 
               16*(-1 + s + t1)*(1 - s2 + t1 - t2)*(-1 + t2) - 179*t2 + 
               151*s2*t2 - 23*Power(s2,2)*t2 + 145*t1*t2 - 
               167*s2*t1*t2 + 102*Power(t1,2)*t2 + 58*Power(t2,2) + 
               50*s2*Power(t2,2) - 148*t1*Power(t2,2) + 32*Power(t2,3) + 
               2*Power(s1,2)*(1 + 7*s2 - 7*t1 + 3*t2) + 
               s1*(-63 + 5*Power(s2,2) + 26*t1 - Power(t1,2) + 
                  s2*(66 - 6*t1 - 73*t2) + 77*t2 + 7*t1*t2 - 
                  28*Power(t2,2)) + 
               s*(-28 - 4*Power(s1,2) + 4*Power(s2,2) - 56*t1 + 
                  Power(t1,2) + s2*(52 - 5*t1 - 60*t2) + 105*t2 + 
                  95*t1*t2 - 112*Power(t2,2) + 
                  s1*(-13 - 7*s2 + 5*t1 + 5*t2)) + 
               ((Power(s,2) + 2*Power(s1,2) - 6*s2 + 2*Power(s2,2) - 
                    5*t1 + s2*t1 - 2*t2 + s2*t2 - 
                    s*(-4 + 3*s1 + 3*s2 + t1 + t2) + 
                    s1*(-7 + 4*s2 + t1 + t2))*
                  (2*Power(t1,2) - s1*(s2 + t1)*(-1 + t2) - 
                    (1 + t2)*(1 + (-1 + s - s2)*t2) + 
                    t1*(-1 + s - t2 + s*t2 - 2*s2*t2)))/(t1 - s2*t2) + 
               8*(-15 + 16*s2 - Power(s2,2) + 4*t1 - 14*s2*t1 + 
                  11*Power(t1,2) + s1*(-3 + 2*s2 + t1)*(-1 + t2) + 
                  26*t2 - 16*s2*t2 + Power(s2,2)*t2 - 17*t1*t2 + 
                  16*s2*t1*t2 - 13*Power(t1,2)*t2 - 10*Power(t2,2) - 
                  2*s2*Power(t2,2) + 15*t1*Power(t2,2) - Power(t2,3) + 
                  Power(s,2)*(-1 + s2 - t1 + t2) + 
                  s*(10 - Power(s2,2) + 10*t1 - 23*t2 - 14*t1*t2 + 
                     15*Power(t2,2) - s1*(-1 + s2 - t1 + t2) + 
                     s2*(-10 + t1 + 12*t2))) - 
               4*(-75 + 86*s2 - 5*Power(s2,2) - Power(s2,3) + 29*t1 - 
                  79*s2*t1 + 2*Power(s2,2)*t1 + 48*Power(t1,2) + 
                  Power(s1,2)*(1 - 3*s2 + 4*t1 - 2*t2) + 124*t2 - 
                  84*s2*t2 + 8*Power(s2,2)*t2 - 92*t1*t2 + 
                  91*s2*t1*t2 - 63*Power(t1,2)*t2 - 40*Power(t2,2) - 
                  20*s2*Power(t2,2) + 81*t1*Power(t2,2) - 
                  10*Power(t2,3) + 
                  Power(s,2)*(-7 + 6*s2 - 6*t1 + 7*t2) - 
                  2*s1*(-13 + Power(s2,2) - Power(t1,2) - 
                     2*t1*(-3 + t2) + 14*t2 - 2*Power(t2,2) - 
                     s2*(-9 + t1 + 11*t2)) + 
                  s*(38 - 5*Power(s2,2) + 42*t1 - 2*Power(t1,2) + 
                     s1*(6 - 3*s2 + 2*t1 - 5*t2) - 97*t2 - 67*t1*t2 + 
                     77*Power(t2,2) + s2*(-44 + 6*t1 + 52*t2))) - 
               2*(154 + 2*Power(s1,3) - 192*s2 + 13*Power(s2,2) + 
                  4*Power(s2,3) - 61*t1 + 174*s2*t1 - 
                  8*Power(s2,2)*t1 - 103*Power(t1,2) + 
                  5*s2*Power(t1,2) + 
                  Power(s,2)*(12 + 2*s1 - 11*s2 + 14*t1 - 15*t2) - 
                  251*t2 + 188*s2*t2 - 23*Power(s2,2)*t2 + 193*t1*t2 - 
                  208*s2*t1*t2 + 133*Power(t1,2)*t2 + 77*Power(t2,2) + 
                  57*s2*Power(t2,2) - 184*t1*Power(t2,2) + 
                  32*Power(t2,3) + 
                  Power(s1,2)*(-8 + 17*s2 - 14*t1 + 9*t2) + 
                  s1*(-71 + 10*Power(s2,2) - 5*Power(t1,2) + 
                     s2*(48 - 3*t1 - 70*t2) + 80*t2 - 22*Power(t2,2) + 
                     t1*(26 + t2)) - 
                  s*(59 + 4*Power(s1,2) - 7*Power(s2,2) + 74*t1 - 
                     5*Power(t1,2) + s1*(4 + 6*s2 - 6*t2) - 169*t2 - 
                     132*t1*t2 + 158*Power(t2,2) + 
                     s2*(-88 + 15*t1 + 95*t2))) + 
               ((3 + s - s1 - s2)*
                  (-2*Power(t1,5) + 
                    (-1 + t2)*Power(1 + t2,2)*
                     Power(1 + (-1 + s)*t2,2) - 
                    2*Power(t1,4)*(s + s1 - 2*t2 + s*t2 - s1*t2) + 
                    Power(s2,3)*
                     (Power(s1,2)*Power(-1 + t2,2) + 
                       2*s1*(t1 - t2)*(-1 + t2)*t2 + 
                       Power(t2,2)*
                        (1 + 2*Power(t1,2) + Power(t2,2) - 
                        2*t1*(1 + t2))) + 
                    Power(t1,3)*
                     (5 - 2*t2 - 5*Power(t2,2) - 
                       Power(s1,2)*(-1 + Power(t2,2)) - 
                       2*s1*(2 - 3*t2 + Power(t2,2)) - 
                       Power(s,2)*(-1 + 2*t2 + Power(t2,2)) + 
                       2*s*(-2 + 3*t2 + 3*Power(t2,2) + 
                        s1*(-2 + t2 + Power(t2,2)))) + 
                    t1*(-(Power(-1 + t2,2)*
                        (-1 + 2*(1 + s1)*t2 + (1 + 2*s1)*Power(t2,2))) \
+ Power(s,2)*t2*(2 + 3*t2 - 4*Power(t2,2) - 3*Power(t2,3)) + 
                       2*s*(-1 + t2)*
                        (t2*(1 + 6*t2 + 3*Power(t2,2)) + 
                        s1*(-1 - 2*t2 + Power(t2,3)))) + 
                    Power(t1,2)*
                     (Power(s,2)*
                        (-1 - 3*t2 + 5*Power(t2,2) + 3*Power(t2,3)) + 
                       (-1 + t2)*
                        (3 + 2*t2 + 3*Power(t2,2) + 
                        2*s1*(-3 + Power(t2,2)) + 
                        Power(s1,2)*(1 + Power(t2,2))) - 
                       2*s*(-3 - 4*t2 + 5*Power(t2,2) + 
                        4*Power(t2,3) + s1*(3 - 5*t2 + 2*Power(t2,3)))) \
+ s2*(1 - 2*t2 + 4*s*t2 + 2*Power(t2,2) + 2*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) + 
                       2*s*Power(t2,3) + Power(t2,4) - 
                       6*s*Power(t2,4) + Power(s,2)*Power(t2,4) - 
                       2*s*Power(t2,5) + Power(t1,4)*(2 + 4*t2) + 
                       2*Power(t1,3)*
                        (-1 - t2 - 4*Power(t2,2) + s*Power(1 + t2,2)) \
+ Power(s1,2)*t1*(-1 + t2)*(-4*t2 + t1*(3 + t2)) + 
                       Power(t1,2)*
                        (-1 - 8*t2 + 7*Power(t2,2) + 8*Power(t2,3) + 
                        Power(s,2)*(1 + Power(t2,2)) - 
                        2*s*(1 - t2 + 5*Power(t2,2) + 3*Power(t2,3))) \
+ 2*t1*t2*(-(Power(s,2)*(1 + Power(t2,2))) - 
                        2*(-2 + t2 + Power(t2,3)) + 
                        s*(-4 - t2 + 6*Power(t2,2) + 3*Power(t2,3))) + 
                       2*s1*(-1 + t2)*
                        (-(Power(t1,3)*t2) + 
                        (1 + t2)*(-1 + t2 + s*Power(t2,2)) - 
                        Power(t1,2)*(-((-4 + t2)*t2) + s*(1 + t2)) + 
                        t1*(-2*s*(-1 + t2) + Power(1 + t2,2)))) + 
                    Power(s2,2)*
                     (Power(s1,2)*(-1 + t2)*
                        (1 + Power(t2,2) - t1*(1 + t2)) - 
                       2*s1*(-1 + t2)*
                        (-1 + Power(t1,2) + t2 + Power(t2,2) + 
                        Power(t2,3) - t1*t2*(1 + 2*t2) + 
                        s*(t1 - Power(t2,2))) + 
                       t2*(-2*Power(t1,3)*(2 + t2) + 
                         4*Power(t1,2)*(1 + t2 + Power(t2,2)) - 
                         3*t1*t2*(-1 + 2*t2 + Power(t2,2)) + 
                         t2*(-5 + 3*t2 + Power(t2,2) + Power(t2,3)) - 
                         2*s*
                         (1 - t2 + Power(t2,2) + Power(t2,3) + 
                         Power(t1,2)*(1 + t2) - 
                         t1*(1 + t2 + 2*Power(t2,2)))))))/
                ((-1 + s2 - t1 + t2)*Power(t1 - s2*t2,2))))/
           ((s - s2 + t1)*(s - s1 + t2))) + 
       ((-4*(-4*Power(t1,4) - 4*Power(t1,3)*(-7 + s1 - 2*t2) + 
               2*Power(s2,4)*t2*(1 + t2) - 
               (1 + (-1 + s)*t2)*
                (-7 - 15*t2 - 11*Power(t2,2) - 2*Power(t2,3) + 
                  Power(t2,4)) + 
               Power(t1,2)*(-49 - 60*t2 - 9*Power(t2,2) + 
                  2*Power(t2,3) - s1*(-13 + Power(t2,2)) + 
                  s*(17 + 10*t2 + Power(t2,2))) + 
               t1*(2 + 44*t2 + 42*Power(t2,2) + 11*Power(t2,3) - 
                  Power(t2,4) - 
                  s1*(7 + 7*t2 - 3*Power(t2,2) - 2*Power(t2,3) + 
                     Power(t2,4)) + 
                  s*(-15 - 30*t2 - 25*Power(t2,2) - Power(t2,3) + 
                     Power(t2,4))) + 
               Power(s2,3)*(-2*t1*(1 + Power(t2,2)) + 
                  t2*(-20 - 2*s*(-1 + t2) - 11*t2 + 5*Power(t2,2) + 
                     s1*(5 + t2))) - 
               Power(s2,2)*(2*Power(t1,2)*(1 + t2) + 
                  s1*(1 + t2)*
                   (1 - 5*t2 - 2*Power(t2,2) + t1*(4 + t2)) + 
                  t1*(-20 - 47*t2 - 3*Power(t2,2) + 4*Power(t2,3) + 
                     s*(2 + t2 - Power(t2,2))) + 
                  t2*(26 + 40*t2 + 18*Power(t2,2) - 4*Power(t2,3) + 
                     s*(-6 - 9*t2 + 5*Power(t2,2)))) + 
               s2*(1 - 8*t2 + 9*s*t2 - 48*Power(t2,2) + 
                  10*s*Power(t2,2) - 36*Power(t2,3) + 5*s*Power(t2,3) - 
                  8*Power(t2,4) - 4*s*Power(t2,4) + Power(t2,5) + 
                  4*Power(t1,3)*(1 + t2) + 
                  4*Power(t1,2)*(-9 + s - 9*t2 - Power(t2,2)) + 
                  s1*(-7 - t2 + 11*Power(t2,2) + 6*Power(t2,3) + 
                     Power(t2,4) + 4*Power(t1,2)*(1 + t2) - 
                     2*t1*(-1 + 8*t2 + 4*Power(t2,2) + Power(t2,3))) + 
                  t1*(25 + 89*t2 + 79*Power(t2,2) + 5*Power(t2,3) - 
                     2*Power(t2,4) + 
                     s*(-7 - 22*t2 + 3*Power(t2,2) + 2*Power(t2,3))))))/
           ((-1 + s2)*(-t1 + s2*t2)) + 
          (4*(-(((-9 + Power(s2,2) - Power(t1,2) - 4*t2 + 
                      s2*(-2 + (-1 + t1)*t2) + 
                      t1*(4 + 2*t2 + Power(t2,2)))*
                    (-2*Power(t1,2) + s1*(s2 + t1)*(-1 + t2) + 
                      (1 + t2)*(1 + (-1 + s - s2)*t2) + 
                      t1*(1 + t2 + 2*s2*t2 - s*(1 + t2))))/(t1 - s2*t2)\
) + 2*(5 + s2 + 6*Power(s2,2) + t1 - s2*t1 - 4*Power(s2,2)*t1 - 
                  3*Power(t1,2) + 3*s2*Power(t1,2) + Power(t1,3) + 
                  2*t2 + 8*s2*t2 + Power(s2,2)*t2 + 3*t1*t2 - 
                  s2*t1*t2 - Power(s2,2)*t1*t2 - Power(t1,2)*t2 + 
                  5*s2*Power(t2,2) + t1*Power(t2,2) - 
                  s2*t1*Power(t2,2) + 3*Power(t2,3) + 
                  s1*(1 + 2*Power(t1,2) + 5*t2 + Power(t2,2) - 
                     t1*(4 + 2*t2 + Power(t2,2)) - 
                     s2*(-5 + t1*(2 + t2))) + 
                  s*(2 - Power(t1,2) + t1*(-2 + 2*t2 + Power(t2,2)) + 
                     s2*(-2 + t2 + t1*(2 + t2))))))/(-1 + t2) + 
          (2*(-1 + t1)*((2*(2*Power(t1,2) - s1*(s2 + t1)*(-1 + t2) - 
                    (1 + t2)*(1 + (-1 + s - s2)*t2) + 
                    t1*(-1 + s - t2 + s*t2 - 2*s2*t2))*
                  (-1 - 6*s2 - 4*Power(s2,2) + Power(s2,3) - 2*t1 + 
                    2*s2*t1 - Power(s2,2)*t1 + 5*Power(t1,2) - 7*t2 - 
                    8*s2*t2 + 3*Power(s2,2)*t2 - 4*t1*t2 - 
                    3*Power(t2,2) + 2*s2*Power(t2,2) + 
                    t1*Power(t2,2) + 
                    s*(2 - Power(s2,2) + s2*(1 + t1 - 3*t2) + 
                       t1*(-1 + t2) + 2*t2 - 2*Power(t2,2)) + 
                    s1*(-5 + Power(s2,2) - s2*(1 + t1 - 3*t2) - 
                       t1*(-4 + t2) - 4*t2 + 2*Power(t2,2))))/
                (t1 - s2*t2) + 
               4*(-3 - 6*s2 - Power(s2,2) - 2*Power(s2,3) + 3*t1 + 
                  14*s2*t1 + 8*Power(s2,2)*t1 - Power(s2,3)*t1 - 
                  9*Power(t1,2) - 7*s2*Power(t1,2) + 
                  Power(s2,2)*Power(t1,2) + Power(t1,3) - s2*t2 - 
                  Power(s2,2)*t2 + 9*t1*t2 + 13*s2*t1*t2 - 
                  3*Power(s2,2)*t1*t2 - 3*Power(t1,2)*t2 + 
                  s2*Power(t1,2)*t2 + 2*Power(t2,2) + s2*Power(t2,2) + 
                  5*t1*Power(t2,2) - 2*s2*t1*Power(t2,2) + 
                  Power(t2,3) + Power(s,2)*(-t1 + t2) + 
                  Power(s1,2)*
                   (1 - Power(s2,2) + s2*(1 + t1 - 2*t2) + 
                     t1*(-1 + t2) + t2 - Power(t2,2)) + 
                  s*(4 - 7*t1 + 2*Power(t1,2) + 
                     Power(s2,2)*(2 + t1) + 2*t2 - 3*t1*t2 - 
                     Power(t1,2)*t2 + 2*Power(t2,2) + 
                     2*t1*Power(t2,2) - 
                     s2*(-1 + t1 + Power(t1,2) - 2*t2 - 3*t1*t2) + 
                     s1*(-2 + Power(s2,2) + 3*t1 - 
                        s2*(3 + t1 - 2*t2) - 4*t2 - t1*t2 + 
                        Power(t2,2))) + 
                  s1*(2 - Power(s2,3) + Power(s2,2)*(5 - 2*t2) + 
                     4*t2 + Power(t1,2)*t2 + Power(t2,2) - 
                     t1*(3 + 2*Power(t2,2)) + 
                     s2*(7 + Power(t1,2) + 7*t2 - Power(t2,2) - 
                        t1*(5 + 2*t2)))) - 
               (2*(-5 + 4*t1 + s2*(-1 + t2) - 4*t2 + Power(t2,2))*
                  (-2*Power(t1,5) + 
                    (-1 + t2)*Power(1 + t2,2)*
                     Power(1 + (-1 + s)*t2,2) - 
                    2*Power(t1,4)*(s + s1 - 2*t2 + s*t2 - s1*t2) + 
                    Power(s2,3)*
                     (Power(s1,2)*Power(-1 + t2,2) + 
                       2*s1*(t1 - t2)*(-1 + t2)*t2 + 
                       Power(t2,2)*
                        (1 + 2*Power(t1,2) + Power(t2,2) - 
                        2*t1*(1 + t2))) + 
                    Power(t1,3)*
                     (5 - 2*t2 - 5*Power(t2,2) - 
                       Power(s1,2)*(-1 + Power(t2,2)) - 
                       2*s1*(2 - 3*t2 + Power(t2,2)) - 
                       Power(s,2)*(-1 + 2*t2 + Power(t2,2)) + 
                       2*s*(-2 + 3*t2 + 3*Power(t2,2) + 
                        s1*(-2 + t2 + Power(t2,2)))) + 
                    t1*(-(Power(-1 + t2,2)*
                        (-1 + 2*(1 + s1)*t2 + (1 + 2*s1)*Power(t2,2))\
) + Power(s,2)*t2*(2 + 3*t2 - 4*Power(t2,2) - 3*Power(t2,3)) + 
                       2*s*(-1 + t2)*
                        (t2*(1 + 6*t2 + 3*Power(t2,2)) + 
                        s1*(-1 - 2*t2 + Power(t2,3)))) + 
                    Power(t1,2)*
                     (Power(s,2)*
                        (-1 - 3*t2 + 5*Power(t2,2) + 3*Power(t2,3)) + 
                       (-1 + t2)*
                        (3 + 2*t2 + 3*Power(t2,2) + 
                        2*s1*(-3 + Power(t2,2)) + 
                        Power(s1,2)*(1 + Power(t2,2))) - 
                       2*s*(-3 - 4*t2 + 5*Power(t2,2) + 
                        4*Power(t2,3) + s1*(3 - 5*t2 + 2*Power(t2,3)))\
) + s2*(1 - 2*t2 + 4*s*t2 + 2*Power(t2,2) + 2*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) + 
                       2*s*Power(t2,3) + Power(t2,4) - 
                       6*s*Power(t2,4) + Power(s,2)*Power(t2,4) - 
                       2*s*Power(t2,5) + Power(t1,4)*(2 + 4*t2) + 
                       2*Power(t1,3)*
                        (-1 - t2 - 4*Power(t2,2) + s*Power(1 + t2,2)) \
+ Power(s1,2)*t1*(-1 + t2)*(-4*t2 + t1*(3 + t2)) + 
                       Power(t1,2)*
                        (-1 - 8*t2 + 7*Power(t2,2) + 8*Power(t2,3) + 
                        Power(s,2)*(1 + Power(t2,2)) - 
                        2*s*(1 - t2 + 5*Power(t2,2) + 3*Power(t2,3))) \
+ 2*t1*t2*(-(Power(s,2)*(1 + Power(t2,2))) - 
                        2*(-2 + t2 + Power(t2,3)) + 
                        s*(-4 - t2 + 6*Power(t2,2) + 3*Power(t2,3))) \
+ 2*s1*(-1 + t2)*(-(Power(t1,3)*t2) + 
                        (1 + t2)*(-1 + t2 + s*Power(t2,2)) - 
                        Power(t1,2)*(-((-4 + t2)*t2) + s*(1 + t2)) + 
                        t1*(-2*s*(-1 + t2) + Power(1 + t2,2)))) + 
                    Power(s2,2)*
                     (Power(s1,2)*(-1 + t2)*
                        (1 + Power(t2,2) - t1*(1 + t2)) - 
                       2*s1*(-1 + t2)*
                        (-1 + Power(t1,2) + t2 + Power(t2,2) + 
                        Power(t2,3) - t1*t2*(1 + 2*t2) + 
                        s*(t1 - Power(t2,2))) + 
                       t2*(-2*Power(t1,3)*(2 + t2) + 
                        4*Power(t1,2)*(1 + t2 + Power(t2,2)) - 
                        3*t1*t2*(-1 + 2*t2 + Power(t2,2)) + 
                        t2*(-5 + 3*t2 + Power(t2,2) + Power(t2,3)) - 
                        2*s*(1 - t2 + Power(t2,2) + Power(t2,3) + 
                        Power(t1,2)*(1 + t2) - 
                        t1*(1 + t2 + 2*Power(t2,2)))))))/
                ((-1 + s2 - t1 + t2)*Power(t1 - s2*t2,2))))/
           ((-1 + s2)*(-s + s2 - t1)) + 
          (4*(-1 + t1)*(-(((2*Power(t1,2) - s1*(s2 + t1)*(-1 + t2) - 
                      (1 + t2)*(1 + (-1 + s - s2)*t2) + 
                      t1*(-1 + s - t2 + s*t2 - 2*s2*t2))*
                    (5 + 2*Power(s2,2) + 12*t1 + s2*t1 - 
                      Power(s2,2)*t1 - 3*Power(t1,2) + 
                      s2*Power(t1,2) - 3*t2 + 4*s2*t2 - t1*t2 - 
                      2*s2*t1*t2 - Power(t1,2)*t2 + 4*Power(t2,2) - 
                      s1*(-7 - Power(t1,2) - t2 + Power(t2,2) + 
                        t1*(3 + t2) + s2*(-4 + t1 + t2)) + 
                      s*(1 + t1 - Power(t1,2) + t1*t2 + Power(t2,2) + 
                        s2*(-1 + t1 + t2))))/(t1 - s2*t2)) + 
               2*(3 + 2*s2 + 2*Power(s2,2) - 2*t1 + 2*s2*t1 + 
                  4*Power(s2,2)*t1 + 4*Power(t1,2) - 
                  2*s2*Power(t1,2) - Power(s2,2)*Power(t1,2) - 
                  Power(t1,3) + s2*Power(t1,3) + 
                  Power(s,2)*(t1 - t2) + 3*t2 + s2*t2 - 
                  4*Power(s2,2)*t2 + 13*s2*t1*t2 - 2*Power(t1,2)*t2 - 
                  2*s2*Power(t1,2)*t2 - 3*Power(t2,2) - 
                  5*s2*Power(t2,2) + 10*t1*Power(t2,2) - 
                  3*Power(t2,3) + 
                  Power(s1,2)*
                   (2 - s2*(-3 + t1) + Power(t1,2) + 3*t2 - 
                     t1*(3 + t2)) + 
                  s*(1 + 2*Power(t1,2) - Power(t1,3) + t1*t2 + 
                     Power(t1,2)*t2 - Power(t2,2) + t1*Power(t2,2) + 
                     s2*(1 + Power(t1,2) + t1*(-1 + t2) + 2*t2) + 
                     s1*(-1 + s2*(-1 + t1) + t1 - Power(t1,2) + t1*t2)\
) + s1*(-3 + Power(t1,3) - 4*t2 + 4*Power(t2,2) + Power(t2,3) - 
                     Power(t1,2)*(4 + t2) + 
                     Power(s2,2)*(3 - t1 + t2) + 
                     t1*(8 + t2 - 2*Power(t2,2)) + 
                     s2*(-2 + t1*(2 - 3*t2) + 5*t2 + 2*Power(t2,2)))) + 
               ((7 + s2 - t2 - t1*(2 + t2))*
                  (-2*Power(t1,5) + 
                    (-1 + t2)*Power(1 + t2,2)*
                     Power(1 + (-1 + s)*t2,2) - 
                    2*Power(t1,4)*(s + s1 - 2*t2 + s*t2 - s1*t2) + 
                    Power(s2,3)*
                     (Power(s1,2)*Power(-1 + t2,2) + 
                       2*s1*(t1 - t2)*(-1 + t2)*t2 + 
                       Power(t2,2)*
                        (1 + 2*Power(t1,2) + Power(t2,2) - 
                        2*t1*(1 + t2))) + 
                    Power(t1,3)*
                     (5 - 2*t2 - 5*Power(t2,2) - 
                       Power(s1,2)*(-1 + Power(t2,2)) - 
                       2*s1*(2 - 3*t2 + Power(t2,2)) - 
                       Power(s,2)*(-1 + 2*t2 + Power(t2,2)) + 
                       2*s*(-2 + 3*t2 + 3*Power(t2,2) + 
                        s1*(-2 + t2 + Power(t2,2)))) + 
                    t1*(-(Power(-1 + t2,2)*
                        (-1 + 2*(1 + s1)*t2 + (1 + 2*s1)*Power(t2,2))\
) + Power(s,2)*t2*(2 + 3*t2 - 4*Power(t2,2) - 3*Power(t2,3)) + 
                       2*s*(-1 + t2)*
                        (t2*(1 + 6*t2 + 3*Power(t2,2)) + 
                        s1*(-1 - 2*t2 + Power(t2,3)))) + 
                    Power(t1,2)*
                     (Power(s,2)*
                        (-1 - 3*t2 + 5*Power(t2,2) + 3*Power(t2,3)) + 
                       (-1 + t2)*
                        (3 + 2*t2 + 3*Power(t2,2) + 
                        2*s1*(-3 + Power(t2,2)) + 
                        Power(s1,2)*(1 + Power(t2,2))) - 
                       2*s*(-3 - 4*t2 + 5*Power(t2,2) + 
                        4*Power(t2,3) + s1*(3 - 5*t2 + 2*Power(t2,3)))\
) + s2*(1 - 2*t2 + 4*s*t2 + 2*Power(t2,2) + 2*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) + 
                       2*s*Power(t2,3) + Power(t2,4) - 
                       6*s*Power(t2,4) + Power(s,2)*Power(t2,4) - 
                       2*s*Power(t2,5) + Power(t1,4)*(2 + 4*t2) + 
                       2*Power(t1,3)*
                        (-1 - t2 - 4*Power(t2,2) + s*Power(1 + t2,2)) \
+ Power(s1,2)*t1*(-1 + t2)*(-4*t2 + t1*(3 + t2)) + 
                       Power(t1,2)*
                        (-1 - 8*t2 + 7*Power(t2,2) + 8*Power(t2,3) + 
                        Power(s,2)*(1 + Power(t2,2)) - 
                        2*s*(1 - t2 + 5*Power(t2,2) + 3*Power(t2,3))) \
+ 2*t1*t2*(-(Power(s,2)*(1 + Power(t2,2))) - 
                        2*(-2 + t2 + Power(t2,3)) + 
                        s*(-4 - t2 + 6*Power(t2,2) + 3*Power(t2,3))) \
+ 2*s1*(-1 + t2)*(-(Power(t1,3)*t2) + 
                        (1 + t2)*(-1 + t2 + s*Power(t2,2)) - 
                        Power(t1,2)*(-((-4 + t2)*t2) + s*(1 + t2)) + 
                        t1*(-2*s*(-1 + t2) + Power(1 + t2,2)))) + 
                    Power(s2,2)*
                     (Power(s1,2)*(-1 + t2)*
                        (1 + Power(t2,2) - t1*(1 + t2)) - 
                       2*s1*(-1 + t2)*
                        (-1 + Power(t1,2) + t2 + Power(t2,2) + 
                        Power(t2,3) - t1*t2*(1 + 2*t2) + 
                        s*(t1 - Power(t2,2))) + 
                       t2*(-2*Power(t1,3)*(2 + t2) + 
                        4*Power(t1,2)*(1 + t2 + Power(t2,2)) - 
                        3*t1*t2*(-1 + 2*t2 + Power(t2,2)) + 
                        t2*(-5 + 3*t2 + Power(t2,2) + Power(t2,3)) - 
                        2*s*(1 - t2 + Power(t2,2) + Power(t2,3) + 
                        Power(t1,2)*(1 + t2) - 
                        t1*(1 + t2 + 2*Power(t2,2)))))))/
                ((-1 + s2 - t1 + t2)*Power(t1 - s2*t2,2))))/
           ((-1 + s1)*(-1 + t2)) + 
          ((-1 + t1)*((-4*(-2*Power(t1,2) + s1*(s2 + t1)*(-1 + t2) + 
                    (1 + t2)*(1 + (-1 + s - s2)*t2) + 
                    t1*(1 + t2 + 2*s2*t2 - s*(1 + t2)))*
                  (4 - 5*s2 - 5*Power(s2,2) + Power(s2,3) - 
                    Power(s,2)*(-1 + t1) - 7*t1 + s2*t1 + 
                    2*Power(s2,2)*t1 - 5*Power(t1,2) - 
                    2*s2*Power(t1,2) - 3*t2 - 9*s2*t2 + 
                    4*Power(s2,2)*t2 - 3*t1*t2 + 2*s2*t1*t2 - 
                    Power(t2,2) + 3*s2*Power(t2,2) + t1*Power(t2,2) + 
                    Power(s1,2)*(2 + s2 - 2*t1 + t2) - 
                    s*(-2 + Power(s2,2) + 5*t1 - Power(t1,2) - 
                       s2*(3 + t1 - 3*t2) - 6*t2 + 2*Power(t2,2) + 
                       s1*(3 + s2 - 3*t1 + t2)) + 
                    s1*(-12 + 4*Power(s2,2) + t1 - Power(t1,2) - 
                       7*t2 + 3*Power(t2,2) + s2*(-4 - 3*t1 + 7*t2))))/
                (t1 - s2*t2) + 
               8*(-1 - 6*s2 - 2*Power(s2,3) + Power(s1,3)*(-1 + t1) + 
                  3*t1 + 5*s2*t1 + 9*Power(s2,2)*t1 - Power(s2,3)*t1 - 
                  3*Power(t1,2) - 4*s2*Power(t1,2) + Power(t1,3) + 
                  s2*Power(t1,3) - 3*t2 + 3*s2*t2 - 2*Power(s2,2)*t2 + 
                  2*t1*t2 + 13*s2*t1*t2 - 4*Power(s2,2)*t1*t2 + 
                  Power(t1,2)*t2 + s2*Power(t1,2)*t2 + Power(t2,2) - 
                  s2*Power(t2,2) + t1*Power(t2,2) - 
                  3*s2*t1*Power(t2,2) + Power(t2,3) + 
                  Power(s,2)*
                   (1 + s1*(-1 + t1) - 3*t1 + Power(t1,2) + t2) + 
                  Power(s1,2)*
                   (5 - 2*Power(s2,2) - 6*t1 + Power(t1,2) + 
                     s2*(2 + 3*t1 - 3*t2) + 5*t2 - Power(t2,2)) + 
                  s*(3 - 2*Power(s1,2)*(-1 + t1) - 7*t1 + 
                     4*Power(t1,2) + Power(s2,2)*(2 + t1) + 6*t2 - 
                     9*t1*t2 - Power(t1,2)*t2 + 3*Power(t2,2) + 
                     2*t1*Power(t2,2) + 
                     s2*(-2*Power(t1,2) + 3*t1*(-1 + t2) + 
                       2*(2 + t2)) + 
                     s1*(Power(s2,2) - 2*Power(t1,2) - 7*(1 + t2) + 
                        t1*(10 + t2) + s2*(-3 - 2*t1 + t2))) - 
                  s1*(4 + Power(s2,3) - 2*Power(t1,2)*(-1 + t2) - 
                     2*t2 - Power(t2,2) + Power(t2,3) + 
                     Power(s2,2)*(-7 + 2*t1 + 2*t2) + 
                     t1*(-2 - 7*t2 + 2*Power(t2,2)) + 
                     s2*(-5 - 4*Power(t1,2) - 12*t2 + 2*Power(t2,2) + 
                        t1*(4 + 6*t2)))) + 
               ((-2 + s2 + t2)*
                  (-2*Power(t1,2) + s1*(s2 + t1)*(-1 + t2) + 
                    (1 + t2)*(1 + (-1 + s - s2)*t2) + 
                    t1*(1 + t2 + 2*s2*t2 - s*(1 + t2)))*
                  (4*(-1 + s2 - t1 + t2)*
                     Power(-2*Power(t1,2) + 
                       s1*(s2 + t1)*(-1 + t2) + 
                       (1 + t2)*(1 + (-1 + s - s2)*t2) + 
                       t1*(1 + t2 + 2*s2*t2 - s*(1 + t2)),2) + 
                    12*(t1 - s2*t2)*
                     (2 - s2 - 2*t1 + s2*t1 - Power(t1,2) + 
                       s2*Power(t1,2) - s2*Power(t1,3) + 
                       Power(t1,4) + Power(s,2)*Power(t1 - t2,2) - 
                       2*t2 + s2*t2 + 4*t1*t2 - s2*t1*t2 - 
                       Power(s2,2)*t1*t2 + s2*Power(t1,2)*t2 + 
                       Power(s2,2)*Power(t1,2)*t2 - 
                       2*Power(t1,3)*t2 - s2*Power(t1,3)*t2 - 
                       Power(t2,2) - s2*Power(t2,2) + 
                       Power(s2,2)*Power(t2,2) - s2*t1*Power(t2,2) - 
                       Power(s2,2)*t1*Power(t2,2) + 
                       2*Power(t1,2)*Power(t2,2) + 
                       2*s2*Power(t1,2)*Power(t2,2) + 
                       s2*Power(t2,3) - 2*t1*Power(t2,3) - 
                       s2*t1*Power(t2,3) + Power(t2,4) + 
                       Power(s1,2)*(-1 + t2)*
                        (s2*(-1 + t1) + t1*(-t1 + t2)) - 
                       s1*(-1 + t2)*
                        (-1 + 2*s2*(-1 + t1) - 
                        Power(s2,2)*(-1 + t1) + Power(t1,3) + 
                        Power(t2,2) - Power(t1,2)*(1 + 2*t2) + 
                        t1*(1 + Power(t2,2))) + 
                       s*(-1 + t1 - Power(t1,2) + Power(t1,3) + t2 + 
                        t1*t2 - Power(t1,2)*t2 + Power(t1,3)*t2 - 
                        Power(t2,2) - t1*Power(t2,2) - 
                        2*Power(t1,2)*Power(t2,2) + Power(t2,3) + 
                        t1*Power(t2,3) - 
                        s1*(-1 + t2)*
                        (1 + s2*(-1 + t1) - Power(t1,2) + 
                        t1*(-2 + t2) + t2) - 
                        s2*(-1 + 2*t2 + Power(t2,2) + 
                        Power(t1,2)*(1 + t2) - t1*t2*(3 + t2))))))/
                ((-1 + s2 - t1 + t2)*Power(-t1 + s2*t2,3)) + 
               (4*(8 + s2 - 2*Power(s2,2) + 6*t1 + 
                    s1*(1 - 2*s2 + t1 - 2*t2) + 2*t2 - 3*s2*t2 - 
                    t1*t2 - Power(t2,2) + s*(1 + s2 - t1 + t2))*
                  (-2*Power(t1,5) + 
                    (-1 + t2)*Power(1 + t2,2)*
                     Power(1 + (-1 + s)*t2,2) - 
                    2*Power(t1,4)*(s + s1 - 2*t2 + s*t2 - s1*t2) + 
                    Power(s2,3)*
                     (Power(s1,2)*Power(-1 + t2,2) + 
                       2*s1*(t1 - t2)*(-1 + t2)*t2 + 
                       Power(t2,2)*
                        (1 + 2*Power(t1,2) + Power(t2,2) - 
                        2*t1*(1 + t2))) + 
                    Power(t1,3)*
                     (5 - 2*t2 - 5*Power(t2,2) - 
                       Power(s1,2)*(-1 + Power(t2,2)) - 
                       2*s1*(2 - 3*t2 + Power(t2,2)) - 
                       Power(s,2)*(-1 + 2*t2 + Power(t2,2)) + 
                       2*s*(-2 + 3*t2 + 3*Power(t2,2) + 
                        s1*(-2 + t2 + Power(t2,2)))) + 
                    t1*(-(Power(-1 + t2,2)*
                        (-1 + 2*(1 + s1)*t2 + (1 + 2*s1)*Power(t2,2))\
) + Power(s,2)*t2*(2 + 3*t2 - 4*Power(t2,2) - 3*Power(t2,3)) + 
                       2*s*(-1 + t2)*
                        (t2*(1 + 6*t2 + 3*Power(t2,2)) + 
                        s1*(-1 - 2*t2 + Power(t2,3)))) + 
                    Power(t1,2)*
                     (Power(s,2)*
                        (-1 - 3*t2 + 5*Power(t2,2) + 3*Power(t2,3)) + 
                       (-1 + t2)*
                        (3 + 2*t2 + 3*Power(t2,2) + 
                        2*s1*(-3 + Power(t2,2)) + 
                        Power(s1,2)*(1 + Power(t2,2))) - 
                       2*s*(-3 - 4*t2 + 5*Power(t2,2) + 
                        4*Power(t2,3) + s1*(3 - 5*t2 + 2*Power(t2,3)))\
) + s2*(1 - 2*t2 + 4*s*t2 + 2*Power(t2,2) + 2*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) + 
                       2*s*Power(t2,3) + Power(t2,4) - 
                       6*s*Power(t2,4) + Power(s,2)*Power(t2,4) - 
                       2*s*Power(t2,5) + Power(t1,4)*(2 + 4*t2) + 
                       2*Power(t1,3)*
                        (-1 - t2 - 4*Power(t2,2) + s*Power(1 + t2,2)) \
+ Power(s1,2)*t1*(-1 + t2)*(-4*t2 + t1*(3 + t2)) + 
                       Power(t1,2)*
                        (-1 - 8*t2 + 7*Power(t2,2) + 8*Power(t2,3) + 
                        Power(s,2)*(1 + Power(t2,2)) - 
                        2*s*(1 - t2 + 5*Power(t2,2) + 3*Power(t2,3))) \
+ 2*t1*t2*(-(Power(s,2)*(1 + Power(t2,2))) - 
                        2*(-2 + t2 + Power(t2,3)) + 
                        s*(-4 - t2 + 6*Power(t2,2) + 3*Power(t2,3))) \
+ 2*s1*(-1 + t2)*(-(Power(t1,3)*t2) + 
                        (1 + t2)*(-1 + t2 + s*Power(t2,2)) - 
                        Power(t1,2)*(-((-4 + t2)*t2) + s*(1 + t2)) + 
                        t1*(-2*s*(-1 + t2) + Power(1 + t2,2)))) + 
                    Power(s2,2)*
                     (Power(s1,2)*(-1 + t2)*
                        (1 + Power(t2,2) - t1*(1 + t2)) - 
                       2*s1*(-1 + t2)*
                        (-1 + Power(t1,2) + t2 + Power(t2,2) + 
                        Power(t2,3) - t1*t2*(1 + 2*t2) + 
                        s*(t1 - Power(t2,2))) + 
                       t2*(-2*Power(t1,3)*(2 + t2) + 
                        4*Power(t1,2)*(1 + t2 + Power(t2,2)) - 
                        3*t1*t2*(-1 + 2*t2 + Power(t2,2)) + 
                        t2*(-5 + 3*t2 + Power(t2,2) + Power(t2,3)) - 
                        2*s*(1 - t2 + Power(t2,2) + Power(t2,3) + 
                        Power(t1,2)*(1 + t2) - 
                        t1*(1 + t2 + 2*Power(t2,2)))))))/
                ((-1 + s2 - t1 + t2)*Power(t1 - s2*t2,2))))/
           ((s - s2 + t1)*(s - s1 + t2)) + 
          ((-1 + t1)*((4*(2*Power(t1,2) - s1*(s2 + t1)*(-1 + t2) - 
                    (1 + t2)*(1 + (-1 + s - s2)*t2) + 
                    t1*(-1 + s - t2 + s*t2 - 2*s2*t2))*
                  (-2 - 7*s2 + Power(s2,2) + Power(s2,3) - 
                    Power(s,2)*(-1 + t1) - 9*t1 - s2*t1 - 
                    2*Power(s2,2)*t1 - 5*Power(t1,2) + 
                    2*s2*Power(t1,2) + t2 - 2*s2*t2 + 
                    4*Power(s2,2)*t2 - 5*t1*t2 - 2*s2*t1*t2 + 
                    3*s2*Power(t2,2) + t1*Power(t2,2) + 
                    Power(s1,2)*(2 + s2 - 2*t1 + t2) + 
                    s*(4 - Power(s2,2) - 5*t1 + Power(t1,2) + 
                       s2*(2 + t1 - 5*t2) + 7*t2 + 2*t1*t2 - 
                       4*Power(t2,2) - s1*(3 + s2 - 3*t1 + t2)) + 
                    s1*(-8 + 3*t1 - Power(t1,2) - 10*t2 + 
                       3*Power(t2,2) + s2*(-3 + t1 + 3*t2))))/
                (t1 - s2*t2) + 
               8*(-4 + 2*s2 + 8*Power(s2,2) + 2*Power(s2,3) + 
                  Power(s1,3)*(-1 + t1) - 3*t1 - 5*s2*t1 - 
                  5*Power(s2,2)*t1 - Power(s2,3)*t1 + 6*Power(t1,2) + 
                  4*s2*Power(t1,2) + 2*Power(s2,2)*Power(t1,2) + 
                  Power(t1,3) - s2*Power(t1,3) + 4*t2 + 11*s2*t2 + 
                  8*Power(s2,2)*t2 - 13*t1*t2 - 7*s2*t1*t2 - 
                  4*Power(s2,2)*t1*t2 + 5*Power(t1,2)*t2 + 
                  3*s2*Power(t1,2)*t2 + 9*Power(t2,2) + 
                  7*s2*Power(t2,2) - 5*t1*Power(t2,2) - 
                  3*s2*t1*Power(t2,2) + 3*Power(t2,3) + 
                  Power(s,2)*
                   (s1*(-1 + t1) - 4*t1 + Power(t1,2) + 3*t2) + 
                  Power(s1,2)*
                   (4 - 9*t1 + Power(t1,2) + s2*(4 + t1 - t2) + 9*t2 - 
                     Power(t2,2)) - 
                  s1*(4 + Power(s2,3) + s2*(-1 + 2*t1*(-1 + t2)) + 
                     Power(t1,2)*(3 - 2*t2) + 9*t2 + 5*Power(t2,2) - 
                     Power(t2,3) + Power(s2,2)*(1 - 2*t1 + 2*t2) + 
                     t1*(-11 - 12*t2 + 4*Power(t2,2))) + 
                  s*(1 - 2*Power(s1,2)*(-1 + t1) - 4*t1 + 
                     Power(s2,2)*t1 + 3*Power(t1,2) + 2*t2 - 7*t1*t2 - 
                     3*Power(t1,2)*t2 + Power(t2,2) + 
                     4*t1*Power(t2,2) - 
                     s2*(1 + 2*Power(t1,2) + 4*t2 - 5*t1*t2) + 
                     s1*(-3 + Power(s2,2) + 12*t1 - 2*Power(t1,2) - 
                        11*t2 - t1*t2 + 2*Power(t2,2) + 
                        s2*(-3 - 2*t1 + 3*t2)))) + 
               ((-2 + s2 + t2)*
                  (-2*Power(t1,2) + s1*(s2 + t1)*(-1 + t2) + 
                    (1 + t2)*(1 + (-1 + s - s2)*t2) + 
                    t1*(1 + t2 + 2*s2*t2 - s*(1 + t2)))*
                  (4*(-1 + s2 - t1 + t2)*
                     Power(-2*Power(t1,2) + s1*(s2 + t1)*(-1 + t2) + 
                       (1 + t2)*(1 + (-1 + s - s2)*t2) + 
                       t1*(1 + t2 + 2*s2*t2 - s*(1 + t2)),2) + 
                    12*(t1 - s2*t2)*
                     (2 - s2 - 2*t1 + s2*t1 - Power(t1,2) + 
                       s2*Power(t1,2) - s2*Power(t1,3) + Power(t1,4) + 
                       Power(s,2)*Power(t1 - t2,2) - 2*t2 + s2*t2 + 
                       4*t1*t2 - s2*t1*t2 - Power(s2,2)*t1*t2 + 
                       s2*Power(t1,2)*t2 + 
                       Power(s2,2)*Power(t1,2)*t2 - 2*Power(t1,3)*t2 - 
                       s2*Power(t1,3)*t2 - Power(t2,2) - 
                       s2*Power(t2,2) + Power(s2,2)*Power(t2,2) - 
                       s2*t1*Power(t2,2) - 
                       Power(s2,2)*t1*Power(t2,2) + 
                       2*Power(t1,2)*Power(t2,2) + 
                       2*s2*Power(t1,2)*Power(t2,2) + s2*Power(t2,3) - 
                       2*t1*Power(t2,3) - s2*t1*Power(t2,3) + 
                       Power(t2,4) + 
                       Power(s1,2)*(-1 + t2)*
                        (s2*(-1 + t1) + t1*(-t1 + t2)) - 
                       s1*(-1 + t2)*
                        (-1 + 2*s2*(-1 + t1) - 
                        Power(s2,2)*(-1 + t1) + Power(t1,3) + 
                        Power(t2,2) - Power(t1,2)*(1 + 2*t2) + 
                        t1*(1 + Power(t2,2))) + 
                       s*(-1 + t1 - Power(t1,2) + Power(t1,3) + t2 + 
                        t1*t2 - Power(t1,2)*t2 + Power(t1,3)*t2 - 
                        Power(t2,2) - t1*Power(t2,2) - 
                        2*Power(t1,2)*Power(t2,2) + Power(t2,3) + 
                        t1*Power(t2,3) - 
                        s1*(-1 + t2)*
                        (1 + s2*(-1 + t1) - Power(t1,2) + 
                        t1*(-2 + t2) + t2) - 
                        s2*(-1 + 2*t2 + Power(t2,2) + 
                        Power(t1,2)*(1 + t2) - t1*t2*(3 + t2))))))/
                ((-1 + s2 - t1 + t2)*Power(-t1 + s2*t2,3)) + 
               (4*(8 - s2 + 6*t1 - 2*s2*t1 + 
                    s1*(1 - 2*s2 + t1 - 2*t2) + 2*t2 - s2*t2 - t1*t2 - 
                    Power(t2,2) + s*(1 + s2 - t1 + t2))*
                  (-2*Power(t1,5) + 
                    (-1 + t2)*Power(1 + t2,2)*
                     Power(1 + (-1 + s)*t2,2) - 
                    2*Power(t1,4)*(s + s1 - 2*t2 + s*t2 - s1*t2) + 
                    Power(s2,3)*
                     (Power(s1,2)*Power(-1 + t2,2) + 
                       2*s1*(t1 - t2)*(-1 + t2)*t2 + 
                       Power(t2,2)*
                        (1 + 2*Power(t1,2) + Power(t2,2) - 
                        2*t1*(1 + t2))) + 
                    Power(t1,3)*
                     (5 - 2*t2 - 5*Power(t2,2) - 
                       Power(s1,2)*(-1 + Power(t2,2)) - 
                       2*s1*(2 - 3*t2 + Power(t2,2)) - 
                       Power(s,2)*(-1 + 2*t2 + Power(t2,2)) + 
                       2*s*(-2 + 3*t2 + 3*Power(t2,2) + 
                        s1*(-2 + t2 + Power(t2,2)))) + 
                    t1*(-(Power(-1 + t2,2)*
                        (-1 + 2*(1 + s1)*t2 + (1 + 2*s1)*Power(t2,2))) \
+ Power(s,2)*t2*(2 + 3*t2 - 4*Power(t2,2) - 3*Power(t2,3)) + 
                       2*s*(-1 + t2)*
                        (t2*(1 + 6*t2 + 3*Power(t2,2)) + 
                        s1*(-1 - 2*t2 + Power(t2,3)))) + 
                    Power(t1,2)*
                     (Power(s,2)*
                        (-1 - 3*t2 + 5*Power(t2,2) + 3*Power(t2,3)) + 
                       (-1 + t2)*
                        (3 + 2*t2 + 3*Power(t2,2) + 
                        2*s1*(-3 + Power(t2,2)) + 
                        Power(s1,2)*(1 + Power(t2,2))) - 
                       2*s*(-3 - 4*t2 + 5*Power(t2,2) + 
                        4*Power(t2,3) + s1*(3 - 5*t2 + 2*Power(t2,3)))) \
+ s2*(1 - 2*t2 + 4*s*t2 + 2*Power(t2,2) + 2*s*Power(t2,2) + 
                       Power(s,2)*Power(t2,2) - 2*Power(t2,3) + 
                       2*s*Power(t2,3) + Power(t2,4) - 
                       6*s*Power(t2,4) + Power(s,2)*Power(t2,4) - 
                       2*s*Power(t2,5) + Power(t1,4)*(2 + 4*t2) + 
                       2*Power(t1,3)*
                        (-1 - t2 - 4*Power(t2,2) + s*Power(1 + t2,2)) \
+ Power(s1,2)*t1*(-1 + t2)*(-4*t2 + t1*(3 + t2)) + 
                       Power(t1,2)*
                        (-1 - 8*t2 + 7*Power(t2,2) + 8*Power(t2,3) + 
                        Power(s,2)*(1 + Power(t2,2)) - 
                        2*s*(1 - t2 + 5*Power(t2,2) + 3*Power(t2,3))) \
+ 2*t1*t2*(-(Power(s,2)*(1 + Power(t2,2))) - 
                        2*(-2 + t2 + Power(t2,3)) + 
                        s*(-4 - t2 + 6*Power(t2,2) + 3*Power(t2,3))) + 
                       2*s1*(-1 + t2)*
                        (-(Power(t1,3)*t2) + 
                        (1 + t2)*(-1 + t2 + s*Power(t2,2)) - 
                        Power(t1,2)*(-((-4 + t2)*t2) + s*(1 + t2)) + 
                        t1*(-2*s*(-1 + t2) + Power(1 + t2,2)))) + 
                    Power(s2,2)*
                     (Power(s1,2)*(-1 + t2)*
                        (1 + Power(t2,2) - t1*(1 + t2)) - 
                       2*s1*(-1 + t2)*
                        (-1 + Power(t1,2) + t2 + Power(t2,2) + 
                        Power(t2,3) - t1*t2*(1 + 2*t2) + 
                        s*(t1 - Power(t2,2))) + 
                       t2*(-2*Power(t1,3)*(2 + t2) + 
                         4*Power(t1,2)*(1 + t2 + Power(t2,2)) - 
                         3*t1*t2*(-1 + 2*t2 + Power(t2,2)) + 
                         t2*(-5 + 3*t2 + Power(t2,2) + Power(t2,3)) - 
                         2*s*
                         (1 - t2 + Power(t2,2) + Power(t2,3) + 
                         Power(t1,2)*(1 + t2) - 
                         t1*(1 + t2 + 2*Power(t2,2)))))))/
                ((-1 + s2 - t1 + t2)*Power(t1 - s2*t2,2))))/
           ((-1 + s1)*(-s + s1 - t2)))/Power(-1 + t1,2))*
     B1(t2,1 - s2 + t1 - t2,t1))/(8.*Power(Pi,2));
    return a;
};