#include <gsl/gsl_math.h>
#include <gsl/gsl_pow_int.h>
#include "nlo_functions.h"

#define Pi M_PI
#define Sqrt sqrt
#define Power gsl_pow_int
#define ln gsl_sf_log_abs
#define Log gsl_sf_log_abs

double Box3(double s,double s1,double s2,double t1,double t2){
    double a=(((-8*s*(-1 + s2 - t1 + t2)*(3 - 3*t1 - s1*(2*s1 + s2 + t1) + 
            s*(-2 + 2*s1 + t1 - t2) + t2 + (2*s1 + s2)*t2)*
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
       (4*(64*s*(s - s1 + t2)*(-1 + s2 - t1 + t2) + 
            (2*(3 - 2*Power(s1,2) - 3*t1 - s1*(s2 + t1 - 2*t2) + 
                 s*(-2 + 2*s1 + t1 - t2) + t2 + s2*t2)*
               (-21 + 2*Power(s,2) + 7*Power(s1,2) + 7*s2 + 
                 Power(s2,2) - 3*t1 + s1*(-3 + 9*s2 - 3*t1 - t2) + 
                 9*t2 - 2*s2*t2 + s*(2 - 9*s1 - 4*s2 + t2)))/
             (-1 + s1 + t1 - t2) + 
            32*(Power(s,3) + (1 + s1)*(-1 + t2)*(-1 + s2 - t1 + t2) - 
               2*Power(s,2)*(-5 + s1 + 6*s2 - 6*t1 + 5*t2) + 
               s*(Power(s1,2) - Power(t1,2) + 
                  (11 - 13*s2 - 12*t2)*t2 + 
                  2*s1*(-6 + 7*s2 - 7*t1 + 6*t2) + t1*(-1 + s2 + 14*t2)\
)) - 16*(10 + 6*Power(s,3) - 7*s2 + Power(s2,2) + 7*t1 - 2*s2*t1 + 
               Power(t1,2) + 
               Power(s,2)*(33 - 16*s1 - 49*s2 + 48*t1 - 35*t2) - 
               20*t2 + 12*s2*t2 - 13*t1*t2 + 10*Power(t2,2) + 
               Power(s1,2)*(-2 + s2 - t1 + t2) + 
               s1*(Power(s2,2) + Power(t1,2) + t1*(7 - 8*t2) - 
                  2*s2*(3 + t1 - 4*t2) + t2*(-6 + 7*t2)) + 
               s*(-10 + 10*Power(s1,2) - 2*Power(s2,2) - 11*t1 - 
                  7*Power(t1,2) + s2*(7 + 9*t1 - 60*t2) + 42*t2 + 
                  64*t1*t2 - 50*Power(t2,2) + 
                  s1*(-41 + 64*s2 - 61*t1 + 46*t2))) - 
            4*(27 + 10*Power(s,3) + 14*s2 + 7*Power(s2,2) - 11*t1 - 
               17*s2*t1 + 12*Power(t1,2) + 
               Power(s,2)*(18 - 28*s1 - 56*s2 + 56*t1 - 33*t2) - 
               7*t2 + 24*s2*t2 + 3*Power(s2,2)*t2 - 42*t1*t2 - 
               4*s2*t1*t2 + 20*Power(t2,2) + 4*s2*Power(t2,2) - 
               Power(s1,2)*(-2 + s2 - 3*t1 + 2*t2) + 
               s1*(1 + 6*Power(t1,2) - 4*s2*(-3 + t1 - 4*t2) + 
                  12*t2 + 18*Power(t2,2) - t1*(5 + 24*t2)) + 
               s*(-22 + 18*Power(s1,2) - 7*Power(s2,2) - 32*t1 - 
                  15*Power(t1,2) + 3*s2*(6 + 7*t1 - 27*t2) + 38*t2 + 
                  90*t1*t2 - 62*Power(t2,2) + 
                  s1*(-31 + 89*s2 - 87*t1 + 59*t2))) + 
            8*(22 + 13*Power(s,3) + 2*Power(s1,3) - 7*s2 + 
               8*Power(s2,2) + Power(s1,2)*(-7 + 5*s2 - t1) + 2*t1 - 
               14*s2*t1 + 8*Power(t1,2) + 
               Power(s,2)*(41 - 36*s1 - 86*s2 + 85*t1 - 55*t2) - 
               33*t2 + 36*s2*t2 + Power(s2,2)*t2 - 43*t1*t2 - 
               2*s2*t1*t2 + 25*Power(t2,2) + 2*s2*Power(t2,2) + 
               s1*(-20 + 3*Power(s2,2) + 5*Power(t1,2) + 
                  t1*(6 - 24*t2) - 6*s2*(t1 - 3*t2) + 11*t2 + 
                  19*Power(t2,2)) + 
               s*(-25 + 21*Power(s1,2) - 7*Power(s2,2) - 31*t1 - 
                  17*Power(t1,2) + s2*(21 + 23*t1 - 115*t2) + 60*t2 + 
                  125*t1*t2 - 91*Power(t2,2) + 
                  s1*(-50 + 121*s2 - 118*t1 + 83*t2))) - 
            (2*(2*s - 3*s1 - 2*s2 + t1)*
               (4*(-1 + s1 + t1 - t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 2*Power(s,3)*(-1 + s2 - t1 + t2)*
                  (2 + 2*Power(s1,2) + Power(t1,2) + 
                    2*s1*(-2 + t1 - t2) + 2*t2 + Power(t2,2) - 
                    2*t1*(1 + t2)) - 
                 4*Power(s,2)*
                  (-3 + 2*t1 + 4*Power(t1,2) - 3*Power(t1,3) + t2 - 
                    6*t1*t2 + 8*Power(t1,2)*t2 + 2*Power(t2,2) - 
                    7*t1*Power(t2,2) + 2*Power(t2,3) + 
                    Power(s2,2)*t2*(2 - t1 + t2) + 
                    2*Power(s1,3)*(-1 + s2 - t1 + t2) + 
                    Power(s1,2)*
                     (2 + Power(s2,2) - 2*Power(t1,2) + 
                       s2*(-3 + t1 - 2*t2) + t2 + 5*t1*t2 - 
                       3*Power(t2,2)) - 
                    s1*(-3 + Power(t1,3) + Power(t1,2)*(5 - 3*t2) + 
                       4*t2 + 3*t1*(-2 + t2)*t2 + Power(t2,2) - 
                       Power(t2,3) + Power(s2,2)*(1 + 2*t2) + 
                       s2*(2 - 4*t1 - Power(t1,2) - 2*t2 + 
                       Power(t2,2))) + 
                    s2*(3 + 2*Power(t2,2) + Power(t2,3) + 
                       Power(t1,2)*(2 + t2) - 
                       t1*(5 + 4*t2 + 2*Power(t2,2)))) + 
                 2*s*(-5 + t1 + 13*Power(t1,2) - 9*Power(t1,3) - t2 - 
                    18*t1*t2 + 19*Power(t1,2)*t2 + Power(t2,2) + 
                    Power(s2,3)*Power(t2,2) - 15*t1*Power(t2,2) + 
                    5*Power(t2,3) + 
                    2*Power(s1,4)*(-1 + s2 - t1 + t2) + 
                    Power(s2,2)*t2*
                     (4 - t2 + Power(t2,2) - t1*(4 + t2)) + 
                    2*Power(s1,3)*
                     (Power(s2,2) - Power(t1,2) - 2*(-1 + t2)*t2 - 
                       s2*(1 + t2) + t1*(-1 + 3*t2)) + 
                    s2*(5 + 6*t2 + 17*Power(t2,2) + 4*Power(t2,3) + 
                       Power(t1,2)*(5 + 8*t2) - 
                       2*t1*(5 + 7*t2 + 6*Power(t2,2))) - 
                    2*s1*(3*Power(t1,3) + Power(t1,2)*(4 - 9*t2) + 
                       t2 + Power(s2,3)*t2 - Power(t2,3) + 
                       t1*(-7 - 6*t2 + 7*Power(t2,2)) - 
                       Power(s2,2)*(-2 + t2 + t1*(2 + t2)) - 
                       s2*(2 + Power(t1,2) - 5*t2 - 4*Power(t2,2) + 
                        Power(t2,3) - t1*(3 + Power(t2,2)))) + 
                    Power(s1,2)*
                     (Power(s2,3) - Power(t1,3) + 
                       Power(t1,2)*(-11 + 3*t2) - 
                       Power(s2,2)*(1 + t1 + 3*t2) - 
                       2*t1*(4 - 7*t2 + 2*Power(t2,2)) + 
                       2*(3 - 2*t2 - 2*Power(t2,2) + Power(t2,3)) + 
                       s2*(Power(t1,2) + 2*t1*(5 + t2) - 
                        2*(3 - 3*t2 + Power(t2,2)))))))/
             (s*Power(-1 + s1 + t1 - t2,2)*(-1 + s2 - t1 + t2)) + 
            ((3 - 2*Power(s1,2) - 3*t1 - s1*(s2 + t1 - 2*t2) + 
                 s*(-2 + 2*s1 + t1 - t2) + t2 + s2*t2)*
               (12*(-1 + s1 + t1 - t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 4*Power(s,3)*(-1 + s2 - t1 + t2)*
                  (-3*(-1 + s1)*(-1 + s1 + t1 - t2) + 
                    Power(2 - 2*s1 - t1 + t2,2)) - 
                 2*Power(s,2)*
                  (-4*(-2 + 2*s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                     (3 - 2*Power(s1,2) - 3*t1 - 
                       s1*(s2 + t1 - 2*t2) + t2 + s2*t2) - 
                    6*(-1 + s1 + t1 - t2)*
                     (3 + s2*(-3 + t1) + 2*t1 - 3*t2 - t1*t2 + 
                       Power(t2,2) + 
                       2*Power(s1,2)*(-1 + s2 - t1 + t2) + 
                       s1*(s2 - t2)*(-1 + s2 - t1 + t2))) + 
                 s*(4*(-1 + s2 - t1 + t2)*
                     Power(3 - 2*Power(s1,2) - 3*t1 - 
                       s1*(s2 + t1 - 2*t2) + t2 + s2*t2,2) + 
                    12*(-1 + s1 + t1 - t2)*
                     (Power(s2,2)*t2 - 
                       Power(s1,3)*(-1 + s2 - t1 + t2) - 
                       Power(s1,2)*(1 + s2 - t2)*(-1 + s2 - t1 + t2) + 
                       s2*(2 + t1*(-2 + t2) - 5*t2 - Power(t2,2)) + 
                       2*(-Power(-1 + t2,2) + t1*(1 + t2)) + 
                       s1*(-2 + Power(s2,2)*(-1 + t2) + t2 + 
                        Power(t2,2) - t1*(4 + t2) + 
                        s2*(3 + t1 + t2 - t1*t2 + Power(t2,2)))))))/
             (s*Power(-1 + s1 + t1 - t2,3)*(-1 + s2 - t1 + t2))))/
        ((s - s2 + t1)*(s - s1 + t2)) + 
       (4*(64*s*(s1 - s2 + t1 - t2)*(-1 + s2 - t1 + t2) + 
            (2*(-17 + 4*Power(s,2) + 9*Power(s1,2) - 2*s2 + 
                 3*Power(s2,2) + s1*(-4 + 13*s2 - 5*t1) + 5*t1 - 
                 2*s2*t1 + s*(2 - 13*s1 - 8*s2 + 3*t1) + t2 - s2*t2)*
               (3 - 2*Power(s1,2) - 3*t1 - s1*(s2 + t1 - 2*t2) + 
                 s*(-2 + 2*s1 + t1 - t2) + t2 + s2*t2))/
             (-1 + s1 + t1 - t2) + 
            32*(Power(s,2)*(-1 + s1) - 
               (-1 + s2 - t1 + t2)*
                (1 - 2*Power(s1,2) + s2 - 2*t1 - t2 + s2*t2 + 
                  s1*(2 - 2*t1 + 2*t2)) - 
               s*(3 + Power(s1,2) - 11*Power(s2,2) - 8*t1 - 
                  13*Power(t1,2) + 6*s2*(1 + 4*t1 - 4*t2) + 9*t2 + 
                  26*t1*t2 - 13*Power(t2,2) + 
                  s1*(-13 + 14*s2 - 14*t1 + 13*t2))) - 
            4*(16 + 3*Power(s,3) + Power(s1,3) + 16*s2 - 
               25*Power(s2,2) + 2*Power(s2,3) - 21*t1 + 80*s2*t1 - 
               4*Power(s2,2)*t1 - 57*Power(t1,2) + 2*s2*Power(t1,2) - 
               27*t2 - 19*s2*t2 - 17*Power(s2,2)*t2 + 64*t1*t2 + 
               18*s2*t1*t2 + Power(t2,2) - 20*s2*Power(t2,2) + 
               Power(s,2)*(-13 + 13*s1 + 9*s2 - 7*t1 + 10*t2) + 
               Power(s1,2)*(-14 + 41*s2 - 42*t1 + 40*t2) + 
               s1*(35 + 9*Power(s2,2) - 26*Power(t1,2) + 
                  s2*(7 + 15*t1 - 30*t2) + 23*t2 - 40*Power(t2,2) + 
                  6*t1*(-5 + 11*t2)) + 
               s*(-43 - 17*Power(s1,2) + 38*Power(s2,2) + 2*t1 + 
                  72*Power(t1,2) + s1*(53 - 114*s2 + 109*t1 - 96*t2) - 
                  5*t2 - 154*t1*t2 + 82*Power(t2,2) + 
                  s2*(11 - 109*t1 + 121*t2))) + 
            8*(14 + Power(s,3) + 5*Power(s1,3) + 10*s2 - 
               34*Power(s2,2) + 3*Power(s2,3) - 27*t1 + 98*s2*t1 - 
               6*Power(s2,2)*t1 - 63*Power(t1,2) + 3*s2*Power(t1,2) - 
               33*t2 - 16*s2*t2 - 19*Power(s2,2)*t2 + 64*t1*t2 + 
               19*s2*t1*t2 + 7*Power(t2,2) - 22*s2*Power(t2,2) + 
               Power(s,2)*(-21 + 22*s1 + 14*s2 - 10*t1 + 11*t2) + 
               Power(s1,2)*(-31 + 54*s2 - 49*t1 + 44*t2) + 
               s1*(9 + 13*Power(s2,2) - 31*Power(t1,2) + 
                  3*s2*(-5 + 6*t1 - 12*t2) + 20*t2 - 44*Power(t2,2) + 
                  t1*(-14 + 75*t2)) + 
               s*(-44 - 28*Power(s1,2) + 62*Power(s2,2) + 11*t1 + 
                  104*Power(t1,2) + 
                  s1*(100 - 160*s2 + 149*t1 - 128*t2) - 24*t2 - 
                  219*t1*t2 + 115*Power(t2,2) + 
                  s2*(13 - 166*t1 + 181*t2))) - 
            16*(7 + Power(s1,3) + 4*s2 - 13*Power(s2,2) + 
               Power(s2,3) - 11*t1 + 35*s2*t1 - 2*Power(s2,2)*t1 - 
               22*Power(t1,2) + s2*Power(t1,2) - 14*t2 - 2*s2*t2 - 
               7*Power(s2,2)*t2 + 19*t1*t2 + 7*s2*t1*t2 + 
               5*Power(t2,2) - 8*s2*Power(t2,2) + 
               Power(s,2)*(-9 + 8*s1 + 4*s2 - 3*t1 + 3*t2) + 
               Power(s1,2)*(-15 + 18*s2 - 17*t1 + 16*t2) + 
               s*(-21 - 9*Power(s1,2) + 40*Power(s2,2) + 17*t1 + 
                  56*Power(t1,2) + s1*(57 - 72*s2 + 70*t1 - 62*t2) - 
                  4*s2*(1 + 24*t1 - 25*t2) - 24*t2 - 115*t1*t2 + 
                  59*Power(t2,2)) + 
               s1*(3*Power(s2,2) - 13*Power(t1,2) + 
                  2*s2*(-6 + 5*t1 - 7*t2) + t1*(-2 + 29*t2) + 
                  4*(2 + t2 - 4*Power(t2,2)))) - 
            ((2 + 5*s - 7*s1 - 5*s2 + 2*t1)*
               (4*(-1 + s1 + t1 - t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 2*Power(s,3)*(-1 + s2 - t1 + t2)*
                  (2 + 2*Power(s1,2) + Power(t1,2) + 
                    2*s1*(-2 + t1 - t2) + 2*t2 + Power(t2,2) - 
                    2*t1*(1 + t2)) - 
                 4*Power(s,2)*
                  (-3 + 2*t1 + 4*Power(t1,2) - 3*Power(t1,3) + t2 - 
                    6*t1*t2 + 8*Power(t1,2)*t2 + 2*Power(t2,2) - 
                    7*t1*Power(t2,2) + 2*Power(t2,3) + 
                    Power(s2,2)*t2*(2 - t1 + t2) + 
                    2*Power(s1,3)*(-1 + s2 - t1 + t2) + 
                    Power(s1,2)*
                     (2 + Power(s2,2) - 2*Power(t1,2) + 
                       s2*(-3 + t1 - 2*t2) + t2 + 5*t1*t2 - 
                       3*Power(t2,2)) - 
                    s1*(-3 + Power(t1,3) + Power(t1,2)*(5 - 3*t2) + 
                       4*t2 + 3*t1*(-2 + t2)*t2 + Power(t2,2) - 
                       Power(t2,3) + Power(s2,2)*(1 + 2*t2) + 
                       s2*(2 - 4*t1 - Power(t1,2) - 2*t2 + 
                       Power(t2,2))) + 
                    s2*(3 + 2*Power(t2,2) + Power(t2,3) + 
                       Power(t1,2)*(2 + t2) - 
                       t1*(5 + 4*t2 + 2*Power(t2,2)))) + 
                 2*s*(-5 + t1 + 13*Power(t1,2) - 9*Power(t1,3) - t2 - 
                    18*t1*t2 + 19*Power(t1,2)*t2 + Power(t2,2) + 
                    Power(s2,3)*Power(t2,2) - 15*t1*Power(t2,2) + 
                    5*Power(t2,3) + 
                    2*Power(s1,4)*(-1 + s2 - t1 + t2) + 
                    Power(s2,2)*t2*
                     (4 - t2 + Power(t2,2) - t1*(4 + t2)) + 
                    2*Power(s1,3)*
                     (Power(s2,2) - Power(t1,2) - 2*(-1 + t2)*t2 - 
                       s2*(1 + t2) + t1*(-1 + 3*t2)) + 
                    s2*(5 + 6*t2 + 17*Power(t2,2) + 4*Power(t2,3) + 
                       Power(t1,2)*(5 + 8*t2) - 
                       2*t1*(5 + 7*t2 + 6*Power(t2,2))) - 
                    2*s1*(3*Power(t1,3) + Power(t1,2)*(4 - 9*t2) + 
                       t2 + Power(s2,3)*t2 - Power(t2,3) + 
                       t1*(-7 - 6*t2 + 7*Power(t2,2)) - 
                       Power(s2,2)*(-2 + t2 + t1*(2 + t2)) - 
                       s2*(2 + Power(t1,2) - 5*t2 - 4*Power(t2,2) + 
                        Power(t2,3) - t1*(3 + Power(t2,2)))) + 
                    Power(s1,2)*
                     (Power(s2,3) - Power(t1,3) + 
                       Power(t1,2)*(-11 + 3*t2) - 
                       Power(s2,2)*(1 + t1 + 3*t2) - 
                       2*t1*(4 - 7*t2 + 2*Power(t2,2)) + 
                       2*(3 - 2*t2 - 2*Power(t2,2) + Power(t2,3)) + 
                       s2*(Power(t1,2) + 2*t1*(5 + t2) - 
                        2*(3 - 3*t2 + Power(t2,2)))))))/
             (s*Power(-1 + s1 + t1 - t2,2)*(-1 + s2 - t1 + t2)) + 
            ((3 - 2*Power(s1,2) - 3*t1 - s1*(s2 + t1 - 2*t2) + 
                 s*(-2 + 2*s1 + t1 - t2) + t2 + s2*t2)*
               (12*(-1 + s1 + t1 - t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 4*Power(s,3)*(-1 + s2 - t1 + t2)*
                  (-3*(-1 + s1)*(-1 + s1 + t1 - t2) + 
                    Power(2 - 2*s1 - t1 + t2,2)) - 
                 2*Power(s,2)*
                  (-4*(-2 + 2*s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                     (3 - 2*Power(s1,2) - 3*t1 - 
                       s1*(s2 + t1 - 2*t2) + t2 + s2*t2) - 
                    6*(-1 + s1 + t1 - t2)*
                     (3 + s2*(-3 + t1) + 2*t1 - 3*t2 - t1*t2 + 
                       Power(t2,2) + 
                       2*Power(s1,2)*(-1 + s2 - t1 + t2) + 
                       s1*(s2 - t2)*(-1 + s2 - t1 + t2))) + 
                 s*(4*(-1 + s2 - t1 + t2)*
                     Power(3 - 2*Power(s1,2) - 3*t1 - 
                       s1*(s2 + t1 - 2*t2) + t2 + s2*t2,2) + 
                    12*(-1 + s1 + t1 - t2)*
                     (Power(s2,2)*t2 - 
                       Power(s1,3)*(-1 + s2 - t1 + t2) - 
                       Power(s1,2)*(1 + s2 - t2)*(-1 + s2 - t1 + t2) + 
                       s2*(2 + t1*(-2 + t2) - 5*t2 - Power(t2,2)) + 
                       2*(-Power(-1 + t2,2) + t1*(1 + t2)) + 
                       s1*(-2 + Power(s2,2)*(-1 + t2) + t2 + 
                        Power(t2,2) - t1*(4 + t2) + 
                        s2*(3 + t1 + t2 - t1*t2 + Power(t2,2)))))))/
             (s*Power(-1 + s1 + t1 - t2,3)*(-1 + s2 - t1 + t2))))/
        ((-1 + s2)*(-s + s2 - t1)) + 
       (4*(64*s*(s1 - t2)*(-1 + s2 - t1 + t2) + 
            (2*(3 - 2*Power(s1,2) - 3*t1 - s1*(s2 + t1 - 2*t2) + 
                 s*(-2 + 2*s1 + t1 - t2) + t2 + s2*t2)*
               (-18 + 9*Power(s,2) + s1 + 9*Power(s1,2) + 7*s2 + 
                 8*s1*s2 + Power(s2,2) - 2*t1 + 2*t2 - 5*s1*t2 - 
                 3*s2*t2 + s*(-2 - 18*s1 - 8*s2 + t1 + 5*t2)))/
             (-1 + s1 + t1 - t2) + 
            32*(Power(s,2)*(-2 + 2*s1 + s2 - t2) + 
               (1 + s1)*(s1 - t2)*(-1 + s2 - t1 + t2) + 
               s*(-1 - 2*Power(s1,2) - s2*t1 + Power(t1,2) + 
                  s1*(15 - 15*s2 + 13*t1 - 11*t2) - 11*t2 + 14*s2*t2 - 
                  14*t1*t2 + 12*Power(t2,2))) + 
            16*(2 + Power(s,3) - Power(s1,3) - 2*s2 + t1 - s2*t1 + 
               Power(t1,2) + Power(s1,2)*(5 - 8*s2 + 7*t1 - 6*t2) - 
               10*t2 + 12*s2*t2 - 13*t1*t2 + 10*Power(t2,2) + 
               Power(s,2)*(15 - 19*s1 - 8*s2 + 7*t2) + 
               s1*(14 + Power(t1,2) + t1*(11 - 8*t2) - 
                  s2*(12 + t1 - 8*t2) - 17*t2 + 7*Power(t2,2)) + 
               s*(5 + 19*Power(s1,2) - t1 - 7*Power(t1,2) + 
                  s2*(3 + 7*t1 - 64*t2) + 37*t2 + 64*t1*t2 - 
                  50*Power(t2,2) + s1*(-68 + 75*s2 - 58*t1 + 42*t2))) - 
            8*(5 + 5*Power(s,3) - 5*Power(s1,3) - 8*s2 - 
               2*Power(s2,2) + 5*t1 - 7*s2*t1 + 8*Power(t1,2) + 
               Power(s1,2)*(-23*s2 + 18*t1 - 13*t2) - 18*t2 + 
               35*s2*t2 + 2*Power(s2,2)*t2 - 43*t1*t2 - 2*s2*t1*t2 + 
               25*Power(t2,2) + 2*s2*Power(t2,2) + 
               Power(s,2)*(30 - 59*s1 - 23*s2 + t1 + 18*t2) - 
               s1*(-36 + Power(s2,2) - 5*Power(t1,2) + 
                  s2*(44 + 4*t1 - 23*t2) + 40*t2 - 19*Power(t2,2) + 
                  t1*(-35 + 24*t2)) + 
               s*(10 + 59*Power(s1,2) + Power(s2,2) - 6*t1 - 
                  18*Power(t1,2) + s2*(20 + 17*t1 - 127*t2) + 48*t2 + 
                  127*t1*t2 - 92*Power(t2,2) + 
                  s1*(-108 + 161*s2 - 112*t1 + 66*t2))) + 
            4*(-22 + Power(s,3) - Power(s1,3) - 15*s2 - Power(s2,2) + 
               6*t1 - 9*s2*t1 + 12*Power(t1,2) + 12*t2 + 27*s2*t2 + 
               3*Power(s2,2)*t2 - 42*t1*t2 - 4*s2*t1*t2 + 
               20*Power(t2,2) + 4*s2*Power(t2,2) + 
               Power(s,2)*(-43*s1 + 3*(6 - 5*s2 + 4*t2)) - 
               Power(s1,2)*(15*s2 + 2*(6 - 8*t1 + 7*t2)) + 
               s*(40 + 43*Power(s1,2) - 4*t1 - 16*Power(t1,2) + 
                  s2*(18 + 15*t1 - 89*t2) + 17*t2 + 93*t1*t2 - 
                  64*Power(t2,2) + s1*(-49 + 112*s2 - 78*t1 + 44*t2)) \
+ s1*(s2*(-39 - 4*t1 + 17*t2) + 
                  3*(-3 + 2*Power(t1,2) + t1*(11 - 8*t2) - 9*t2 + 
                     6*Power(t2,2)))) + 
            ((-7*s + 7*s1 + 3*s2 - 2*t2)*
               (4*(-1 + s1 + t1 - t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 2*Power(s,3)*(-1 + s2 - t1 + t2)*
                  (2 + 2*Power(s1,2) + Power(t1,2) + 
                    2*s1*(-2 + t1 - t2) + 2*t2 + Power(t2,2) - 
                    2*t1*(1 + t2)) - 
                 4*Power(s,2)*
                  (-3 + 2*t1 + 4*Power(t1,2) - 3*Power(t1,3) + t2 - 
                    6*t1*t2 + 8*Power(t1,2)*t2 + 2*Power(t2,2) - 
                    7*t1*Power(t2,2) + 2*Power(t2,3) + 
                    Power(s2,2)*t2*(2 - t1 + t2) + 
                    2*Power(s1,3)*(-1 + s2 - t1 + t2) + 
                    Power(s1,2)*
                     (2 + Power(s2,2) - 2*Power(t1,2) + 
                       s2*(-3 + t1 - 2*t2) + t2 + 5*t1*t2 - 
                       3*Power(t2,2)) - 
                    s1*(-3 + Power(t1,3) + Power(t1,2)*(5 - 3*t2) + 
                       4*t2 + 3*t1*(-2 + t2)*t2 + Power(t2,2) - 
                       Power(t2,3) + Power(s2,2)*(1 + 2*t2) + 
                       s2*(2 - 4*t1 - Power(t1,2) - 2*t2 + 
                       Power(t2,2))) + 
                    s2*(3 + 2*Power(t2,2) + Power(t2,3) + 
                       Power(t1,2)*(2 + t2) - 
                       t1*(5 + 4*t2 + 2*Power(t2,2)))) + 
                 2*s*(-5 + t1 + 13*Power(t1,2) - 9*Power(t1,3) - t2 - 
                    18*t1*t2 + 19*Power(t1,2)*t2 + Power(t2,2) + 
                    Power(s2,3)*Power(t2,2) - 15*t1*Power(t2,2) + 
                    5*Power(t2,3) + 
                    2*Power(s1,4)*(-1 + s2 - t1 + t2) + 
                    Power(s2,2)*t2*
                     (4 - t2 + Power(t2,2) - t1*(4 + t2)) + 
                    2*Power(s1,3)*
                     (Power(s2,2) - Power(t1,2) - 2*(-1 + t2)*t2 - 
                       s2*(1 + t2) + t1*(-1 + 3*t2)) + 
                    s2*(5 + 6*t2 + 17*Power(t2,2) + 4*Power(t2,3) + 
                       Power(t1,2)*(5 + 8*t2) - 
                       2*t1*(5 + 7*t2 + 6*Power(t2,2))) - 
                    2*s1*(3*Power(t1,3) + Power(t1,2)*(4 - 9*t2) + 
                       t2 + Power(s2,3)*t2 - Power(t2,3) + 
                       t1*(-7 - 6*t2 + 7*Power(t2,2)) - 
                       Power(s2,2)*(-2 + t2 + t1*(2 + t2)) - 
                       s2*(2 + Power(t1,2) - 5*t2 - 4*Power(t2,2) + 
                        Power(t2,3) - t1*(3 + Power(t2,2)))) + 
                    Power(s1,2)*
                     (Power(s2,3) - Power(t1,3) + 
                       Power(t1,2)*(-11 + 3*t2) - 
                       Power(s2,2)*(1 + t1 + 3*t2) - 
                       2*t1*(4 - 7*t2 + 2*Power(t2,2)) + 
                       2*(3 - 2*t2 - 2*Power(t2,2) + Power(t2,3)) + 
                       s2*(Power(t1,2) + 2*t1*(5 + t2) - 
                        2*(3 - 3*t2 + Power(t2,2)))))))/
             (s*Power(-1 + s1 + t1 - t2,2)*(-1 + s2 - t1 + t2)) + 
            ((3 - 2*Power(s1,2) - 3*t1 - s1*(s2 + t1 - 2*t2) + 
                 s*(-2 + 2*s1 + t1 - t2) + t2 + s2*t2)*
               (12*(-1 + s1 + t1 - t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 4*Power(s,3)*(-1 + s2 - t1 + t2)*
                  (-3*(-1 + s1)*(-1 + s1 + t1 - t2) + 
                    Power(2 - 2*s1 - t1 + t2,2)) - 
                 2*Power(s,2)*
                  (-4*(-2 + 2*s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                     (3 - 2*Power(s1,2) - 3*t1 - 
                       s1*(s2 + t1 - 2*t2) + t2 + s2*t2) - 
                    6*(-1 + s1 + t1 - t2)*
                     (3 + s2*(-3 + t1) + 2*t1 - 3*t2 - t1*t2 + 
                       Power(t2,2) + 
                       2*Power(s1,2)*(-1 + s2 - t1 + t2) + 
                       s1*(s2 - t2)*(-1 + s2 - t1 + t2))) + 
                 s*(4*(-1 + s2 - t1 + t2)*
                     Power(3 - 2*Power(s1,2) - 3*t1 - 
                       s1*(s2 + t1 - 2*t2) + t2 + s2*t2,2) + 
                    12*(-1 + s1 + t1 - t2)*
                     (Power(s2,2)*t2 - 
                       Power(s1,3)*(-1 + s2 - t1 + t2) - 
                       Power(s1,2)*(1 + s2 - t2)*(-1 + s2 - t1 + t2) + 
                       s2*(2 + t1*(-2 + t2) - 5*t2 - Power(t2,2)) + 
                       2*(-Power(-1 + t2,2) + t1*(1 + t2)) + 
                       s1*(-2 + Power(s2,2)*(-1 + t2) + t2 + 
                        Power(t2,2) - t1*(4 + t2) + 
                        s2*(3 + t1 + t2 - t1*t2 + Power(t2,2)))))))/
             (s*Power(-1 + s1 + t1 - t2,3)*(-1 + s2 - t1 + t2))))/
        ((-1 + t1)*(-1 + t2)) + 
       (4*(-64*s*(-1 + s1)*(-1 + s2 - t1 + t2) + 
            (2*(-13 + 6*Power(s,2) + 7*Power(s1,2) + 5*s2 + 
                 Power(s2,2) - t1 - s*(10 + 13*s1 + 4*s2 - 5*t2) + 
                 s1*(1 + 5*s2 + t1 - 5*t2) + t2 - 2*s2*t2)*
               (3 - 2*Power(s1,2) - 3*t1 - s1*(s2 + t1 - 2*t2) + 
                 s*(-2 + 2*s1 + t1 - t2) + t2 + s2*t2))/
             (-1 + s1 + t1 - t2) - 
            32*(Power(s,3) - Power(s,2)*(3 + s1 - 2*t2) + 
               (-1 + Power(s1,2))*(-1 + s2 - t1 + t2) + 
               s*(-11 - 12*t1 + Power(t1,2) + 
                  s1*(14 - 13*s2 + 13*t1 - 14*t2) + 9*t2 - 2*t1*t2 + 
                  2*Power(t2,2) + s2*(13 - t1 + t2))) + 
            16*(15 + 8*Power(s,3) - 15*s2 + Power(s2,2) + 13*t1 + 
               s2*t1 - 2*Power(t1,2) - 13*t2 + 5*s2*t2 - 
               2*Power(s2,2)*t2 + 2*s2*t1*t2 - 2*Power(t2,2) - 
               2*s2*Power(t2,2) + 
               Power(s1,2)*(-9 + 10*s2 - 10*t1 + 10*t2) + 
               Power(s,2)*(-32 - 9*s1 + 15*t2) + 
               s1*(-8 + Power(s2,2) - 2*Power(t1,2) + 
                  s2*(1 + t1 - 3*t2) + 11*t2 - 4*Power(t2,2) + 
                  t1*(-5 + 6*t2)) + 
               s*(-30 + Power(s1,2) - 51*t1 + 7*Power(t1,2) + 
                  s1*(74 - 60*s2 + 60*t1 - 68*t2) + 15*t2 - 16*t1*t2 + 
                  16*Power(t2,2) + s2*(58 - 7*t1 + 9*t2))) + 
            4*(-38 + 16*Power(s,3) - 23*s2 + 36*t1 + 4*s2*t1 - 
               6*Power(t1,2) + 
               Power(s1,2)*(17*s2 - 19*t1 + 18*(-1 + t2)) + 11*t2 + 
               s2*t2 - 3*Power(s2,2)*t2 + t1*t2 + 4*s2*t1*t2 - 
               13*Power(t2,2) - 4*s2*Power(t2,2) + 
               Power(s,2)*(-97 - 14*s1 + 2*s2 + 29*t2) + 
               s1*(-67 - 4*Power(t1,2) + s2*(15 + 2*t1 - 6*t2) + 
                  44*t2 - 8*Power(t2,2) + 4*t1*(-5 + 3*t2)) + 
               s*(86 - 2*Power(s1,2) + Power(s2,2) - 58*t1 + 
                  9*Power(t1,2) + s1*(139 - 81*s2 + 77*t1 - 91*t2) - 
                  54*t2 - 28*t1*t2 + 32*Power(t2,2) + 
                  s2*(54 - 9*t1 + 21*t2))) - 
            8*(20 + 21*Power(s,3) - 2*Power(s1,3) - 36*s2 + 
               Power(s2,2) + 40*t1 + 5*s2*t1 - 8*Power(t1,2) - 22*t2 + 
               10*s2*t2 - 5*Power(s2,2)*t2 + 2*t1*t2 + 6*s2*t1*t2 - 
               12*Power(t2,2) - 6*s2*Power(t2,2) + 
               Power(s1,2)*(-20 + 22*s2 - 26*t1 + 27*t2) + 
               Power(s,2)*(-105 - 25*s1 + 38*t2) + 
               s1*(-41 + Power(s2,2) - 6*Power(t1,2) + 
                  s2*(9 + 3*t1 - 7*t2) + 43*t2 - 12*Power(t2,2) + 
                  2*t1*(-10 + 9*t2)) + 
               s*(13 + 6*Power(s1,2) + Power(s2,2) - 92*t1 + 
                  15*Power(t1,2) + 
                  s1*(172 - 113*s2 + 115*t1 - 137*t2) - 26*t2 - 
                  38*t1*t2 + 40*Power(t2,2) + s2*(104 - 15*t1 + 23*t2))\
) + (2*(1 - 3*s + 3*s1 + s2 - t2)*
               (4*(-1 + s1 + t1 - t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 2*Power(s,3)*(-1 + s2 - t1 + t2)*
                  (2 + 2*Power(s1,2) + Power(t1,2) + 
                    2*s1*(-2 + t1 - t2) + 2*t2 + Power(t2,2) - 
                    2*t1*(1 + t2)) - 
                 4*Power(s,2)*
                  (-3 + 2*t1 + 4*Power(t1,2) - 3*Power(t1,3) + t2 - 
                    6*t1*t2 + 8*Power(t1,2)*t2 + 2*Power(t2,2) - 
                    7*t1*Power(t2,2) + 2*Power(t2,3) + 
                    Power(s2,2)*t2*(2 - t1 + t2) + 
                    2*Power(s1,3)*(-1 + s2 - t1 + t2) + 
                    Power(s1,2)*
                     (2 + Power(s2,2) - 2*Power(t1,2) + 
                       s2*(-3 + t1 - 2*t2) + t2 + 5*t1*t2 - 
                       3*Power(t2,2)) - 
                    s1*(-3 + Power(t1,3) + Power(t1,2)*(5 - 3*t2) + 
                       4*t2 + 3*t1*(-2 + t2)*t2 + Power(t2,2) - 
                       Power(t2,3) + Power(s2,2)*(1 + 2*t2) + 
                       s2*(2 - 4*t1 - Power(t1,2) - 2*t2 + 
                       Power(t2,2))) + 
                    s2*(3 + 2*Power(t2,2) + Power(t2,3) + 
                       Power(t1,2)*(2 + t2) - 
                       t1*(5 + 4*t2 + 2*Power(t2,2)))) + 
                 2*s*(-5 + t1 + 13*Power(t1,2) - 9*Power(t1,3) - t2 - 
                    18*t1*t2 + 19*Power(t1,2)*t2 + Power(t2,2) + 
                    Power(s2,3)*Power(t2,2) - 15*t1*Power(t2,2) + 
                    5*Power(t2,3) + 
                    2*Power(s1,4)*(-1 + s2 - t1 + t2) + 
                    Power(s2,2)*t2*
                     (4 - t2 + Power(t2,2) - t1*(4 + t2)) + 
                    2*Power(s1,3)*
                     (Power(s2,2) - Power(t1,2) - 2*(-1 + t2)*t2 - 
                       s2*(1 + t2) + t1*(-1 + 3*t2)) + 
                    s2*(5 + 6*t2 + 17*Power(t2,2) + 4*Power(t2,3) + 
                       Power(t1,2)*(5 + 8*t2) - 
                       2*t1*(5 + 7*t2 + 6*Power(t2,2))) - 
                    2*s1*(3*Power(t1,3) + Power(t1,2)*(4 - 9*t2) + 
                       t2 + Power(s2,3)*t2 - Power(t2,3) + 
                       t1*(-7 - 6*t2 + 7*Power(t2,2)) - 
                       Power(s2,2)*(-2 + t2 + t1*(2 + t2)) - 
                       s2*(2 + Power(t1,2) - 5*t2 - 4*Power(t2,2) + 
                        Power(t2,3) - t1*(3 + Power(t2,2)))) + 
                    Power(s1,2)*
                     (Power(s2,3) - Power(t1,3) + 
                       Power(t1,2)*(-11 + 3*t2) - 
                       Power(s2,2)*(1 + t1 + 3*t2) - 
                       2*t1*(4 - 7*t2 + 2*Power(t2,2)) + 
                       2*(3 - 2*t2 - 2*Power(t2,2) + Power(t2,3)) + 
                       s2*(Power(t1,2) + 2*t1*(5 + t2) - 
                        2*(3 - 3*t2 + Power(t2,2)))))))/
             (s*Power(-1 + s1 + t1 - t2,2)*(-1 + s2 - t1 + t2)) + 
            ((3 - 2*Power(s1,2) - 3*t1 - s1*(s2 + t1 - 2*t2) + 
                 s*(-2 + 2*s1 + t1 - t2) + t2 + s2*t2)*
               (12*(-1 + s1 + t1 - t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 4*Power(s,3)*(-1 + s2 - t1 + t2)*
                  (-3*(-1 + s1)*(-1 + s1 + t1 - t2) + 
                    Power(2 - 2*s1 - t1 + t2,2)) - 
                 2*Power(s,2)*
                  (-4*(-2 + 2*s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                     (3 - 2*Power(s1,2) - 3*t1 - 
                       s1*(s2 + t1 - 2*t2) + t2 + s2*t2) - 
                    6*(-1 + s1 + t1 - t2)*
                     (3 + s2*(-3 + t1) + 2*t1 - 3*t2 - t1*t2 + 
                       Power(t2,2) + 
                       2*Power(s1,2)*(-1 + s2 - t1 + t2) + 
                       s1*(s2 - t2)*(-1 + s2 - t1 + t2))) + 
                 s*(4*(-1 + s2 - t1 + t2)*
                     Power(3 - 2*Power(s1,2) - 3*t1 - 
                       s1*(s2 + t1 - 2*t2) + t2 + s2*t2,2) + 
                    12*(-1 + s1 + t1 - t2)*
                     (Power(s2,2)*t2 - 
                       Power(s1,3)*(-1 + s2 - t1 + t2) - 
                       Power(s1,2)*(1 + s2 - t2)*(-1 + s2 - t1 + t2) + 
                       s2*(2 + t1*(-2 + t2) - 5*t2 - Power(t2,2)) + 
                       2*(-Power(-1 + t2,2) + t1*(1 + t2)) + 
                       s1*(-2 + Power(s2,2)*(-1 + t2) + t2 + 
                        Power(t2,2) - t1*(4 + t2) + 
                        s2*(3 + t1 + t2 - t1*t2 + Power(t2,2)))))))/
             (s*Power(-1 + s1 + t1 - t2,3)*(-1 + s2 - t1 + t2))))/
        ((-1 + s1)*(-1 + t2)) + 
       (4*(64*s*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2) + 
            (4*(3 - 2*Power(s1,2) - 3*t1 - s1*(s2 + t1 - 2*t2) + 
                 s*(-2 + 2*s1 + t1 - t2) + t2 + s2*t2)*
               (-8 + 3*Power(s,2) + 6*Power(s1,2) + Power(s2,2) + 
                 t1 + s1*(-2 + 6*s2 - 3*t2) + t2 - 2*s2*t2 + 
                 s*(2 - 9*s1 - 4*s2 + 2*t2)))/(-1 + s1 + t1 - t2) - 
            32*(-2*Power(s1,3) + Power(s,2)*(-1 + s2) - s2 + 2*t1 + 
               s2*t1 - 2*Power(t1,2) - t2 - s2*t2 + t1*t2 + s2*t1*t2 + 
               Power(t2,2) - s2*Power(t2,2) - 
               Power(s1,2)*(s2 + 3*t1 - 4*(1 + t2)) - 
               2*s1*(Power(t1,2) + t2*(3 - s2 + t2) - t1*(1 + 2*t2)) + 
               s*(-8 + 2*Power(s1,2) - 3*t1 + 13*Power(t1,2) + 
                  s1*(9 - 13*s2 + 15*t1 - 15*t2) + 2*t2 - 26*t1*t2 + 
                  13*Power(t2,2) + s2*(10 - 12*t1 + 13*t2))) - 
            16*(-7 + 18*Power(s1,3) + 12*s2 - Power(s2,2) - 15*t1 - 
               12*s2*t1 + Power(s2,2)*t1 + 22*Power(t1,2) - 
               s2*Power(t1,2) + 
               2*Power(s1,2)*(-12 + 5*s2 + 12*t1 - 17*t2) + 
               Power(s,2)*(6 + 2*s1 - 9*s2 + 4*t1 - 4*t2) + 20*t2 + 
               8*s2*t2 - 19*t1*t2 - 7*s2*t1*t2 - 5*Power(t2,2) + 
               8*s2*Power(t2,2) + 
               s1*(-7 + 13*Power(t1,2) + 3*s2*(-1 + t1 - 6*t2) + 
                  33*t2 + 16*Power(t2,2) - t1*(6 + 29*t2)) + 
               s*(26 - 20*Power(s1,2) + Power(s2,2) + 23*t1 - 
                  56*Power(t1,2) + s2*(-32 + 46*t1 - 55*t2) - 19*t2 + 
                  115*t1*t2 - 59*Power(t2,2) + 
                  s1*(-30 + 58*s2 - 79*t1 + 81*t2))) - 
            8*(28 + 2*Power(s,3) - 52*Power(s1,3) - 28*s2 + 
               2*Power(s2,2) + 27*t1 + 33*s2*t1 - 3*Power(s2,2)*t1 - 
               63*Power(t1,2) + 3*s2*Power(t1,2) - 59*t2 - 17*s2*t2 - 
               Power(s2,2)*t2 + 64*t1*t2 + 19*s2*t1*t2 + 
               7*Power(t2,2) - 22*s2*Power(t2,2) + 
               Power(s,2)*(-11 - 12*s1 + 18*s2 - 10*t1 + 12*t2) + 
               Power(s1,2)*(34 - 31*s2 - 63*t1 + 94*t2) - 
               s1*(-22 + Power(s2,2) + 31*Power(t1,2) + 
                  t1*(8 - 75*t2) + s2*(-6 + 10*t1 - 49*t2) + 46*t2 + 
                  44*Power(t2,2)) + 
               s*(-52 + 62*Power(s1,2) - 53*t1 + 104*Power(t1,2) + 
                  s1*(55 - 102*s2 + 166*t1 - 177*t2) + 50*t2 - 
                  219*t1*t2 + 115*Power(t2,2) + 
                  s2*(45 - 80*t1 + 102*t2))) + 
            4*(16 + 2*Power(s,3) - 40*Power(s1,3) - 34*s2 + 
               Power(s2,2) + 25*t1 + 30*s2*t1 - 2*Power(s2,2)*t1 - 
               57*Power(t1,2) + 2*s2*Power(t1,2) - 
               6*Power(s1,2)*(3*s2 + 9*t1 - 13*t2) - 49*t2 - 9*s2*t2 - 
               4*Power(s2,2)*t2 + 64*t1*t2 + 18*s2*t1*t2 + 
               Power(t2,2) - 20*s2*Power(t2,2) + 
               Power(s,2)*(-11 - 4*s1 + 12*s2 - 4*t1 + 6*t2) + 
               s*(-19 + 42*Power(s1,2) + 2*Power(s2,2) - 46*t1 + 
                  72*Power(t1,2) + 42*t2 - 154*t1*t2 + 
                  82*Power(t2,2) - 
                  2*s1*(-27 + 38*s2 - 60*t1 + 63*t2) + 
                  s2*(31 - 56*t1 + 76*t2)) + 
               s1*(2*Power(s2,2) - 26*Power(t1,2) + 
                  3*t1*(-7 + 22*t2) + s2*(-1 - 8*t1 + 34*t2) - 
                  2*(7 + 2*t2 + 20*Power(t2,2)))) - 
            (2*(1 + 3*s - 4*s1 - 2*s2 + t2)*
               (4*(-1 + s1 + t1 - t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 2*Power(s,3)*(-1 + s2 - t1 + t2)*
                  (2 + 2*Power(s1,2) + Power(t1,2) + 
                    2*s1*(-2 + t1 - t2) + 2*t2 + Power(t2,2) - 
                    2*t1*(1 + t2)) - 
                 4*Power(s,2)*
                  (-3 + 2*t1 + 4*Power(t1,2) - 3*Power(t1,3) + t2 - 
                    6*t1*t2 + 8*Power(t1,2)*t2 + 2*Power(t2,2) - 
                    7*t1*Power(t2,2) + 2*Power(t2,3) + 
                    Power(s2,2)*t2*(2 - t1 + t2) + 
                    2*Power(s1,3)*(-1 + s2 - t1 + t2) + 
                    Power(s1,2)*
                     (2 + Power(s2,2) - 2*Power(t1,2) + 
                       s2*(-3 + t1 - 2*t2) + t2 + 5*t1*t2 - 
                       3*Power(t2,2)) - 
                    s1*(-3 + Power(t1,3) + Power(t1,2)*(5 - 3*t2) + 
                       4*t2 + 3*t1*(-2 + t2)*t2 + Power(t2,2) - 
                       Power(t2,3) + Power(s2,2)*(1 + 2*t2) + 
                       s2*(2 - 4*t1 - Power(t1,2) - 2*t2 + 
                       Power(t2,2))) + 
                    s2*(3 + 2*Power(t2,2) + Power(t2,3) + 
                       Power(t1,2)*(2 + t2) - 
                       t1*(5 + 4*t2 + 2*Power(t2,2)))) + 
                 2*s*(-5 + t1 + 13*Power(t1,2) - 9*Power(t1,3) - t2 - 
                    18*t1*t2 + 19*Power(t1,2)*t2 + Power(t2,2) + 
                    Power(s2,3)*Power(t2,2) - 15*t1*Power(t2,2) + 
                    5*Power(t2,3) + 
                    2*Power(s1,4)*(-1 + s2 - t1 + t2) + 
                    Power(s2,2)*t2*
                     (4 - t2 + Power(t2,2) - t1*(4 + t2)) + 
                    2*Power(s1,3)*
                     (Power(s2,2) - Power(t1,2) - 2*(-1 + t2)*t2 - 
                       s2*(1 + t2) + t1*(-1 + 3*t2)) + 
                    s2*(5 + 6*t2 + 17*Power(t2,2) + 4*Power(t2,3) + 
                       Power(t1,2)*(5 + 8*t2) - 
                       2*t1*(5 + 7*t2 + 6*Power(t2,2))) - 
                    2*s1*(3*Power(t1,3) + Power(t1,2)*(4 - 9*t2) + 
                       t2 + Power(s2,3)*t2 - Power(t2,3) + 
                       t1*(-7 - 6*t2 + 7*Power(t2,2)) - 
                       Power(s2,2)*(-2 + t2 + t1*(2 + t2)) - 
                       s2*(2 + Power(t1,2) - 5*t2 - 4*Power(t2,2) + 
                        Power(t2,3) - t1*(3 + Power(t2,2)))) + 
                    Power(s1,2)*
                     (Power(s2,3) - Power(t1,3) + 
                       Power(t1,2)*(-11 + 3*t2) - 
                       Power(s2,2)*(1 + t1 + 3*t2) - 
                       2*t1*(4 - 7*t2 + 2*Power(t2,2)) + 
                       2*(3 - 2*t2 - 2*Power(t2,2) + Power(t2,3)) + 
                       s2*(Power(t1,2) + 2*t1*(5 + t2) - 
                        2*(3 - 3*t2 + Power(t2,2)))))))/
             (s*Power(-1 + s1 + t1 - t2,2)*(-1 + s2 - t1 + t2)) + 
            ((3 - 2*Power(s1,2) - 3*t1 - s1*(s2 + t1 - 2*t2) + 
                 s*(-2 + 2*s1 + t1 - t2) + t2 + s2*t2)*
               (12*(-1 + s1 + t1 - t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 4*Power(s,3)*(-1 + s2 - t1 + t2)*
                  (-3*(-1 + s1)*(-1 + s1 + t1 - t2) + 
                    Power(2 - 2*s1 - t1 + t2,2)) - 
                 2*Power(s,2)*
                  (-4*(-2 + 2*s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                     (3 - 2*Power(s1,2) - 3*t1 - 
                       s1*(s2 + t1 - 2*t2) + t2 + s2*t2) - 
                    6*(-1 + s1 + t1 - t2)*
                     (3 + s2*(-3 + t1) + 2*t1 - 3*t2 - t1*t2 + 
                       Power(t2,2) + 
                       2*Power(s1,2)*(-1 + s2 - t1 + t2) + 
                       s1*(s2 - t2)*(-1 + s2 - t1 + t2))) + 
                 s*(4*(-1 + s2 - t1 + t2)*
                     Power(3 - 2*Power(s1,2) - 3*t1 - 
                       s1*(s2 + t1 - 2*t2) + t2 + s2*t2,2) + 
                    12*(-1 + s1 + t1 - t2)*
                     (Power(s2,2)*t2 - 
                       Power(s1,3)*(-1 + s2 - t1 + t2) - 
                       Power(s1,2)*(1 + s2 - t2)*(-1 + s2 - t1 + t2) + 
                       s2*(2 + t1*(-2 + t2) - 5*t2 - Power(t2,2)) + 
                       2*(-Power(-1 + t2,2) + t1*(1 + t2)) + 
                       s1*(-2 + Power(s2,2)*(-1 + t2) + t2 + 
                        Power(t2,2) - t1*(4 + t2) + 
                        s2*(3 + t1 + t2 - t1*t2 + Power(t2,2)))))))/
             (s*Power(-1 + s1 + t1 - t2,3)*(-1 + s2 - t1 + t2))))/
        ((-1 + s2)*(-1 + t1)) + 
       ((8*(-15 + 2*Power(s,2) + 7*Power(s1,2) + 5*s2 + Power(s2,2) - 
               t1 + s1*(-1 + 5*s2 + t1 - 5*t2) - 
               s*(4 + 9*s1 + 4*t1 - 5*t2) + 5*t2 - 2*s2*t2)*
             (3 - 2*Power(s1,2) - 3*t1 - s1*(s2 + t1 - 2*t2) + 
               s*(-2 + 2*s1 + t1 - t2) + t2 + s2*t2))/
           (-1 + s1 + t1 - t2) + 
          16*(-22 + 4*Power(s1,3) - 5*s2 + Power(s2,2) + 2*s2*t1 - 
             2*Power(t1,2) + Power(s1,2)*(-10 + 5*s2 + t1 - 4*t2) + 
             Power(s,2)*(-1 + 4*s1 - 2*s2 + 4*t1 - 3*t2) + 5*t2 + 
             2*s2*t2 - Power(s2,2)*t2 + 3*t1*t2 - 3*Power(t2,2) + 
             s1*(-17 + 2*Power(s2,2) - 2*t1 + 13*t2 - 4*s2*t2) + 
             s*(10 - 8*Power(s1,2) + Power(s2,2) - t1 + Power(t1,2) - 
                5*t2 + s2*(-1 - 3*t1 + 3*t2) + 
                s1*(8 - 3*s2 - 5*t1 + 7*t2))) - 
          (8*(-1 + 2*s - 3*s1 - s2 + t2)*
             (4*(-1 + s1 + t1 - t2)*
                Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
               2*Power(s,3)*(-1 + s2 - t1 + t2)*
                (2 + 2*Power(s1,2) + Power(t1,2) + 
                  2*s1*(-2 + t1 - t2) + 2*t2 + Power(t2,2) - 
                  2*t1*(1 + t2)) - 
               4*Power(s,2)*(-3 + 2*t1 + 4*Power(t1,2) - 
                  3*Power(t1,3) + t2 - 6*t1*t2 + 8*Power(t1,2)*t2 + 
                  2*Power(t2,2) - 7*t1*Power(t2,2) + 2*Power(t2,3) + 
                  Power(s2,2)*t2*(2 - t1 + t2) + 
                  2*Power(s1,3)*(-1 + s2 - t1 + t2) + 
                  Power(s1,2)*
                   (2 + Power(s2,2) - 2*Power(t1,2) + 
                     s2*(-3 + t1 - 2*t2) + t2 + 5*t1*t2 - 
                     3*Power(t2,2)) - 
                  s1*(-3 + Power(t1,3) + Power(t1,2)*(5 - 3*t2) + 
                     4*t2 + 3*t1*(-2 + t2)*t2 + Power(t2,2) - 
                     Power(t2,3) + Power(s2,2)*(1 + 2*t2) + 
                     s2*(2 - 4*t1 - Power(t1,2) - 2*t2 + Power(t2,2))\
) + s2*(3 + 2*Power(t2,2) + Power(t2,3) + Power(t1,2)*(2 + t2) - 
                     t1*(5 + 4*t2 + 2*Power(t2,2)))) + 
               2*s*(-5 + t1 + 13*Power(t1,2) - 9*Power(t1,3) - t2 - 
                  18*t1*t2 + 19*Power(t1,2)*t2 + Power(t2,2) + 
                  Power(s2,3)*Power(t2,2) - 15*t1*Power(t2,2) + 
                  5*Power(t2,3) + 2*Power(s1,4)*(-1 + s2 - t1 + t2) + 
                  Power(s2,2)*t2*
                   (4 - t2 + Power(t2,2) - t1*(4 + t2)) + 
                  2*Power(s1,3)*
                   (Power(s2,2) - Power(t1,2) - 2*(-1 + t2)*t2 - 
                     s2*(1 + t2) + t1*(-1 + 3*t2)) + 
                  s2*(5 + 6*t2 + 17*Power(t2,2) + 4*Power(t2,3) + 
                     Power(t1,2)*(5 + 8*t2) - 
                     2*t1*(5 + 7*t2 + 6*Power(t2,2))) - 
                  2*s1*(3*Power(t1,3) + Power(t1,2)*(4 - 9*t2) + t2 + 
                     Power(s2,3)*t2 - Power(t2,3) + 
                     t1*(-7 - 6*t2 + 7*Power(t2,2)) - 
                     Power(s2,2)*(-2 + t2 + t1*(2 + t2)) - 
                     s2*(2 + Power(t1,2) - 5*t2 - 4*Power(t2,2) + 
                        Power(t2,3) - t1*(3 + Power(t2,2)))) + 
                  Power(s1,2)*
                   (Power(s2,3) - Power(t1,3) + 
                     Power(t1,2)*(-11 + 3*t2) - 
                     Power(s2,2)*(1 + t1 + 3*t2) - 
                     2*t1*(4 - 7*t2 + 2*Power(t2,2)) + 
                     2*(3 - 2*t2 - 2*Power(t2,2) + Power(t2,3)) + 
                     s2*(Power(t1,2) + 2*t1*(5 + t2) - 
                        2*(3 - 3*t2 + Power(t2,2)))))))/
           (s*Power(-1 + s1 + t1 - t2,2)*(-1 + s2 - t1 + t2)) + 
          (4*(3 - 2*Power(s1,2) - 3*t1 - s1*(s2 + t1 - 2*t2) + 
               s*(-2 + 2*s1 + t1 - t2) + t2 + s2*t2)*
             (12*(-1 + s1 + t1 - t2)*
                Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
               4*Power(s,3)*(-1 + s2 - t1 + t2)*
                (-3*(-1 + s1)*(-1 + s1 + t1 - t2) + 
                  Power(2 - 2*s1 - t1 + t2,2)) - 
               2*Power(s,2)*(-4*(-2 + 2*s1 + t1 - t2)*
                   (-1 + s2 - t1 + t2)*
                   (3 - 2*Power(s1,2) - 3*t1 - s1*(s2 + t1 - 2*t2) + 
                     t2 + s2*t2) - 
                  6*(-1 + s1 + t1 - t2)*
                   (3 + s2*(-3 + t1) + 2*t1 - 3*t2 - t1*t2 + 
                     Power(t2,2) + 2*Power(s1,2)*(-1 + s2 - t1 + t2) + 
                     s1*(s2 - t2)*(-1 + s2 - t1 + t2))) + 
               s*(4*(-1 + s2 - t1 + t2)*
                   Power(3 - 2*Power(s1,2) - 3*t1 - 
                     s1*(s2 + t1 - 2*t2) + t2 + s2*t2,2) + 
                  12*(-1 + s1 + t1 - t2)*
                   (Power(s2,2)*t2 - Power(s1,3)*(-1 + s2 - t1 + t2) - 
                     Power(s1,2)*(1 + s2 - t2)*(-1 + s2 - t1 + t2) + 
                     s2*(2 + t1*(-2 + t2) - 5*t2 - Power(t2,2)) + 
                     2*(-Power(-1 + t2,2) + t1*(1 + t2)) + 
                     s1*(-2 + Power(s2,2)*(-1 + t2) + t2 + 
                        Power(t2,2) - t1*(4 + t2) + 
                        s2*(3 + t1 + t2 - t1*t2 + Power(t2,2)))))))/
           (s*Power(-1 + s1 + t1 - t2,3)*(-1 + s2 - t1 + t2)))/
        ((-1 + s1)*(-s + s1 - t2)))*B3(s,1 - s2 + t1 - t2,2 + s - s1 - s2))/
   (16.*Power(Pi,2)) + (((-8*s*(-1 + s1 + t1 - t2)*
          (-1 - 3*t1 + s*t1 + t2 - s*t2 + s2*(2 + t2) + 
            s1*(-2 + s2 - t1 + 2*t2))*
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
       (4*(-64*s*(-1 + s1 + t1 - t2)*(-1 + s + t2) - 
            (2*(-1 - 3*t1 + s*t1 + t2 - s*t2 + s2*(2 + t2) + 
                 s1*(-2 + s2 - t1 + 2*t2))*
               (-16 + 3*Power(s,2) + 4*Power(s1,2) + 3*s2 + 
                 Power(s2,2) - 2*t1 - t2 + 5*s2*t2 - t1*t2 + 
                 Power(t2,2) - 2*s*(-3 + 4*s1 + 2*s2 + 3*t2) + 
                 s1*(-6 + 6*s2 - t1 + 7*t2)))/(-1 + s2 - t1 + t2) + 
            32*(Power(s,2)*(-11 + 12*s1 + 13*t1 - 12*t2) + 
               (-1 + t2)*(1 - t1 - t2 + s2*t2 + 
                  s1*(-2 + s2 - t1 + 2*t2)) + 
               s*(13 + s2*(-1 + t1) - 13*t1 - Power(t1,2) - t2 + 
                  16*t1*t2 - 14*Power(t2,2) + 
                  s1*(-15 + 2*s2 - 2*t1 + 15*t2))) - 
            16*(-2 + Power(s,3) + s2 + 5*t1 - s2*t1 + Power(t1,2) + 
               Power(s,2)*(-36 + 47*s1 - s2 + 56*t1 - 48*t2) + 7*t2 - 
               3*s2*t2 + Power(s2,2)*t2 - 12*t1*t2 - s2*t1*t2 - 
               5*Power(t2,2) + 8*s2*Power(t2,2) + 
               Power(s1,2)*(-3 + 2*t2) + 
               s1*(11 + Power(s2,2) + Power(t1,2) + t1*(9 - 10*t2) - 
                  2*s2*(5 + t1 - 6*t2) - 26*t2 + 16*Power(t2,2)) + 
               s*(48 + 4*Power(s1,2) - 57*t1 - 8*Power(t1,2) + 
                  s2*(-4 + 9*t1 - 3*t2) - 3*t2 + 83*t1*t2 - 
                  68*Power(t2,2) + s1*(-74 + 16*s2 - 11*t1 + 68*t2))) - 
            4*(8 + 4*Power(s,3) + 3*Power(s1,3) + 18*s2 + 
               4*Power(s2,2) - 12*t1 - 12*s2*t1 + 12*Power(t1,2) + 
               Power(s,2)*(-22 + 48*s1 - 4*s2 + 73*t1 - 59*t2) + 
               Power(s1,2)*(-10 + 6*t1 - t2) + 36*t2 + 17*s2*t2 + 
               3*Power(s2,2)*t2 - 48*t1*t2 - 6*s2*t1*t2 + 
               2*Power(t2,2) + 19*s2*Power(t2,2) + 
               s1*(3 + 3*Power(s2,2) + 6*Power(t1,2) + 
                  t1*(17 - 31*t2) - 19*t2 + 38*Power(t2,2) + 
                  s2*(-13 - 6*t1 + 27*t2)) + 
               s*(31 + 11*Power(s1,2) - 65*t1 - 18*Power(t1,2) + 
                  s2*(2 + 15*t1 - 2*t2) - 15*t2 + 126*t1*t2 - 
                  99*Power(t2,2) + s1*(-107 + 33*s2 - 25*t1 + 97*t2))) \
+ 8*(4 + 4*Power(s,3) + 2*Power(s1,3) + 8*s2 + 4*Power(s2,2) - 4*t1 - 
               10*s2*t1 + 8*Power(t1,2) + 
               Power(s,2)*(-49 + 79*s1 - 4*s2 + 105*t1 - 84*t2) + 
               14*t2 + 11*s2*t2 + 4*Power(s2,2)*t2 - 44*t1*t2 - 
               5*s2*t1*t2 - 6*Power(t2,2) + 22*s2*Power(t2,2) + 
               Power(s1,2)*(-16 + s2 + 3*t1 + 6*t2) + 
               s1*(10 + 4*Power(s2,2) + 5*Power(t1,2) + 
                  t1*(21 - 32*t2) - 43*t2 + 44*Power(t2,2) + 
                  s2*(-25 - 7*t1 + 37*t2)) + 
               s*(71 + 12*Power(s1,2) - 98*t1 - 21*Power(t1,2) + 
                  3*s2*(-1 + 7*t1 - 3*t2) - 12*t2 + 173*t1*t2 - 
                  138*Power(t2,2) + s1*(-140 + 39*s2 - 27*t1 + 128*t2))\
) + ((-1 - 3*t1 + s*t1 + t2 - s*t2 + s2*(2 + t2) + 
                 s1*(-2 + s2 - t1 + 2*t2))*
               (4*s*(-1 + s1 + t1 - t2)*
                  Power(-1 - 3*t1 + s*t1 + t2 - s*t2 + 
                    s2*(2 + t2) + s1*(-2 + s2 - t1 + 2*t2),2) + 
                 12*(-1 + s2 - t1 + t2)*
                  (Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                    Power(s,2)*(-1 + s1)*
                     (t1*(-1 + t2) + (s1 - t2)*t2) - 
                    s*(-1 + s1)*
                     (Power(s1,2)*(-1 + t2) + 2*(1 + t1 - t2)*t2 + 
                       s2*(1 + t1*(-1 + t2) - 2*t2 - Power(t2,2)) + 
                       s1*(t1*(-3 + t2) - Power(-1 + t2,2) + 
                        s2*(1 + t2))))))/
             (s*(-1 + s1 + t1 - t2)*Power(1 - s2 + t1 - t2,3)) - 
            (2*(3 + 5*s - 5*s1 - 4*s2 + t1 - 3*t2)*
               (Power(s,3)*Power(t1 - t2,2)*(-1 + s1 + t1 - t2) + 
                 2*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                 2*Power(s,2)*
                  (-Power(t1,2) + 3*Power(t1,3) + t2 - 
                    8*Power(t1,2)*t2 + Power(t2,2) + 
                    7*t1*Power(t2,2) - 2*Power(t2,3) + 
                    Power(s1,2)*
                     (Power(t1,2) + (-1 + t2)*t2 - 
                       t1*(-2 + s2 + 2*t2)) - 
                    s2*(Power(t1,2)*(2 + t2) + 
                       t2*(2 + 4*t2 + Power(t2,2)) - 
                       t1*(1 + 6*t2 + 2*Power(t2,2))) + 
                    s1*(Power(t1,3) - Power(t1,2)*(-3 + s2 + 3*t2) + 
                       t2*(2*s2 + t2 + s2*t2 - Power(t2,2)) + 
                       t1*(-2 - 4*t2 + 3*Power(t2,2)))) + 
                 s*(-1 - 5*t1 - 3*Power(t1,2) + 9*Power(t1,3) + 
                    Power(s1,3)*
                     (Power(s2,2) + Power(t1,2) - 
                       2*s2*(1 + t1 - t2) - 2*t1*(-1 + t2) + 
                       2*Power(-1 + t2,2)) - 3*t2 - 10*t1*t2 - 
                    19*Power(t1,2)*t2 + 9*Power(t2,2) + 
                    15*t1*Power(t2,2) - 5*Power(t2,3) + 
                    Power(s1,2)*
                     (Power(t1,3) - 
                       2*s2*
                       (2 + Power(t1,2) - 2*t1*(-1 + t2) - 2*t2) + 
                       Power(s2,2)*(1 + t1 - t2) - 
                       3*Power(t1,2)*(-1 + t2) - 
                       2*Power(-1 + t2,2)*t2 + 
                       4*t1*(2 - 3*t2 + Power(t2,2))) + 
                    Power(s2,2)*
                     (-2 - 12*t2 - 7*Power(t2,2) - Power(t2,3) + 
                       t1*(2 + 6*t2 + Power(t2,2))) - 
                    2*s2*(-1 - 6*t2 + 5*Power(t2,2) + 2*Power(t2,3) + 
                       Power(t1,2)*(5 + 4*t2) - 
                       2*t1*(2 + 6*t2 + 3*Power(t2,2))) + 
                    s1*(6*Power(t1,3) + Power(t1,2)*(23 - 18*t2) - 
                       Power(-1 + t2,2)*(1 + 2*t2) + 
                       2*t1*(2 - 9*t2 + 7*Power(t2,2)) + 
                       Power(s2,2)*(6*t1 + t2*(4 + t2)) - 
                       2*s2*(-2 + 6*Power(t1,2) + Power(t2,2) + 
                        Power(t2,3) - t1*(-10 + 4*t2 + Power(t2,2))))))\
)/(s*(-1 + s1 + t1 - t2)*Power(-1 + s2 - t1 + t2,2))))/
        ((s - s2 + t1)*(s - s1 + t2)) + 
       (4*(64*s*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) - 
            (2*(-1 - 3*t1 + s*t1 + t2 - s*t2 + s2*(2 + t2) + 
                 s1*(-2 + s2 - t1 + 2*t2))*
               (-9 + 3*Power(s,2) + 9*Power(s1,2) - 5*s2 + 
                 Power(s2,2) - 2*t1 + 4*s2*t1 - Power(t1,2) - 3*t2 + 
                 s2*t2 + t1*t2 - s*(-7 + 13*s1 + 4*s2 + 5*t1 + t2) + 
                 s1*(10*s2 + 3*(-4 + t1 + t2))))/(-1 + s2 - t1 + t2) - 
            32*(1 + s2 - 3*t1 - s2*t1 + 2*Power(t1,2) + 
               Power(s1,2)*(2 - s2 + t1 - 2*t2) + Power(s2,2)*t2 + 
               t1*t2 - 2*s2*t1*t2 - Power(t2,2) + s2*Power(t2,2) + 
               s1*(-3 + Power(s2,2) + 2*Power(t1,2) + t1*(4 - 5*t2) + 
                  t2 + 2*Power(t2,2) + s2*(-1 - 3*t1 + 2*t2)) + 
               s*(12 + t1 - 15*Power(t1,2) + 
                  s2*(-12 + 14*t1 - 13*t2) - 2*t2 + 29*t1*t2 - 
                  14*Power(t2,2) + s1*(-13 + 13*s2 - 14*t1 + 14*t2))) - 
            16*(-6 + Power(s1,3) + Power(s,2)*(s1 - s2) - 4*s2 - 
               4*Power(s2,2) + 17*t1 + 18*s2*t1 - 21*Power(t1,2) + 
               Power(s1,2)*(9*s2 - 7*t1 + 16*(-1 + t2)) + t2 - 
               10*s2*t2 - 7*Power(s2,2)*t2 + 4*t1*t2 + 15*s2*t1*t2 + 
               5*Power(t2,2) - 8*s2*Power(t2,2) - 
               s1*(-18 + 7*Power(s2,2) + 15*Power(t1,2) + 
                  t1*(36 - 38*t2) + 5*t2 + 16*Power(t2,2) + 
                  s2*(-11 - 23*t1 + 14*t2)) + 
               s*(-46 - 3*Power(s1,2) + Power(s2,2) - 5*t1 + 
                  73*Power(t1,2) + s1*(60 - 58*s2 + 64*t1 - 68*t2) + 
                  16*t2 - 141*t1*t2 + 67*Power(t2,2) + 
                  s2*(44 - 65*t1 + 57*t2))) - 
            4*(-21 + Power(s1,3) + s2 - 16*Power(s2,2) + 38*t1 + 
               62*s2*t1 + Power(s2,2)*t1 - 63*Power(t1,2) + 
               Power(s,2)*(s1 - 4*(-1 + s2 + t1)) + 2*t2 - 39*s2*t2 - 
               16*Power(s2,2)*t2 + 46*t1*t2 + 33*s2*t1*t2 + 
               Power(t2,2) - 20*s2*Power(t2,2) + 
               Power(s1,2)*(-17 + 20*s2 - 17*t1 + 39*t2) + 
               s*(-49 - 11*Power(s1,2) + 4*Power(s2,2) - 3*t1 + 
                  102*Power(t1,2) + 
                  s1*(65 - 71*s2 + 81*t1 - 100*t2) + 22*t2 - 
                  205*t1*t2 + 98*Power(t2,2) + s2*(28 - 85*t1 + 74*t2)\
) + s1*(45 - 14*Power(s2,2) - 33*Power(t1,2) + 
                  s2*(43 + 50*t1 - 30*t2) + 9*t2 - 40*Power(t2,2) + 
                  t1*(-87 + 86*t2))) + 
            8*(-13 + 5*Power(s1,3) + Power(s,2)*(5*s1 - 4*s2) - 7*s2 - 
               16*Power(s2,2) + 35*t1 + 62*s2*t1 + Power(s2,2)*t1 - 
               64*Power(t1,2) - 39*s2*t2 - 18*Power(s2,2)*t2 + 
               33*t1*t2 + 39*s2*t1*t2 + 7*Power(t2,2) - 
               22*s2*Power(t2,2) + 
               Power(s1,2)*(-38 + 28*s2 - 17*t1 + 45*t2) + 
               s*(-74 - 17*Power(s1,2) + 4*Power(s2,2) - 4*t1 + 
                  145*Power(t1,2) + 
                  s1*(117 - 111*s2 + 116*t1 - 139*t2) + 35*t2 - 
                  285*t1*t2 + 135*Power(t2,2) + 
                  s2*(64 - 127*t1 + 107*t2)) - 
               s1*(16*Power(s2,2) + 39*Power(t1,2) + 
                  t1*(101 - 100*t2) + s2*(-35 - 62*t1 + 34*t2) + 
                  4*(-9 + t2 + 11*Power(t2,2)))) + 
            ((-1 - 3*t1 + s*t1 + t2 - s*t2 + s2*(2 + t2) + 
                 s1*(-2 + s2 - t1 + 2*t2))*
               (4*s*(-1 + s1 + t1 - t2)*
                  Power(-1 - 3*t1 + s*t1 + t2 - s*t2 + 
                    s2*(2 + t2) + s1*(-2 + s2 - t1 + 2*t2),2) + 
                 12*(-1 + s2 - t1 + t2)*
                  (Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                    Power(s,2)*(-1 + s1)*
                     (t1*(-1 + t2) + (s1 - t2)*t2) - 
                    s*(-1 + s1)*
                     (Power(s1,2)*(-1 + t2) + 2*(1 + t1 - t2)*t2 + 
                       s2*(1 + t1*(-1 + t2) - 2*t2 - Power(t2,2)) + 
                       s1*(t1*(-3 + t2) - Power(-1 + t2,2) + 
                        s2*(1 + t2))))))/
             (s*(-1 + s1 + t1 - t2)*Power(1 - s2 + t1 - t2,3)) - 
            (2*(5 + 5*s - 7*s1 - 4*s2 - t1 - t2)*
               (Power(s,3)*Power(t1 - t2,2)*(-1 + s1 + t1 - t2) + 
                 2*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                 2*Power(s,2)*
                  (-Power(t1,2) + 3*Power(t1,3) + t2 - 
                    8*Power(t1,2)*t2 + Power(t2,2) + 
                    7*t1*Power(t2,2) - 2*Power(t2,3) + 
                    Power(s1,2)*
                     (Power(t1,2) + (-1 + t2)*t2 - 
                       t1*(-2 + s2 + 2*t2)) - 
                    s2*(Power(t1,2)*(2 + t2) + 
                       t2*(2 + 4*t2 + Power(t2,2)) - 
                       t1*(1 + 6*t2 + 2*Power(t2,2))) + 
                    s1*(Power(t1,3) - Power(t1,2)*(-3 + s2 + 3*t2) + 
                       t2*(2*s2 + t2 + s2*t2 - Power(t2,2)) + 
                       t1*(-2 - 4*t2 + 3*Power(t2,2)))) + 
                 s*(-1 - 5*t1 - 3*Power(t1,2) + 9*Power(t1,3) + 
                    Power(s1,3)*
                     (Power(s2,2) + Power(t1,2) - 
                       2*s2*(1 + t1 - t2) - 2*t1*(-1 + t2) + 
                       2*Power(-1 + t2,2)) - 3*t2 - 10*t1*t2 - 
                    19*Power(t1,2)*t2 + 9*Power(t2,2) + 
                    15*t1*Power(t2,2) - 5*Power(t2,3) + 
                    Power(s1,2)*
                     (Power(t1,3) - 
                       2*s2*
                       (2 + Power(t1,2) - 2*t1*(-1 + t2) - 2*t2) + 
                       Power(s2,2)*(1 + t1 - t2) - 
                       3*Power(t1,2)*(-1 + t2) - 
                       2*Power(-1 + t2,2)*t2 + 
                       4*t1*(2 - 3*t2 + Power(t2,2))) + 
                    Power(s2,2)*
                     (-2 - 12*t2 - 7*Power(t2,2) - Power(t2,3) + 
                       t1*(2 + 6*t2 + Power(t2,2))) - 
                    2*s2*(-1 - 6*t2 + 5*Power(t2,2) + 2*Power(t2,3) + 
                       Power(t1,2)*(5 + 4*t2) - 
                       2*t1*(2 + 6*t2 + 3*Power(t2,2))) + 
                    s1*(6*Power(t1,3) + Power(t1,2)*(23 - 18*t2) - 
                       Power(-1 + t2,2)*(1 + 2*t2) + 
                       2*t1*(2 - 9*t2 + 7*Power(t2,2)) + 
                       Power(s2,2)*(6*t1 + t2*(4 + t2)) - 
                       2*s2*(-2 + 6*Power(t1,2) + Power(t2,2) + 
                        Power(t2,3) - t1*(-10 + 4*t2 + Power(t2,2))))))\
)/(s*(-1 + s1 + t1 - t2)*Power(-1 + s2 - t1 + t2,2))))/
        ((-1 + s2)*(-s + s2 - t1)) + 
       (4*(64*s*(-1 + s1 + t1 - t2)*(-1 + t2) - 
            (2*(-1 - 3*t1 + s*t1 + t2 - s*t2 + s2*(2 + t2) + 
                 s1*(-2 + s2 - t1 + 2*t2))*
               (-14 + Power(s,2) + 6*Power(s1,2) + 4*s2 + 
                 Power(s2,2) - 2*t1 - 5*t2 + 2*s2*t2 + t1*t2 - 
                 Power(t2,2) - s*(3 + 6*s1 + 3*s2 + 2*t2) + 
                 s1*(-3 + 6*s2 + t1 + 4*t2)))/(-1 + s2 - t1 + t2) - 
            32*(-((s1 - t2)*(1 - t1 - t2 + s2*t2 + 
                    s1*(-2 + s2 - t1 + 2*t2))) + 
               s*(12 + Power(s1,2) - 12*t1 - Power(t1,2) + 
                  s2*(-1 + t1 - t2) - t2 + 16*t1*t2 - 14*Power(t2,2) + 
                  s1*(-14 + s2 - t1 + 14*t2))) + 
            16*(1 - Power(s1,3) - 2*s2 + Power(s,2)*(-1 + t1) + 
               Power(t1,2) + Power(s1,2)*(15 - 8*s2 + 7*t1 - 16*t2) + 
               6*t2 + 4*s2*t2 - 12*t1*t2 - s2*t1*t2 - 5*Power(t2,2) + 
               8*s2*Power(t2,2) - 
               s1*(4 - Power(t1,2) + 8*t2 - 16*Power(t2,2) + 
                  s2*(3 + t1 + t2) + t1*(-13 + 10*t2)) + 
               s*(45 + 8*Power(s1,2) - 48*t1 - 8*Power(t1,2) + 
                  s2*(-6 + 7*t1 - 7*t2) - 7*t2 + 82*t1*t2 - 
                  67*Power(t2,2) + s1*(-63 + 8*s2 - 8*t1 + 68*t2))) + 
            4*(-7 - 4*Power(s1,3) - 16*s2 - t1 - 6*s2*t1 + 
               12*Power(t1,2) + 
               Power(s1,2)*(20 - 17*s2 + 14*t1 - 35*t2) + 
               Power(s,2)*(-5 + s1 + 7*t1 - t2) + 9*t2 + 21*s2*t2 + 
               Power(s2,2)*t2 - 48*t1*t2 - 6*s2*t1*t2 + 
               2*Power(t2,2) + 19*s2*Power(t2,2) + 
               s1*(-15 + Power(s2,2) + 44*t1 + 6*Power(t1,2) - 
                  32*t2 - 31*t1*t2 + 38*Power(t2,2) - 
                  s2*(26 + 6*t1 + t2)) + 
               s*(34 + 17*Power(s1,2) + s2 - 42*t1 + 11*s2*t1 - 
                  18*Power(t1,2) - 16*t2 - 11*s2*t2 + 120*t1*t2 - 
                  93*Power(t2,2) + s1*(-59 + 16*s2 - 19*t1 + 93*t2))) - 
            8*(4 - 5*Power(s1,3) - 9*s2 - Power(s2,2) - 3*s2*t1 + 
               8*Power(t1,2) + 
               Power(s1,2)*(34 - 22*s2 + 17*t1 - 44*t2) + 
               Power(s,2)*(-6 + 6*t1 - t2) + 16*t2 + 19*s2*t2 - 
               44*t1*t2 - 5*s2*t1*t2 - 6*Power(t2,2) + 
               22*s2*Power(t2,2) - 
               s1*(5 - 5*Power(t1,2) + 25*t2 - 44*Power(t2,2) + 
                  5*s2*(4 + t1 + t2) + t1*(-46 + 32*t2)) + 
               s*(63 + 22*Power(s1,2) - 77*t1 - 21*Power(t1,2) + 
                  s2*(-5 + 15*t1 - 14*t2) - 18*t2 + 168*t1*t2 - 
                  133*Power(t2,2) + s1*(-106 + 21*s2 - 22*t1 + 137*t2))\
) + ((-1 - 3*t1 + s*t1 + t2 - s*t2 + s2*(2 + t2) + 
                 s1*(-2 + s2 - t1 + 2*t2))*
               (4*s*(-1 + s1 + t1 - t2)*
                  Power(-1 - 3*t1 + s*t1 + t2 - s*t2 + 
                    s2*(2 + t2) + s1*(-2 + s2 - t1 + 2*t2),2) + 
                 12*(-1 + s2 - t1 + t2)*
                  (Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                    Power(s,2)*(-1 + s1)*
                     (t1*(-1 + t2) + (s1 - t2)*t2) - 
                    s*(-1 + s1)*
                     (Power(s1,2)*(-1 + t2) + 2*(1 + t1 - t2)*t2 + 
                       s2*(1 + t1*(-1 + t2) - 2*t2 - Power(t2,2)) + 
                       s1*(t1*(-3 + t2) - Power(-1 + t2,2) + 
                        s2*(1 + t2))))))/
             (s*(-1 + s1 + t1 - t2)*Power(1 - s2 + t1 - t2,3)) - 
            (2*(3 + 3*s - 6*s1 - 3*s2 - t1 - t2)*
               (Power(s,3)*Power(t1 - t2,2)*(-1 + s1 + t1 - t2) + 
                 2*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                 2*Power(s,2)*
                  (-Power(t1,2) + 3*Power(t1,3) + t2 - 
                    8*Power(t1,2)*t2 + Power(t2,2) + 
                    7*t1*Power(t2,2) - 2*Power(t2,3) + 
                    Power(s1,2)*
                     (Power(t1,2) + (-1 + t2)*t2 - 
                       t1*(-2 + s2 + 2*t2)) - 
                    s2*(Power(t1,2)*(2 + t2) + 
                       t2*(2 + 4*t2 + Power(t2,2)) - 
                       t1*(1 + 6*t2 + 2*Power(t2,2))) + 
                    s1*(Power(t1,3) - Power(t1,2)*(-3 + s2 + 3*t2) + 
                       t2*(2*s2 + t2 + s2*t2 - Power(t2,2)) + 
                       t1*(-2 - 4*t2 + 3*Power(t2,2)))) + 
                 s*(-1 - 5*t1 - 3*Power(t1,2) + 9*Power(t1,3) + 
                    Power(s1,3)*
                     (Power(s2,2) + Power(t1,2) - 
                       2*s2*(1 + t1 - t2) - 2*t1*(-1 + t2) + 
                       2*Power(-1 + t2,2)) - 3*t2 - 10*t1*t2 - 
                    19*Power(t1,2)*t2 + 9*Power(t2,2) + 
                    15*t1*Power(t2,2) - 5*Power(t2,3) + 
                    Power(s1,2)*
                     (Power(t1,3) - 
                       2*s2*
                       (2 + Power(t1,2) - 2*t1*(-1 + t2) - 2*t2) + 
                       Power(s2,2)*(1 + t1 - t2) - 
                       3*Power(t1,2)*(-1 + t2) - 
                       2*Power(-1 + t2,2)*t2 + 
                       4*t1*(2 - 3*t2 + Power(t2,2))) + 
                    Power(s2,2)*
                     (-2 - 12*t2 - 7*Power(t2,2) - Power(t2,3) + 
                       t1*(2 + 6*t2 + Power(t2,2))) - 
                    2*s2*(-1 - 6*t2 + 5*Power(t2,2) + 2*Power(t2,3) + 
                       Power(t1,2)*(5 + 4*t2) - 
                       2*t1*(2 + 6*t2 + 3*Power(t2,2))) + 
                    s1*(6*Power(t1,3) + Power(t1,2)*(23 - 18*t2) - 
                       Power(-1 + t2,2)*(1 + 2*t2) + 
                       2*t1*(2 - 9*t2 + 7*Power(t2,2)) + 
                       Power(s2,2)*(6*t1 + t2*(4 + t2)) - 
                       2*s2*(-2 + 6*Power(t1,2) + Power(t2,2) + 
                        Power(t2,3) - t1*(-10 + 4*t2 + Power(t2,2))))))\
)/(s*(-1 + s1 + t1 - t2)*Power(-1 + s2 - t1 + t2,2))))/
        ((-1 + t1)*(-1 + t2)) + 
       (4*(64*s*(t1 - t2)*(-1 + s1 + t1 - t2) - 
            4*(5 + 40*Power(s1,3) + 8*s2 + Power(s2,2) - 48*t1 - 
               18*s2*t1 + Power(s2,2)*t1 + 63*Power(t1,2) + 
               2*Power(s,2)*(1 + t1) + 
               Power(s1,2)*(-3 + 17*s2 + 59*t1 - 79*t2) + 20*t2 + 
               46*s2*t2 - 46*t1*t2 - 33*s2*t1*t2 - Power(t2,2) + 
               20*s2*Power(t2,2) + 
               s1*(-16 + 2*Power(s2,2) + 33*Power(t1,2) + 
                  t1*(57 - 86*t2) + 2*s2*(-8 + t1 - 17*t2) + 30*t2 + 
                  40*Power(t2,2)) - 
               s*(-23 + 21*Power(s1,2) - 34*t1 + 102*Power(t1,2) + 
                  s1*(21 - 19*s2 + 111*t1 - 140*t2) - 
                  2*s2*(-7 + 9*t1 - 7*t2) + 66*t2 - 205*t1*t2 + 
                  98*Power(t2,2))) - 
            (2*(-1 - 3*t1 + s*t1 + t2 - s*t2 + s2*(2 + t2) + 
                 s1*(-2 + s2 - t1 + 2*t2))*
               (-9 + Power(s,2) + 12*Power(s1,2) - 3*s2 + 
                 Power(s2,2) - 4*t1 + 3*s2*t1 + Power(t1,2) - t2 - 
                 s2*t2 - t1*t2 + s*(2 - 9*s1 - 3*s2 - 3*t1 + t2) + 
                 s1*(9*s2 + 8*t1 - 3*(4 + t2))))/(-1 + s2 - t1 + t2) - 
            16*(11 + 18*Power(s1,3) + 2*s2 - 32*t1 - 3*s2*t1 + 
               21*Power(t1,2) + 
               Power(s1,2)*(-25 + 9*s2 + 26*t1 - 34*t2) + 
               Power(s,2)*(t1 - t2) + 4*t2 + 17*s2*t2 - 4*t1*t2 - 
               15*s2*t1*t2 - 5*Power(t2,2) + 8*s2*Power(t2,2) + 
               s1*(-14 + 15*Power(t1,2) + s2*(-3 + 2*t1 - 17*t2) + 
                  41*t2 - 38*t1*t2 + 16*Power(t2,2)) - 
               s*(-12 + 10*Power(s1,2) - 48*t1 + 73*Power(t1,2) + 
                  s1*(4 - 9*s2 + 77*t1 - 88*t2) + 62*t2 - 141*t1*t2 + 
                  67*Power(t2,2) + s2*(7 - 8*t1 + 6*t2))) + 
            8*(20 + 52*Power(s1,3) + s2 - 72*t1 - 13*s2*t1 + 
               Power(s2,2)*t1 + 64*Power(t1,2) + 
               Power(s1,2)*(-41 + 27*s2 + 73*t1 - 93*t2) + 
               Power(s,2)*(1 + s1 + 3*t1 - 2*t2) + 17*t2 + 49*s2*t2 - 
               33*t1*t2 - 39*s2*t1*t2 - 7*Power(t2,2) + 
               22*s2*Power(t2,2) + 
               s1*(-42 + 2*Power(s2,2) + 39*Power(t1,2) + 
                  t1*(29 - 100*t2) + s2*(-17 + 8*t1 - 44*t2) + 
                  69*t2 + 44*Power(t2,2)) - 
               s*(-39 + 31*Power(s1,2) - 66*t1 + 145*Power(t1,2) + 
                  s1*(17 - 21*s2 + 160*t1 - 191*t2) + 107*t2 - 
                  285*t1*t2 + 135*Power(t2,2) + s2*(16 - 20*t1 + 15*t2)\
)) - 32*(-2*Power(s1,3) - Power(s1,2)*(-4 + s2 + 3*t1 - 4*t2) - 
               (-2 + 2*t1 - t2)*(-1 + t1 + t2 - s2*t2) + 
               s1*(1 + t1 - 2*Power(t1,2) - 7*t2 + 2*s2*t2 + 
                  5*t1*t2 - 2*Power(t2,2)) + 
               s*(-1 + Power(s1,2) + s2 - 13*t1 - s2*t1 + 
                  15*Power(t1,2) + 14*t2 + s2*t2 - 29*t1*t2 + 
                  14*Power(t2,2) - s1*(s2 - 15*t1 + 16*t2))) + 
            ((-1 - 3*t1 + s*t1 + t2 - s*t2 + s2*(2 + t2) + 
                 s1*(-2 + s2 - t1 + 2*t2))*
               (4*s*(-1 + s1 + t1 - t2)*
                  Power(-1 - 3*t1 + s*t1 + t2 - s*t2 + 
                    s2*(2 + t2) + s1*(-2 + s2 - t1 + 2*t2),2) + 
                 12*(-1 + s2 - t1 + t2)*
                  (Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                    Power(s,2)*(-1 + s1)*
                     (t1*(-1 + t2) + (s1 - t2)*t2) - 
                    s*(-1 + s1)*
                     (Power(s1,2)*(-1 + t2) + 2*(1 + t1 - t2)*t2 + 
                       s2*(1 + t1*(-1 + t2) - 2*t2 - Power(t2,2)) + 
                       s1*(t1*(-3 + t2) - Power(-1 + t2,2) + 
                        s2*(1 + t2))))))/
             (s*(-1 + s1 + t1 - t2)*Power(1 - s2 + t1 - t2,3)) - 
            (2*(5 + 3*s - 8*s1 - 3*s2 - 3*t1 + t2)*
               (Power(s,3)*Power(t1 - t2,2)*(-1 + s1 + t1 - t2) + 
                 2*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                 2*Power(s,2)*
                  (-Power(t1,2) + 3*Power(t1,3) + t2 - 
                    8*Power(t1,2)*t2 + Power(t2,2) + 
                    7*t1*Power(t2,2) - 2*Power(t2,3) + 
                    Power(s1,2)*
                     (Power(t1,2) + (-1 + t2)*t2 - 
                       t1*(-2 + s2 + 2*t2)) - 
                    s2*(Power(t1,2)*(2 + t2) + 
                       t2*(2 + 4*t2 + Power(t2,2)) - 
                       t1*(1 + 6*t2 + 2*Power(t2,2))) + 
                    s1*(Power(t1,3) - Power(t1,2)*(-3 + s2 + 3*t2) + 
                       t2*(2*s2 + t2 + s2*t2 - Power(t2,2)) + 
                       t1*(-2 - 4*t2 + 3*Power(t2,2)))) + 
                 s*(-1 - 5*t1 - 3*Power(t1,2) + 9*Power(t1,3) + 
                    Power(s1,3)*
                     (Power(s2,2) + Power(t1,2) - 
                       2*s2*(1 + t1 - t2) - 2*t1*(-1 + t2) + 
                       2*Power(-1 + t2,2)) - 3*t2 - 10*t1*t2 - 
                    19*Power(t1,2)*t2 + 9*Power(t2,2) + 
                    15*t1*Power(t2,2) - 5*Power(t2,3) + 
                    Power(s1,2)*
                     (Power(t1,3) - 
                       2*s2*
                       (2 + Power(t1,2) - 2*t1*(-1 + t2) - 2*t2) + 
                       Power(s2,2)*(1 + t1 - t2) - 
                       3*Power(t1,2)*(-1 + t2) - 
                       2*Power(-1 + t2,2)*t2 + 
                       4*t1*(2 - 3*t2 + Power(t2,2))) + 
                    Power(s2,2)*
                     (-2 - 12*t2 - 7*Power(t2,2) - Power(t2,3) + 
                       t1*(2 + 6*t2 + Power(t2,2))) - 
                    2*s2*(-1 - 6*t2 + 5*Power(t2,2) + 2*Power(t2,3) + 
                       Power(t1,2)*(5 + 4*t2) - 
                       2*t1*(2 + 6*t2 + 3*Power(t2,2))) + 
                    s1*(6*Power(t1,3) + Power(t1,2)*(23 - 18*t2) - 
                       Power(-1 + t2,2)*(1 + 2*t2) + 
                       2*t1*(2 - 9*t2 + 7*Power(t2,2)) + 
                       Power(s2,2)*(6*t1 + t2*(4 + t2)) - 
                       2*s2*(-2 + 6*Power(t1,2) + Power(t2,2) + 
                        Power(t2,3) - t1*(-10 + 4*t2 + Power(t2,2))))))\
)/(s*(-1 + s1 + t1 - t2)*Power(-1 + s2 - t1 + t2,2))))/
        ((-1 + s2)*(-1 + t1)) + 
       (4*(64*s*(-1 + s1)*(-1 + s1 + t1 - t2) - 
            (2*(-1 - 3*t1 + s*t1 + t2 - s*t2 + s2*(2 + t2) + 
                 s1*(-2 + s2 - t1 + 2*t2))*
               (-14 + 2*Power(s,2) + 2*Power(s1,2) + 3*s2 + 
                 Power(s2,2) + 2*s2*t2 - 2*s*(1 + 2*s1 + s2 + 2*t2) + 
                 s1*(-1 + 3*s2 + 4*t2)))/(-1 + s2 - t1 + t2) - 
            32*(s*(13 + 12*Power(s1,2) - 13*t1 - Power(t1,2) + 
                  s1*(-24 + s2 + 12*t1 - 11*t2) + s2*(-1 + t1 - t2) + 
                  12*t2 + 2*t1*t2 - 2*Power(t2,2)) + 
               (-1 + s1)*(1 - t1 - t2 + s2*t2 + 
                  s1*(-2 + s2 - t1 + 2*t2))) + 
            16*(3*Power(s1,3) - s2 + 6*t1 - s2*t1 + 2*Power(t1,2) + 
               Power(s1,2)*(8*s2 - 3*(5 + t1 - 3*t2)) + 
               Power(s,2)*(-1 + 2*s1 + 2*t1 - t2) + 4*t2 - 8*s2*t2 + 
               Power(s2,2)*t2 - 2*s2*t1*t2 + 2*Power(t2,2) + 
               2*s2*Power(t2,2) + 
               s*(52 + 46*Power(s1,2) - 54*t1 - 8*Power(t1,2) + 
                  s1*(-88 + 5*s2 + 46*t1 - 37*t2) + 
                  s2*(-5 + 6*t1 - 6*t2) + 43*t2 + 18*t1*t2 - 
                  18*Power(t2,2)) + 
               s1*(17 + Power(s2,2) + 2*Power(t1,2) - 
                  s2*(9 + t1 - 8*t2) - 22*t2 + 4*Power(t2,2) - 
                  t1*(1 + 6*t2))) + 
            4*(-12 + 8*Power(s1,3) + 2*t1 - 4*s2*t1 + 6*Power(t1,2) + 
               Power(s1,2)*(13*s2 - 2*(8 + t1 - 7*t2)) - 37*t2 - 
               3*s2*t2 + Power(s2,2)*t2 - t1*t2 - 4*s2*t1*t2 + 
               13*Power(t2,2) + 4*s2*Power(t2,2) + 
               s*(57 + 42*Power(s1,2) - 59*t1 - 10*Power(t1,2) + 
                  s1*(-70 + 7*s2 + 44*t1 - 22*t2) + 
                  s2*(-2 + 6*t1 - 4*t2) + 27*t2 + 32*t1*t2 - 
                  36*Power(t2,2)) + 
               2*Power(s,2)*(5*s1 + 6*t1 - 4*(1 + t2)) + 
               s1*(11 + Power(s2,2) + 4*Power(t1,2) + t1*(3 - 12*t2) - 
                  55*t2 + 8*Power(t2,2) + s2*(-19 - 2*t1 + 15*t2))) - 
            8*(3 + 10*Power(s1,3) + s2 + 9*t1 - 5*s2*t1 + 
               8*Power(t1,2) + Power(s,2)*(-6 + 9*s1 + 10*t1 - 7*t2) - 
               5*t2 - 15*s2*t2 + 2*Power(s2,2)*t2 - 2*t1*t2 - 
               6*s2*t1*t2 + 12*Power(t2,2) + 6*s2*Power(t2,2) + 
               Power(s1,2)*(-27 + 18*s2 - 3*t1 + 16*t2) + 
               s*(81 + 73*Power(s1,2) - 93*t1 - 17*Power(t1,2) + 
                  s1*(-133 + 10*s2 + 74*t1 - 47*t2) + 
                  s2*(-6 + 11*t1 - 9*t2) + 62*t2 + 44*t1*t2 - 
                  46*Power(t2,2)) + 
               s1*(29 + 2*Power(s2,2) + 6*Power(t1,2) - 62*t2 - 
                  18*t1*t2 + 12*Power(t2,2) + s2*(-22 - 3*t1 + 18*t2))) \
+ ((-1 - 3*t1 + s*t1 + t2 - s*t2 + s2*(2 + t2) + 
                 s1*(-2 + s2 - t1 + 2*t2))*
               (4*s*(-1 + s1 + t1 - t2)*
                  Power(-1 - 3*t1 + s*t1 + t2 - s*t2 + 
                    s2*(2 + t2) + s1*(-2 + s2 - t1 + 2*t2),2) + 
                 12*(-1 + s2 - t1 + t2)*
                  (Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                    Power(s,2)*(-1 + s1)*
                     (t1*(-1 + t2) + (s1 - t2)*t2) - 
                    s*(-1 + s1)*
                     (Power(s1,2)*(-1 + t2) + 2*(1 + t1 - t2)*t2 + 
                       s2*(1 + t1*(-1 + t2) - 2*t2 - Power(t2,2)) + 
                       s1*(t1*(-3 + t2) - Power(-1 + t2,2) + 
                        s2*(1 + t2))))))/
             (s*(-1 + s1 + t1 - t2)*Power(1 - s2 + t1 - t2,3)) + 
            (2*(-2*s + 2*s1 + s2 + t2)*
               (2*Power(s,3)*Power(t1 - t2,2)*(-1 + s1 + t1 - t2) + 
                 4*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                 4*Power(s,2)*
                  (-Power(t1,2) + 3*Power(t1,3) + t2 - 
                    8*Power(t1,2)*t2 + Power(t2,2) + 
                    7*t1*Power(t2,2) - 2*Power(t2,3) + 
                    Power(s1,2)*
                     (Power(t1,2) + (-1 + t2)*t2 - 
                       t1*(-2 + s2 + 2*t2)) - 
                    s2*(Power(t1,2)*(2 + t2) + 
                       t2*(2 + 4*t2 + Power(t2,2)) - 
                       t1*(1 + 6*t2 + 2*Power(t2,2))) + 
                    s1*(Power(t1,3) - Power(t1,2)*(-3 + s2 + 3*t2) + 
                       t2*(2*s2 + t2 + s2*t2 - Power(t2,2)) + 
                       t1*(-2 - 4*t2 + 3*Power(t2,2)))) + 
                 2*s*(-1 - 5*t1 - 3*Power(t1,2) + 9*Power(t1,3) + 
                    Power(s1,3)*
                     (Power(s2,2) + Power(t1,2) - 
                       2*s2*(1 + t1 - t2) - 2*t1*(-1 + t2) + 
                       2*Power(-1 + t2,2)) - 3*t2 - 10*t1*t2 - 
                    19*Power(t1,2)*t2 + 9*Power(t2,2) + 
                    15*t1*Power(t2,2) - 5*Power(t2,3) + 
                    Power(s1,2)*
                     (Power(t1,3) - 
                       2*s2*
                       (2 + Power(t1,2) - 2*t1*(-1 + t2) - 2*t2) + 
                       Power(s2,2)*(1 + t1 - t2) - 
                       3*Power(t1,2)*(-1 + t2) - 
                       2*Power(-1 + t2,2)*t2 + 
                       4*t1*(2 - 3*t2 + Power(t2,2))) + 
                    Power(s2,2)*
                     (-2 - 12*t2 - 7*Power(t2,2) - Power(t2,3) + 
                       t1*(2 + 6*t2 + Power(t2,2))) - 
                    2*s2*(-1 - 6*t2 + 5*Power(t2,2) + 2*Power(t2,3) + 
                       Power(t1,2)*(5 + 4*t2) - 
                       2*t1*(2 + 6*t2 + 3*Power(t2,2))) + 
                    s1*(6*Power(t1,3) + Power(t1,2)*(23 - 18*t2) - 
                       Power(-1 + t2,2)*(1 + 2*t2) + 
                       2*t1*(2 - 9*t2 + 7*Power(t2,2)) + 
                       Power(s2,2)*(6*t1 + t2*(4 + t2)) - 
                       2*s2*(-2 + 6*Power(t1,2) + Power(t2,2) + 
                        Power(t2,3) - t1*(-10 + 4*t2 + Power(t2,2))))))\
)/(s*(-1 + s1 + t1 - t2)*Power(-1 + s2 - t1 + t2,2))))/
        ((-1 + s1)*(-s + s1 - t2)) + 
       ((-8*(-12 + 2*Power(s,2) + s1 + 2*Power(s1,2) + 3*s2 + 
               3*s1*s2 + Power(s2,2) - 4*t2 + 4*s1*t2 + 2*s2*t2 - 
               2*s*(4 + 2*s1 + s2 + 2*t2))*
             (-1 - 3*t1 + s*t1 + t2 - s*t2 + s2*(2 + t2) + 
               s1*(-2 + s2 - t1 + 2*t2)))/(-1 + s2 - t1 + t2) + 
          16*(-14 - 4*s2 + Power(s2,2) - 2*s2*t1 + 2*Power(t1,2) - 
             13*t2 + 2*s2*t2 + Power(s2,2)*t2 - 3*t1*t2 + 
             3*Power(t2,2) + 2*Power(s,2)*(2 + t2) + 
             Power(s1,2)*(-2 + s2 + 2*t2) + 
             s1*(-3 - 2*s2 + Power(s2,2) + t1 - 6*t2 + 3*s2*t2) - 
             s*(-5 + 5*s2 + 2*t1 + 2*t2 + 2*s2*t2 + s1*(5 + s2 + 4*t2))) \
+ (4*(-1 - 3*t1 + s*t1 + t2 - s*t2 + s2*(2 + t2) + 
               s1*(-2 + s2 - t1 + 2*t2))*
             (4*s*(-1 + s1 + t1 - t2)*
                Power(-1 - 3*t1 + s*t1 + t2 - s*t2 + s2*(2 + t2) + 
                  s1*(-2 + s2 - t1 + 2*t2),2) + 
               12*(-1 + s2 - t1 + t2)*
                (Power(-1 + s1*s2 + t1 - s1*t1 + t2 - s2*t2,2) + 
                  Power(s,2)*(-1 + s1)*(t1*(-1 + t2) + (s1 - t2)*t2) - 
                  s*(-1 + s1)*
                   (Power(s1,2)*(-1 + t2) + 2*(1 + t1 - t2)*t2 + 
                     s2*(1 + t1*(-1 + t2) - 2*t2 - Power(t2,2)) + 
                     s1*(t1*(-3 + t2) - Power(-1 + t2,2) + s2*(1 + t2))\
))))/(s*(-1 + s1 + t1 - t2)*Power(1 - s2 + t1 - t2,3)) + 
          (8*(-2*s + 2*s1 + s2 + t2)*
             (2*Power(s,3)*Power(t1 - t2,2)*(-1 + s1 + t1 - t2) + 
               4*(-1 + s2 - t1 + t2)*
                Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
               4*Power(s,2)*(-Power(t1,2) + 3*Power(t1,3) + t2 - 
                  8*Power(t1,2)*t2 + Power(t2,2) + 7*t1*Power(t2,2) - 
                  2*Power(t2,3) + 
                  Power(s1,2)*
                   (Power(t1,2) + (-1 + t2)*t2 - t1*(-2 + s2 + 2*t2)) \
- s2*(Power(t1,2)*(2 + t2) + t2*(2 + 4*t2 + Power(t2,2)) - 
                     t1*(1 + 6*t2 + 2*Power(t2,2))) + 
                  s1*(Power(t1,3) - Power(t1,2)*(-3 + s2 + 3*t2) + 
                     t2*(2*s2 + t2 + s2*t2 - Power(t2,2)) + 
                     t1*(-2 - 4*t2 + 3*Power(t2,2)))) + 
               2*s*(-1 - 5*t1 - 3*Power(t1,2) + 9*Power(t1,3) + 
                  Power(s1,3)*
                   (Power(s2,2) + Power(t1,2) - 2*s2*(1 + t1 - t2) - 
                     2*t1*(-1 + t2) + 2*Power(-1 + t2,2)) - 3*t2 - 
                  10*t1*t2 - 19*Power(t1,2)*t2 + 9*Power(t2,2) + 
                  15*t1*Power(t2,2) - 5*Power(t2,3) + 
                  Power(s1,2)*
                   (Power(t1,3) - 
                     2*s2*(2 + Power(t1,2) - 2*t1*(-1 + t2) - 2*t2) + 
                     Power(s2,2)*(1 + t1 - t2) - 
                     3*Power(t1,2)*(-1 + t2) - 2*Power(-1 + t2,2)*t2 + 
                     4*t1*(2 - 3*t2 + Power(t2,2))) + 
                  Power(s2,2)*
                   (-2 - 12*t2 - 7*Power(t2,2) - Power(t2,3) + 
                     t1*(2 + 6*t2 + Power(t2,2))) - 
                  2*s2*(-1 - 6*t2 + 5*Power(t2,2) + 2*Power(t2,3) + 
                     Power(t1,2)*(5 + 4*t2) - 
                     2*t1*(2 + 6*t2 + 3*Power(t2,2))) + 
                  s1*(6*Power(t1,3) + Power(t1,2)*(23 - 18*t2) - 
                     Power(-1 + t2,2)*(1 + 2*t2) + 
                     2*t1*(2 - 9*t2 + 7*Power(t2,2)) + 
                     Power(s2,2)*(6*t1 + t2*(4 + t2)) - 
                     2*s2*(-2 + 6*Power(t1,2) + Power(t2,2) + 
                        Power(t2,3) - t1*(-10 + 4*t2 + Power(t2,2)))))))/
           (s*(-1 + s1 + t1 - t2)*Power(-1 + s2 - t1 + t2,2)))/
        ((-1 + s1)*(-1 + t2)))*B3(s,1 - s1 - t1 + t2,2 + s - s1 - s2))/
   (16.*Power(Pi,2)) + (((8*s*(-1 + s2 - t1 + t2)*
          (1 - t1 + s*t1 - s1*(2 + s2 + t1) + 3*t2 - s*t2 + 
            s2*(2 - 2*t1 + t2))*
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
       (16*(-6 - 17*s2 - 10*Power(s2,2) + 2*Power(s2,3) + 6*t1 - 
             2*s2*t1 + 6*Power(s2,2)*t1 - 4*Power(t1,2) + 
             2*Power(s1,2)*(1 + s2 + t1) + 
             Power(s,2)*(-1 + 2*s2 + 2*t1) - 
             s*(-9 - 10*s2 + 4*Power(s2,2) + t1 + 8*s2*t1 + 
                s1*(2 + 4*s2 + 4*t1) - 8*t2) + 
             s1*(-3*s2 + 4*Power(s2,2) - 3*t1 + 8*s2*t1 - 2*t2) - 
             13*t2 - 6*s2*t2 + 3*t1*t2 - Power(t2,2)) + 
          (16*(-5 + Power(s,2) + s1 + Power(s1,2) - 5*s2 + 4*s1*s2 + 
               3*Power(s2,2) - 2*t1 + 2*s1*t1 + 3*s2*t1 - 
               s*(-1 + 2*s1 + 4*s2 + 2*t1) - 2*t2)*
             (1 - t1 + s*t1 - s1*(2 + s2 + t1) + 3*t2 - s*t2 + 
               s2*(2 - 2*t1 + t2)))/(-1 + s1 + t1 - t2) - 
          (2*(8*s - 4*(-2 + 2*s1 + 3*s2 + t1))*
             (-2*Power(s,3)*Power(t1 - t2,2)*(1 - s2 + t1 - t2) + 
               4*(-1 + s1 + t1 - t2)*
                Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
               4*Power(s,2)*
                ((2 + s1 + s2)*Power(t1,3) - 
                  Power(t1,2)*
                   (1 + s2 + Power(s2,2) + 7*t2 + 3*s2*t2 + 
                     s1*(-4 + s2 + 2*t2)) + 
                  t1*(-1 + 8*Power(t2,2) + Power(s2,2)*(1 + 2*t2) + 
                     s2*t2*(4 + 3*t2) + 
                     s1*(2 - 2*s2 - 6*t2 + Power(t2,2))) + 
                  t2*(t2 - 3*Power(t2,2) - Power(s2,2)*(2 + t2) + 
                     s1*(-1 + Power(s2,2) + 2*t2 + s2*t2) - 
                     s2*(-2 + 3*t2 + Power(t2,2)))) + 
               2*s*(-1 - 3*t1 + 9*Power(t1,2) - 5*Power(t1,3) - 
                  5*t2 - 10*t1*t2 + 15*Power(t1,2)*t2 - 
                  3*Power(t2,2) - 19*t1*Power(t2,2) + 9*Power(t2,3) + 
                  Power(s2,3)*
                   (2 + 2*Power(t1,2) + 2*t2 + Power(t2,2) - 
                     2*t1*(2 + t2)) + 
                  Power(s1,2)*
                   (Power(s2,3) - Power(t1,3) + 
                     Power(t1,2)*(-7 + t2) + 6*t1*(-2 + t2) + 
                     2*(-1 + t2) + Power(s2,2)*(1 - t1 + t2) + 
                     s2*(4*t1 + Power(t1,2) + 6*t2)) + 
                  s2*(-1 - 2*Power(t1,3) + 4*t2 + 23*Power(t2,2) + 
                     6*Power(t2,3) - 18*t1*t2*(1 + t2) + 
                     Power(t1,2)*(3 + 14*t2)) + 
                  Power(s2,2)*
                   (-2*Power(t1,3) + 4*Power(t1,2)*(1 + t2) + 
                     t2*(8 + 3*t2 + Power(t2,2)) - 
                     t1*(2 + 12*t2 + 3*Power(t2,2))) + 
                  2*s1*(1 - 2*Power(t1,3) + 
                     Power(s2,3)*(-1 + t1 - t2) + 4*t2 - 
                     5*Power(t2,2) + Power(t1,2)*(-5 + 6*t2) + 
                     t1*(6 + 12*t2 - 4*Power(t2,2)) - 
                     s2*(-2 + Power(t1,3) - Power(t1,2)*(-1 + t2) + 
                        10*t2 - 4*t1*t2 + 6*Power(t2,2)) - 
                     Power(s2,2)*
                      (2 + 2*t2 + Power(t2,2) - 2*t1*(1 + t2))))))/
           (s*Power(-1 + s1 + t1 - t2,2)*(-1 + s2 - t1 + t2)) - 
          (4*(-1 + t1 - s*t1 + s1*(2 + s2 + t1) + 
               s2*(-2 + 2*t1 - t2) - 3*t2 + s*t2)*
             (4*s*(-1 + s2 - t1 + t2)*
                Power(-1 + t1 - s*t1 + s1*(2 + s2 + t1) + 
                  s2*(-2 + 2*t1 - t2) - 3*t2 + s*t2,2) + 
               12*(-1 + s1 + t1 - t2)*
                (Power(-1 + s1*s2 + t1 - s1*t1 + t2 - s2*t2,2) + 
                  Power(s,2)*(-1 + s2)*
                   (s2*t1 - Power(t1,2) - t2 + t1*t2) - 
                  s*(-1 + s2)*
                   (Power(s2,2)*(-1 + t1) + 2*t1*(1 - t1 + t2) + 
                     s1*(1 + s2 - 2*t1 + s2*t1 - Power(t1,2) - t2 + 
                        t1*t2) - 
                     s2*(1 + Power(t1,2) + 3*t2 - t1*(2 + t2))))))/
           (s*Power(-1 + s1 + t1 - t2,3)*(-1 + s2 - t1 + t2)))/
        ((-1 + s2)*(-1 + t1)) + 
       (2*(128*s*(-1 + s + t1)*(1 - s2 + t1 - t2) - 
            (8*(-1 + t1 - s*t1 + s1*(2 + s2 + t1) + 
                 s2*(-2 + 2*t1 - t2) - 3*t2 + s*t2)*
               (-4 + 3*Power(s,2) + Power(s1,2) - 9*s2 + 
                 6*Power(s2,2) - t1 + 3*s2*t2 + 
                 2*s1*(-2 + 3*s2 + t2) - s*(-9 + 4*s1 + 9*s2 + 3*t2)))/
             (-1 + s1 + t1 - t2) + 
            64*((-1 + s1 + 2*s2)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
               Power(s,2)*(-13 + 14*s2 - 15*t1 + 16*t2) + 
               s*(17 - 2*Power(s2,2) - 3*t1 - 16*Power(t1,2) + 
                  s2*(-17 + 22*t1 - 6*t2) - 14*t2 + 19*t1*t2 - 
                  2*Power(t2,2) + s1*(-2 + 3*s2 + t2))) - 
            32*(-2 + Power(s,3) + 6*s2 - 23*Power(s2,2) + 
               20*Power(s2,3) - 4*t1 + 26*s2*t1 - 38*Power(s2,2)*t1 - 
               6*Power(t1,2) + 18*s2*Power(t1,2) + 2*t2 - 8*s2*t2 + 
               23*Power(s2,2)*t2 + 4*t1*t2 - 21*s2*t1*t2 + 
               3*s2*Power(t2,2) + 
               Power(s,2)*(-52 - 2*s1 + 65*s2 - 76*t1 + 85*t2) + 
               s*(89 + Power(s1,2) - 19*Power(s2,2) - 16*t1 - 
                  85*Power(t1,2) + s2*(-97 + 140*t1 - 59*t2) - 
                  71*t2 + 114*t1*t2 - 21*Power(t2,2) + 
                  s1*(-13 + 28*s2 + t1 + 8*t2)) + 
               s1*(-1 + 10*Power(s2,2) + 9*Power(t1,2) - t2 + 
                  Power(t2,2) - 2*t1*(-6 + 5*t2) + 
                  s2*(-12 - 19*t1 + 11*t2))) - 
            8*(-22 + 4*Power(s,3) + 9*s2 - 4*Power(s2,2) + 
               44*Power(s2,3) + 12*t1 + 5*s2*t1 - 88*Power(s2,2)*t1 + 
               48*s2*Power(t1,2) + Power(s1,2)*(7 + 2*s2 + 2*t1) - 
               2*Power(s,2)*(17 + 4*s1 - 38*s2 + 50*t1 - 61*t2) + 
               21*t2 + 34*s2*t2 + 60*Power(s2,2)*t2 - 29*t1*t2 - 
               66*s2*t1*t2 + 21*Power(t2,2) + 18*s2*Power(t2,2) + 
               2*s1*(-6 + 9*Power(s2,2) + 19*t1 + 12*Power(t1,2) - 
                  6*t2 - 14*t1*t2 + 2*Power(t2,2) + 
                  s2*(-13 - 17*t1 + 9*t2)) + 
               s*(138 + 4*Power(s1,2) - 32*Power(s2,2) + t1 - 
                  144*Power(t1,2) + s2*(-199 + 246*t1 - 134*t2) - 
                  124*t2 + 220*t1*t2 - 62*Power(t2,2) + 
                  2*s1*(-15 + 36*s2 - 6*t1 + 14*t2))) + 
            16*(-7 + 3*Power(s,3) + 6*s2 - 44*Power(s2,2) + 
               58*Power(s2,3) - 5*t1 + 35*s2*t1 - 104*Power(s2,2)*t1 - 
               8*Power(t1,2) + 52*s2*Power(t1,2) + 9*t2 + 6*s2*t2 + 
               71*Power(s2,2)*t2 - 10*t1*t2 - 67*s2*t1*t2 + 
               10*Power(t2,2) + 15*s2*Power(t2,2) + 
               Power(s1,2)*(3*s2 + t1 + t2) + 
               Power(s,2)*(-85 - 6*s1 + 126*s2 - 150*t1 + 178*t2) + 
               s*(187 + 3*Power(s1,2) - 55*Power(s2,2) - 15*t1 - 
                  186*Power(t1,2) + s2*(-216 + 323*t1 - 167*t2) - 
                  157*t2 + 270*t1*t2 - 65*Power(t2,2) + 
                  s1*(-20 + 72*s2 - 4*t1 + 22*t2)) + 
               s1*(-11 + 31*Power(s2,2) + 26*Power(t1,2) + 
                  t1*(39 - 30*t2) - 5*t2 + 4*Power(t2,2) + 
                  s2*(-44 - 47*t1 + 31*t2))) - 
            (4*(3 + 3*s - 2*s1 - 4*s2 - t2)*
               (-2*Power(s,3)*Power(t1 - t2,2)*(1 - s2 + t1 - t2) + 
                 4*(-1 + s1 + t1 - t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 4*Power(s,2)*
                  ((2 + s1 + s2)*Power(t1,3) - 
                    Power(t1,2)*
                     (1 + s2 + Power(s2,2) + 7*t2 + 3*s2*t2 + 
                       s1*(-4 + s2 + 2*t2)) + 
                    t1*(-1 + 8*Power(t2,2) + 
                       Power(s2,2)*(1 + 2*t2) + s2*t2*(4 + 3*t2) + 
                       s1*(2 - 2*s2 - 6*t2 + Power(t2,2))) + 
                    t2*(t2 - 3*Power(t2,2) - Power(s2,2)*(2 + t2) + 
                       s1*(-1 + Power(s2,2) + 2*t2 + s2*t2) - 
                       s2*(-2 + 3*t2 + Power(t2,2)))) + 
                 2*s*(-1 - 3*t1 + 9*Power(t1,2) - 5*Power(t1,3) - 
                    5*t2 - 10*t1*t2 + 15*Power(t1,2)*t2 - 
                    3*Power(t2,2) - 19*t1*Power(t2,2) + 
                    9*Power(t2,3) + 
                    Power(s2,3)*
                     (2 + 2*Power(t1,2) + 2*t2 + Power(t2,2) - 
                       2*t1*(2 + t2)) + 
                    Power(s1,2)*
                     (Power(s2,3) - Power(t1,3) + 
                       Power(t1,2)*(-7 + t2) + 6*t1*(-2 + t2) + 
                       2*(-1 + t2) + Power(s2,2)*(1 - t1 + t2) + 
                       s2*(4*t1 + Power(t1,2) + 6*t2)) + 
                    s2*(-1 - 2*Power(t1,3) + 4*t2 + 23*Power(t2,2) + 
                       6*Power(t2,3) - 18*t1*t2*(1 + t2) + 
                       Power(t1,2)*(3 + 14*t2)) + 
                    Power(s2,2)*
                     (-2*Power(t1,3) + 4*Power(t1,2)*(1 + t2) + 
                       t2*(8 + 3*t2 + Power(t2,2)) - 
                       t1*(2 + 12*t2 + 3*Power(t2,2))) + 
                    2*s1*(1 - 2*Power(t1,3) + 
                       Power(s2,3)*(-1 + t1 - t2) + 4*t2 - 
                       5*Power(t2,2) + Power(t1,2)*(-5 + 6*t2) + 
                       t1*(6 + 12*t2 - 4*Power(t2,2)) - 
                       s2*(-2 + Power(t1,3) - 
                        Power(t1,2)*(-1 + t2) + 10*t2 - 4*t1*t2 + 
                        6*Power(t2,2)) - 
                       Power(s2,2)*
                        (2 + 2*t2 + Power(t2,2) - 2*t1*(1 + t2))))))/
             (s*Power(-1 + s1 + t1 - t2,2)*(-1 + s2 - t1 + t2)) - 
            (2*(-1 + t1 - s*t1 + s1*(2 + s2 + t1) + 
                 s2*(-2 + 2*t1 - t2) - 3*t2 + s*t2)*
               (4*s*(-1 + s2 - t1 + t2)*
                  Power(-1 + t1 - s*t1 + s1*(2 + s2 + t1) + 
                    s2*(-2 + 2*t1 - t2) - 3*t2 + s*t2,2) + 
                 12*(-1 + s1 + t1 - t2)*
                  (Power(-1 + s1*s2 + t1 - s1*t1 + t2 - s2*t2,2) + 
                    Power(s,2)*(-1 + s2)*
                     (s2*t1 - Power(t1,2) - t2 + t1*t2) - 
                    s*(-1 + s2)*
                     (Power(s2,2)*(-1 + t1) + 2*t1*(1 - t1 + t2) + 
                       s1*(1 + s2 - 2*t1 + s2*t1 - Power(t1,2) - t2 + 
                        t1*t2) - 
                       s2*(1 + Power(t1,2) + 3*t2 - t1*(2 + t2))))))/
             (s*Power(-1 + s1 + t1 - t2,3)*(-1 + s2 - t1 + t2))))/
        ((s - s2 + t1)*(s - s1 + t2)) + 
       (2*(128*s*(t1 - t2)*(1 - s2 + t1 - t2) - 
            (8*(-1 + t1 - s*t1 + s1*(2 + s2 + t1) + 
                 s2*(-2 + 2*t1 - t2) - 3*t2 + s*t2)*
               (3*Power(s,2) + Power(s1,2) - 10*s2 + 6*Power(s2,2) - 
                 4*t2 + 3*s2*t2 + 2*s1*(-2 + 3*s2 + t2) - 
                 s*(-5 + 4*s1 + 9*s2 + 3*t2)))/(-1 + s1 + t1 - t2) + 
            64*(Power(s,2)*(-1 + t2) - 
               (3 - 2*t1 + s1*(-1 + s2 + 2*t1) + 
                  s2*(-3 + 4*t1 - 2*t2) - 2*t2)*(-1 + s2 - t1 + t2) + 
               s*(2 - 14*t1 - 18*Power(t1,2) + 
                  2*s2*(-1 + 9*t1 - 9*t2) + 15*t2 + 34*t1*t2 - 
                  17*Power(t2,2) + s1*(-1 + 2*s2 - 2*t1 + t2))) + 
            32*(-24 + Power(s,3) + 50*s2 - 28*Power(s2,2) - 12*t1 - 
               15*s2*t1 + 36*Power(s2,2)*t1 + 8*Power(t1,2) - 
               36*s2*Power(t1,2) + 52*t2 - 37*s2*t2 - 
               18*Power(s2,2)*t2 + 16*t1*t2 + 54*s2*t1*t2 - 
               28*Power(t2,2) - 18*s2*Power(t2,2) + 
               Power(s1,2)*(-1 + s2 - t1 + t2) - 
               Power(s,2)*(-10 + s1 + 4*s2 + 11*t2) + 
               s1*(6 + 11*Power(s2,2) - 18*Power(t1,2) - 5*t2 + 
                  s2*(-19 + 7*t1 + 11*t2) + t1*(-11 + 18*t2)) + 
               s*(-23 + 2*Power(s2,2) + 66*t1 + 103*Power(t1,2) + 
                  s1*(11 - 18*s2 + 19*t1 - 10*t2) - 69*t2 - 
                  188*t1*t2 + 93*Power(t2,2) + 
                  s2*(22 - 103*t1 + 103*t2))) + 
            8*(-45 + 4*Power(s,3) + 61*s2 - 65*Power(s2,2) - 
               4*Power(s2,3) - 35*t1 + 2*s2*t1 + 80*Power(s2,2)*t1 - 
               10*Power(t1,2) - 80*s2*Power(t1,2) + 118*t2 - 
               115*s2*t2 - 42*Power(s2,2)*t2 + 85*t1*t2 + 
               120*s2*t1*t2 - 93*Power(t2,2) - 40*s2*Power(t2,2) + 
               Power(s1,2)*(-11 + 2*s2 - 6*t1 + 4*t2) - 
               2*Power(s,2)*(-1 + 2*s1 + 6*s2 + 13*t2) + 
               2*s1*(4 + 12*Power(s2,2) - 20*Power(t1,2) + 6*t2 + 
                  2*s2*(-9 + 2*t1 + 7*t2) + t1*(-17 + 20*t2)) + 
               s*(-47 + 2*Power(s2,2) + 86*t1 + 176*Power(t1,2) + 
                  s1*(46 - 52*s2 + 50*t1 - 36*t2) - 60*t2 - 
                  312*t1*t2 + 150*Power(t2,2) + 
                  s2*(101 - 176*t1 + 168*t2))) - 
            16*(-61 + 5*Power(s,3) + 111*s2 - 65*Power(s2,2) - 
               6*Power(s2,3) - 39*t1 - 16*s2*t1 + 96*Power(s2,2)*t1 + 
               2*Power(t1,2) - 96*s2*Power(t1,2) + 145*t2 - 
               112*s2*t2 - 52*Power(s2,2)*t2 + 76*t1*t2 + 
               144*s2*t1*t2 - 96*Power(t2,2) - 48*s2*Power(t2,2) + 
               Power(s1,2)*(-4 + s2 - 5*t1 + 3*t2) - 
               2*Power(s,2)*(-10 + 3*s1 + 10*s2 + 17*t2) + 
               2*s1*(4 + 12*Power(s2,2) - 17*t1 - 24*Power(t1,2) + 
                  2*t2 + 24*t1*t2 + s2*(-19 + 7*t1 + 14*t2)) + 
               s*(-62 + Power(s1,2) + 14*Power(s2,2) + 127*t1 + 
                  230*Power(t1,2) + s1*(35 - 46*s2 + 55*t1 - 30*t2) - 
                  122*t2 - 412*t1*t2 + 201*Power(t2,2) + 
                  s2*(75 - 230*t1 + 232*t2))) + 
            (4*(-3 - 3*s + 2*s1 + 4*s2 + t2)*
               (-2*Power(s,3)*Power(t1 - t2,2)*(1 - s2 + t1 - t2) + 
                 4*(-1 + s1 + t1 - t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 4*Power(s,2)*
                  ((2 + s1 + s2)*Power(t1,3) - 
                    Power(t1,2)*
                     (1 + s2 + Power(s2,2) + 7*t2 + 3*s2*t2 + 
                       s1*(-4 + s2 + 2*t2)) + 
                    t1*(-1 + 8*Power(t2,2) + 
                       Power(s2,2)*(1 + 2*t2) + s2*t2*(4 + 3*t2) + 
                       s1*(2 - 2*s2 - 6*t2 + Power(t2,2))) + 
                    t2*(t2 - 3*Power(t2,2) - Power(s2,2)*(2 + t2) + 
                       s1*(-1 + Power(s2,2) + 2*t2 + s2*t2) - 
                       s2*(-2 + 3*t2 + Power(t2,2)))) + 
                 2*s*(-1 - 3*t1 + 9*Power(t1,2) - 5*Power(t1,3) - 
                    5*t2 - 10*t1*t2 + 15*Power(t1,2)*t2 - 
                    3*Power(t2,2) - 19*t1*Power(t2,2) + 
                    9*Power(t2,3) + 
                    Power(s2,3)*
                     (2 + 2*Power(t1,2) + 2*t2 + Power(t2,2) - 
                       2*t1*(2 + t2)) + 
                    Power(s1,2)*
                     (Power(s2,3) - Power(t1,3) + 
                       Power(t1,2)*(-7 + t2) + 6*t1*(-2 + t2) + 
                       2*(-1 + t2) + Power(s2,2)*(1 - t1 + t2) + 
                       s2*(4*t1 + Power(t1,2) + 6*t2)) + 
                    s2*(-1 - 2*Power(t1,3) + 4*t2 + 23*Power(t2,2) + 
                       6*Power(t2,3) - 18*t1*t2*(1 + t2) + 
                       Power(t1,2)*(3 + 14*t2)) + 
                    Power(s2,2)*
                     (-2*Power(t1,3) + 4*Power(t1,2)*(1 + t2) + 
                       t2*(8 + 3*t2 + Power(t2,2)) - 
                       t1*(2 + 12*t2 + 3*Power(t2,2))) + 
                    2*s1*(1 - 2*Power(t1,3) + 
                       Power(s2,3)*(-1 + t1 - t2) + 4*t2 - 
                       5*Power(t2,2) + Power(t1,2)*(-5 + 6*t2) + 
                       t1*(6 + 12*t2 - 4*Power(t2,2)) - 
                       s2*(-2 + Power(t1,3) - 
                        Power(t1,2)*(-1 + t2) + 10*t2 - 4*t1*t2 + 
                        6*Power(t2,2)) - 
                       Power(s2,2)*
                        (2 + 2*t2 + Power(t2,2) - 2*t1*(1 + t2))))))/
             (s*Power(-1 + s1 + t1 - t2,2)*(-1 + s2 - t1 + t2)) - 
            (2*(-1 + t1 - s*t1 + s1*(2 + s2 + t1) + 
                 s2*(-2 + 2*t1 - t2) - 3*t2 + s*t2)*
               (4*s*(-1 + s2 - t1 + t2)*
                  Power(-1 + t1 - s*t1 + s1*(2 + s2 + t1) + 
                    s2*(-2 + 2*t1 - t2) - 3*t2 + s*t2,2) + 
                 12*(-1 + s1 + t1 - t2)*
                  (Power(-1 + s1*s2 + t1 - s1*t1 + t2 - s2*t2,2) + 
                    Power(s,2)*(-1 + s2)*
                     (s2*t1 - Power(t1,2) - t2 + t1*t2) - 
                    s*(-1 + s2)*
                     (Power(s2,2)*(-1 + t1) + 2*t1*(1 - t1 + t2) + 
                       s1*(1 + s2 - 2*t1 + s2*t1 - Power(t1,2) - t2 + 
                        t1*t2) - 
                       s2*(1 + Power(t1,2) + 3*t2 - t1*(2 + t2))))))/
             (s*Power(-1 + s1 + t1 - t2,3)*(-1 + s2 - t1 + t2))))/
        ((-1 + s1)*(-1 + t2)) + 
       (4*(64*s*(-1 + s2)*(-1 + s2 - t1 + t2) + 
            (4*(-6 + 3*Power(s,2) + Power(s1,2) - 5*s2 + 
                 3*Power(s2,2) - 2*t1 + 3*s2*t1 + 2*s1*(2*s2 + t1) - 
                 s*(-4 + 4*s1 + 6*s2 + 3*t1))*
               (1 - t1 + s*t1 - s1*(2 + s2 + t1) + 3*t2 - s*t2 + 
                 s2*(2 - 2*t1 + t2)))/(-1 + s1 + t1 - t2) + 
            32*(Power(s,2)*(t1 - t2) + 
               (1 + s2)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
               s*(-15 - 15*Power(s2,2) - 11*t1 + 2*Power(t1,2) + 
                  s2*(30 + 11*t1 - 13*t2) + s1*(2 - 2*s2 + t1 - t2) + 
                  13*t2 - 3*t1*t2 + Power(t2,2))) - 
            4*(-14 - 14*s2 + 29*Power(s2,2) + 8*Power(s2,3) + 
               Power(s1,2)*(1 - 2*s2 - 2*t1) + 30*t1 - 60*s2*t1 + 
               2*Power(s2,2)*t1 + 44*Power(t1,2) - 10*s2*Power(t1,2) + 
               Power(s,2)*(23 + 4*s1 - 16*s2 + 38*t1 - 38*t2) + 5*t2 + 
               17*s2*t2 + 2*Power(s2,2)*t2 - 19*t1*t2 + 18*s2*t1*t2 - 
               15*Power(t2,2) - 8*s2*Power(t2,2) + 
               s*(-83 - 4*Power(s1,2) - 86*Power(s2,2) + 8*t1 + 
                  68*Power(t1,2) + s2*(153 - 28*t1 - 36*t2) + 
                  s1*(11 - 42*s2 + 6*t1 - 18*t2) + 57*t2 - 98*t1*t2 + 
                  30*Power(t2,2)) + 
               s1*(23 - 12*Power(t1,2) + 7*t2 + 
                  s2*(11 + 8*t1 + 4*t2) + t1*(-1 + 12*t2))) - 
            16*(-1 - 10*s2 + 6*Power(s2,2) + 7*Power(s2,3) + 8*t1 - 
               25*s2*t1 - 10*Power(s2,2)*t1 + 17*Power(t1,2) + 
               4*s2*Power(t1,2) + 
               Power(s,2)*(3 + s1 - 3*s2 + 13*t1 - 12*t2) + 13*s2*t2 + 
               6*Power(s2,2)*t2 - 14*t1*t2 - 3*s2*t1*t2 - 
               Power(t2,2) - s2*Power(t2,2) + 
               s1*(2 + t1 - 2*Power(t1,2) + s2*(-2 + 3*t1) + 
                  2*t1*t2) - 
               s*(70 + Power(s1,2) + 71*Power(s2,2) + 28*t1 - 
                  22*Power(t1,2) - 53*t2 + 32*t1*t2 - 10*Power(t2,2) + 
                  s1*(-17 + 18*s2 - 5*t1 + 8*t2) + 
                  s2*(-140 - 29*t1 + 52*t2))) - 
            8*(9 + Power(s,3) - Power(s1,2) + 26*s2 - 23*Power(s2,2) - 
               15*Power(s2,3) - 17*t1 + 82*s2*t1 + 8*Power(s2,2)*t1 - 
               54*Power(t1,2) + s2*Power(t1,2) - 2*t2 - 32*s2*t2 - 
               9*Power(s2,2)*t2 + 35*t1*t2 - 7*s2*t1*t2 + 
               9*Power(t2,2) + 6*s2*Power(t2,2) + 
               Power(s,2)*(-15 - 6*s1 + 11*s2 - 44*t1 + 39*t2) - 
               s1*(11 + 2*Power(s2,2) - 10*Power(t1,2) + 3*t2 + 
                  s2*(-3 + 14*t1 + 2*t2) + t1*(3 + 10*t2)) + 
               s*(124 + 5*Power(s1,2) + 136*Power(s2,2) + 11*t1 - 
                  70*Power(t1,2) - 90*t2 + 101*t1*t2 - 
                  31*Power(t2,2) + s1*(-34 + 53*s2 - 5*t1 + 21*t2) + 
                  s2*(-259 - 2*t1 + 78*t2))) - 
            (2*(2 + 3*s - 2*s1 - 3*s2 - t1)*
               (-2*Power(s,3)*Power(t1 - t2,2)*(1 - s2 + t1 - t2) + 
                 4*(-1 + s1 + t1 - t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 4*Power(s,2)*
                  ((2 + s1 + s2)*Power(t1,3) - 
                    Power(t1,2)*
                     (1 + s2 + Power(s2,2) + 7*t2 + 3*s2*t2 + 
                       s1*(-4 + s2 + 2*t2)) + 
                    t1*(-1 + 8*Power(t2,2) + 
                       Power(s2,2)*(1 + 2*t2) + s2*t2*(4 + 3*t2) + 
                       s1*(2 - 2*s2 - 6*t2 + Power(t2,2))) + 
                    t2*(t2 - 3*Power(t2,2) - Power(s2,2)*(2 + t2) + 
                       s1*(-1 + Power(s2,2) + 2*t2 + s2*t2) - 
                       s2*(-2 + 3*t2 + Power(t2,2)))) + 
                 2*s*(-1 - 3*t1 + 9*Power(t1,2) - 5*Power(t1,3) - 
                    5*t2 - 10*t1*t2 + 15*Power(t1,2)*t2 - 
                    3*Power(t2,2) - 19*t1*Power(t2,2) + 
                    9*Power(t2,3) + 
                    Power(s2,3)*
                     (2 + 2*Power(t1,2) + 2*t2 + Power(t2,2) - 
                       2*t1*(2 + t2)) + 
                    Power(s1,2)*
                     (Power(s2,3) - Power(t1,3) + 
                       Power(t1,2)*(-7 + t2) + 6*t1*(-2 + t2) + 
                       2*(-1 + t2) + Power(s2,2)*(1 - t1 + t2) + 
                       s2*(4*t1 + Power(t1,2) + 6*t2)) + 
                    s2*(-1 - 2*Power(t1,3) + 4*t2 + 23*Power(t2,2) + 
                       6*Power(t2,3) - 18*t1*t2*(1 + t2) + 
                       Power(t1,2)*(3 + 14*t2)) + 
                    Power(s2,2)*
                     (-2*Power(t1,3) + 4*Power(t1,2)*(1 + t2) + 
                       t2*(8 + 3*t2 + Power(t2,2)) - 
                       t1*(2 + 12*t2 + 3*Power(t2,2))) + 
                    2*s1*(1 - 2*Power(t1,3) + 
                       Power(s2,3)*(-1 + t1 - t2) + 4*t2 - 
                       5*Power(t2,2) + Power(t1,2)*(-5 + 6*t2) + 
                       t1*(6 + 12*t2 - 4*Power(t2,2)) - 
                       s2*(-2 + Power(t1,3) - 
                        Power(t1,2)*(-1 + t2) + 10*t2 - 4*t1*t2 + 
                        6*Power(t2,2)) - 
                       Power(s2,2)*
                        (2 + 2*t2 + Power(t2,2) - 2*t1*(1 + t2))))))/
             (s*Power(-1 + s1 + t1 - t2,2)*(-1 + s2 - t1 + t2)) - 
            ((-1 + t1 - s*t1 + s1*(2 + s2 + t1) + 
                 s2*(-2 + 2*t1 - t2) - 3*t2 + s*t2)*
               (4*s*(-1 + s2 - t1 + t2)*
                  Power(-1 + t1 - s*t1 + s1*(2 + s2 + t1) + 
                    s2*(-2 + 2*t1 - t2) - 3*t2 + s*t2,2) + 
                 12*(-1 + s1 + t1 - t2)*
                  (Power(-1 + s1*s2 + t1 - s1*t1 + t2 - s2*t2,2) + 
                    Power(s,2)*(-1 + s2)*
                     (s2*t1 - Power(t1,2) - t2 + t1*t2) - 
                    s*(-1 + s2)*
                     (Power(s2,2)*(-1 + t1) + 2*t1*(1 - t1 + t2) + 
                       s1*(1 + s2 - 2*t1 + s2*t1 - Power(t1,2) - t2 + 
                        t1*t2) - 
                       s2*(1 + Power(t1,2) + 3*t2 - t1*(2 + t2))))))/
             (s*Power(-1 + s1 + t1 - t2,3)*(-1 + s2 - t1 + t2))))/
        ((-1 + s2)*(-s + s2 - t1)) + 
       (4*(-64*s*(-1 + t1)*(1 - s2 + t1 - t2) - 
            (4*(-1 + t1 - s*t1 + s1*(2 + s2 + t1) + 
                 s2*(-2 + 2*t1 - t2) - 3*t2 + s*t2)*
               (-2 + Power(s,2) + Power(s1,2) - 9*s2 + 
                 6*Power(s2,2) - t1 - 4*t2 + 3*s2*t2 + 
                 2*s1*(-1 + 3*s2 + t2) - 2*s*(-1 + s1 + 3*s2 + t2)))/
             (-1 + s1 + t1 - t2) - 
            32*(-((-1 + s1 + 2*s2)*(-1 + t1)*(-1 + s2 - t1 + t2)) + 
               s*(14 - 2*t1 - 16*Power(t1,2) + 
                  s2*(-15 + 16*t1 - 2*t2) - 12*t2 + 19*t1*t2 - 
                  2*Power(t2,2) + s1*(-1 + s2 - t1 + t2))) - 
            16*(-7 + 21*s2 - 18*Power(s2,2) + Power(s2,3) + t1 - 
               6*s2*t1 + 18*Power(s2,2)*t1 + 6*Power(t1,2) - 
               18*s2*Power(t1,2) + Power(s,2)*(s2 - t2) + 5*t2 - 
               19*s2*t2 - 2*Power(s2,2)*t2 - 4*t1*t2 + 21*s2*t1*t2 - 
               3*s2*Power(t2,2) + Power(s1,2)*(-1 + s2 - t1 + t2) + 
               s1*(10 + 4*Power(s2,2) - 9*Power(t1,2) - 7*t2 - 
                  Power(t2,2) + 3*s2*(-3 + 2*t1 + t2) + 
                  t1*(-3 + 10*t2)) + 
               s*(-61 - 4*Power(s2,2) + 21*t1 + 85*Power(t1,2) + 
                  40*t2 - 114*t1*t2 + 21*Power(t2,2) - 
                  2*s1*(-5 + 6*s2 - 5*t1 + 5*t2) + 
                  s2*(76 - 85*t1 + 19*t2))) - 
            4*(-11 + 46*s2 - 16*Power(s2,2) - 5*t1 - s2*t1 + 
               48*Power(s2,2)*t1 - 48*s2*Power(t1,2) + 
               Power(s,2)*(1 + 2*s2 + 2*t1 - 8*t2) + 8*t2 - 46*s2*t2 - 
               16*Power(s2,2)*t2 + 29*t1*t2 + 66*s2*t1*t2 - 
               21*Power(t2,2) - 18*s2*Power(t2,2) + 
               Power(s1,2)*(-3 + 2*s2 - 6*t1 + 4*t2) - 
               s*(64 + 16*Power(s2,2) - 65*t1 - 140*Power(t1,2) + 
                  s2*(-89 + 140*t1 - 52*t2) + t2 + 214*t1*t2 - 
                  60*Power(t2,2) + s1*(-17 + 34*s2 - 26*t1 + 26*t2)) + 
               s1*(21 + 16*Power(s2,2) - 24*Power(t1,2) + 13*t2 - 
                  4*Power(t2,2) + s2*(9 + 4*t1 + 16*t2) + 
                  t1*(-27 + 28*t2))) + 
            8*(-15 + 45*s2 - 45*Power(s2,2) + 6*Power(s2,3) + 3*t1 - 
               10*s2*t1 + 52*Power(s2,2)*t1 + 8*Power(t1,2) - 
               52*s2*Power(t1,2) + Power(s,2)*(5*s2 + t1 - 5*t2) + 
               9*t2 - 61*s2*t2 - 9*Power(s2,2)*t2 + 10*t1*t2 + 
               67*s2*t1*t2 - 10*Power(t2,2) - 15*s2*Power(t2,2) + 
               Power(s1,2)*(-4 + 5*s2 - 5*t1 + 5*t2) + 
               s1*(25 + 22*Power(s2,2) - 26*Power(t1,2) - 7*t2 - 
                  4*Power(t2,2) + 3*t1*(-7 + 10*t2) + 
                  s2*(-15 + 10*t1 + 18*t2)) + 
               s*(-104 - 22*Power(s2,2) + 68*t1 + 184*Power(t1,2) + 
                  s1*(26 - 41*s2 + 29*t1 - 31*t2) + 42*t2 - 
                  267*t1*t2 + 64*Power(t2,2) + 
                  s2*(150 - 184*t1 + 52*t2))) - 
            (2*(3 + 2*s - 2*s1 - 4*s2 - t2)*
               (-2*Power(s,3)*Power(t1 - t2,2)*(1 - s2 + t1 - t2) + 
                 4*(-1 + s1 + t1 - t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 4*Power(s,2)*
                  ((2 + s1 + s2)*Power(t1,3) - 
                    Power(t1,2)*
                     (1 + s2 + Power(s2,2) + 7*t2 + 3*s2*t2 + 
                       s1*(-4 + s2 + 2*t2)) + 
                    t1*(-1 + 8*Power(t2,2) + 
                       Power(s2,2)*(1 + 2*t2) + s2*t2*(4 + 3*t2) + 
                       s1*(2 - 2*s2 - 6*t2 + Power(t2,2))) + 
                    t2*(t2 - 3*Power(t2,2) - Power(s2,2)*(2 + t2) + 
                       s1*(-1 + Power(s2,2) + 2*t2 + s2*t2) - 
                       s2*(-2 + 3*t2 + Power(t2,2)))) + 
                 2*s*(-1 - 3*t1 + 9*Power(t1,2) - 5*Power(t1,3) - 
                    5*t2 - 10*t1*t2 + 15*Power(t1,2)*t2 - 
                    3*Power(t2,2) - 19*t1*Power(t2,2) + 
                    9*Power(t2,3) + 
                    Power(s2,3)*
                     (2 + 2*Power(t1,2) + 2*t2 + Power(t2,2) - 
                       2*t1*(2 + t2)) + 
                    Power(s1,2)*
                     (Power(s2,3) - Power(t1,3) + 
                       Power(t1,2)*(-7 + t2) + 6*t1*(-2 + t2) + 
                       2*(-1 + t2) + Power(s2,2)*(1 - t1 + t2) + 
                       s2*(4*t1 + Power(t1,2) + 6*t2)) + 
                    s2*(-1 - 2*Power(t1,3) + 4*t2 + 23*Power(t2,2) + 
                       6*Power(t2,3) - 18*t1*t2*(1 + t2) + 
                       Power(t1,2)*(3 + 14*t2)) + 
                    Power(s2,2)*
                     (-2*Power(t1,3) + 4*Power(t1,2)*(1 + t2) + 
                       t2*(8 + 3*t2 + Power(t2,2)) - 
                       t1*(2 + 12*t2 + 3*Power(t2,2))) + 
                    2*s1*(1 - 2*Power(t1,3) + 
                       Power(s2,3)*(-1 + t1 - t2) + 4*t2 - 
                       5*Power(t2,2) + Power(t1,2)*(-5 + 6*t2) + 
                       t1*(6 + 12*t2 - 4*Power(t2,2)) - 
                       s2*(-2 + Power(t1,3) - 
                        Power(t1,2)*(-1 + t2) + 10*t2 - 4*t1*t2 + 
                        6*Power(t2,2)) - 
                       Power(s2,2)*
                        (2 + 2*t2 + Power(t2,2) - 2*t1*(1 + t2))))))/
             (s*Power(-1 + s1 + t1 - t2,2)*(-1 + s2 - t1 + t2)) - 
            ((-1 + t1 - s*t1 + s1*(2 + s2 + t1) + 
                 s2*(-2 + 2*t1 - t2) - 3*t2 + s*t2)*
               (4*s*(-1 + s2 - t1 + t2)*
                  Power(-1 + t1 - s*t1 + s1*(2 + s2 + t1) + 
                    s2*(-2 + 2*t1 - t2) - 3*t2 + s*t2,2) + 
                 12*(-1 + s1 + t1 - t2)*
                  (Power(-1 + s1*s2 + t1 - s1*t1 + t2 - s2*t2,2) + 
                    Power(s,2)*(-1 + s2)*
                     (s2*t1 - Power(t1,2) - t2 + t1*t2) - 
                    s*(-1 + s2)*
                     (Power(s2,2)*(-1 + t1) + 2*t1*(1 - t1 + t2) + 
                       s1*(1 + s2 - 2*t1 + s2*t1 - Power(t1,2) - t2 + 
                        t1*t2) - 
                       s2*(1 + Power(t1,2) + 3*t2 - t1*(2 + t2))))))/
             (s*Power(-1 + s1 + t1 - t2,3)*(-1 + s2 - t1 + t2))))/
        ((-1 + t1)*(-1 + t2)) + 
       (4*(64*s*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) - 
            (4*(-1 + t1 - s*t1 + s1*(2 + s2 + t1) + 
                 s2*(-2 + 2*t1 - t2) - 3*t2 + s*t2)*
               (-1 + 3*Power(s,2) + Power(s1,2) - 10*s2 + 
                 6*Power(s2,2) - 2*t2 + 3*s2*t2 + 
                 s1*(-5 + 6*s2 + 2*t2) - s*(-8 + 4*s1 + 9*s2 + 3*t2)))/
             (-1 + s1 + t1 - t2) - 
            32*(-1 + 3*s2 - 3*Power(s2,2) + 2*Power(s2,3) + t1 - 
               4*Power(s2,2)*t1 - 2*Power(t1,2) + 4*s2*Power(t1,2) + 
               Power(s1,2)*(s2 + t1) + Power(s,2)*(-1 + t2) - t2 + 
               2*s2*t2 + 3*Power(s2,2)*t2 - 6*s2*t1*t2 + 
               2*Power(t2,2) + 2*s2*Power(t2,2) + 
               2*s1*(1 + Power(t1,2) + s2*(-2 + t1 - t2) - t2 - 
                  t1*t2) + s*
                (14 - Power(s2,2) - 18*Power(t1,2) + 
                  2*s2*(-7 + 9*t1 - 9*t2) + 34*t1*t2 - 17*Power(t2,2) + 
                  2*s1*(-7 + 8*s2 - 8*t1 + 8*t2))) - 
            16*(4 + Power(s,3) - 26*s2 + 23*Power(s2,2) - 
               18*Power(s2,3) + 8*s2*t1 + 36*Power(s2,2)*t1 + 
               8*Power(t1,2) - 36*s2*Power(t1,2) - 
               2*s1*(8 + Power(s2,2) + 9*Power(t1,2) + t1*(8 - 9*t2) + 
                  s2*(-23 + 8*t1 - 8*t2) - 18*t2) + 4*t2 - 24*s2*t2 - 
               27*Power(s2,2)*t2 + 16*t1*t2 + 54*s2*t1*t2 - 
               28*Power(t2,2) - 18*s2*Power(t2,2) - 
               Power(s1,2)*(5 + 10*s2 + 8*t1 + t2) - 
               Power(s,2)*(-11 + 2*s1 + 2*s2 + 9*t2) + 
               s*(-63 + Power(s1,2) + 9*Power(s2,2) + 2*t1 + 
                  103*Power(t1,2) + s1*(59 - 81*s2 + 83*t1 - 83*t2) + 
                  2*t2 - 188*t1*t2 + 93*Power(t2,2) + 
                  s2*(65 - 103*t1 + 101*t2))) - 
            4*(-13 + 4*Power(s,3) - 14*s2 + 40*Power(s2,2) - 
               36*Power(s2,3) + 21*t1 + 38*s2*t1 + 80*Power(s2,2)*t1 - 
               10*Power(t1,2) - 80*s2*Power(t1,2) + 30*t2 - 70*s2*t2 - 
               58*Power(s2,2)*t2 + 85*t1*t2 + 120*s2*t1*t2 - 
               93*Power(t2,2) - 40*s2*Power(t2,2) - 
               2*Power(s,2)*(-16 + 4*s1 + 6*s2 - 4*t1 + 13*t2) + 
               s1*(-41 - 4*Power(s2,2) - 82*t1 - 40*Power(t1,2) - 
                  4*s2*(-31 + 7*t1 - 8*t2) + 123*t2 + 40*t1*t2) - 
               2*Power(s1,2)*(11*s2 + 7*t1 + 2*(5 + t2)) + 
               s*(-95 + 4*Power(s1,2) + 18*Power(s2,2) - t1 + 
                  176*Power(t1,2) + 
                  s1*(45 - 116*s2 + 116*t1 - 120*t2) + 19*t2 - 
                  312*t1*t2 + 150*Power(t2,2) + 
                  s2*(96 - 184*t1 + 168*t2))) + 
            8*(3 + 3*Power(s,3) - 53*s2 + 42*Power(s2,2) - 
               42*Power(s2,3) + 15*t1 + 36*s2*t1 + 96*Power(s2,2)*t1 + 
               2*Power(t1,2) - 96*s2*Power(t1,2) - 
               2*s1*(20 + 24*Power(t1,2) + t1*(37 - 24*t2) + 
                  s2*(-63 + 19*t1 - 22*t2) - 63*t2) + 
               Power(s,2)*(26 - 6*s1 - 4*s2 + 4*t1 - 24*t2) + 13*t2 - 
               80*s2*t2 - 68*Power(s2,2)*t2 + 76*t1*t2 + 144*s2*t1*t2 - 
               96*Power(t2,2) - 48*s2*Power(t2,2) - 
               Power(s1,2)*(25 + 25*s2 + 19*t1 + 3*t2) + 
               s*(-124 + 3*Power(s1,2) + 14*Power(s2,2) + 4*t1 + 
                  230*Power(t1,2) + 16*t2 - 412*t1*t2 + 
                  201*Power(t2,2) - 
                  2*s1*(-52 + 86*s2 - 84*t1 + 87*t2) + 
                  s2*(148 - 234*t1 + 216*t2))) + 
            ((-6 - 6*s + 4*s1 + 8*s2 + 2*t2)*
               (-2*Power(s,3)*Power(t1 - t2,2)*(1 - s2 + t1 - t2) + 
                 4*(-1 + s1 + t1 - t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 4*Power(s,2)*
                  ((2 + s1 + s2)*Power(t1,3) - 
                    Power(t1,2)*
                     (1 + s2 + Power(s2,2) + 7*t2 + 3*s2*t2 + 
                       s1*(-4 + s2 + 2*t2)) + 
                    t1*(-1 + 8*Power(t2,2) + 
                       Power(s2,2)*(1 + 2*t2) + s2*t2*(4 + 3*t2) + 
                       s1*(2 - 2*s2 - 6*t2 + Power(t2,2))) + 
                    t2*(t2 - 3*Power(t2,2) - Power(s2,2)*(2 + t2) + 
                       s1*(-1 + Power(s2,2) + 2*t2 + s2*t2) - 
                       s2*(-2 + 3*t2 + Power(t2,2)))) + 
                 2*s*(-1 - 3*t1 + 9*Power(t1,2) - 5*Power(t1,3) - 
                    5*t2 - 10*t1*t2 + 15*Power(t1,2)*t2 - 
                    3*Power(t2,2) - 19*t1*Power(t2,2) + 
                    9*Power(t2,3) + 
                    Power(s2,3)*
                     (2 + 2*Power(t1,2) + 2*t2 + Power(t2,2) - 
                       2*t1*(2 + t2)) + 
                    Power(s1,2)*
                     (Power(s2,3) - Power(t1,3) + 
                       Power(t1,2)*(-7 + t2) + 6*t1*(-2 + t2) + 
                       2*(-1 + t2) + Power(s2,2)*(1 - t1 + t2) + 
                       s2*(4*t1 + Power(t1,2) + 6*t2)) + 
                    s2*(-1 - 2*Power(t1,3) + 4*t2 + 23*Power(t2,2) + 
                       6*Power(t2,3) - 18*t1*t2*(1 + t2) + 
                       Power(t1,2)*(3 + 14*t2)) + 
                    Power(s2,2)*
                     (-2*Power(t1,3) + 4*Power(t1,2)*(1 + t2) + 
                       t2*(8 + 3*t2 + Power(t2,2)) - 
                       t1*(2 + 12*t2 + 3*Power(t2,2))) + 
                    2*s1*(1 - 2*Power(t1,3) + 
                       Power(s2,3)*(-1 + t1 - t2) + 4*t2 - 
                       5*Power(t2,2) + Power(t1,2)*(-5 + 6*t2) + 
                       t1*(6 + 12*t2 - 4*Power(t2,2)) - 
                       s2*(-2 + Power(t1,3) - Power(t1,2)*(-1 + t2) + 
                        10*t2 - 4*t1*t2 + 6*Power(t2,2)) - 
                       Power(s2,2)*
                        (2 + 2*t2 + Power(t2,2) - 2*t1*(1 + t2))))))/
             (s*Power(-1 + s1 + t1 - t2,2)*(-1 + s2 - t1 + t2)) - 
            ((-1 + t1 - s*t1 + s1*(2 + s2 + t1) + 
                 s2*(-2 + 2*t1 - t2) - 3*t2 + s*t2)*
               (4*s*(-1 + s2 - t1 + t2)*
                  Power(-1 + t1 - s*t1 + s1*(2 + s2 + t1) + 
                    s2*(-2 + 2*t1 - t2) - 3*t2 + s*t2,2) + 
                 12*(-1 + s1 + t1 - t2)*
                  (Power(-1 + s1*s2 + t1 - s1*t1 + t2 - s2*t2,2) + 
                    Power(s,2)*(-1 + s2)*
                     (s2*t1 - Power(t1,2) - t2 + t1*t2) - 
                    s*(-1 + s2)*
                     (Power(s2,2)*(-1 + t1) + 2*t1*(1 - t1 + t2) + 
                       s1*(1 + s2 - 2*t1 + s2*t1 - Power(t1,2) - t2 + 
                        t1*t2) - 
                       s2*(1 + Power(t1,2) + 3*t2 - t1*(2 + t2))))))/
             (s*Power(-1 + s1 + t1 - t2,3)*(-1 + s2 - t1 + t2))))/
        ((-1 + s1)*(-s + s1 - t2)))*B3(1 - s2 + t1 - t2,s,2 + s - s1 - s2))/
   (16.*Power(Pi,2)) + (((8*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)*
          (1 - t1 + s1*(-s2 + t1) - t2 + s2*t2 - s*(2 + t1 + t2))*
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
                t1*(-5 - 2*t2 + Power(t2,2))))) + 
       (4*(64*(s1 - t2)*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) - 
            (2*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
                 s*(2 + t1 + t2))*
               (-2 + 2*Power(s,2) + 6*Power(s1,2) - 7*s2 + 
                 3*Power(s2,2) - 13*t1 + 6*s2*t1 + 3*Power(t1,2) + 
                 2*s1*(-8 + 4*s2 + 5*t1) + 2*t2 + s2*t2 - t1*t2 - 
                 s*(-9 + 8*s1 + 5*s2 + 7*t1 + t2)))/s + 
            16*(10 + 17*s2 - 11*Power(s2,2) + Power(s2,3) - 9*t1 + 
               2*s2*t1 - Power(s2,2)*t1 - 21*Power(t1,2) + 
               8*s2*Power(t1,2) + 10*Power(t1,3) - 34*t2 + 50*s2*t2 + 
               2*Power(s2,2)*t2 + 26*t1*t2 - 78*s2*t1*t2 + 
               41*Power(t1,2)*t2 - 19*Power(t2,2) + 
               79*s2*Power(t2,2) - 112*t1*Power(t2,2) + 
               61*Power(t2,3) + Power(s,2)*(t1 + t2) + 
               Power(s1,2)*(-53 + 75*s2 - 58*t1 + 58*t2) + 
               s1*(40 + 9*Power(s2,2) + 6*t1 - 32*Power(t1,2) + 
                  2*s2*(-42 + 29*t1 - 64*t2) + 46*t2 + 144*t1*t2 - 
                  103*Power(t2,2)) - 
               s*(-4 + 7*s2 + Power(s2,2) - 19*s2*t1 + 
                  26*Power(t1,2) + s1*(5 - 11*s2 + 11*t1 - 27*t2) + 
                  2*t2 + 19*s2*t2 - 44*t1*t2 + 34*Power(t2,2))) - 
            4*(-66 - s2 + 13*Power(s2,2) - Power(s2,3) + 59*t1 + 
               3*s2*t1 + Power(s2,2)*t1 + 36*Power(t1,2) - 
               15*s2*Power(t1,2) - 21*Power(t1,3) + 
               Power(s1,2)*(59 - 110*s2 + 78*t1 - 74*t2) + 
               Power(s,2)*(1 + s2 - 5*t1 - 3*t2) + 58*t2 - 84*s2*t2 - 
               3*Power(s2,2)*t2 - 31*t1*t2 + 124*s2*t1*t2 - 
               47*Power(t1,2)*t2 + 24*Power(t2,2) - 
               128*s2*Power(t2,2) + 158*t1*Power(t2,2) - 
               90*Power(t2,3) - 
               s1*(12 + 16*Power(s2,2) - 28*Power(t1,2) + 
                  2*s2*(-67 + 40*t1 - 95*t2) + 51*t2 - 
                  134*Power(t2,2) + 5*t1*(5 + 36*t2)) + 
               s*(4 - 8*t1 + 54*Power(t1,2) + 26*t2 - 100*t1*t2 + 
                  74*Power(t2,2) - 2*s1*(5 + 17*s2 - 15*t1 + 33*t2) + 
                  s2*(20 - 44*t1 + 42*t2))) + 
            8*(-47 + Power(s1,3) - 32*s2 + 23*Power(s2,2) - 
               2*Power(s2,3) + 47*t1 - 7*s2*t1 + 4*Power(s2,2)*t1 + 
               42*Power(t1,2) - 18*s2*Power(t1,2) - 26*Power(t1,3) + 
               Power(s1,2)*(84 - 151*s2 + 115*t1 - 109*t2) + 
               Power(s,2)*(-2 + s1 + s2 - 3*t1 - 3*t2) + 61*t2 - 
               102*s2*t2 - 5*Power(s2,2)*t2 - 48*t1*t2 + 
               166*s2*t1*t2 - 70*Power(t1,2)*t2 + 37*Power(t2,2) - 
               170*s2*Power(t2,2) + 219*t1*Power(t2,2) - 
               123*Power(t2,3) + 
               s*(-2*Power(s1,2) + Power(s2,2) + 3*t1 + 
                  65*Power(t1,2) + s1*(10 - 38*s2 + 28*t1 - 76*t2) + 
                  16*t2 - 118*t1*t2 + 89*Power(t2,2) + 
                  s2*(23 - 55*t1 + 50*t2)) - 
               s1*(44 + 20*Power(s2,2) - 49*Power(t1,2) + 
                  4*s2*(-43 + 27*t1 - 65*t2) + 75*t2 - 
                  193*Power(t2,2) + t1*(37 + 263*t2))) - 
            32*(2*s2 - Power(s2,2) - 3*Power(t1,2) + s2*Power(t1,2) + 
               Power(t1,3) + Power(s1,2)*(15*s2 - 13*(1 + t1 - t2)) - 
               10*t2 + 12*s2*t2 + 4*t1*t2 - 15*s2*t1*t2 + 
               11*Power(t1,2)*t2 - 3*Power(t2,2) + 15*s2*Power(t2,2) - 
               25*t1*Power(t2,2) + 13*Power(t2,3) + 
               s*(1 - 3*Power(t1,2) + s2*(-1 + 2*t1 - 2*t2) + 
                  5*t1*t2 - 4*Power(t2,2) + s1*(-1 + s2 - t1 + 3*t2)) \
+ s1*(Power(s2,2) - 10*Power(t1,2) + s2*(-16 + 13*t1 - 27*t2) + 
                  35*t1*t2 + 12*(1 + t2 - 2*Power(t2,2)))) - 
            ((-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + s*(2 + t1 + t2))*
               (-4*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
                    s*(2 + t1 + t2),2) - 
                 12*s*(2 - s2 - 2*t1 + s2*t1 - Power(t1,2) + 
                    s2*Power(t1,2) - s2*Power(t1,3) + Power(t1,4) + 
                    Power(s,2)*Power(t1 - t2,2) - 2*t2 + s2*t2 + 
                    4*t1*t2 - s2*t1*t2 - Power(s2,2)*t1*t2 + 
                    s2*Power(t1,2)*t2 + Power(s2,2)*Power(t1,2)*t2 - 
                    2*Power(t1,3)*t2 - s2*Power(t1,3)*t2 - 
                    Power(t2,2) - s2*Power(t2,2) + 
                    Power(s2,2)*Power(t2,2) - s2*t1*Power(t2,2) - 
                    Power(s2,2)*t1*Power(t2,2) + 
                    2*Power(t1,2)*Power(t2,2) + 
                    2*s2*Power(t1,2)*Power(t2,2) + s2*Power(t2,3) - 
                    2*t1*Power(t2,3) - s2*t1*Power(t2,3) + 
                    Power(t2,4) + 
                    Power(s1,2)*(-1 + t2)*
                     (s2*(-1 + t1) + t1*(-t1 + t2)) - 
                    s1*(-1 + t2)*
                     (-1 + 2*s2*(-1 + t1) - Power(s2,2)*(-1 + t1) + 
                       Power(t1,3) + Power(t2,2) - 
                       Power(t1,2)*(1 + 2*t2) + t1*(1 + Power(t2,2))) \
+ s*(-1 + t1 - Power(t1,2) + Power(t1,3) + t2 + t1*t2 - 
                       Power(t1,2)*t2 + Power(t1,3)*t2 - 
                       Power(t2,2) - t1*Power(t2,2) - 
                       2*Power(t1,2)*Power(t2,2) + Power(t2,3) + 
                       t1*Power(t2,3) - 
                       s1*(-1 + t2)*
                        (1 + s2*(-1 + t1) - Power(t1,2) + 
                        t1*(-2 + t2) + t2) - 
                       s2*(-1 + 2*t2 + Power(t2,2) + 
                        Power(t1,2)*(1 + t2) - t1*t2*(3 + t2))))))/
             (Power(s,3)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)) + 
            (2*(-6 - 4*s + 6*s1 + 4*s2 + 4*t1)*
               (-2*Power(s,3)*Power(t1 - t2,2) + 
                 (-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                 2*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
                    s*(2 + t1 + t2),2) - 
                 Power(s,2)*(-6 - 2*t1 + Power(t1,2) + 
                    6*Power(t1,3) + Power(t1,4) - 2*t2 - 8*t1*t2 - 
                    6*Power(t1,2)*t2 + 2*Power(t1,3)*t2 + 
                    Power(t2,2) - 6*t1*Power(t2,2) - 
                    6*Power(t1,2)*Power(t2,2) + 6*Power(t2,3) + 
                    2*t1*Power(t2,3) + Power(t2,4) + 
                    s2*(6 - Power(t1,3) + 4*t2 + 3*Power(t2,2) + 
                       Power(t2,3) - Power(t1,2)*(5 + 3*t2) + 
                       t1*t2*(8 + 3*t2)) - 
                    s1*(-6 - Power(t1,3) + 5*Power(t2,2) + 
                       Power(t2,3) - 3*Power(t1,2)*(1 + t2) + 
                       t1*(-4 - 8*t2 + 3*Power(t2,2)) + 
                       s2*(6 + Power(t1,2) + 2*t2 + Power(t2,2) + 
                        t1*(2 + 4*t2)))) + 
                 2*s*(-4 + 3*t1 + 4*Power(t1,2) - Power(t1,3) - 
                    2*Power(t1,4) + 3*t2 - 6*t1*t2 + Power(t1,2)*t2 + 
                    2*Power(t1,3)*t2 + 4*Power(t2,2) + 
                    t1*Power(t2,2) - Power(t2,3) + 2*t1*Power(t2,3) - 
                    2*Power(t2,4) + 
                    Power(s2,2)*t2*
                     (2 - 2*Power(t1,2) + 2*t2 + t1*t2 + Power(t2,2)) \
+ s2*(3 - 2*t2 - 2*Power(t2,2) + Power(t2,4) + 
                       2*Power(t1,3)*(1 + t2) + 
                       Power(t1,2)*(-1 + 2*t2 - 3*Power(t2,2)) - 
                       2*t1*(2 + t2 + 2*Power(t2,2))) + 
                    Power(s1,2)*
                     (Power(s2,2)*(2 + t1 + t2) + 
                       s2*(-3 - 2*Power(t1,2) + 2*t2 + Power(t2,2) - 
                        2*t1*(2 + t2)) + 
                       t1*(2 + Power(t1,2) - 2*Power(t2,2) + 
                        t1*(2 + t2))) + 
                    s1*(Power(t1,4) + 2*t1*(-1 + t2)*Power(1 + t2,2) + 
                       Power(-1 + t2,2)*(3 + 2*t2) + 
                       Power(s2,2)*
                        (-3 + Power(t1,2) - 2*t1*(-1 + t2) - 4*t2 - 
                        2*Power(t2,2)) - 
                       Power(t1,2)*(2 + 4*t2 + 3*Power(t2,2)) - 
                       2*s2*(Power(t1,3) - Power(t1,2)*(-1 + t2) + 
                        (-1 + t2)*Power(1 + t2,2) - 
                        t1*(1 + 7*t2 + Power(t2,2)))))))/
             (Power(s,2)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2))))/
        ((-1 + s1)*(-s + s1 - t2)) + 
       (4*(64*(-1 + s1 + t1 - t2)*(-1 + t2)*(-1 + s2 - t1 + t2) - 
            (2*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
                 s*(2 + t1 + t2))*
               (-6 + 2*Power(s,2) + 4*Power(s1,2) + Power(s2,2) - 
                 12*t1 + 6*s2*t1 + 2*Power(t1,2) + 
                 s1*(-7 + 4*s2 + 9*t1) - t2 + t1*t2 - Power(t2,2) - 
                 s*(1 + 6*s1 + 3*s2 + 7*t1 + t2)))/s + 
            32*(10 - 11*s2 - t1 + 13*s2*t1 - 12*Power(t1,2) - 
               s2*Power(t1,2) + Power(t1,3) - 11*t2 + s2*t2 + 
               22*t1*t2 - 12*s2*t1*t2 + 11*Power(t1,2)*t2 - 
               12*Power(t2,2) + 14*s2*Power(t2,2) - 
               25*t1*Power(t2,2) + 13*Power(t2,3) + 
               s*(1 + t1 - Power(t1,2) + s2*(-1 + t1 - t2) + t2 + 
                  t1*t2 - 2*Power(t2,2) + s1*(-1 + s2 - t1 + t2)) + 
               s1*(-11 + Power(t1,2) + 25*t2 - 14*Power(t2,2) + 
                  t1*(-11 + 14*t2) - s2*(-12 + t1 + 15*t2))) - 
            16*(31 - 35*s2 - 16*t1 + 57*s2*t1 - Power(s2,2)*t1 - 
               49*Power(t1,2) - 9*s2*Power(t1,2) + 10*Power(t1,3) - 
               42*t2 + 11*s2*t2 + Power(s2,2)*t2 + 82*t1*t2 - 
               50*s2*t1*t2 + 40*Power(t1,2)*t2 - 49*Power(t2,2) + 
               69*s2*Power(t2,2) - 110*t1*Power(t2,2) + 
               60*Power(t2,3) + Power(s1,2)*(-1 + s2 - t1 + t2) + 
               Power(s,2)*(t1 + t2) + 
               s*(7 + 6*t1 - 9*Power(t1,2) + s2*(-8 + 9*t1 - 10*t2) + 
                  6*t2 + 8*t1*t2 - 17*Power(t2,2) + 
                  s1*(-8 + 9*s2 - 9*t1 + 8*t2)) + 
               s1*(-32 + 10*Power(t1,2) + 97*t2 - 66*Power(t2,2) - 
                  2*s2*(-20 + 5*t1 + 37*t2) + t1*(-31 + 66*t2))) + 
            4*(-41 + Power(s1,3) + 7*s2 + 2*Power(s2,2) + 33*t1 - 
               57*s2*t1 + 34*Power(t1,2) + 21*s2*Power(t1,2) - 
               22*Power(t1,3) + 
               Power(s1,2)*(2 - 17*s2 + 18*t1 - 14*t2) + 
               Power(s,2)*(-7 + s1 - 3*t1 - 5*t2) + 59*t2 - 50*s2*t2 - 
               63*t1*t2 + 61*s2*t1*t2 - 41*Power(t1,2)*t2 + 
               48*Power(t2,2) - 103*s2*Power(t2,2) + 
               147*t1*Power(t2,2) - 84*Power(t2,3) + 
               s*(12 - 2*Power(s1,2) + 7*t1 + 17*Power(t1,2) + 
                  s1*(23 - 21*s2 + 13*t1 - 19*t2) + 7*t2 - 18*t1*t2 + 
                  37*Power(t2,2) + s2*(16 - 23*t1 + 23*t2)) - 
               s1*(1 + 6*Power(s2,2) + 11*Power(t1,2) + 82*t2 - 
                  100*Power(t2,2) + t1*(32 + 111*t2) - 
                  s2*(5 + 16*t1 + 117*t2))) + 
            8*(50 - 40*s2 - Power(s2,2) - 49*t1 + 95*s2*t1 - 
               Power(s2,2)*t1 - 78*Power(t1,2) - 23*s2*Power(t1,2) + 
               27*Power(t1,3) - 79*t2 + 45*s2*t2 + 2*Power(s2,2)*t2 + 
               120*t1*t2 - 90*s2*t1*t2 + 65*Power(t1,2)*t2 - 
               79*Power(t2,2) + 141*s2*Power(t2,2) - 
               210*t1*Power(t2,2) + 118*Power(t2,3) + 
               Power(s,2)*(4 + 5*t1 + 5*t2) + 
               Power(s1,2)*(-7 + 11*s2 - 8*t1 + 9*t2) + 
               s*(5 + 6*t1 - 24*Power(t1,2) + 
                  s2*(-21 + 23*t1 - 28*t2) + 5*t2 + 20*t1*t2 - 
                  44*Power(t2,2) + s1*(-22 + 24*s2 - 25*t1 + 21*t2)) + 
               s1*(-26 + 3*Power(s2,2) + 23*Power(t1,2) + 
                  s2*(34 - 20*t1 - 155*t2) + 150*t2 - 
                  135*Power(t2,2) + t1*(-17 + 141*t2))) - 
            ((-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + s*(2 + t1 + t2))*
               (-4*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
                    s*(2 + t1 + t2),2) - 
                 12*s*(2 - s2 - 2*t1 + s2*t1 - Power(t1,2) + 
                    s2*Power(t1,2) - s2*Power(t1,3) + Power(t1,4) + 
                    Power(s,2)*Power(t1 - t2,2) - 2*t2 + s2*t2 + 
                    4*t1*t2 - s2*t1*t2 - Power(s2,2)*t1*t2 + 
                    s2*Power(t1,2)*t2 + Power(s2,2)*Power(t1,2)*t2 - 
                    2*Power(t1,3)*t2 - s2*Power(t1,3)*t2 - 
                    Power(t2,2) - s2*Power(t2,2) + 
                    Power(s2,2)*Power(t2,2) - s2*t1*Power(t2,2) - 
                    Power(s2,2)*t1*Power(t2,2) + 
                    2*Power(t1,2)*Power(t2,2) + 
                    2*s2*Power(t1,2)*Power(t2,2) + s2*Power(t2,3) - 
                    2*t1*Power(t2,3) - s2*t1*Power(t2,3) + 
                    Power(t2,4) + 
                    Power(s1,2)*(-1 + t2)*
                     (s2*(-1 + t1) + t1*(-t1 + t2)) - 
                    s1*(-1 + t2)*
                     (-1 + 2*s2*(-1 + t1) - Power(s2,2)*(-1 + t1) + 
                       Power(t1,3) + Power(t2,2) - 
                       Power(t1,2)*(1 + 2*t2) + t1*(1 + Power(t2,2))) \
+ s*(-1 + t1 - Power(t1,2) + Power(t1,3) + t2 + t1*t2 - 
                       Power(t1,2)*t2 + Power(t1,3)*t2 - 
                       Power(t2,2) - t1*Power(t2,2) - 
                       2*Power(t1,2)*Power(t2,2) + Power(t2,3) + 
                       t1*Power(t2,3) - 
                       s1*(-1 + t2)*
                        (1 + s2*(-1 + t1) - Power(t1,2) + 
                        t1*(-2 + t2) + t2) - 
                       s2*(-1 + 2*t2 + Power(t2,2) + 
                        Power(t1,2)*(1 + t2) - t1*t2*(3 + t2))))))/
             (Power(s,3)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)) + 
            ((-8 - 8*s + 10*s1 + 6*s2 + 8*t1)*
               (-2*Power(s,3)*Power(t1 - t2,2) + 
                 (-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                 2*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
                    s*(2 + t1 + t2),2) - 
                 Power(s,2)*(-6 - 2*t1 + Power(t1,2) + 
                    6*Power(t1,3) + Power(t1,4) - 2*t2 - 8*t1*t2 - 
                    6*Power(t1,2)*t2 + 2*Power(t1,3)*t2 + 
                    Power(t2,2) - 6*t1*Power(t2,2) - 
                    6*Power(t1,2)*Power(t2,2) + 6*Power(t2,3) + 
                    2*t1*Power(t2,3) + Power(t2,4) + 
                    s2*(6 - Power(t1,3) + 4*t2 + 3*Power(t2,2) + 
                       Power(t2,3) - Power(t1,2)*(5 + 3*t2) + 
                       t1*t2*(8 + 3*t2)) - 
                    s1*(-6 - Power(t1,3) + 5*Power(t2,2) + 
                       Power(t2,3) - 3*Power(t1,2)*(1 + t2) + 
                       t1*(-4 - 8*t2 + 3*Power(t2,2)) + 
                       s2*(6 + Power(t1,2) + 2*t2 + Power(t2,2) + 
                        t1*(2 + 4*t2)))) + 
                 2*s*(-4 + 3*t1 + 4*Power(t1,2) - Power(t1,3) - 
                    2*Power(t1,4) + 3*t2 - 6*t1*t2 + Power(t1,2)*t2 + 
                    2*Power(t1,3)*t2 + 4*Power(t2,2) + 
                    t1*Power(t2,2) - Power(t2,3) + 2*t1*Power(t2,3) - 
                    2*Power(t2,4) + 
                    Power(s2,2)*t2*
                     (2 - 2*Power(t1,2) + 2*t2 + t1*t2 + Power(t2,2)) \
+ s2*(3 - 2*t2 - 2*Power(t2,2) + Power(t2,4) + 
                       2*Power(t1,3)*(1 + t2) + 
                       Power(t1,2)*(-1 + 2*t2 - 3*Power(t2,2)) - 
                       2*t1*(2 + t2 + 2*Power(t2,2))) + 
                    Power(s1,2)*
                     (Power(s2,2)*(2 + t1 + t2) + 
                       s2*(-3 - 2*Power(t1,2) + 2*t2 + Power(t2,2) - 
                        2*t1*(2 + t2)) + 
                       t1*(2 + Power(t1,2) - 2*Power(t2,2) + 
                        t1*(2 + t2))) + 
                    s1*(Power(t1,4) + 2*t1*(-1 + t2)*Power(1 + t2,2) + 
                       Power(-1 + t2,2)*(3 + 2*t2) + 
                       Power(s2,2)*
                        (-3 + Power(t1,2) - 2*t1*(-1 + t2) - 4*t2 - 
                        2*Power(t2,2)) - 
                       Power(t1,2)*(2 + 4*t2 + 3*Power(t2,2)) - 
                       2*s2*(Power(t1,3) - Power(t1,2)*(-1 + t2) + 
                        (-1 + t2)*Power(1 + t2,2) - 
                        t1*(1 + 7*t2 + Power(t2,2)))))))/
             (Power(s,2)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2))))/
        ((-1 + s1)*(-1 + t2)) + 
       (4*(64*(-1 + t1)*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) - 
            (2*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
                 s*(2 + t1 + t2))*
               (-8 + Power(s,2) + 2*Power(s1,2) - 9*s2 + 
                 4*Power(s2,2) + 2*s2*t1 - 13*t2 + 6*s2*t2 + t1*t2 + 
                 Power(t2,2) - s*(-2 + 3*s1 + 5*s2 + t1 + 5*t2) + 
                 s1*(7*s2 + t1 + 6*t2)))/s + 
            32*(13 - 13*s2 - 12*t1 + 26*s2*t1 - 14*Power(t1,2) - 
               13*s2*Power(t1,2) + 13*Power(t1,3) - 3*t2 - 10*s2*t2 - 
               Power(s2,2)*t2 + 25*t1*t2 + 14*s2*t1*t2 - 
               25*Power(t1,2)*t2 - 10*Power(t2,2) - 2*s2*Power(t2,2) + 
               12*t1*Power(t2,2) + 
               s1*(Power(s2,2) + t1 + 15*Power(t1,2) + 13*(-1 + t2) - 
                  15*t1*t2 + 2*s2*(6 - 8*t1 + t2)) + 
               s*(2 - 3*Power(t1,2) + s2*(-2 + 2*t1 - t2) - t2 + 
                  3*t1*t2 + s1*(-1 + s2 - t1 + t2))) + 
            16*(-51 + 55*s2 + 2*Power(s2,2) + 46*t1 - 111*s2*t1 - 
               3*Power(s2,2)*t1 + 63*Power(t1,2) + 60*s2*Power(t1,2) - 
               60*Power(t1,3) + 24*t2 + 31*s2*t2 + 10*Power(s2,2)*t2 - 
               105*t1*t2 - 69*s2*t1*t2 + 111*Power(t1,2)*t2 + 
               38*Power(t2,2) + 16*s2*Power(t2,2) - 
               50*t1*Power(t2,2) - Power(t2,3) - 
               Power(s1,2)*(-1 + 2*s2 + t2) + 
               Power(s,2)*(-1 + s2 - t1 + t2) - 
               s*(16 + Power(s2,2) - 3*t1 - 24*Power(t1,2) + 
                  2*s2*(-8 + 7*t1 - 2*t2) - 8*t2 + 22*t1*t2 + 
                  2*Power(t2,2) + s1*(-8 + 8*s2 - 9*t1 + 9*t2)) - 
               s1*(-51 + 11*Power(s2,2) + 76*Power(t1,2) + 
                  t1*(7 - 73*t2) + 60*t2 - 2*Power(t2,2) + 
                  s2*(43 - 82*t1 + 15*t2))) + 
            4*(-38 + 51*s2 + 6*Power(s2,2) - Power(s2,3) + 
               Power(s,2)*(-3 + s2 - 2*t1) + 47*t1 - 127*s2*t1 - 
               9*Power(s2,2)*t1 + 61*Power(t1,2) + 80*s2*Power(t1,2) - 
               82*Power(t1,3) + 45*t2 - 7*s2*t2 + 22*Power(s2,2)*t2 - 
               102*t1*t2 - 87*s2*t1*t2 + 138*Power(t1,2)*t2 + 
               33*Power(t2,2) + 20*s2*Power(t2,2) - 
               46*t1*Power(t2,2) - 10*Power(t2,3) - 
               Power(s1,2)*(3 + 10*s2 - 2*t1 + 3*t2) + 
               s*(-20 + 21*t1 + 48*Power(t1,2) + 
                  s1*(17 - 17*s2 + 18*t1 - 23*t2) + 11*t2 - 
                  51*t1*t2 + 2*Power(t2,2) + s2*(25 - 33*t1 + 10*t2)) \
- s1*(-27 + 29*Power(s2,2) + 40*t1 + 112*Power(t1,2) + 78*t2 - 
                  97*t1*t2 - 12*Power(t2,2) + s2*(25 - 121*t1 + 17*t2))\
) - 8*(-80 + 103*s2 + 8*Power(s2,2) - Power(s2,3) + 70*t1 - 
               198*s2*t1 - 11*Power(s2,2)*t1 + 109*Power(t1,2) + 
               116*s2*Power(t1,2) - 117*Power(t1,3) + 
               Power(s1,2)*(-10*s2 + t1 - 5*t2) + 69*t2 + 34*s2*t2 + 
               25*Power(s2,2)*t2 - 179*t1*t2 - 134*s2*t1*t2 + 
               207*Power(t1,2)*t2 + 66*Power(t2,2) + 
               33*s2*Power(t2,2) - 83*t1*Power(t2,2) - 7*Power(t2,3) + 
               Power(s,2)*(-3 + 2*s2 - 3*t1 + t2) - 
               s*(40 + Power(s2,2) - 18*t1 - 60*Power(t1,2) + 
                  s2*(-35 + 36*t1 - 12*t2) - 13*t2 + 57*t1*t2 + 
                  2*Power(t2,2) + s1*(-21 + 19*s2 - 23*t1 + 23*t2)) - 
               s1*(-75 + 34*Power(s2,2) + 156*Power(t1,2) + 
                  t1*(29 - 142*t2) + 114*t2 - 9*Power(t2,2) + 
                  s2*(59 - 168*t1 + 34*t2))) - 
            ((-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + s*(2 + t1 + t2))*
               (-4*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
                    s*(2 + t1 + t2),2) - 
                 12*s*(2 - s2 - 2*t1 + s2*t1 - Power(t1,2) + 
                    s2*Power(t1,2) - s2*Power(t1,3) + Power(t1,4) + 
                    Power(s,2)*Power(t1 - t2,2) - 2*t2 + s2*t2 + 
                    4*t1*t2 - s2*t1*t2 - Power(s2,2)*t1*t2 + 
                    s2*Power(t1,2)*t2 + Power(s2,2)*Power(t1,2)*t2 - 
                    2*Power(t1,3)*t2 - s2*Power(t1,3)*t2 - 
                    Power(t2,2) - s2*Power(t2,2) + 
                    Power(s2,2)*Power(t2,2) - s2*t1*Power(t2,2) - 
                    Power(s2,2)*t1*Power(t2,2) + 
                    2*Power(t1,2)*Power(t2,2) + 
                    2*s2*Power(t1,2)*Power(t2,2) + s2*Power(t2,3) - 
                    2*t1*Power(t2,3) - s2*t1*Power(t2,3) + 
                    Power(t2,4) + 
                    Power(s1,2)*(-1 + t2)*
                     (s2*(-1 + t1) + t1*(-t1 + t2)) - 
                    s1*(-1 + t2)*
                     (-1 + 2*s2*(-1 + t1) - Power(s2,2)*(-1 + t1) + 
                       Power(t1,3) + Power(t2,2) - 
                       Power(t1,2)*(1 + 2*t2) + t1*(1 + Power(t2,2))) \
+ s*(-1 + t1 - Power(t1,2) + Power(t1,3) + t2 + t1*t2 - 
                       Power(t1,2)*t2 + Power(t1,3)*t2 - 
                       Power(t2,2) - t1*Power(t2,2) - 
                       2*Power(t1,2)*Power(t2,2) + Power(t2,3) + 
                       t1*Power(t2,3) - 
                       s1*(-1 + t2)*
                        (1 + s2*(-1 + t1) - Power(t1,2) + 
                        t1*(-2 + t2) + t2) - 
                       s2*(-1 + 2*t2 + Power(t2,2) + 
                        Power(t1,2)*(1 + t2) - t1*t2*(3 + t2))))))/
             (Power(s,3)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)) - 
            (2*(5 + 3*s - 4*s1 - 5*s2 - t1 - 3*t2)*
               (-2*Power(s,3)*Power(t1 - t2,2) + 
                 (-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                 2*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
                    s*(2 + t1 + t2),2) - 
                 Power(s,2)*(-6 - 2*t1 + Power(t1,2) + 
                    6*Power(t1,3) + Power(t1,4) - 2*t2 - 8*t1*t2 - 
                    6*Power(t1,2)*t2 + 2*Power(t1,3)*t2 + 
                    Power(t2,2) - 6*t1*Power(t2,2) - 
                    6*Power(t1,2)*Power(t2,2) + 6*Power(t2,3) + 
                    2*t1*Power(t2,3) + Power(t2,4) + 
                    s2*(6 - Power(t1,3) + 4*t2 + 3*Power(t2,2) + 
                       Power(t2,3) - Power(t1,2)*(5 + 3*t2) + 
                       t1*t2*(8 + 3*t2)) - 
                    s1*(-6 - Power(t1,3) + 5*Power(t2,2) + 
                       Power(t2,3) - 3*Power(t1,2)*(1 + t2) + 
                       t1*(-4 - 8*t2 + 3*Power(t2,2)) + 
                       s2*(6 + Power(t1,2) + 2*t2 + Power(t2,2) + 
                        t1*(2 + 4*t2)))) + 
                 2*s*(-4 + 3*t1 + 4*Power(t1,2) - Power(t1,3) - 
                    2*Power(t1,4) + 3*t2 - 6*t1*t2 + Power(t1,2)*t2 + 
                    2*Power(t1,3)*t2 + 4*Power(t2,2) + 
                    t1*Power(t2,2) - Power(t2,3) + 2*t1*Power(t2,3) - 
                    2*Power(t2,4) + 
                    Power(s2,2)*t2*
                     (2 - 2*Power(t1,2) + 2*t2 + t1*t2 + Power(t2,2)) \
+ s2*(3 - 2*t2 - 2*Power(t2,2) + Power(t2,4) + 
                       2*Power(t1,3)*(1 + t2) + 
                       Power(t1,2)*(-1 + 2*t2 - 3*Power(t2,2)) - 
                       2*t1*(2 + t2 + 2*Power(t2,2))) + 
                    Power(s1,2)*
                     (Power(s2,2)*(2 + t1 + t2) + 
                       s2*(-3 - 2*Power(t1,2) + 2*t2 + Power(t2,2) - 
                        2*t1*(2 + t2)) + 
                       t1*(2 + Power(t1,2) - 2*Power(t2,2) + 
                        t1*(2 + t2))) + 
                    s1*(Power(t1,4) + 2*t1*(-1 + t2)*Power(1 + t2,2) + 
                       Power(-1 + t2,2)*(3 + 2*t2) + 
                       Power(s2,2)*
                        (-3 + Power(t1,2) - 2*t1*(-1 + t2) - 4*t2 - 
                        2*Power(t2,2)) - 
                       Power(t1,2)*(2 + 4*t2 + 3*Power(t2,2)) - 
                       2*s2*(Power(t1,3) - Power(t1,2)*(-1 + t2) + 
                        (-1 + t2)*Power(1 + t2,2) - 
                        t1*(1 + 7*t2 + Power(t2,2)))))))/
             (Power(s,2)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2))))/
        ((-1 + s2)*(-1 + t1)) + 
       (4*(64*s*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2) - 
            (2*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
                 s*(2 + t1 + t2))*
               (2*Power(s,2) + 6*Power(s1,2) - 15*s2 + 
                 5*Power(s2,2) - 11*t1 + 8*s2*t1 - Power(t1,2) - 
                 2*t2 + 3*s2*t2 + 3*t1*t2 + 
                 2*s1*(-10 + 7*s2 + 2*t1 + 3*t2) - 
                 s*(-13 + 8*s1 + 7*s2 + 5*t1 + 3*t2)))/s + 
            32*(2 + (-1 + s1)*Power(s2,2) - 5*t1 + 5*s1*t1 - 
               Power(s1,2)*t1 + Power(t1,2) - s1*Power(t1,2) - s1*t2 + 
               3*t1*t2 - s1*t1*t2 - 3*Power(t2,2) + s1*Power(t2,2) - 
               t1*Power(t2,2) + Power(t2,3) + Power(s,2)*(t1 + t2) + 
               s2*(1 - 4*s1 + Power(s1,2) - 2*t1*(-1 + t2) + 
                  2*Power(t2,2)) + 
               s*(12 - 3*t1 - 13*Power(t1,2) + 
                  s2*(-13 + 15*t1 - 15*t2) + 30*t1*t2 - 
                  15*Power(t2,2) + s1*(-14 + 15*s2 - 14*t1 + 15*t2))) - 
            16*(21 + 8*s2 - 11*Power(s2,2) + Power(s2,3) - 40*t1 + 
               18*s2*t1 + Power(s2,2)*t1 + 8*Power(t1,2) - 
               3*s2*Power(t1,2) + Power(t1,3) - 10*t2 + 6*s2*t2 + 
               Power(s2,2)*t2 + 19*t1*t2 - 18*s2*t1*t2 - 
               19*Power(t2,2) + 20*s2*Power(t2,2) - 
               11*t1*Power(t2,2) + 10*Power(t2,3) + 
               Power(s1,2)*(-3 + 14*s2 - 13*t1 + 4*t2) + 
               Power(s,2)*(3 + s2 + 7*t1 + 9*t2) + 
               s1*(-2 + 12*Power(s2,2) - 11*Power(t1,2) + 
                  t1*(47 - 2*t2) - 6*t2 + 4*Power(t2,2) - 
                  2*s2*(21 + 2*t2)) + 
               s*(37 - 2*Power(s2,2) - 23*t1 - 56*Power(t1,2) + 
                  s2*(-51 + 75*t1 - 77*t2) + t2 + 149*t1*t2 - 
                  77*Power(t2,2) + s1*(-58 + 71*s2 - 62*t1 + 73*t2))) - 
            4*(39 + 2*s2 - 9*Power(s2,2) + 3*Power(s2,3) - 66*t1 + 
               42*s2*t1 + Power(s2,2)*t1 + 6*Power(t1,2) - 
               9*s2*Power(t1,2) + 9*Power(t1,3) - 33*t2 + 50*s2*t2 + 
               3*Power(s2,2)*t2 + 31*t1*t2 - 48*s2*t1*t2 - 
               3*Power(t1,2)*t2 - 16*Power(t2,2) + 56*s2*Power(t2,2) - 
               38*t1*Power(t2,2) + 32*Power(t2,3) + 
               8*Power(s1,2)*(5*s2 - 4*t1 + 2*t2) + 
               Power(s,2)*(21 + 3*s2 + 13*t1 + 19*t2) + 
               s1*(-29 + 30*Power(s2,2) - 14*Power(t1,2) + 6*t2 - 
                  8*Power(t2,2) + t1*(103 + 4*t2) - 
                  2*s2*(37 + 2*t1 + 13*t2)) - 
               s*(-23 + 6*Power(s2,2) + 34*t1 + 80*Power(t1,2) + 
                  s1*(71 - 100*s2 + 84*t1 - 108*t2) + 3*t2 - 
                  226*t1*t2 + 118*Power(t2,2) + 
                  s2*(65 - 116*t1 + 116*t2))) + 
            8*(61 + Power(s1,3) + 13*s2 - 25*Power(s2,2) + 
               4*Power(s2,3) - 97*t1 + 38*s2*t1 + 6*Power(s2,2)*t1 + 
               14*Power(t1,2) - 10*s2*Power(t1,2) + 6*Power(t1,3) - 
               39*t2 + 38*s2*t2 + 4*Power(s2,2)*t2 + 39*t1*t2 - 
               50*s2*t1*t2 - Power(t1,2)*t2 - 34*Power(t2,2) + 
               60*s2*Power(t2,2) - 37*t1*Power(t2,2) + 
               32*Power(t2,3) + 
               Power(s1,2)*(-15 + 49*s2 - 37*t1 + 19*t2) + 
               Power(s,2)*(14 + s1 + 4*s2 + 18*t1 + 24*t2) + 
               s1*(-17 + 39*Power(s2,2) - 26*Power(t1,2) + 
                  2*s2*(-60 + 2*t1 - 9*t2) - 7*t2 + t1*(117 + 5*t2)) - 
               s*(-43 + 2*Power(s1,2) + 8*Power(s2,2) + 41*t1 + 
                  108*Power(t1,2) + 
                  s1*(91 - 135*s2 + 121*t1 - 145*t2) + t2 - 
                  305*t1*t2 + 161*Power(t2,2) + 
                  s2*(85 - 151*t1 + 161*t2))) - 
            ((-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + s*(2 + t1 + t2))*
               (-4*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
                    s*(2 + t1 + t2),2) - 
                 12*s*(2 - s2 - 2*t1 + s2*t1 - Power(t1,2) + 
                    s2*Power(t1,2) - s2*Power(t1,3) + Power(t1,4) + 
                    Power(s,2)*Power(t1 - t2,2) - 2*t2 + s2*t2 + 
                    4*t1*t2 - s2*t1*t2 - Power(s2,2)*t1*t2 + 
                    s2*Power(t1,2)*t2 + Power(s2,2)*Power(t1,2)*t2 - 
                    2*Power(t1,3)*t2 - s2*Power(t1,3)*t2 - 
                    Power(t2,2) - s2*Power(t2,2) + 
                    Power(s2,2)*Power(t2,2) - s2*t1*Power(t2,2) - 
                    Power(s2,2)*t1*Power(t2,2) + 
                    2*Power(t1,2)*Power(t2,2) + 
                    2*s2*Power(t1,2)*Power(t2,2) + s2*Power(t2,3) - 
                    2*t1*Power(t2,3) - s2*t1*Power(t2,3) + 
                    Power(t2,4) + 
                    Power(s1,2)*(-1 + t2)*
                     (s2*(-1 + t1) + t1*(-t1 + t2)) - 
                    s1*(-1 + t2)*
                     (-1 + 2*s2*(-1 + t1) - Power(s2,2)*(-1 + t1) + 
                       Power(t1,3) + Power(t2,2) - 
                       Power(t1,2)*(1 + 2*t2) + t1*(1 + Power(t2,2))) \
+ s*(-1 + t1 - Power(t1,2) + Power(t1,3) + t2 + t1*t2 - 
                       Power(t1,2)*t2 + Power(t1,3)*t2 - 
                       Power(t2,2) - t1*Power(t2,2) - 
                       2*Power(t1,2)*Power(t2,2) + Power(t2,3) + 
                       t1*Power(t2,3) - 
                       s1*(-1 + t2)*
                        (1 + s2*(-1 + t1) - Power(t1,2) + 
                        t1*(-2 + t2) + t2) - 
                       s2*(-1 + 2*t2 + Power(t2,2) + 
                        Power(t1,2)*(1 + t2) - t1*t2*(3 + t2))))))/
             (Power(s,3)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)) - 
            (4*(4 + 2*s - 3*s1 - 3*s2 - t1 - t2)*
               (-2*Power(s,3)*Power(t1 - t2,2) + 
                 (-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                 2*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
                    s*(2 + t1 + t2),2) - 
                 Power(s,2)*(-6 - 2*t1 + Power(t1,2) + 
                    6*Power(t1,3) + Power(t1,4) - 2*t2 - 8*t1*t2 - 
                    6*Power(t1,2)*t2 + 2*Power(t1,3)*t2 + 
                    Power(t2,2) - 6*t1*Power(t2,2) - 
                    6*Power(t1,2)*Power(t2,2) + 6*Power(t2,3) + 
                    2*t1*Power(t2,3) + Power(t2,4) + 
                    s2*(6 - Power(t1,3) + 4*t2 + 3*Power(t2,2) + 
                       Power(t2,3) - Power(t1,2)*(5 + 3*t2) + 
                       t1*t2*(8 + 3*t2)) - 
                    s1*(-6 - Power(t1,3) + 5*Power(t2,2) + 
                       Power(t2,3) - 3*Power(t1,2)*(1 + t2) + 
                       t1*(-4 - 8*t2 + 3*Power(t2,2)) + 
                       s2*(6 + Power(t1,2) + 2*t2 + Power(t2,2) + 
                        t1*(2 + 4*t2)))) + 
                 2*s*(-4 + 3*t1 + 4*Power(t1,2) - Power(t1,3) - 
                    2*Power(t1,4) + 3*t2 - 6*t1*t2 + Power(t1,2)*t2 + 
                    2*Power(t1,3)*t2 + 4*Power(t2,2) + 
                    t1*Power(t2,2) - Power(t2,3) + 2*t1*Power(t2,3) - 
                    2*Power(t2,4) + 
                    Power(s2,2)*t2*
                     (2 - 2*Power(t1,2) + 2*t2 + t1*t2 + Power(t2,2)) \
+ s2*(3 - 2*t2 - 2*Power(t2,2) + Power(t2,4) + 
                       2*Power(t1,3)*(1 + t2) + 
                       Power(t1,2)*(-1 + 2*t2 - 3*Power(t2,2)) - 
                       2*t1*(2 + t2 + 2*Power(t2,2))) + 
                    Power(s1,2)*
                     (Power(s2,2)*(2 + t1 + t2) + 
                       s2*(-3 - 2*Power(t1,2) + 2*t2 + Power(t2,2) - 
                        2*t1*(2 + t2)) + 
                       t1*(2 + Power(t1,2) - 2*Power(t2,2) + 
                        t1*(2 + t2))) + 
                    s1*(Power(t1,4) + 2*t1*(-1 + t2)*Power(1 + t2,2) + 
                       Power(-1 + t2,2)*(3 + 2*t2) + 
                       Power(s2,2)*
                        (-3 + Power(t1,2) - 2*t1*(-1 + t2) - 4*t2 - 
                        2*Power(t2,2)) - 
                       Power(t1,2)*(2 + 4*t2 + 3*Power(t2,2)) - 
                       2*s2*(Power(t1,3) - Power(t1,2)*(-1 + t2) + 
                        (-1 + t2)*Power(1 + t2,2) - 
                        t1*(1 + 7*t2 + Power(t2,2)))))))/
             (Power(s,2)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2))))/
        ((s - s2 + t1)*(s - s1 + t2)) + 
       (4*(64*(s2 - t1)*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) - 
            (2*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
                 s*(2 + t1 + t2))*
               (-2 + 2*Power(s,2) + 4*Power(s1,2) - 16*s2 + 
                 5*Power(s2,2) - 3*t1 + 5*s2*t1 - Power(t1,2) - 
                 11*t2 + 6*s2*t2 + 2*t1*t2 + Power(t2,2) - 
                 s*(-7 + 6*s1 + 7*s2 + 3*t1 + 5*t2) + 
                 s1*(-10 + 11*s2 + 2*t1 + 7*t2)))/s - 
            32*(2 - 13*t1 + 15*s1*t1 + 14*s1*Power(t1,2) + 
               13*Power(t1,3) - 3*t2 - 14*s1*t1*t2 - 
               25*Power(t1,2)*t2 + Power(t2,2) + 12*t1*Power(t2,2) - 
               s*(-2 + s2*(2 - 3*t1) + t1 + 2*Power(t1,2) + 2*t2 - 
                  3*t1*t2 + Power(t2,2) - s1*(-1 + s2 - t1 + t2)) + 
               Power(s2,2)*(15*s1 + 13*t1 - 2*(6 + 7*t2)) + 
               s2*(10 - 26*Power(t1,2) + 4*t2 - 13*Power(t2,2) + 
                  s1*(-15 - 29*t1 + 14*t2) + t1*(13 + 38*t2))) + 
            4*(64 + Power(s1,3) - 11*s2 - 44*Power(s2,2) - 
               3*Power(s2,3) - 72*t1 + 68*s2*t1 + 71*Power(s2,2)*t1 - 
               10*Power(t1,2) - 149*s2*Power(t1,2) + 84*Power(t1,3) + 
               Power(s,2)*(s1 - 3*s2 + 6*t1 - 2*t2) - 61*t2 + 
               90*s2*t2 - 93*Power(s2,2)*t2 + 4*t1*t2 + 189*s2*t1*t2 - 
               149*Power(t1,2)*t2 + 9*Power(t2,2) - 
               53*s2*Power(t2,2) + 51*t1*Power(t2,2) + 
               14*Power(t2,3) + Power(s1,2)*(8*s2 + t1 - 3*(1 + t2)) + 
               s*(13 - 2*Power(s1,2) + 6*Power(s2,2) - 28*t1 - 
                  45*Power(t1,2) - 26*t2 + 76*t1*t2 - 
                  27*Power(t2,2) + s2*(-9 + 53*t1 + t2) + 
                  s1*(-15 + 33*s2 - 35*t1 + 43*t2)) + 
               s1*(-22 + 96*Power(s2,2) + 101*Power(t1,2) + 
                  t1*(113 - 97*t2) + t2 - 12*Power(t2,2) + 
                  s2*(-112 - 186*t1 + 51*t2))) - 
            16*(-19 - 24*s2 + 45*Power(s2,2) + Power(s2,3) + 54*t1 - 
               58*s2*t1 - 59*Power(s2,2)*t1 + Power(t1,2) + 
               119*s2*Power(t1,2) - 60*Power(t1,3) + 
               Power(s1,2)*(-1 + t2) + 26*t2 - 39*s2*t2 + 
               69*Power(s2,2)*t2 - t1*t2 - 169*s2*t1*t2 + 
               113*Power(t1,2)*t2 - 7*Power(t2,2) + 
               57*s2*Power(t2,2) - 51*t1*Power(t2,2) - 2*Power(t2,3) + 
               Power(s,2)*(-1 + s2 - t1 + t2) + 
               s*(-15 - 2*Power(s2,2) + 11*t1 + 16*Power(t1,2) + 
                  s1*(11 - 11*s2 + 10*t1 - 12*t2) + 
                  s2*(15 - 23*t1 - 3*t2) + 16*t2 - 25*t1*t2 + 
                  8*Power(t2,2)) + 
               s1*(3 - 72*Power(s2,2) - 67*Power(t1,2) + 
                  s2*(69 + 140*t1 - 60*t2) - 3*t2 + 2*Power(t2,2) + 
                  t1*(-75 + 67*t2))) + 
            8*(-57 - 11*s2 + 69*Power(s2,2) + 4*Power(s2,3) + 98*t1 - 
               115*s2*t1 - 107*Power(s2,2)*t1 + 11*Power(t1,2) + 
               224*s2*Power(t1,2) - 118*Power(t1,3) - 
               Power(s1,2)*(4 + s2 - 6*t2) + 66*t2 - 110*s2*t2 + 
               138*Power(s2,2)*t2 - 10*t1*t2 - 302*s2*t1*t2 + 
               216*Power(t1,2)*t2 - 16*Power(t2,2) + 
               97*s2*Power(t2,2) - 87*t1*Power(t2,2) - 
               11*Power(t2,3) + Power(s,2)*(-3 + 4*s2 - 4*t1 + 4*t2) + 
               s*(-29 - 8*Power(s2,2) + 37*t1 + 47*Power(t1,2) + 
                  s1*(32 - 38*s2 + 32*t1 - 45*t2) + 
                  s2*(30 - 64*t1 - 9*t2) + 38*t2 - 78*t1*t2 + 
                  25*Power(t2,2)) + 
               s1*(16 - 137*Power(s2,2) - 137*Power(t1,2) + 
                  s2*(131 + 276*t1 - 93*t2) - 12*t2 + 11*Power(t2,2) + 
                  t1*(-155 + 136*t2))) - 
            ((-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + s*(2 + t1 + t2))*
               (-4*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
                    s*(2 + t1 + t2),2) - 
                 12*s*(2 - s2 - 2*t1 + s2*t1 - Power(t1,2) + 
                    s2*Power(t1,2) - s2*Power(t1,3) + Power(t1,4) + 
                    Power(s,2)*Power(t1 - t2,2) - 2*t2 + s2*t2 + 
                    4*t1*t2 - s2*t1*t2 - Power(s2,2)*t1*t2 + 
                    s2*Power(t1,2)*t2 + Power(s2,2)*Power(t1,2)*t2 - 
                    2*Power(t1,3)*t2 - s2*Power(t1,3)*t2 - 
                    Power(t2,2) - s2*Power(t2,2) + 
                    Power(s2,2)*Power(t2,2) - s2*t1*Power(t2,2) - 
                    Power(s2,2)*t1*Power(t2,2) + 
                    2*Power(t1,2)*Power(t2,2) + 
                    2*s2*Power(t1,2)*Power(t2,2) + s2*Power(t2,3) - 
                    2*t1*Power(t2,3) - s2*t1*Power(t2,3) + 
                    Power(t2,4) + 
                    Power(s1,2)*(-1 + t2)*
                     (s2*(-1 + t1) + t1*(-t1 + t2)) - 
                    s1*(-1 + t2)*
                     (-1 + 2*s2*(-1 + t1) - Power(s2,2)*(-1 + t1) + 
                       Power(t1,3) + Power(t2,2) - 
                       Power(t1,2)*(1 + 2*t2) + t1*(1 + Power(t2,2))) \
+ s*(-1 + t1 - Power(t1,2) + Power(t1,3) + t2 + t1*t2 - 
                       Power(t1,2)*t2 + Power(t1,3)*t2 - 
                       Power(t2,2) - t1*Power(t2,2) - 
                       2*Power(t1,2)*Power(t2,2) + Power(t2,3) + 
                       t1*Power(t2,3) - 
                       s1*(-1 + t2)*
                        (1 + s2*(-1 + t1) - Power(t1,2) + 
                        t1*(-2 + t2) + t2) - 
                       s2*(-1 + 2*t2 + Power(t2,2) + 
                        Power(t1,2)*(1 + t2) - t1*t2*(3 + t2))))))/
             (Power(s,3)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)) + 
            (2*(-7 - 4*s + 5*s1 + 6*s2 + t1 + 3*t2)*
               (-2*Power(s,3)*Power(t1 - t2,2) + 
                 (-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
                 2*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
                    s*(2 + t1 + t2),2) - 
                 Power(s,2)*(-6 - 2*t1 + Power(t1,2) + 
                    6*Power(t1,3) + Power(t1,4) - 2*t2 - 8*t1*t2 - 
                    6*Power(t1,2)*t2 + 2*Power(t1,3)*t2 + 
                    Power(t2,2) - 6*t1*Power(t2,2) - 
                    6*Power(t1,2)*Power(t2,2) + 6*Power(t2,3) + 
                    2*t1*Power(t2,3) + Power(t2,4) + 
                    s2*(6 - Power(t1,3) + 4*t2 + 3*Power(t2,2) + 
                       Power(t2,3) - Power(t1,2)*(5 + 3*t2) + 
                       t1*t2*(8 + 3*t2)) - 
                    s1*(-6 - Power(t1,3) + 5*Power(t2,2) + 
                       Power(t2,3) - 3*Power(t1,2)*(1 + t2) + 
                       t1*(-4 - 8*t2 + 3*Power(t2,2)) + 
                       s2*(6 + Power(t1,2) + 2*t2 + Power(t2,2) + 
                        t1*(2 + 4*t2)))) + 
                 2*s*(-4 + 3*t1 + 4*Power(t1,2) - Power(t1,3) - 
                    2*Power(t1,4) + 3*t2 - 6*t1*t2 + Power(t1,2)*t2 + 
                    2*Power(t1,3)*t2 + 4*Power(t2,2) + 
                    t1*Power(t2,2) - Power(t2,3) + 2*t1*Power(t2,3) - 
                    2*Power(t2,4) + 
                    Power(s2,2)*t2*
                     (2 - 2*Power(t1,2) + 2*t2 + t1*t2 + Power(t2,2)) \
+ s2*(3 - 2*t2 - 2*Power(t2,2) + Power(t2,4) + 
                       2*Power(t1,3)*(1 + t2) + 
                       Power(t1,2)*(-1 + 2*t2 - 3*Power(t2,2)) - 
                       2*t1*(2 + t2 + 2*Power(t2,2))) + 
                    Power(s1,2)*
                     (Power(s2,2)*(2 + t1 + t2) + 
                       s2*(-3 - 2*Power(t1,2) + 2*t2 + Power(t2,2) - 
                        2*t1*(2 + t2)) + 
                       t1*(2 + Power(t1,2) - 2*Power(t2,2) + 
                        t1*(2 + t2))) + 
                    s1*(Power(t1,4) + 2*t1*(-1 + t2)*Power(1 + t2,2) + 
                       Power(-1 + t2,2)*(3 + 2*t2) + 
                       Power(s2,2)*
                        (-3 + Power(t1,2) - 2*t1*(-1 + t2) - 4*t2 - 
                        2*Power(t2,2)) - 
                       Power(t1,2)*(2 + 4*t2 + 3*Power(t2,2)) - 
                       2*s2*(Power(t1,3) - Power(t1,2)*(-1 + t2) + 
                        (-1 + t2)*Power(1 + t2,2) - 
                        t1*(1 + 7*t2 + Power(t2,2)))))))/
             (Power(s,2)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2))))/
        ((-1 + s2)*(-s + s2 - t1)) + 
       ((-8*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + s*(2 + t1 + t2))*
             (-11 + Power(s,2) + 2*Power(s1,2) + 2*Power(s2,2) - 
               7*t1 + 4*s2*t1 - 7*t2 + 3*s2*t2 + 2*t1*t2 - 
               s*(2 + 3*s1 + 3*s2 + 3*t1 + 3*t2) + 
               s1*(5*s2 + 3*t1 + 4*t2)))/s + 
          16*(-4 - 11*s2 + Power(s2,2) - 10*t1 + s2*t1 + 
             2*Power(s2,2)*t1 - 2*Power(t1,2) - 3*t2 - 5*s2*t2 + 
             Power(s2,2)*t2 - 7*t1*t2 + 3*s2*t1*t2 + Power(t2,2) + 
             Power(s,2)*(2 + t1 + t2) + 
             Power(s1,2)*(1 + s2 + t1 + 2*t2) + 
             s1*(-8 + s2 + Power(s2,2) + 3*s2*t1 - 7*t2 + 3*s2*t2 + 
                3*t1*t2) - s*
              (-6 - 2*t2 + 3*t1*t2 + s2*(2 + 3*t1 + 2*t2) + 
                s1*(2 + s2 + 2*t1 + 3*t2))) - 
          (4*(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + s*(2 + t1 + t2))*
             (-4*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
                  s*(2 + t1 + t2),2) - 
               12*s*(2 - s2 - 2*t1 + s2*t1 - Power(t1,2) + 
                  s2*Power(t1,2) - s2*Power(t1,3) + Power(t1,4) + 
                  Power(s,2)*Power(t1 - t2,2) - 2*t2 + s2*t2 + 
                  4*t1*t2 - s2*t1*t2 - Power(s2,2)*t1*t2 + 
                  s2*Power(t1,2)*t2 + Power(s2,2)*Power(t1,2)*t2 - 
                  2*Power(t1,3)*t2 - s2*Power(t1,3)*t2 - Power(t2,2) - 
                  s2*Power(t2,2) + Power(s2,2)*Power(t2,2) - 
                  s2*t1*Power(t2,2) - Power(s2,2)*t1*Power(t2,2) + 
                  2*Power(t1,2)*Power(t2,2) + 
                  2*s2*Power(t1,2)*Power(t2,2) + s2*Power(t2,3) - 
                  2*t1*Power(t2,3) - s2*t1*Power(t2,3) + Power(t2,4) + 
                  Power(s1,2)*(-1 + t2)*
                   (s2*(-1 + t1) + t1*(-t1 + t2)) - 
                  s1*(-1 + t2)*
                   (-1 + 2*s2*(-1 + t1) - Power(s2,2)*(-1 + t1) + 
                     Power(t1,3) + Power(t2,2) - 
                     Power(t1,2)*(1 + 2*t2) + t1*(1 + Power(t2,2))) + 
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
           (Power(s,3)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)) + 
          (8*(-3*s + 2*(-2 + 2*s1 + 2*s2 + t1 + t2))*
             (-2*Power(s,3)*Power(t1 - t2,2) + 
               (-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
               2*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
                  s*(2 + t1 + t2),2) - 
               Power(s,2)*(-6 - 2*t1 + Power(t1,2) + 6*Power(t1,3) + 
                  Power(t1,4) - 2*t2 - 8*t1*t2 - 6*Power(t1,2)*t2 + 
                  2*Power(t1,3)*t2 + Power(t2,2) - 6*t1*Power(t2,2) - 
                  6*Power(t1,2)*Power(t2,2) + 6*Power(t2,3) + 
                  2*t1*Power(t2,3) + Power(t2,4) + 
                  s2*(6 - Power(t1,3) + 4*t2 + 3*Power(t2,2) + 
                     Power(t2,3) - Power(t1,2)*(5 + 3*t2) + 
                     t1*t2*(8 + 3*t2)) - 
                  s1*(-6 - Power(t1,3) + 5*Power(t2,2) + Power(t2,3) - 
                     3*Power(t1,2)*(1 + t2) + 
                     t1*(-4 - 8*t2 + 3*Power(t2,2)) + 
                     s2*(6 + Power(t1,2) + 2*t2 + Power(t2,2) + 
                        t1*(2 + 4*t2)))) + 
               2*s*(-4 + 3*t1 + 4*Power(t1,2) - Power(t1,3) - 
                  2*Power(t1,4) + 3*t2 - 6*t1*t2 + Power(t1,2)*t2 + 
                  2*Power(t1,3)*t2 + 4*Power(t2,2) + t1*Power(t2,2) - 
                  Power(t2,3) + 2*t1*Power(t2,3) - 2*Power(t2,4) + 
                  Power(s2,2)*t2*
                   (2 - 2*Power(t1,2) + 2*t2 + t1*t2 + Power(t2,2)) + 
                  s2*(3 - 2*t2 - 2*Power(t2,2) + Power(t2,4) + 
                     2*Power(t1,3)*(1 + t2) + 
                     Power(t1,2)*(-1 + 2*t2 - 3*Power(t2,2)) - 
                     2*t1*(2 + t2 + 2*Power(t2,2))) + 
                  Power(s1,2)*
                   (Power(s2,2)*(2 + t1 + t2) + 
                     s2*(-3 - 2*Power(t1,2) + 2*t2 + Power(t2,2) - 
                        2*t1*(2 + t2)) + 
                     t1*(2 + Power(t1,2) - 2*Power(t2,2) + t1*(2 + t2))\
) + s1*(Power(t1,4) + 2*t1*(-1 + t2)*Power(1 + t2,2) + 
                     Power(-1 + t2,2)*(3 + 2*t2) + 
                     Power(s2,2)*
                      (-3 + Power(t1,2) - 2*t1*(-1 + t2) - 4*t2 - 
                        2*Power(t2,2)) - 
                     Power(t1,2)*(2 + 4*t2 + 3*Power(t2,2)) - 
                     2*s2*(Power(t1,3) - Power(t1,2)*(-1 + t2) + 
                        (-1 + t2)*Power(1 + t2,2) - 
                        t1*(1 + 7*t2 + Power(t2,2)))))))/
           (Power(s,2)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)))/
        ((-1 + t1)*(-1 + t2)))*
     B3(1 - s2 + t1 - t2,1 - s1 - t1 + t2,2 + s - s1 - s2))/
   (16.*Power(Pi,2)) + (((8*s*(-1 + s1 + t1 - t2)*
          (-3 + 2*Power(s2,2) + s1*(s2 - t1) - t1 - 2*s2*t1 + 
            s*(2 - 2*s2 + t1 - t2) + 3*t2 + s2*t2)*
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
       ((8*(3 - 2*Power(s2,2) + t1 + 2*s2*t1 + s1*(-s2 + t1) - 
               3*t2 - s2*t2 + s*(-2 + 2*s2 - t1 + t2))*
             (-19 + 6*Power(s,2) + 2*Power(s1,2) - 5*s2 + 
               9*Power(s2,2) + 5*t1 - 3*s2*t1 + 
               s1*(5 + 10*s2 - 2*t1 - t2) + t2 - 2*s2*t2 + 
               s*(1 - 8*s1 - 15*s2 + 2*t1 + 2*t2)))/(-1 + s2 - t1 + t2) \
- 16*(12 + 2*Power(s,3) + 32*s2 + 3*Power(s2,2) - 5*Power(s2,3) + t1 - 
             3*s2*t1 + 3*Power(s2,2)*t1 - 5*Power(t1,2) + 
             Power(s1,2)*(-3 - 3*s2 + t1) - 5*t2 - 6*s2*t2 + 
             Power(s2,2)*t2 + 6*t1*t2 - Power(t2,2) + 
             Power(s,2)*(2 - 4*s1 - 9*s2 + t1 + t2) + 
             s*(-20 + 2*Power(s1,2) + 12*Power(s2,2) + t1 + 
                s1*(12*s2 - 2*t1 - t2) + 4*t2 - 2*s2*(2 + 2*t1 + t2)) \
+ s1*(14 - 8*Power(s2,2) - 6*t1 + t2 + s2*(-2 + 4*t1 + t2))) - 
          (8*(1 + 6*s - 4*s1 - 7*s2 + t1 + t2)*
             (2*(-1 + s2 - t1 + t2)*
                Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
               Power(s,3)*(-1 + s1 + t1 - t2)*
                (2 + 2*Power(s2,2) + Power(t1,2) - 
                  2*s2*(2 + t1 - t2) - 2*t1*(-1 + t2) - 2*t2 + 
                  Power(t2,2)) - 
               2*Power(s,2)*
                (-3 + t1 + 2*Power(t1,2) + 2*Power(t1,3) + 
                  Power(s1,2)*
                   (Power(s2,2) - s2*(1 + 2*t1) + t1*(2 + t1 - t2)) \
+ 2*Power(s2,3)*(-1 + t1 - t2) + 2*t2 - 6*t1*t2 - 7*Power(t1,2)*t2 + 
                  4*Power(t2,2) + 8*t1*Power(t2,2) - 3*Power(t2,3) + 
                  Power(s2,2)*
                   (2 + t1 - 3*Power(t1,2) + 5*t1*t2 - 
                     2*Power(t2,2)) + 
                  s1*(3 + 2*Power(s2,3) + Power(t1,3) - 
                     2*Power(t1,2)*(-1 + t2) - 5*t2 + 
                     t1*(-4 + t2)*t2 + 2*Power(t2,2) + 
                     Power(s2,2)*(-3 - 2*t1 + t2) + 
                     s2*(-2 + 2*t1 - Power(t1,2) + 4*t2 + 
                       Power(t2,2))) + 
                  s2*(3 + Power(t1,3) - 5*Power(t2,2) - 
                     Power(t2,3) - Power(t1,2)*(1 + 3*t2) + 
                     t1*(-4 + 6*t2 + 3*Power(t2,2)))) + 
               s*(-5 + Power(s1,3)*Power(s2 - t1,2) - t1 + 
                  Power(t1,2) + 5*Power(t1,3) + 
                  2*Power(s2,4)*(-1 + t1 - t2) + t2 - 18*t1*t2 - 
                  15*Power(t1,2)*t2 + 13*Power(t2,2) + 
                  19*t1*Power(t2,2) - 9*Power(t2,3) - 
                  2*Power(s2,3)*
                   (-2*t1 + 2*Power(t1,2) + t2 - 3*t1*t2 + 
                     Power(t2,2)) + 
                  Power(s1,2)*(s2 - t1)*
                   (-4 + 2*Power(s2,2) + t1 - Power(t1,2) + 4*t2 + 
                     t1*t2 - s2*(1 + t1 + t2)) + 
                  Power(s2,2)*
                   (6 + 2*Power(t1,3) - 8*t2 - 11*Power(t2,2) - 
                     Power(t2,3) - 4*Power(t1,2)*(1 + t2) + 
                     t1*(-4 + 14*t2 + 3*Power(t2,2))) + 
                  2*s2*(Power(t1,3) - 7*Power(t1,2)*t2 + 
                     t2*(7 - 4*t2 - 3*Power(t2,2)) + 
                     t1*(-1 + 6*t2 + 9*Power(t2,2))) + 
                  s1*(2*Power(s2,4) + 4*Power(t1,3) - 
                     2*Power(s2,3)*(1 + t1) + 
                     Power(t1,2)*(17 - 12*t2) + 5*Power(-1 + t2,2) + 
                     2*t1*(3 - 7*t2 + 4*Power(t2,2)) + 
                     Power(s2,2)*
                      (-6 - 2*Power(t1,2) + 10*t2 + Power(t2,2) + 
                        2*t1*(3 + t2)) + 
                     2*s2*(2 - 5*t1 + Power(t1,3) - 3*t2 + 
                        Power(t2,2) - Power(t1,2)*(4 + t2))))))/
           (s*(-1 + s1 + t1 - t2)*Power(-1 + s2 - t1 + t2,2)) - 
          (4*(3 - 2*Power(s2,2) + t1 + 2*s2*t1 + s1*(-s2 + t1) - 
               3*t2 - s2*t2 + s*(-2 + 2*s2 - t1 + t2))*
             (12*(-1 + s2 - t1 + t2)*
                Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
               4*Power(s,3)*(-1 + s1 + t1 - t2)*
                (Power(2 - 2*s2 + t1 - t2,2) - 
                  3*(-1 + s2)*(-1 + s2 - t1 + t2)) - 
               2*Power(s,2)*(-6*(-1 + s2 - t1 + t2)*
                   (3 + Power(s1,2)*s2 - 3*t1 + Power(t1,2) + 
                     2*Power(s2,2)*(-1 + t1 - t2) + 2*t2 - t1*t2 + 
                     s2*t1*(1 - t1 + t2) + 
                     s1*(-3 + 2*Power(s2,2) + t2 - s2*(1 + t2))) + 
                  4*(-2 + 2*s2 - t1 + t2)*
                   (Power(s1,2)*(s2 - t1) + 
                     s1*(-3 + 2*Power(s2,2) - Power(t1,2) - 
                        s2*(1 + t1) + 3*t2 + t1*t2) + 
                     (-1 + t1 - t2)*
                      (-3 + 2*Power(s2,2) - t1 + 3*t2 + 
                        s2*(-2*t1 + t2)))) + 
               s*(4*(-1 + s1 + t1 - t2)*
                   Power(3 - 2*Power(s2,2) + t1 + 2*s2*t1 + 
                     s1*(-s2 + t1) - 3*t2 - s2*t2,2) - 
                  12*(-1 + s2 - t1 + t2)*
                   (Power(s1,2)*(1 + s2)*(s2 - t1) + 
                     Power(s2,3)*(-1 + t1 - t2) - 
                     Power(s2,2)*(-1 + t1)*(-1 + t1 - t2) - 
                     s2*(-2 + t1 + Power(t1,2) - 4*t2 - t1*t2) + 
                     2*(1 + Power(t1,2) - t2 - t1*(2 + t2)) + 
                     s1*(-2 + Power(s2,3) + 5*t1 + Power(t1,2) + 
                        2*t2 - Power(s2,2)*t2 - t1*t2 - 
                        s2*(3 + t1 + Power(t1,2) + t2 - t1*t2))))))/
           (s*(-1 + s1 + t1 - t2)*Power(1 - s2 + t1 - t2,3)))/
        ((-1 + s2)*(-s + s2 - t1)) + 
       (4*(-64*s*(-1 + s1 + t1 - t2)*(s1 - s2 + t1 - t2) + 
            (2*(3 - 2*Power(s2,2) + t1 + 2*s2*t1 + s1*(-s2 + t1) - 
                 3*t2 - s2*t2 + s*(-2 + 2*s2 - t1 + t2))*
               (6*Power(s,2) + s*(11 - 12*s1 - 18*s2 + 3*t2) + 
                 2*(-4 + 3*Power(s1,2) - 7*s2 + 6*Power(s2,2) + 
                    s1*(-6 + 9*s2 - 2*t2) + 4*t2 - 3*s2*t2)))/
             (-1 + s2 - t1 + t2) + 
            32*(1 + 2*Power(s,3) - s2 + 2*Power(s2,2) - 
               2*Power(s2,3) - 3*t1 - 2*Power(s1,2)*t1 + s2*t1 + 
               4*Power(s2,2)*t1 + 2*Power(t1,2) - 4*s2*Power(t1,2) + 
               t2 - 2*s2*t2 - 3*Power(s2,2)*t2 + 6*s2*t1*t2 - 
               2*Power(t2,2) - 2*s2*Power(t2,2) + 
               Power(s,2)*(-1 - 2*s1 - 2*s2 + 2*t1 + t2) + 
               s1*(-2 + s2 + Power(s2,2) + 3*t1 - 3*s2*t1 - 
                  2*Power(t1,2) + 2*t2 + 2*s2*t2 + 2*t1*t2) + 
               s*(-3 + 13*Power(s1,2) + 2*Power(s2,2) - 11*t1 + 
                  17*Power(t1,2) + 10*t2 - 32*t1*t2 + 16*Power(t2,2) + 
                  s2*(14 - 19*t1 + 17*t2) - 
                  s1*(7 + 16*s2 - 28*t1 + 29*t2))) - 
            4*(8 + 30*Power(s,3) + 36*s2 + 22*Power(s2,2) - 
               36*Power(s2,3) + Power(s1,2)*(-15 + 4*s2 - 36*t1) - 
               30*t1 + 56*s2*t1 + 80*Power(s2,2)*t1 - 10*Power(t1,2) - 
               80*s2*Power(t1,2) - 
               2*Power(s,2)*(25 + 7*s1 + 14*s2 - 18*t1 - 3*t2) + 
               11*t2 - 78*s2*t2 - 62*Power(s2,2)*t2 + 85*t1*t2 + 
               120*s2*t1*t2 - 93*Power(t2,2) - 40*s2*Power(t2,2) + 
               s1*(-35 + 18*Power(s2,2) - 40*Power(t1,2) + 106*t2 + 
                  s2*(60 - 46*t1 + 28*t2) + t1*(-9 + 40*t2)) + 
               s*(-49 + 74*Power(s1,2) + 34*Power(s2,2) - 61*t1 + 
                  166*Power(t1,2) + 
                  s1*(27 - 162*s2 + 216*t1 - 214*t2) - t2 - 
                  296*t1*t2 + 143*Power(t2,2) + 
                  s2*(113 - 206*t1 + 180*t2))) + 
            8*(11 + 38*Power(s,3) + Power(s1,3) - 3*s2 + 
               22*Power(s2,2) - 42*Power(s2,3) - 41*t1 + 57*s2*t1 + 
               96*Power(s2,2)*t1 + 2*Power(t1,2) - 96*s2*Power(t1,2) + 
               Power(s1,2)*(-24 + 8*s2 - 44*t1 - t2) + 13*t2 - 
               70*s2*t2 - 76*Power(s2,2)*t2 + 76*t1*t2 + 
               144*s2*t1*t2 - 96*Power(t2,2) - 48*s2*Power(t2,2) + 
               Power(s,2)*(-57 - 22*s1 - 32*s2 + 42*t1 + 11*t2) + 
               s1*(-42 + 28*Power(s2,2) - 48*Power(t1,2) + 116*t2 + 
                  s2*(43 - 59*t1 + 34*t2) + t1*(7 + 48*t2)) + 
               s*(-41 + 106*Power(s1,2) + 36*Power(s2,2) - 87*t1 + 
                  213*Power(t1,2) + 
                  s1*(31 - 215*s2 + 290*t1 - 292*t2) + 13*t2 - 
                  384*t1*t2 + 188*Power(t2,2) + 
                  s2*(160 - 261*t1 + 231*t2))) - 
            16*(3 + 16*Power(s,3) - 6*s2 + 14*Power(s2,2) - 
               18*Power(s2,3) - 19*t1 + 16*s2*t1 + 36*Power(s2,2)*t1 + 
               8*Power(t1,2) - 36*s2*Power(t1,2) - 
               Power(s1,2)*(5 + 17*t1) + 7*t2 - 24*s2*t2 - 
               27*Power(s2,2)*t2 + 16*t1*t2 + 54*s2*t1*t2 - 
               28*Power(t2,2) - 18*s2*Power(t2,2) + 
               Power(s,2)*(-15 - 13*s1 - 16*s2 + 16*t1 + 6*t2) + 
               s*(-18 + 58*Power(s1,2) + 18*Power(s2,2) - 46*t1 + 
                  95*Power(t1,2) + 26*t2 - 174*t1*t2 + 
                  86*Power(t2,2) + s2*(69 - 113*t1 + 99*t2) - 
                  s1*(7 + 89*s2 - 139*t1 + 143*t2)) + 
               s1*(-15 + 7*Power(s2,2) - 18*Power(t1,2) + 32*t2 + 
                  6*t1*(2 + 3*t2) - 8*s2*(3*t1 - 2*(1 + t2)))) - 
            (4*(2 + 3*s - 3*s1 - 4*s2 + t2)*
               (2*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 Power(s,3)*(-1 + s1 + t1 - t2)*
                  (2 + 2*Power(s2,2) + Power(t1,2) - 
                    2*s2*(2 + t1 - t2) - 2*t1*(-1 + t2) - 2*t2 + 
                    Power(t2,2)) - 
                 2*Power(s,2)*
                  (-3 + t1 + 2*Power(t1,2) + 2*Power(t1,3) + 
                    Power(s1,2)*
                     (Power(s2,2) - s2*(1 + 2*t1) + 
                       t1*(2 + t1 - t2)) + 
                    2*Power(s2,3)*(-1 + t1 - t2) + 2*t2 - 6*t1*t2 - 
                    7*Power(t1,2)*t2 + 4*Power(t2,2) + 
                    8*t1*Power(t2,2) - 3*Power(t2,3) + 
                    Power(s2,2)*
                     (2 + t1 - 3*Power(t1,2) + 5*t1*t2 - 
                       2*Power(t2,2)) + 
                    s1*(3 + 2*Power(s2,3) + Power(t1,3) - 
                       2*Power(t1,2)*(-1 + t2) - 5*t2 + 
                       t1*(-4 + t2)*t2 + 2*Power(t2,2) + 
                       Power(s2,2)*(-3 - 2*t1 + t2) + 
                       s2*(-2 + 2*t1 - Power(t1,2) + 4*t2 + 
                       Power(t2,2))) + 
                    s2*(3 + Power(t1,3) - 5*Power(t2,2) - 
                       Power(t2,3) - Power(t1,2)*(1 + 3*t2) + 
                       t1*(-4 + 6*t2 + 3*Power(t2,2)))) + 
                 s*(-5 + Power(s1,3)*Power(s2 - t1,2) - t1 + 
                    Power(t1,2) + 5*Power(t1,3) + 
                    2*Power(s2,4)*(-1 + t1 - t2) + t2 - 18*t1*t2 - 
                    15*Power(t1,2)*t2 + 13*Power(t2,2) + 
                    19*t1*Power(t2,2) - 9*Power(t2,3) - 
                    2*Power(s2,3)*
                     (-2*t1 + 2*Power(t1,2) + t2 - 3*t1*t2 + 
                       Power(t2,2)) + 
                    Power(s1,2)*(s2 - t1)*
                     (-4 + 2*Power(s2,2) + t1 - Power(t1,2) + 4*t2 + 
                       t1*t2 - s2*(1 + t1 + t2)) + 
                    Power(s2,2)*
                     (6 + 2*Power(t1,3) - 8*t2 - 11*Power(t2,2) - 
                       Power(t2,3) - 4*Power(t1,2)*(1 + t2) + 
                       t1*(-4 + 14*t2 + 3*Power(t2,2))) + 
                    2*s2*(Power(t1,3) - 7*Power(t1,2)*t2 + 
                       t2*(7 - 4*t2 - 3*Power(t2,2)) + 
                       t1*(-1 + 6*t2 + 9*Power(t2,2))) + 
                    s1*(2*Power(s2,4) + 4*Power(t1,3) - 
                       2*Power(s2,3)*(1 + t1) + 
                       Power(t1,2)*(17 - 12*t2) + 
                       5*Power(-1 + t2,2) + 
                       2*t1*(3 - 7*t2 + 4*Power(t2,2)) + 
                       Power(s2,2)*
                        (-6 - 2*Power(t1,2) + 10*t2 + Power(t2,2) + 
                        2*t1*(3 + t2)) + 
                       2*s2*(2 - 5*t1 + Power(t1,3) - 3*t2 + 
                        Power(t2,2) - Power(t1,2)*(4 + t2))))))/
             (s*(-1 + s1 + t1 - t2)*Power(-1 + s2 - t1 + t2,2)) - 
            ((3 - 2*Power(s2,2) + t1 + 2*s2*t1 + s1*(-s2 + t1) - 
                 3*t2 - s2*t2 + s*(-2 + 2*s2 - t1 + t2))*
               (12*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 4*Power(s,3)*(-1 + s1 + t1 - t2)*
                  (Power(2 - 2*s2 + t1 - t2,2) - 
                    3*(-1 + s2)*(-1 + s2 - t1 + t2)) - 
                 2*Power(s,2)*
                  (-6*(-1 + s2 - t1 + t2)*
                     (3 + Power(s1,2)*s2 - 3*t1 + Power(t1,2) + 
                       2*Power(s2,2)*(-1 + t1 - t2) + 2*t2 - t1*t2 + 
                       s2*t1*(1 - t1 + t2) + 
                       s1*(-3 + 2*Power(s2,2) + t2 - s2*(1 + t2))) + 
                    4*(-2 + 2*s2 - t1 + t2)*
                     (Power(s1,2)*(s2 - t1) + 
                       s1*(-3 + 2*Power(s2,2) - Power(t1,2) - 
                        s2*(1 + t1) + 3*t2 + t1*t2) + 
                       (-1 + t1 - t2)*
                        (-3 + 2*Power(s2,2) - t1 + 3*t2 + 
                        s2*(-2*t1 + t2)))) + 
                 s*(4*(-1 + s1 + t1 - t2)*
                     Power(3 - 2*Power(s2,2) + t1 + 2*s2*t1 + 
                       s1*(-s2 + t1) - 3*t2 - s2*t2,2) - 
                    12*(-1 + s2 - t1 + t2)*
                     (Power(s1,2)*(1 + s2)*(s2 - t1) + 
                       Power(s2,3)*(-1 + t1 - t2) - 
                       Power(s2,2)*(-1 + t1)*(-1 + t1 - t2) - 
                       s2*(-2 + t1 + Power(t1,2) - 4*t2 - t1*t2) + 
                       2*(1 + Power(t1,2) - t2 - t1*(2 + t2)) + 
                       s1*(-2 + Power(s2,3) + 5*t1 + Power(t1,2) + 
                        2*t2 - Power(s2,2)*t2 - t1*t2 - 
                        s2*(3 + t1 + Power(t1,2) + t2 - t1*t2))))))/
             (s*(-1 + s1 + t1 - t2)*Power(1 - s2 + t1 - t2,3))))/
        ((-1 + s1)*(-s + s1 - t2)) + 
       (4*(64*s*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2) + 
            (2*(3 - 2*Power(s2,2) + t1 + 2*s2*t1 + s1*(-s2 + t1) - 
                 3*t2 - s2*t2 + s*(-2 + 2*s2 - t1 + t2))*
               (12*Power(s,2) + s*(9 - 18*s1 - 24*s2 + 7*t2) + 
                 2*(-3 + 3*Power(s1,2) - 7*s2 + 6*Power(s2,2) + 
                    s1*(-5 + 9*s2 - 2*t2) + 2*t2 - 3*s2*t2)))/
             (-1 + s2 - t1 + t2) - 
            32*(2*Power(s,3) + Power(s1,2)*(s2 - t1) + 
               Power(s,2)*(-2 - s1 - 2*s2 + 2*t1 + t2) + 
               s1*(-1 + 2*Power(s2,2) - 2*Power(t1,2) + t2 + 
                  2*t1*(1 + t2) - s2*(3 + t2)) + 
               (-1 + 2*t1 - 2*t2)*
                (-1 + 2*Power(s2,2) + t1 + t2 + s2*(-2 - 2*t1 + t2)) + 
               s*(-15 + 2*t1 + 17*Power(t1,2) + 
                  s1*(16 - 16*s2 + 15*t1 - 15*t2) - 4*t2 - 32*t1*t2 + 
                  16*Power(t2,2) + s2*(18 - 19*t1 + 17*t2))) + 
            16*(13 + 18*Power(s,3) + 12*s2 - 19*Power(s2,2) + 
               Power(s1,2)*(-1 + 9*s2 - 10*t1) - 29*t1 - 5*s2*t1 + 
               36*Power(s2,2)*t1 + 8*Power(t1,2) - 36*s2*Power(t1,2) - 
               2*Power(s,2)*(11 + 6*s1 + 10*s2 - 8*t1 - 4*t2) + 
               13*t2 + 17*s2*t2 - 36*Power(s2,2)*t2 + 16*t1*t2 + 
               54*s2*t1*t2 - 28*Power(t2,2) - 18*s2*Power(t2,2) + 
               s1*(-13 + 20*Power(s2,2) - 18*Power(t1,2) + 17*t2 + 
                  3*t1*(5 + 6*t2) - s2*(25 + 3*t1 + 7*t2)) + 
               s*(-70 + 4*Power(s1,2) + 2*Power(s2,2) + 12*t1 + 
                  95*Power(t1,2) + s1*(85 - 85*s2 + 81*t1 - 81*t2) - 
                  38*t2 - 174*t1*t2 + 86*Power(t2,2) + 
                  s2*(107 - 113*t1 + 97*t2))) + 
            4*(14 + 40*Power(s,3) - 29*s2 - 47*Power(s2,2) - 
               4*Power(s2,3) + 8*Power(s1,2)*(-2 + 2*s2 - 3*t1) - 
               64*t1 + 24*s2*t1 + 80*Power(s2,2)*t1 - 10*Power(t1,2) - 
               80*s2*Power(t1,2) - 
               2*Power(s,2)*(50 + 12*s1 + 20*s2 - 18*t1 - 7*t2) + 
               57*t2 + 13*s2*t2 - 78*Power(s2,2)*t2 + 85*t1*t2 + 
               120*s2*t1*t2 - 93*Power(t2,2) - 40*s2*Power(t2,2) + 
               s1*(-38 + 42*Power(s2,2) - 40*Power(t1,2) + 69*t2 - 
                  2*s2*(26 + 7*t1 + 4*t2) + t1*(13 + 40*t2)) + 
               s*(-37 + 16*Power(s1,2) + 4*Power(s2,2) + 6*t1 + 
                  166*Power(t1,2) - 106*t2 - 296*t1*t2 + 
                  143*Power(t2,2) - 
                  2*s1*(-81 + 76*s2 - 73*t1 + 70*t2) + 
                  s2*(238 - 206*t1 + 176*t2))) - 
            8*(27 + 52*Power(s,3) - Power(s1,3) + 13*s2 - 
               45*Power(s2,2) - 6*Power(s2,3) - 79*t1 + 11*s2*t1 + 
               96*Power(s2,2)*t1 + 2*Power(t1,2) - 96*s2*Power(t1,2) + 
               49*t2 + 22*s2*t2 - 92*Power(s2,2)*t2 + 76*t1*t2 + 
               144*s2*t1*t2 - 96*Power(t2,2) - 48*s2*Power(t2,2) + 
               Power(s1,2)*(-5 + 16*s2 - 28*t1 + t2) + 
               Power(s,2)*(-82 - 41*s1 - 64*s2 + 42*t1 + 23*t2) + 
               s1*(-38 + 44*Power(s2,2) - 48*Power(t1,2) + 62*t2 - 
                  s2*(51 + 13*t1 + 10*t2) + t1*(29 + 48*t2)) + 
               s*(-115 + 22*Power(s1,2) + 18*Power(s2,2) + 18*t1 + 
                  213*Power(t1,2) + 
                  s1*(186 - 172*s2 + 183*t1 - 185*t2) - 108*t2 - 
                  384*t1*t2 + 188*Power(t2,2) + 
                  s2*(258 - 261*t1 + 215*t2))) - 
            (4*(2 + 4*s - 3*s1 - 4*s2 + t2)*
               (2*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 Power(s,3)*(-1 + s1 + t1 - t2)*
                  (2 + 2*Power(s2,2) + Power(t1,2) - 
                    2*s2*(2 + t1 - t2) - 2*t1*(-1 + t2) - 2*t2 + 
                    Power(t2,2)) - 
                 2*Power(s,2)*
                  (-3 + t1 + 2*Power(t1,2) + 2*Power(t1,3) + 
                    Power(s1,2)*
                     (Power(s2,2) - s2*(1 + 2*t1) + 
                       t1*(2 + t1 - t2)) + 
                    2*Power(s2,3)*(-1 + t1 - t2) + 2*t2 - 6*t1*t2 - 
                    7*Power(t1,2)*t2 + 4*Power(t2,2) + 
                    8*t1*Power(t2,2) - 3*Power(t2,3) + 
                    Power(s2,2)*
                     (2 + t1 - 3*Power(t1,2) + 5*t1*t2 - 
                       2*Power(t2,2)) + 
                    s1*(3 + 2*Power(s2,3) + Power(t1,3) - 
                       2*Power(t1,2)*(-1 + t2) - 5*t2 + 
                       t1*(-4 + t2)*t2 + 2*Power(t2,2) + 
                       Power(s2,2)*(-3 - 2*t1 + t2) + 
                       s2*(-2 + 2*t1 - Power(t1,2) + 4*t2 + 
                       Power(t2,2))) + 
                    s2*(3 + Power(t1,3) - 5*Power(t2,2) - 
                       Power(t2,3) - Power(t1,2)*(1 + 3*t2) + 
                       t1*(-4 + 6*t2 + 3*Power(t2,2)))) + 
                 s*(-5 + Power(s1,3)*Power(s2 - t1,2) - t1 + 
                    Power(t1,2) + 5*Power(t1,3) + 
                    2*Power(s2,4)*(-1 + t1 - t2) + t2 - 18*t1*t2 - 
                    15*Power(t1,2)*t2 + 13*Power(t2,2) + 
                    19*t1*Power(t2,2) - 9*Power(t2,3) - 
                    2*Power(s2,3)*
                     (-2*t1 + 2*Power(t1,2) + t2 - 3*t1*t2 + 
                       Power(t2,2)) + 
                    Power(s1,2)*(s2 - t1)*
                     (-4 + 2*Power(s2,2) + t1 - Power(t1,2) + 4*t2 + 
                       t1*t2 - s2*(1 + t1 + t2)) + 
                    Power(s2,2)*
                     (6 + 2*Power(t1,3) - 8*t2 - 11*Power(t2,2) - 
                       Power(t2,3) - 4*Power(t1,2)*(1 + t2) + 
                       t1*(-4 + 14*t2 + 3*Power(t2,2))) + 
                    2*s2*(Power(t1,3) - 7*Power(t1,2)*t2 + 
                       t2*(7 - 4*t2 - 3*Power(t2,2)) + 
                       t1*(-1 + 6*t2 + 9*Power(t2,2))) + 
                    s1*(2*Power(s2,4) + 4*Power(t1,3) - 
                       2*Power(s2,3)*(1 + t1) + 
                       Power(t1,2)*(17 - 12*t2) + 
                       5*Power(-1 + t2,2) + 
                       2*t1*(3 - 7*t2 + 4*Power(t2,2)) + 
                       Power(s2,2)*
                        (-6 - 2*Power(t1,2) + 10*t2 + Power(t2,2) + 
                        2*t1*(3 + t2)) + 
                       2*s2*(2 - 5*t1 + Power(t1,3) - 3*t2 + 
                        Power(t2,2) - Power(t1,2)*(4 + t2))))))/
             (s*(-1 + s1 + t1 - t2)*Power(-1 + s2 - t1 + t2,2)) - 
            ((3 - 2*Power(s2,2) + t1 + 2*s2*t1 + s1*(-s2 + t1) - 
                 3*t2 - s2*t2 + s*(-2 + 2*s2 - t1 + t2))*
               (12*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 4*Power(s,3)*(-1 + s1 + t1 - t2)*
                  (Power(2 - 2*s2 + t1 - t2,2) - 
                    3*(-1 + s2)*(-1 + s2 - t1 + t2)) - 
                 2*Power(s,2)*
                  (-6*(-1 + s2 - t1 + t2)*
                     (3 + Power(s1,2)*s2 - 3*t1 + Power(t1,2) + 
                       2*Power(s2,2)*(-1 + t1 - t2) + 2*t2 - t1*t2 + 
                       s2*t1*(1 - t1 + t2) + 
                       s1*(-3 + 2*Power(s2,2) + t2 - s2*(1 + t2))) + 
                    4*(-2 + 2*s2 - t1 + t2)*
                     (Power(s1,2)*(s2 - t1) + 
                       s1*(-3 + 2*Power(s2,2) - Power(t1,2) - 
                        s2*(1 + t1) + 3*t2 + t1*t2) + 
                       (-1 + t1 - t2)*
                        (-3 + 2*Power(s2,2) - t1 + 3*t2 + 
                        s2*(-2*t1 + t2)))) + 
                 s*(4*(-1 + s1 + t1 - t2)*
                     Power(3 - 2*Power(s2,2) + t1 + 2*s2*t1 + 
                       s1*(-s2 + t1) - 3*t2 - s2*t2,2) - 
                    12*(-1 + s2 - t1 + t2)*
                     (Power(s1,2)*(1 + s2)*(s2 - t1) + 
                       Power(s2,3)*(-1 + t1 - t2) - 
                       Power(s2,2)*(-1 + t1)*(-1 + t1 - t2) - 
                       s2*(-2 + t1 + Power(t1,2) - 4*t2 - t1*t2) + 
                       2*(1 + Power(t1,2) - t2 - t1*(2 + t2)) + 
                       s1*(-2 + Power(s2,3) + 5*t1 + Power(t1,2) + 
                        2*t2 - Power(s2,2)*t2 - t1*t2 - 
                        s2*(3 + t1 + Power(t1,2) + t2 - t1*t2))))))/
             (s*(-1 + s1 + t1 - t2)*Power(1 - s2 + t1 - t2,3))))/
        ((-1 + s1)*(-1 + t2)) + 
       (4*(64*s*(s - s2 + t1)*(-1 + s1 + t1 - t2) + 
            (2*(3 - 2*Power(s2,2) + t1 + 2*s2*t1 + s1*(-s2 + t1) - 
                 3*t2 - s2*t2 + s*(-2 + 2*s2 - t1 + t2))*
               (-17 + 3*Power(s,2) + 4*Power(s1,2) - 9*s2 + 
                 12*Power(s2,2) + t1 - 3*s2*t1 + 
                 s1*(-5 + 15*s2 - 2*t1 - 2*t2) + 9*t2 - 3*s2*t2 + 
                 s*(9 - 8*s1 - 15*s2 + t1 + t2)))/(-1 + s2 - t1 + t2) + 
            32*(Power(s,3) + 
               Power(s,2)*(9 - 14*s1 + s2 - 13*t1 + 15*t2) + 
               s*(1 - 4*Power(s2,2) + 11*t1 - 15*Power(t1,2) + 
                  s2*(-11 + 20*t1 - 18*t2) - 4*t2 + 17*t1*t2 - 
                  Power(t2,2) + s1*(2 + 15*s2 - 14*t1 + t2)) + 
               (s2 - t1)*(-1 + 2*Power(s2,2) + t1 - s1*t1 + t2 + 
                  s2*(-2 + s1 - 2*t1 + t2))) + 
            8*(-25 + 8*Power(s,3) - 51*s2 - 7*Power(s2,2) + 
               58*Power(s2,3) + 45*t1 + 13*s2*t1 - 
               107*Power(s2,2)*t1 - 8*Power(t1,2) + 
               52*s2*Power(t1,2) + Power(s1,2)*(2 + s2 + 2*t1) + 
               25*t2 + 56*s2*t2 + 35*Power(s2,2)*t2 - 50*t1*t2 - 
               40*s2*t1*t2 + 10*Power(t2,2) + 7*s2*Power(t2,2) + 
               Power(s,2)*(17 - 117*s1 + 39*s2 - 119*t1 + 140*t2) + 
               s1*(-8 + 33*Power(s2,2) + 26*Power(t1,2) + 
                  t1*(5 - 7*t2) - 3*t2 + s2*(3 - 53*t1 + 2*t2)) + 
               s*(59 - 9*Power(s1,2) - 105*Power(s2,2) + 95*t1 - 
                  165*Power(t1,2) + s2*(-108 + 287*t1 - 230*t2) - 
                  130*t2 + 213*t1*t2 - 35*Power(t2,2) + 
                  s1*(77 + 148*s2 - 140*t1 + 39*t2))) - 
            4*(-35 + 3*Power(s,3) - Power(s1,3) - 19*s2 + 
               25*Power(s2,2) + 44*Power(s2,3) + 51*t1 - 17*s2*t1 - 
               91*Power(s2,2)*t1 + 48*s2*Power(t1,2) + 33*t2 + 
               45*s2*t2 + 37*Power(s2,2)*t2 - 58*t1*t2 - 44*s2*t1*t2 + 
               20*Power(t2,2) + 10*s2*Power(t2,2) + 
               Power(s1,2)*(11 - 6*s2 + 5*t1 + t2) + 
               Power(s,2)*(7 - 75*s1 + 28*s2 - 79*t1 + 91*t2) + 
               s1*(3 + 15*Power(s2,2) + 24*Power(t1,2) + 
                  t1*(2 - 10*t2) - 20*t2 + s2*(25 - 38*t1 + 4*t2)) + 
               s*(61 - 11*Power(s1,2) - 75*Power(s2,2) + 68*t1 - 
                  126*Power(t1,2) + s2*(-109 + 220*t1 - 178*t2) - 
                  103*t2 + 173*t1*t2 - 38*Power(t2,2) + 
                  s1*(49 + 120*s2 - 112*t1 + 45*t2))) - 
            16*(-4 + 5*Power(s,3) - 15*s2 - 12*Power(s2,2) + 
               20*Power(s2,3) + 14*t1 + 17*s2*t1 - 38*Power(s2,2)*t1 - 
               6*Power(t1,2) + 18*s2*Power(t1,2) + 3*t2 + 14*s2*t2 + 
               11*Power(s2,2)*t2 - 13*t1*t2 - 11*s2*t1*t2 + 
               Power(t2,2) + s2*Power(t2,2) + 
               Power(s,2)*(23 - 64*s1 + 13*s2 - 62*t1 + 73*t2) + 
               s1*(-1 + 11*Power(s2,2) + 9*Power(t1,2) - 
                  t1*(-2 + t2) + s2*(-1 - 20*t1 + t2)) - 
               s*(-14 + Power(s1,2) + 38*Power(s2,2) - 47*t1 + 
                  77*Power(t1,2) + 43*t2 - 93*t1*t2 + 10*Power(t2,2) - 
                  s1*(25 + 74*s2 - 66*t1 + 10*t2) + 
                  s2*(49 - 124*t1 + 103*t2))) - 
            (2*(3 + 5*s - 5*s1 - 8*s2 + t1 + t2)*
               (2*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 Power(s,3)*(-1 + s1 + t1 - t2)*
                  (2 + 2*Power(s2,2) + Power(t1,2) - 
                    2*s2*(2 + t1 - t2) - 2*t1*(-1 + t2) - 2*t2 + 
                    Power(t2,2)) - 
                 2*Power(s,2)*
                  (-3 + t1 + 2*Power(t1,2) + 2*Power(t1,3) + 
                    Power(s1,2)*
                     (Power(s2,2) - s2*(1 + 2*t1) + 
                       t1*(2 + t1 - t2)) + 
                    2*Power(s2,3)*(-1 + t1 - t2) + 2*t2 - 6*t1*t2 - 
                    7*Power(t1,2)*t2 + 4*Power(t2,2) + 
                    8*t1*Power(t2,2) - 3*Power(t2,3) + 
                    Power(s2,2)*
                     (2 + t1 - 3*Power(t1,2) + 5*t1*t2 - 
                       2*Power(t2,2)) + 
                    s1*(3 + 2*Power(s2,3) + Power(t1,3) - 
                       2*Power(t1,2)*(-1 + t2) - 5*t2 + 
                       t1*(-4 + t2)*t2 + 2*Power(t2,2) + 
                       Power(s2,2)*(-3 - 2*t1 + t2) + 
                       s2*(-2 + 2*t1 - Power(t1,2) + 4*t2 + 
                       Power(t2,2))) + 
                    s2*(3 + Power(t1,3) - 5*Power(t2,2) - 
                       Power(t2,3) - Power(t1,2)*(1 + 3*t2) + 
                       t1*(-4 + 6*t2 + 3*Power(t2,2)))) + 
                 s*(-5 + Power(s1,3)*Power(s2 - t1,2) - t1 + 
                    Power(t1,2) + 5*Power(t1,3) + 
                    2*Power(s2,4)*(-1 + t1 - t2) + t2 - 18*t1*t2 - 
                    15*Power(t1,2)*t2 + 13*Power(t2,2) + 
                    19*t1*Power(t2,2) - 9*Power(t2,3) - 
                    2*Power(s2,3)*
                     (-2*t1 + 2*Power(t1,2) + t2 - 3*t1*t2 + 
                       Power(t2,2)) + 
                    Power(s1,2)*(s2 - t1)*
                     (-4 + 2*Power(s2,2) + t1 - Power(t1,2) + 4*t2 + 
                       t1*t2 - s2*(1 + t1 + t2)) + 
                    Power(s2,2)*
                     (6 + 2*Power(t1,3) - 8*t2 - 11*Power(t2,2) - 
                       Power(t2,3) - 4*Power(t1,2)*(1 + t2) + 
                       t1*(-4 + 14*t2 + 3*Power(t2,2))) + 
                    2*s2*(Power(t1,3) - 7*Power(t1,2)*t2 + 
                       t2*(7 - 4*t2 - 3*Power(t2,2)) + 
                       t1*(-1 + 6*t2 + 9*Power(t2,2))) + 
                    s1*(2*Power(s2,4) + 4*Power(t1,3) - 
                       2*Power(s2,3)*(1 + t1) + 
                       Power(t1,2)*(17 - 12*t2) + 
                       5*Power(-1 + t2,2) + 
                       2*t1*(3 - 7*t2 + 4*Power(t2,2)) + 
                       Power(s2,2)*
                        (-6 - 2*Power(t1,2) + 10*t2 + Power(t2,2) + 
                        2*t1*(3 + t2)) + 
                       2*s2*(2 - 5*t1 + Power(t1,3) - 3*t2 + 
                        Power(t2,2) - Power(t1,2)*(4 + t2))))))/
             (s*(-1 + s1 + t1 - t2)*Power(-1 + s2 - t1 + t2,2)) - 
            ((3 - 2*Power(s2,2) + t1 + 2*s2*t1 + s1*(-s2 + t1) - 
                 3*t2 - s2*t2 + s*(-2 + 2*s2 - t1 + t2))*
               (12*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 4*Power(s,3)*(-1 + s1 + t1 - t2)*
                  (Power(2 - 2*s2 + t1 - t2,2) - 
                    3*(-1 + s2)*(-1 + s2 - t1 + t2)) - 
                 2*Power(s,2)*
                  (-6*(-1 + s2 - t1 + t2)*
                     (3 + Power(s1,2)*s2 - 3*t1 + Power(t1,2) + 
                       2*Power(s2,2)*(-1 + t1 - t2) + 2*t2 - t1*t2 + 
                       s2*t1*(1 - t1 + t2) + 
                       s1*(-3 + 2*Power(s2,2) + t2 - s2*(1 + t2))) + 
                    4*(-2 + 2*s2 - t1 + t2)*
                     (Power(s1,2)*(s2 - t1) + 
                       s1*(-3 + 2*Power(s2,2) - Power(t1,2) - 
                        s2*(1 + t1) + 3*t2 + t1*t2) + 
                       (-1 + t1 - t2)*
                        (-3 + 2*Power(s2,2) - t1 + 3*t2 + 
                        s2*(-2*t1 + t2)))) + 
                 s*(4*(-1 + s1 + t1 - t2)*
                     Power(3 - 2*Power(s2,2) + t1 + 2*s2*t1 + 
                       s1*(-s2 + t1) - 3*t2 - s2*t2,2) - 
                    12*(-1 + s2 - t1 + t2)*
                     (Power(s1,2)*(1 + s2)*(s2 - t1) + 
                       Power(s2,3)*(-1 + t1 - t2) - 
                       Power(s2,2)*(-1 + t1)*(-1 + t1 - t2) - 
                       s2*(-2 + t1 + Power(t1,2) - 4*t2 - t1*t2) + 
                       2*(1 + Power(t1,2) - t2 - t1*(2 + t2)) + 
                       s1*(-2 + Power(s2,3) + 5*t1 + Power(t1,2) + 
                        2*t2 - Power(s2,2)*t2 - t1*t2 - 
                        s2*(3 + t1 + Power(t1,2) + t2 - t1*t2))))))/
             (s*(-1 + s1 + t1 - t2)*Power(1 - s2 + t1 - t2,3))))/
        ((s - s2 + t1)*(s - s1 + t2)) + 
       (4*(64*s*(s2 - t1)*(-1 + s1 + t1 - t2) + 
            (2*(3 - 2*Power(s2,2) + t1 + 2*s2*t1 + s1*(-s2 + t1) - 
                 3*t2 - s2*t2 + s*(-2 + 2*s2 - t1 + t2))*
               (-13 + 12*Power(s,2) + 4*Power(s1,2) - 9*s2 + 
                 12*Power(s2,2) + t1 - 3*s2*t1 + 
                 s1*(-1 + 15*s2 - 2*t1 - 2*t2) + t2 - 3*s2*t2 + 
                 s*(8 - 15*s1 - 24*s2 + 3*t1 + 4*t2)))/
             (-1 + s2 - t1 + t2) + 
            16*(-8 + Power(s,3) - 10*s2 + 17*Power(s2,2) - 
               Power(s2,3) + 16*t1 - 4*s2*t1 - 18*Power(s2,2)*t1 - 
               6*Power(t1,2) + 18*s2*Power(t1,2) + 
               Power(s1,2)*(-s2 + t1) + 
               Power(s,2)*(-4*s1 - 19*s2 + 6*(2 + t1)) + 9*t2 + 
               10*s2*t2 + 3*Power(s2,2)*t2 - 13*t1*t2 - 11*s2*t1*t2 + 
               Power(t2,2) + s2*Power(t2,2) - 
               s1*(1 + 4*Power(s2,2) + 6*s2*(-1 + t1) - 
                  9*Power(t1,2) + t2 + t1*(6 + t2)) + 
               s*(30 + Power(s1,2) + 19*Power(s2,2) + 46*t1 - 
                  78*Power(t1,2) + s2*(-95 + 73*t1 - 72*t2) - 12*t2 + 
                  94*t1*t2 - 10*Power(t2,2) + 
                  s1*(-9 + 85*s2 - 78*t1 + 10*t2))) + 
            32*(Power(s,2)*(-1 + 2*s2 - t1) + 
               (-1 + t1)*(-1 + 2*Power(s2,2) + t1 - s1*t1 + t2 + 
                  s2*(-2 + s1 - 2*t1 + t2)) + 
               s*(-4 - 2*Power(s2,2) - 12*t1 + 15*Power(t1,2) + 2*t2 - 
                  17*t1*t2 + Power(t2,2) - 
                  s1*(-1 + 15*s2 - 15*t1 + t2) + 
                  s2*(17 - 14*t1 + 14*t2))) - 
            8*(-21 + 6*Power(s,3) - 10*s2 + 35*Power(s2,2) - 
               6*Power(s2,3) + 41*t1 - 21*s2*t1 - 49*Power(s2,2)*t1 - 
               8*Power(t1,2) + 52*s2*Power(t1,2) + 
               Power(s1,2)*(-2 - 6*s2 + 7*t1) + 29*t2 + 35*s2*t2 + 
               17*Power(s2,2)*t2 - 50*t1*t2 - 40*s2*t1*t2 + 
               10*Power(t2,2) + 7*s2*Power(t2,2) + 
               Power(s,2)*(34 - 23*s1 - 64*s2 + 12*t1 + 3*t2) - 
               s1*(6 + 23*Power(s2,2) - 26*Power(t1,2) + 9*t2 + 
                  s2*(3 + 5*t1 + t2) + t1*(8 + 7*t2)) + 
               s*(70 + 6*Power(s1,2) + 64*Power(s2,2) + 67*t1 - 
                  170*Power(t1,2) + 2*s2*(-100 + 80*t1 - 83*t2) - 
                  19*t2 + 218*t1*t2 - 35*Power(t2,2) + 
                  s1*(-19 + 215*s2 - 174*t1 + 33*t2))) - 
            4*(27 - Power(s1,3) + 31*s2 - 9*Power(s2,2) - 45*t1 + 
               23*s2*t1 + 45*Power(s2,2)*t1 - 48*s2*Power(t1,2) + 
               Power(s,2)*(-21 + 15*s1 + 44*s2 - 5*t1 - t2) - 33*t2 - 
               43*s2*t2 - 19*Power(s2,2)*t2 + 58*t1*t2 + 44*s2*t1*t2 - 
               20*Power(t2,2) - 10*s2*Power(t2,2) + 
               Power(s1,2)*(5 - 9*t1 + t2) + 
               s1*(17 + 15*Power(s2,2) - 24*Power(t1,2) + 16*t2 + 
                  10*t1*t2 + s2*(29 - 2*t1 + 6*t2)) + 
               s*(-83 - 44*Power(s2,2) - 29*t1 + 132*Power(t1,2) + 
                  s1*(4 - 164*s2 + 136*t1 - 39*t2) + 11*t2 - 
                  179*t1*t2 + 38*Power(t2,2) + 
                  s2*(121 - 130*t1 + 132*t2))) - 
            (2*(3 + 8*s - 5*s1 - 8*s2 + t1 + t2)*
               (2*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 Power(s,3)*(-1 + s1 + t1 - t2)*
                  (2 + 2*Power(s2,2) + Power(t1,2) - 
                    2*s2*(2 + t1 - t2) - 2*t1*(-1 + t2) - 2*t2 + 
                    Power(t2,2)) - 
                 2*Power(s,2)*
                  (-3 + t1 + 2*Power(t1,2) + 2*Power(t1,3) + 
                    Power(s1,2)*
                     (Power(s2,2) - s2*(1 + 2*t1) + 
                       t1*(2 + t1 - t2)) + 
                    2*Power(s2,3)*(-1 + t1 - t2) + 2*t2 - 6*t1*t2 - 
                    7*Power(t1,2)*t2 + 4*Power(t2,2) + 
                    8*t1*Power(t2,2) - 3*Power(t2,3) + 
                    Power(s2,2)*
                     (2 + t1 - 3*Power(t1,2) + 5*t1*t2 - 
                       2*Power(t2,2)) + 
                    s1*(3 + 2*Power(s2,3) + Power(t1,3) - 
                       2*Power(t1,2)*(-1 + t2) - 5*t2 + 
                       t1*(-4 + t2)*t2 + 2*Power(t2,2) + 
                       Power(s2,2)*(-3 - 2*t1 + t2) + 
                       s2*(-2 + 2*t1 - Power(t1,2) + 4*t2 + 
                       Power(t2,2))) + 
                    s2*(3 + Power(t1,3) - 5*Power(t2,2) - 
                       Power(t2,3) - Power(t1,2)*(1 + 3*t2) + 
                       t1*(-4 + 6*t2 + 3*Power(t2,2)))) + 
                 s*(-5 + Power(s1,3)*Power(s2 - t1,2) - t1 + 
                    Power(t1,2) + 5*Power(t1,3) + 
                    2*Power(s2,4)*(-1 + t1 - t2) + t2 - 18*t1*t2 - 
                    15*Power(t1,2)*t2 + 13*Power(t2,2) + 
                    19*t1*Power(t2,2) - 9*Power(t2,3) - 
                    2*Power(s2,3)*
                     (-2*t1 + 2*Power(t1,2) + t2 - 3*t1*t2 + 
                       Power(t2,2)) + 
                    Power(s1,2)*(s2 - t1)*
                     (-4 + 2*Power(s2,2) + t1 - Power(t1,2) + 4*t2 + 
                       t1*t2 - s2*(1 + t1 + t2)) + 
                    Power(s2,2)*
                     (6 + 2*Power(t1,3) - 8*t2 - 11*Power(t2,2) - 
                       Power(t2,3) - 4*Power(t1,2)*(1 + t2) + 
                       t1*(-4 + 14*t2 + 3*Power(t2,2))) + 
                    2*s2*(Power(t1,3) - 7*Power(t1,2)*t2 + 
                       t2*(7 - 4*t2 - 3*Power(t2,2)) + 
                       t1*(-1 + 6*t2 + 9*Power(t2,2))) + 
                    s1*(2*Power(s2,4) + 4*Power(t1,3) - 
                       2*Power(s2,3)*(1 + t1) + 
                       Power(t1,2)*(17 - 12*t2) + 
                       5*Power(-1 + t2,2) + 
                       2*t1*(3 - 7*t2 + 4*Power(t2,2)) + 
                       Power(s2,2)*
                        (-6 - 2*Power(t1,2) + 10*t2 + Power(t2,2) + 
                        2*t1*(3 + t2)) + 
                       2*s2*(2 - 5*t1 + Power(t1,3) - 3*t2 + 
                        Power(t2,2) - Power(t1,2)*(4 + t2))))))/
             (s*(-1 + s1 + t1 - t2)*Power(-1 + s2 - t1 + t2,2)) - 
            ((3 - 2*Power(s2,2) + t1 + 2*s2*t1 + s1*(-s2 + t1) - 
                 3*t2 - s2*t2 + s*(-2 + 2*s2 - t1 + t2))*
               (12*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 4*Power(s,3)*(-1 + s1 + t1 - t2)*
                  (Power(2 - 2*s2 + t1 - t2,2) - 
                    3*(-1 + s2)*(-1 + s2 - t1 + t2)) - 
                 2*Power(s,2)*
                  (-6*(-1 + s2 - t1 + t2)*
                     (3 + Power(s1,2)*s2 - 3*t1 + Power(t1,2) + 
                       2*Power(s2,2)*(-1 + t1 - t2) + 2*t2 - t1*t2 + 
                       s2*t1*(1 - t1 + t2) + 
                       s1*(-3 + 2*Power(s2,2) + t2 - s2*(1 + t2))) + 
                    4*(-2 + 2*s2 - t1 + t2)*
                     (Power(s1,2)*(s2 - t1) + 
                       s1*(-3 + 2*Power(s2,2) - Power(t1,2) - 
                        s2*(1 + t1) + 3*t2 + t1*t2) + 
                       (-1 + t1 - t2)*
                        (-3 + 2*Power(s2,2) - t1 + 3*t2 + 
                        s2*(-2*t1 + t2)))) + 
                 s*(4*(-1 + s1 + t1 - t2)*
                     Power(3 - 2*Power(s2,2) + t1 + 2*s2*t1 + 
                       s1*(-s2 + t1) - 3*t2 - s2*t2,2) - 
                    12*(-1 + s2 - t1 + t2)*
                     (Power(s1,2)*(1 + s2)*(s2 - t1) + 
                       Power(s2,3)*(-1 + t1 - t2) - 
                       Power(s2,2)*(-1 + t1)*(-1 + t1 - t2) - 
                       s2*(-2 + t1 + Power(t1,2) - 4*t2 - t1*t2) + 
                       2*(1 + Power(t1,2) - t2 - t1*(2 + t2)) + 
                       s1*(-2 + Power(s2,3) + 5*t1 + Power(t1,2) + 
                        2*t2 - Power(s2,2)*t2 - t1*t2 - 
                        s2*(3 + t1 + Power(t1,2) + t2 - t1*t2))))))/
             (s*(-1 + s1 + t1 - t2)*Power(1 - s2 + t1 - t2,3))))/
        ((-1 + t1)*(-1 + t2)) + 
       (4*(-64*s*(-1 + s2)*(-1 + s1 + t1 - t2) + 
            (2*(3 - 2*Power(s2,2) + t1 + 2*s2*t1 + s1*(-s2 + t1) - 
                 3*t2 - s2*t2 + s*(-2 + 2*s2 - t1 + t2))*
               (-17 + 9*Power(s,2) + 2*Power(s1,2) - 5*s2 + 
                 9*Power(s2,2) + 5*t1 - 3*s2*t1 + 
                 s1*(7 + 10*s2 - 2*t1 - t2) - 3*t2 - 2*s2*t2 + 
                 s*(3 - 10*s1 - 18*s2 + 3*t1 + 2*t2)))/
             (-1 + s2 - t1 + t2) - 
            32*(Power(s,2)*(-4 + 4*s2 - t1 + t2) - 
               (-1 + t1)*(-1 + 2*Power(s2,2) + t1 - s1*t1 + t2 + 
                  s2*(-2 + s1 - 2*t1 + t2)) + 
               s*(-10 - 4*Power(s2,2) + 8*t1 + Power(t1,2) - 14*t2 - 
                  2*t1*t2 + Power(t2,2) - 
                  s1*(-14 + 16*s2 - 3*t1 + t2) + s2*(16 - 9*t1 + 13*t2))\
) + 16*(-8 - 7*s2 + 18*Power(s2,2) + 8*t1 - 8*s2*t1 - 
               17*Power(s2,2)*t1 + 14*s2*Power(t1,2) - 
               Power(s1,2)*(s2 + t1) + 19*t2 + 6*s2*t2 + 
               Power(s2,2)*t2 - 13*t1*t2 - 5*s2*t1*t2 - Power(t2,2) - 
               s2*Power(t2,2) + 
               Power(s,2)*(-32 + 34*s2 - 12*t1 + 11*t2) + 
               s1*(-8 - Power(s2,2) + 7*Power(t1,2) + t1*(-4 + t2) + 
                  t2 + s2*(9 - 11*t1 + 2*t2)) + 
               s*(-28 + Power(s1,2) - 34*Power(s2,2) + 25*t1 + 
                  9*Power(t1,2) + s1*(67 - 84*s2 + 30*t1 - 11*t2) - 
                  71*t2 - 19*t1*t2 + 9*Power(t2,2) + 
                  s2*(72 - 22*t1 + 56*t2))) + 
            4*(-31 + Power(s,3) - 21*s2 + 25*Power(s2,2) - 
               Power(s2,3) - 4*t1 - 32*s2*t1 - 23*Power(s2,2)*t1 + 
               15*Power(t1,2) + 12*s2*Power(t1,2) - 
               Power(s1,2)*(6 + s2 + 7*t1) + 76*t2 + 16*s2*t2 - 
               5*Power(s2,2)*t2 - 20*t1*t2 + 10*s2*t1*t2 - 
               15*Power(t2,2) - 8*s2*Power(t2,2) + 
               Power(s,2)*(-72 + 6*s1 + 65*s2 - 19*t1 + 15*t2) + 
               s1*(-57 + 2*Power(s2,2) + 6*Power(t1,2) + 18*t2 + 
                  s2*(2 - 28*t1 + 7*t2) + t1*(-2 + 8*t2)) + 
               s*(41 + 3*Power(s1,2) - 65*Power(s2,2) + 18*t1 + 
                  36*Power(t1,2) + s1*(95 - 142*s2 + 78*t1 - 29*t2) - 
                  109*t2 - 66*t1*t2 + 26*Power(t2,2) + 
                  s2*(104 - 24*t1 + 90*t2))) - 
            8*(-20 + 3*Power(s,3) + 4*s2 + 44*Power(s2,2) - 
               3*Power(s2,3) - 4*t1 - 40*s2*t1 - 36*Power(s2,2)*t1 + 
               14*Power(t1,2) + 26*s2*Power(t1,2) - 
               Power(s1,2)*(5 + 4*s2 + 5*t1) - 
               Power(s,2)*(82 + s1 - 77*s2 + 28*t1 - 26*t2) + 77*t2 + 
               19*s2*t2 - 35*t1*t2 - s2*t1*t2 - 9*Power(t2,2) - 
               6*s2*Power(t2,2) + 
               s1*(-42 - 5*Power(s2,2) - 5*t1 + 13*Power(t1,2) + 
                  12*t2 + 6*t1*t2 + s2*(11 - 30*t1 + 8*t2)) + 
               s*(-24 + 5*Power(s1,2) - 77*Power(s2,2) + 43*t1 + 
                  32*Power(t1,2) + s1*(139 - 175*s2 + 85*t1 - 33*t2) - 
                  153*t2 - 63*t1*t2 + 27*Power(t2,2) + 
                  s2*(133 - 31*t1 + 112*t2))) - 
            (2*(1 + 7*s - 4*s1 - 7*s2 + t1 + t2)*
               (2*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 Power(s,3)*(-1 + s1 + t1 - t2)*
                  (2 + 2*Power(s2,2) + Power(t1,2) - 
                    2*s2*(2 + t1 - t2) - 2*t1*(-1 + t2) - 2*t2 + 
                    Power(t2,2)) - 
                 2*Power(s,2)*
                  (-3 + t1 + 2*Power(t1,2) + 2*Power(t1,3) + 
                    Power(s1,2)*
                     (Power(s2,2) - s2*(1 + 2*t1) + t1*(2 + t1 - t2)) \
+ 2*Power(s2,3)*(-1 + t1 - t2) + 2*t2 - 6*t1*t2 - 7*Power(t1,2)*t2 + 
                    4*Power(t2,2) + 8*t1*Power(t2,2) - 
                    3*Power(t2,3) + 
                    Power(s2,2)*
                     (2 + t1 - 3*Power(t1,2) + 5*t1*t2 - 
                       2*Power(t2,2)) + 
                    s1*(3 + 2*Power(s2,3) + Power(t1,3) - 
                       2*Power(t1,2)*(-1 + t2) - 5*t2 + 
                       t1*(-4 + t2)*t2 + 2*Power(t2,2) + 
                       Power(s2,2)*(-3 - 2*t1 + t2) + 
                       s2*(-2 + 2*t1 - Power(t1,2) + 4*t2 + 
                        Power(t2,2))) + 
                    s2*(3 + Power(t1,3) - 5*Power(t2,2) - 
                       Power(t2,3) - Power(t1,2)*(1 + 3*t2) + 
                       t1*(-4 + 6*t2 + 3*Power(t2,2)))) + 
                 s*(-5 + Power(s1,3)*Power(s2 - t1,2) - t1 + 
                    Power(t1,2) + 5*Power(t1,3) + 
                    2*Power(s2,4)*(-1 + t1 - t2) + t2 - 18*t1*t2 - 
                    15*Power(t1,2)*t2 + 13*Power(t2,2) + 
                    19*t1*Power(t2,2) - 9*Power(t2,3) - 
                    2*Power(s2,3)*
                     (-2*t1 + 2*Power(t1,2) + t2 - 3*t1*t2 + 
                       Power(t2,2)) + 
                    Power(s1,2)*(s2 - t1)*
                     (-4 + 2*Power(s2,2) + t1 - Power(t1,2) + 4*t2 + 
                       t1*t2 - s2*(1 + t1 + t2)) + 
                    Power(s2,2)*
                     (6 + 2*Power(t1,3) - 8*t2 - 11*Power(t2,2) - 
                       Power(t2,3) - 4*Power(t1,2)*(1 + t2) + 
                       t1*(-4 + 14*t2 + 3*Power(t2,2))) + 
                    2*s2*(Power(t1,3) - 7*Power(t1,2)*t2 + 
                       t2*(7 - 4*t2 - 3*Power(t2,2)) + 
                       t1*(-1 + 6*t2 + 9*Power(t2,2))) + 
                    s1*(2*Power(s2,4) + 4*Power(t1,3) - 
                       2*Power(s2,3)*(1 + t1) + 
                       Power(t1,2)*(17 - 12*t2) + 5*Power(-1 + t2,2) + 
                       2*t1*(3 - 7*t2 + 4*Power(t2,2)) + 
                       Power(s2,2)*
                        (-6 - 2*Power(t1,2) + 10*t2 + Power(t2,2) + 
                        2*t1*(3 + t2)) + 
                       2*s2*(2 - 5*t1 + Power(t1,3) - 3*t2 + 
                        Power(t2,2) - Power(t1,2)*(4 + t2))))))/
             (s*(-1 + s1 + t1 - t2)*Power(-1 + s2 - t1 + t2,2)) - 
            ((3 - 2*Power(s2,2) + t1 + 2*s2*t1 + s1*(-s2 + t1) - 
                 3*t2 - s2*t2 + s*(-2 + 2*s2 - t1 + t2))*
               (12*(-1 + s2 - t1 + t2)*
                  Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) + 
                 4*Power(s,3)*(-1 + s1 + t1 - t2)*
                  (Power(2 - 2*s2 + t1 - t2,2) - 
                    3*(-1 + s2)*(-1 + s2 - t1 + t2)) - 
                 2*Power(s,2)*
                  (-6*(-1 + s2 - t1 + t2)*
                     (3 + Power(s1,2)*s2 - 3*t1 + Power(t1,2) + 
                       2*Power(s2,2)*(-1 + t1 - t2) + 2*t2 - t1*t2 + 
                       s2*t1*(1 - t1 + t2) + 
                       s1*(-3 + 2*Power(s2,2) + t2 - s2*(1 + t2))) + 
                    4*(-2 + 2*s2 - t1 + t2)*
                     (Power(s1,2)*(s2 - t1) + 
                       s1*(-3 + 2*Power(s2,2) - Power(t1,2) - 
                        s2*(1 + t1) + 3*t2 + t1*t2) + 
                       (-1 + t1 - t2)*
                        (-3 + 2*Power(s2,2) - t1 + 3*t2 + 
                        s2*(-2*t1 + t2)))) + 
                 s*(4*(-1 + s1 + t1 - t2)*
                     Power(3 - 2*Power(s2,2) + t1 + 2*s2*t1 + 
                       s1*(-s2 + t1) - 3*t2 - s2*t2,2) - 
                    12*(-1 + s2 - t1 + t2)*
                     (Power(s1,2)*(1 + s2)*(s2 - t1) + 
                       Power(s2,3)*(-1 + t1 - t2) - 
                       Power(s2,2)*(-1 + t1)*(-1 + t1 - t2) - 
                       s2*(-2 + t1 + Power(t1,2) - 4*t2 - t1*t2) + 
                       2*(1 + Power(t1,2) - t2 - t1*(2 + t2)) + 
                       s1*(-2 + Power(s2,3) + 5*t1 + Power(t1,2) + 
                         2*t2 - Power(s2,2)*t2 - t1*t2 - 
                         s2*(3 + t1 + Power(t1,2) + t2 - t1*t2))))))/
             (s*(-1 + s1 + t1 - t2)*Power(1 - s2 + t1 - t2,3))))/
        ((-1 + s2)*(-1 + t1)))*B3(1 - s1 - t1 + t2,s,2 + s - s1 - s2))/
   (16.*Power(Pi,2)) + (((-2*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)*
          (1 - 2*Power(s,2) - t1 + s1*(-s2 + t1) + 
            s*(2 + 2*s1 + 2*s2 - t1 - t2) - t2 + s2*t2)*
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
               Power(s1,2)*(13 - 2*Power(s2,2) + Power(t1,2) - 3*t2 - 
                  t1*(2 + 3*t2) + s2*(2 + t1 + 3*t2)) + 
               Power(s,2)*(16 - 5*Power(s1,2) + 4*t1 + Power(t1,2) - 
                  8*t2 - t1*t2 - Power(t2,2) + 
                  s2*(-2 - 2*t1 + 5*t2) + s1*(3 - 7*s2 + 2*t1 + 5*t2)) \
+ s1*(21 + Power(s2,3) - 18*t2 + 5*Power(t2,2) - 
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
               3*Power(t1,2)*t2 - s2*Power(t1,2)*t2 - 12*Power(t2,2) - 
               2*s2*Power(t2,2) + 5*t1*Power(t2,2) - 2*Power(t2,3) + 
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
               Power(s,2)*(9 - 2*s1*(-4 + t1) + Power(t1,2) - 10*t2 + 
                  t1*(5 - 2*s2 + t2)) + 
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
                  3*t1 - 2*Power(t1,2) - t2 - 3*t1*t2 + Power(t2,2) + 
                  s1*(1 - 14*s2 + 9*t1 + t2) + 2*s2*(-1 + 5*t1 + 2*t2)\
) + s1*(3 - 4*Power(s2,3) + 3*t2 - 6*Power(t2,2) + 
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
               Power(s,2)*(21 - 9*Power(s1,2) - 9*Power(s2,2) - 2*t1 - 
                  4*t2 - 4*t1*t2 + s1*(5 - 19*s2 + 6*t1 + 7*t2) + 
                  s2*(4 + 5*t1 + 9*t2)) + 
               s*(-23 + 5*Power(s1,3) + 5*Power(s2,3) + 16*t1 + 
                  Power(s1,2)*(-6 + 17*s2 - 6*t1 - 5*t2) + 20*t2 - 
                  5*t1*t2 - 2*t1*Power(t2,2) + Power(t2,3) - 
                  Power(s2,2)*(4 + 4*t1 + 9*t2) + 
                  s2*(-37 + 6*t1 + 6*t2 + 7*t1*t2 + Power(t2,2)) + 
                  s1*(-36 + 17*Power(s2,2) + 5*t1 + 7*t2 + 8*t1*t2 - 
                     Power(t2,2) - 5*s2*(2 + 2*t1 + 3*t2))) - 
               s1*(-17 + 5*Power(s2,3) + 13*t2 + Power(t2,2) + 
                  Power(t2,3) - Power(s2,2)*(5 + 4*t1 + 8*t2) - 
                  2*t1*(-9 + 3*t2 + Power(t2,2)) + 
                  s2*(-33 + 5*t2 + t1*(5 + 7*t2))))/((-1 + t1)*(-1 + t2))\
))/((s - s2 + t1)*Power(-1 + s1 + t1 - t2,3)*(-1 + s2 - t1 + t2) - 
          Power(-1 + s2 + s*t1 - s2*t1 + Power(t1,2) + t2 - s*t2 - t1*t2,
           2) - Power(-1 + s1 + t1 - t2,2)*
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
       (64*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)*(-1 + s + t2) + 
          (2*(-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
               s*(-2 - 2*s1 - 2*s2 + t1 + t2))*
             (-16 + 12*Power(s,2) + s1 + 4*Power(s1,2) + 14*s1*s2 + 
               9*Power(s2,2) + 6*t1 - 9*s1*t1 - 10*s2*t1 + 
               2*Power(t1,2) + t2 - 2*s2*t2 + t1*t2 - Power(t2,2) + 
               s*(-3 - 15*s1 - 21*s2 + 12*t1 + t2)))/s - 
          32*(13 + 2*Power(s,3) - 12*s2 - 3*t1 + 13*s2*t1 + 
             Power(s2,2)*t1 - 13*Power(t1,2) - 2*s2*Power(t1,2) + 
             Power(t1,3) - s1*(11*s2 - 11*t1 + 12*(-1 + t2))*
              (-1 + t2) - 12*t2 - Power(s2,2)*t2 + 28*t1*t2 - 
             9*s2*t1*t2 + 10*Power(t1,2)*t2 - 13*Power(t2,2) + 
             10*s2*Power(t2,2) - 23*t1*Power(t2,2) + 12*Power(t2,3) + 
             Power(s,2)*(-2 - 2*s1 - 2*s2 + t1 + 3*t2) + 
             s*(-14 + t1 + 14*Power(t1,2) - t2 - 27*t1*t2 + 
                15*Power(t2,2) - 3*s1*(-5 + 4*s2 - 4*t1 + 5*t2) + 
                s2*(15 - 14*t1 + 11*t2))) + 
          16*(61 + 18*Power(s,3) - 45*s2 - 4*Power(s2,2) - 34*t1 + 
             66*s2*t1 + 11*Power(s2,2)*t1 - 63*Power(t1,2) - 
             21*s2*Power(t1,2) + 10*Power(t1,3) - 
             2*Power(s,2)*(12 + 10*s1 + 10*s2 - 5*t1 - 12*t2) - 48*t2 - 
             2*s2*t2 - 9*Power(s2,2)*t2 + 142*t1*t2 - 22*s2*t1*t2 + 
             32*Power(t1,2)*t2 - 65*Power(t2,2) + 34*s2*Power(t2,2) - 
             94*t1*Power(t2,2) + 52*Power(t2,3) + 
             Power(s1,2)*(-3 + 2*s2 - 2*t1 + 4*t2) + 
             s*(-66 + 2*Power(s1,2) + 2*Power(s2,2) + 5*t1 + 
                70*Power(t1,2) + s1*(84 - 50*s2 + 52*t1 - 80*t2) - 
                13*t2 - 131*t1*t2 + 77*Power(t2,2) + 
                s2*(84 - 72*t1 + 44*t2)) + 
             s1*(-48 + 3*Power(s2,2) + s2*(30 - 3*t1 - 41*t2) + 
                102*t2 - 55*Power(t2,2) + t1*(-36 + 46*t2))) + 
          4*(51 + 40*Power(s,3) - 3*Power(s1,3) - 98*s2 - 
             19*Power(s2,2) - Power(s2,3) - 40*t1 + 127*s2*t1 + 
             24*Power(s2,2)*t1 - 109*Power(t1,2) - 45*s2*Power(t1,2) + 
             22*Power(t1,3) - 67*t2 + 7*s2*t2 - 18*Power(s2,2)*t2 + 
             220*t1*t2 - 5*s2*t1*t2 + 23*Power(t1,2)*t2 - 
             96*Power(t2,2) + 28*s2*Power(t2,2) - 111*t1*Power(t2,2) + 
             66*Power(t2,3) + Power(s1,2)*(-9 + 6*s2 - 12*t1 + 20*t2) + 
             Power(s,2)*(-102 - 45*s1 - 41*s2 + 20*t1 + 54*t2) + 
             s*(-10 + 8*Power(s1,2) + 2*Power(s2,2) - 33*t1 + 
                100*Power(t1,2) + s1*(170 - 64*s2 + 72*t1 - 135*t2) - 
                6*s2*(-30 + 17*t1 - 7*t2) - 44*t2 - 182*t1*t2 + 
                116*Power(t2,2)) + 
             s1*(-69 + 11*Power(s2,2) + 2*Power(t1,2) + 126*t2 - 
                79*Power(t2,2) - s2*(17 + 16*t1 + 36*t2) + 
                t1*(-3 + 56*t2))) - 
          8*(114 + 52*Power(s,3) - 2*Power(s1,3) - 81*s2 - 
             17*Power(s2,2) - 3*Power(s2,3) - 93*t1 + 148*s2*t1 + 
             34*Power(s2,2)*t1 - 134*Power(t1,2) - 58*s2*Power(t1,2) + 
             27*Power(t1,3) - 87*t2 - 3*s2*t2 - 21*Power(s2,2)*t2 + 
             301*t1*t2 - 21*s2*t1*t2 + 44*Power(t1,2)*t2 - 
             136*Power(t2,2) + 51*s2*Power(t2,2) - 168*t1*Power(t2,2) + 
             97*Power(t2,3) + Power(s1,2)*(-8 + 3*s2 - 7*t1 + 18*t2) + 
             Power(s,2)*(-93 - 62*s1 - 63*s2 + 32*t1 + 64*t2) + 
             s*(-103 + 12*Power(s1,2) + 14*Power(s2,2) - 7*t1 + 
                144*Power(t1,2) + s1*(201 - 78*s2 + 91*t1 - 176*t2) - 
                43*t2 - 259*t1*t2 + 159*Power(t2,2) + 
                s2*(208 - 157*t1 + 70*t2)) + 
             s1*(-88 + 5*Power(s2,2) - Power(t1,2) + 
                s2*(17 - 6*t1 - 63*t2) + 186*t2 - 109*Power(t2,2) + 
                t1*(-38 + 83*t2))) + 
          (2*(-8*s + 5*s1 + 7*s2 - 4*t1)*
             (6*Power(s,4)*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) + 
               (-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
               2*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - 
                  s2*t2 + s*(-2 - 2*s1 - 2*s2 + t1 + t2),2) - 
               2*Power(s,3)*(6 - 3*t1 - 5*Power(t1,2) + 
                  3*Power(t1,3) + 6*Power(s2,2)*(-1 + t1 - t2) - 
                  3*t2 + 10*t1*t2 - 3*Power(t1,2)*t2 - 
                  5*Power(t2,2) - 3*t1*Power(t2,2) + 3*Power(t2,3) + 
                  6*Power(s1,2)*(-1 + s2 - t1 + t2) - 
                  3*s2*(-3*t1 + 3*Power(t1,2) + t2 - 4*t1*t2 + 
                     Power(t2,2)) + 
                  3*s1*(2*Power(s2,2) - Power(t1,2) - 
                     3*(-1 + t2)*t2 - s2*(2 + t1 + t2) + 
                     t1*(-1 + 4*t2))) + 
               Power(s,2)*(-4*t1 + 5*Power(t1,2) - Power(t1,4) + 
                  6*Power(s2,3)*(-1 + t1 - t2) - 4*t2 - 4*t1*t2 - 
                  2*Power(t1,3)*t2 + 5*Power(t2,2) + 
                  6*Power(t1,2)*Power(t2,2) - 2*t1*Power(t2,3) - 
                  Power(t2,4) + 6*Power(s1,3)*(-1 + s2 - t1 + t2) + 
                  6*Power(s2,2)*
                   (-1 - 2*Power(t1,2) + Power(t2,2) + t1*(3 + t2)) + 
                  6*Power(s1,2)*
                   (-1 + 3*Power(s2,2) + Power(t1,2) + 3*t2 + 
                     t1*t2 - 2*Power(t2,2) + s2*(-2 - 4*t1 + t2)) + 
                  s2*(12 + 7*Power(t1,3) - 4*t2 - 9*Power(t2,2) + 
                     11*Power(t2,3) + Power(t1,2)*(-13 + 3*t2) + 
                     t1*(-6 + 16*t2 - 21*Power(t2,2))) + 
                  s1*(12 + 6*Power(s2,3) + 11*Power(t1,3) + 
                     6*Power(s2,2)*(-2 + t1 - 4*t2) - 6*t2 - 
                     13*Power(t2,2) + 7*Power(t2,3) - 
                     3*Power(t1,2)*(3 + 7*t2) + 
                     t1*(-4 + 16*t2 + 3*Power(t2,2)) + 
                     s2*(-6 - 23*Power(t1,2) + 14*t2 - 
                        23*Power(t2,2) + 2*t1*(7 + 26*t2)))) - 
               2*s*(-2 + 3*t1 + 2*Power(t1,2) - 5*Power(t1,3) + 
                  2*Power(t1,4) + 3*t2 - 6*t1*t2 + 5*Power(t1,2)*t2 - 
                  2*Power(t1,3)*t2 + 2*Power(t2,2) + 
                  5*t1*Power(t2,2) - 5*Power(t2,3) - 
                  2*t1*Power(t2,3) + 2*Power(t2,4) + 
                  3*Power(s2,3)*t2*(1 - t1 + t2) + 
                  3*Power(s1,3)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                  Power(s1,2)*
                   (3*Power(s2,3) + 2*Power(t1,3) + 
                     Power(t1,2)*(1 - 7*t2) + 3*Power(-1 + t2,2) - 
                     2*Power(s2,2)*(1 + 2*t1 + 2*t2) + 
                     s2*(-3 + t1 - Power(t1,2) + 10*t2 + 11*t1*t2 - 
                        7*Power(t2,2)) + t1*(1 - 6*t2 + 5*Power(t2,2))\
) + Power(s2,2)*(3 + t2 + Power(t2,2) + 2*Power(t2,3) + 
                     Power(t1,2)*(3 + 5*t2) - 
                     t1*(6 + 6*t2 + 7*Power(t2,2))) - 
                  s2*(Power(t1,3)*(5 + 2*t2) - 
                     Power(t1,2)*(10 + 7*t2 + 3*Power(t2,2)) + 
                     t1*(5 + 4*t2 + 5*Power(t2,2)) + 
                     t2*(1 + t2 - 3*Power(t2,2) + Power(t2,3))) + 
                  s1*(3*Power(t1,3) - Power(t1,4) + 
                     3*Power(s2,3)*(-1 + t1 - 2*t2) - 
                     5*Power(-1 + t2,2)*t2 + 
                     Power(t1,2)*(-1 - 5*t2 + 3*Power(t2,2)) - 
                     t1*(1 + 4*t2 - 7*Power(t2,2) + 2*Power(t2,3)) + 
                     Power(s2,2)*
                      (-3 - 7*Power(t1,2) + t2 - Power(t2,2) + 
                        t1*(10 + 11*t2)) + 
                     s2*(4 + 5*Power(t1,3) + t2 - 10*Power(t2,2) + 
                        5*Power(t2,3) - 5*Power(t1,2)*(2 + t2) + 
                        t1*(1 + 10*t2 - 5*Power(t2,2)))))))/
           (Power(s,2)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)) + 
          ((-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
               s*(-2 - 2*s1 - 2*s2 + t1 + t2))*
             (-4*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - 
                  s2*t2 + s*(-2 - 2*s1 - 2*s2 + t1 + t2),2) + 
               12*s*(-Power(s2,2) + 2*s2*t1 + 2*Power(s2,2)*t1 - 
                  Power(t1,2) - 4*s2*Power(t1,2) - 
                  Power(s2,2)*Power(t1,2) + 2*Power(t1,3) + 
                  2*s2*Power(t1,3) - Power(t1,4) - Power(s2,2)*t2 - 
                  Power(s2,3)*t2 + 3*s2*t1*t2 + 3*Power(s2,2)*t1*t2 + 
                  Power(s2,3)*t1*t2 - 2*Power(t1,2)*t2 - 
                  4*s2*Power(t1,2)*t2 - 2*Power(s2,2)*Power(t1,2)*t2 + 
                  2*Power(t1,3)*t2 + s2*Power(t1,3)*t2 - Power(t2,2) + 
                  2*s2*Power(t2,2) - 2*Power(s2,2)*Power(t2,2) - 
                  Power(s2,3)*Power(t2,2) - 2*t1*Power(t2,2) + 
                  4*s2*t1*Power(t2,2) + 3*Power(s2,2)*t1*Power(t2,2) - 
                  2*Power(t1,2)*Power(t2,2) - 
                  2*s2*Power(t1,2)*Power(t2,2) + 2*Power(t2,3) - 
                  2*s2*Power(t2,3) - Power(s2,2)*Power(t2,3) + 
                  2*t1*Power(t2,3) + s2*t1*Power(t2,3) - Power(t2,4) - 
                  Power(s1,3)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                  Power(s,3)*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) - 
                  Power(s1,2)*
                   (Power(s2,3) + Power(t1,3) + 
                     Power(t1,2)*(2 - 3*t2) + Power(-1 + t2,2) - 
                     Power(s2,2)*(t1 + t2) + 
                     t1*(1 - 3*t2 + 2*Power(t2,2)) - 
                     s2*(1 + Power(t1,2) + t1*(2 - 4*t2) - 3*t2 + 
                        2*Power(t2,2))) + 
                  Power(s,2)*
                   (-2 + t1 + Power(t1,2) - Power(t1,3) + t2 - 
                     2*t1*t2 + Power(t1,2)*t2 + Power(t2,2) + 
                     t1*Power(t2,2) - Power(t2,3) - 
                     2*Power(s1,2)*(-1 + s2 - t1 + t2) + 
                     Power(s2,2)*(2 - 2*t1 + 2*t2) + 
                     s2*(-3*t1 + 3*Power(t1,2) + t2 - 4*t1*t2 + 
                        Power(t2,2)) + 
                     s1*(-2*Power(s2,2) + t1 + Power(t1,2) - 4*t1*t2 + 
                        3*(-1 + t2)*t2 + s2*(2 + t1 + t2))) + 
                  s1*(Power(t1,3)*(-2 + t2) + 2*Power(-1 + t2,2)*t2 + 
                     Power(s2,3)*(1 - t1 + 2*t2) + 
                     Power(t1,2)*(2 + 4*t2 - 2*Power(t2,2)) + 
                     t1*t2*(3 - 4*t2 + Power(t2,2)) + 
                     Power(s2,2)*
                      (2*Power(t1,2) + Power(1 + t2,2) - 
                        t1*(3 + 4*t2)) - 
                     s2*(Power(t1,3) - Power(t1,2)*(4 + t2) + 
                        t1*(3 + 6*t2 - Power(t2,2)) + 
                        t2*(3 - 4*t2 + Power(t2,2)))) + 
                  s*(-2*t1 + 2*Power(t1,2) + 
                     Power(s2,3)*(-1 + t1 - t2) - 2*t2 - 3*t1*t2 - 
                     Power(t1,3)*t2 + 2*Power(t2,2) + 
                     2*Power(t1,2)*Power(t2,2) - t1*Power(t2,3) + 
                     Power(s1,3)*(-1 + s2 - t1 + t2) + 
                     Power(s2,2)*
                      (-1 - 2*Power(t1,2) + Power(t2,2) + t1*(3 + t2)) \
+ Power(s1,2)*(-1 + 3*Power(s2,2) + Power(t1,2) + 3*t2 + t1*t2 - 
                        2*Power(t2,2) + s2*(-2 - 4*t1 + t2)) + 
                     s1*(2 + Power(s2,3) + 2*Power(t1,3) + 
                        Power(s2,2)*(-2 + t1 - 4*t2) - t2 - 
                        4*Power(t1,2)*t2 - 2*Power(t2,2) + 
                        Power(t2,3) + 
                        s2*(-1 + t1 - 4*Power(t1,2) + t2 + 9*t1*t2 - 
                        4*Power(t2,2)) + t1*(2 + t2 + Power(t2,2))) + 
                     s2*(Power(t1,3) + Power(t1,2)*(-2 + t2) + 
                        t1*(-1 + t2 - 4*Power(t2,2)) + 
                        2*(1 + t2 + Power(t2,3)))))))/
           (Power(s,3)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)))/
        ((-1 + s1)*(-1 + t2)) + 
       ((2*(-23 + 6*Power(s,2) + 2*Power(s1,2) + 8*s2 + 
               6*Power(s2,2) - s*(1 + 8*s1 + 12*s2 - 7*t1) + 
               s1*(5 + 8*s2 - 5*t1) - 2*t1 - 8*s2*t1 + 
               3*Power(t1,2) + 7*t2 - t1*t2)*
             (-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
               s*(-2 - 2*s1 - 2*s2 + t1 + t2)))/s - 
          4*(33 + 2*Power(s,3) + 25*s2 - 12*Power(s2,2) - 
             2*Power(s2,3) - 13*t1 + 12*s2*t1 + 4*Power(s2,2)*t1 - 
             3*Power(t1,2) - 3*s2*Power(t1,2) + Power(t1,3) + 
             Power(s1,2)*(-2*s2 + t1) + 
             Power(s,2)*(3 - 4*s1 - 6*s2 + 2*t1) - 12*t2 - 12*s2*t2 + 
             5*t1*t2 + 2*s2*t1*t2 - Power(t1,2)*t2 + 3*Power(t2,2) + 
             s*(-33 + 2*Power(s1,2) + 9*s2 + 6*Power(s2,2) + 
                s1*(-6 + 8*s2 - 3*t1) - 6*t1 - 6*s2*t1 + 
                2*Power(t1,2) + 17*t2 - 2*t1*t2) + 
             s1*(18 - 4*Power(s2,2) + 4*s2*(-1 + t1) + t1 - 
                Power(t1,2) - 9*t2 + t1*t2)) - 
          (4*(-1 + 3*s - 2*s1 - 3*s2 + 2*t1)*
             (6*Power(s,4)*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) + 
               (-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
               2*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - 
                  s2*t2 + s*(-2 - 2*s1 - 2*s2 + t1 + t2),2) - 
               2*Power(s,3)*(6 - 3*t1 - 5*Power(t1,2) + 
                  3*Power(t1,3) + 6*Power(s2,2)*(-1 + t1 - t2) - 
                  3*t2 + 10*t1*t2 - 3*Power(t1,2)*t2 - 
                  5*Power(t2,2) - 3*t1*Power(t2,2) + 3*Power(t2,3) + 
                  6*Power(s1,2)*(-1 + s2 - t1 + t2) - 
                  3*s2*(-3*t1 + 3*Power(t1,2) + t2 - 4*t1*t2 + 
                     Power(t2,2)) + 
                  3*s1*(2*Power(s2,2) - Power(t1,2) - 
                     3*(-1 + t2)*t2 - s2*(2 + t1 + t2) + 
                     t1*(-1 + 4*t2))) + 
               Power(s,2)*(-4*t1 + 5*Power(t1,2) - Power(t1,4) + 
                  6*Power(s2,3)*(-1 + t1 - t2) - 4*t2 - 4*t1*t2 - 
                  2*Power(t1,3)*t2 + 5*Power(t2,2) + 
                  6*Power(t1,2)*Power(t2,2) - 2*t1*Power(t2,3) - 
                  Power(t2,4) + 6*Power(s1,3)*(-1 + s2 - t1 + t2) + 
                  6*Power(s2,2)*
                   (-1 - 2*Power(t1,2) + Power(t2,2) + t1*(3 + t2)) + 
                  6*Power(s1,2)*
                   (-1 + 3*Power(s2,2) + Power(t1,2) + 3*t2 + 
                     t1*t2 - 2*Power(t2,2) + s2*(-2 - 4*t1 + t2)) + 
                  s2*(12 + 7*Power(t1,3) - 4*t2 - 9*Power(t2,2) + 
                     11*Power(t2,3) + Power(t1,2)*(-13 + 3*t2) + 
                     t1*(-6 + 16*t2 - 21*Power(t2,2))) + 
                  s1*(12 + 6*Power(s2,3) + 11*Power(t1,3) + 
                     6*Power(s2,2)*(-2 + t1 - 4*t2) - 6*t2 - 
                     13*Power(t2,2) + 7*Power(t2,3) - 
                     3*Power(t1,2)*(3 + 7*t2) + 
                     t1*(-4 + 16*t2 + 3*Power(t2,2)) + 
                     s2*(-6 - 23*Power(t1,2) + 14*t2 - 
                        23*Power(t2,2) + 2*t1*(7 + 26*t2)))) - 
               2*s*(-2 + 3*t1 + 2*Power(t1,2) - 5*Power(t1,3) + 
                  2*Power(t1,4) + 3*t2 - 6*t1*t2 + 5*Power(t1,2)*t2 - 
                  2*Power(t1,3)*t2 + 2*Power(t2,2) + 
                  5*t1*Power(t2,2) - 5*Power(t2,3) - 
                  2*t1*Power(t2,3) + 2*Power(t2,4) + 
                  3*Power(s2,3)*t2*(1 - t1 + t2) + 
                  3*Power(s1,3)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                  Power(s1,2)*
                   (3*Power(s2,3) + 2*Power(t1,3) + 
                     Power(t1,2)*(1 - 7*t2) + 3*Power(-1 + t2,2) - 
                     2*Power(s2,2)*(1 + 2*t1 + 2*t2) + 
                     s2*(-3 + t1 - Power(t1,2) + 10*t2 + 11*t1*t2 - 
                        7*Power(t2,2)) + t1*(1 - 6*t2 + 5*Power(t2,2))\
) + Power(s2,2)*(3 + t2 + Power(t2,2) + 2*Power(t2,3) + 
                     Power(t1,2)*(3 + 5*t2) - 
                     t1*(6 + 6*t2 + 7*Power(t2,2))) - 
                  s2*(Power(t1,3)*(5 + 2*t2) - 
                     Power(t1,2)*(10 + 7*t2 + 3*Power(t2,2)) + 
                     t1*(5 + 4*t2 + 5*Power(t2,2)) + 
                     t2*(1 + t2 - 3*Power(t2,2) + Power(t2,3))) + 
                  s1*(3*Power(t1,3) - Power(t1,4) + 
                     3*Power(s2,3)*(-1 + t1 - 2*t2) - 
                     5*Power(-1 + t2,2)*t2 + 
                     Power(t1,2)*(-1 - 5*t2 + 3*Power(t2,2)) - 
                     t1*(1 + 4*t2 - 7*Power(t2,2) + 2*Power(t2,3)) + 
                     Power(s2,2)*
                      (-3 - 7*Power(t1,2) + t2 - Power(t2,2) + 
                        t1*(10 + 11*t2)) + 
                     s2*(4 + 5*Power(t1,3) + t2 - 10*Power(t2,2) + 
                        5*Power(t2,3) - 5*Power(t1,2)*(2 + t2) + 
                        t1*(1 + 10*t2 - 5*Power(t2,2)))))))/
           (Power(s,2)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)) + 
          ((-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
               s*(-2 - 2*s1 - 2*s2 + t1 + t2))*
             (-4*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - 
                  s2*t2 + s*(-2 - 2*s1 - 2*s2 + t1 + t2),2) + 
               12*s*(-Power(s2,2) + 2*s2*t1 + 2*Power(s2,2)*t1 - 
                  Power(t1,2) - 4*s2*Power(t1,2) - 
                  Power(s2,2)*Power(t1,2) + 2*Power(t1,3) + 
                  2*s2*Power(t1,3) - Power(t1,4) - Power(s2,2)*t2 - 
                  Power(s2,3)*t2 + 3*s2*t1*t2 + 3*Power(s2,2)*t1*t2 + 
                  Power(s2,3)*t1*t2 - 2*Power(t1,2)*t2 - 
                  4*s2*Power(t1,2)*t2 - 2*Power(s2,2)*Power(t1,2)*t2 + 
                  2*Power(t1,3)*t2 + s2*Power(t1,3)*t2 - Power(t2,2) + 
                  2*s2*Power(t2,2) - 2*Power(s2,2)*Power(t2,2) - 
                  Power(s2,3)*Power(t2,2) - 2*t1*Power(t2,2) + 
                  4*s2*t1*Power(t2,2) + 3*Power(s2,2)*t1*Power(t2,2) - 
                  2*Power(t1,2)*Power(t2,2) - 
                  2*s2*Power(t1,2)*Power(t2,2) + 2*Power(t2,3) - 
                  2*s2*Power(t2,3) - Power(s2,2)*Power(t2,3) + 
                  2*t1*Power(t2,3) + s2*t1*Power(t2,3) - Power(t2,4) - 
                  Power(s1,3)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                  Power(s,3)*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) - 
                  Power(s1,2)*
                   (Power(s2,3) + Power(t1,3) + 
                     Power(t1,2)*(2 - 3*t2) + Power(-1 + t2,2) - 
                     Power(s2,2)*(t1 + t2) + 
                     t1*(1 - 3*t2 + 2*Power(t2,2)) - 
                     s2*(1 + Power(t1,2) + t1*(2 - 4*t2) - 3*t2 + 
                        2*Power(t2,2))) + 
                  Power(s,2)*
                   (-2 + t1 + Power(t1,2) - Power(t1,3) + t2 - 
                     2*t1*t2 + Power(t1,2)*t2 + Power(t2,2) + 
                     t1*Power(t2,2) - Power(t2,3) - 
                     2*Power(s1,2)*(-1 + s2 - t1 + t2) + 
                     Power(s2,2)*(2 - 2*t1 + 2*t2) + 
                     s2*(-3*t1 + 3*Power(t1,2) + t2 - 4*t1*t2 + 
                        Power(t2,2)) + 
                     s1*(-2*Power(s2,2) + t1 + Power(t1,2) - 4*t1*t2 + 
                        3*(-1 + t2)*t2 + s2*(2 + t1 + t2))) + 
                  s1*(Power(t1,3)*(-2 + t2) + 2*Power(-1 + t2,2)*t2 + 
                     Power(s2,3)*(1 - t1 + 2*t2) + 
                     Power(t1,2)*(2 + 4*t2 - 2*Power(t2,2)) + 
                     t1*t2*(3 - 4*t2 + Power(t2,2)) + 
                     Power(s2,2)*
                      (2*Power(t1,2) + Power(1 + t2,2) - 
                        t1*(3 + 4*t2)) - 
                     s2*(Power(t1,3) - Power(t1,2)*(4 + t2) + 
                        t1*(3 + 6*t2 - Power(t2,2)) + 
                        t2*(3 - 4*t2 + Power(t2,2)))) + 
                  s*(-2*t1 + 2*Power(t1,2) + 
                     Power(s2,3)*(-1 + t1 - t2) - 2*t2 - 3*t1*t2 - 
                     Power(t1,3)*t2 + 2*Power(t2,2) + 
                     2*Power(t1,2)*Power(t2,2) - t1*Power(t2,3) + 
                     Power(s1,3)*(-1 + s2 - t1 + t2) + 
                     Power(s2,2)*
                      (-1 - 2*Power(t1,2) + Power(t2,2) + t1*(3 + t2)) \
+ Power(s1,2)*(-1 + 3*Power(s2,2) + Power(t1,2) + 3*t2 + t1*t2 - 
                        2*Power(t2,2) + s2*(-2 - 4*t1 + t2)) + 
                     s1*(2 + Power(s2,3) + 2*Power(t1,3) + 
                        Power(s2,2)*(-2 + t1 - 4*t2) - t2 - 
                        4*Power(t1,2)*t2 - 2*Power(t2,2) + 
                        Power(t2,3) + 
                        s2*(-1 + t1 - 4*Power(t1,2) + t2 + 9*t1*t2 - 
                        4*Power(t2,2)) + t1*(2 + t2 + Power(t2,2))) + 
                     s2*(Power(t1,3) + Power(t1,2)*(-2 + t2) + 
                        t1*(-1 + t2 - 4*Power(t2,2)) + 
                        2*(1 + t2 + Power(t2,3)))))))/
           (Power(s,3)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)))/
        ((s - s2 + t1)*(s - s1 + t2)) + 
       (64*(-1 + s1 + t1 - t2)*(s - s1 + t2)*(-1 + s2 - t1 + t2) + 
          (2*(-17 + 6*Power(s,2) + 2*Power(s1,2) + 6*s2 + 
               6*Power(s2,2) - s*(3 + 8*s1 + 12*s2 - 7*t1) + 
               s1*(3 + 8*s2 - 5*t1) - 8*s2*t1 + 3*Power(t1,2) + 
               3*t2 - t1*t2)*
             (-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
               s*(-2 - 2*s1 - 2*s2 + t1 + t2)))/s + 
          32*(2*Power(s,3) + 2*s2 - Power(s2,2) - Power(s2,3) - 2*t1 + 
             2*s2*t1 + 3*Power(s2,2)*t1 - Power(t1,2) - 
             3*s2*Power(t1,2) + Power(t1,3) + 
             Power(s1,2)*(11*s2 - 11*t1 + 12*(-1 + t2)) - 11*t2 + 
             9*s2*t2 - 2*Power(s2,2)*t2 + 4*t1*t2 - 8*s2*t1*t2 + 
             10*Power(t1,2)*t2 - Power(t2,2) + 10*s2*Power(t2,2) - 
             23*t1*Power(t2,2) + 12*Power(t2,3) + 
             Power(s,2)*(-4*s1 - 2*s2 + t1 + 3*t2) + 
             s*(-14 + 2*Power(s1,2) + Power(s2,2) + 2*t1 + 
                14*Power(t1,2) + s1*(14 - 11*s2 + 12*t1 - 17*t2) - 
                3*s2*(-4 + 5*t1 - 4*t2) - 27*t1*t2 + 15*Power(t2,2)) + 
             s1*(11 + Power(s2,2) - 11*Power(t1,2) + 
                s2*(-11 + 10*t1 - 21*t2) + 13*t2 - 24*Power(t2,2) + 
                t1*(-2 + 34*t2))) - 
          16*(16*Power(s,3) - 3*Power(s1,3) + 18*s2 - 14*Power(s2,2) - 
             9*Power(s2,3) - 20*t1 + 32*s2*t1 + 27*Power(s2,2)*t1 - 
             18*Power(t1,2) - 28*s2*Power(t1,2) + 10*Power(t1,3) - 
             37*t2 + 11*s2*t2 - 17*Power(s2,2)*t2 + 50*t1*t2 - 
             14*s2*t1*t2 + 31*Power(t1,2)*t2 - 20*Power(t2,2) + 
             34*s2*Power(t2,2) - 92*t1*Power(t2,2) + 51*Power(t2,3) + 
             Power(s,2)*(-4 - 31*s1 - 17*s2 + 8*t1 + 22*t2) + 
             Power(s1,2)*(-49 + 40*s2 - 45*t1 + 55*t2) + 
             s1*(35 + 7*Power(s2,2) - 41*Power(t1,2) + 
                s2*(-29 + 32*t1 - 75*t2) + 69*t2 - 103*Power(t2,2) + 
                8*t1*(-3 + 17*t2)) + 
             s*(-63 + 18*Power(s1,2) + 10*Power(s2,2) + 14*t1 + 
                67*Power(t1,2) + s1*(65 - 41*s2 + 53*t1 - 91*t2) - 
                4*t2 - 127*t1*t2 + 74*Power(t2,2) + 
                s2*(49 - 76*t1 + 52*t2))) + 
          8*(26 + 38*Power(s,3) - 10*Power(s1,3) + 48*s2 - 
             41*Power(s2,2) - 22*Power(s2,3) - 70*t1 + 101*s2*t1 + 
             68*Power(s2,2)*t1 - 58*Power(t1,2) - 72*s2*Power(t1,2) + 
             26*Power(t1,3) - 56*t2 - 25*s2*t2 - 43*Power(s2,2)*t2 + 
             157*t1*t2 + 3*s2*t1*t2 + 39*Power(t1,2)*t2 - 
             72*Power(t2,2) + 48*s2*Power(t2,2) - 157*t1*Power(t2,2) + 
             92*Power(t2,3) + 
             Power(s,2)*(-17 - 74*s1 - 39*s2 + 18*t1 + 54*t2) + 
             Power(s1,2)*(-92 + 60*s2 - 75*t1 + 103*t2) + 
             s*(-125 + 46*Power(s1,2) + 23*Power(s2,2) + 38*t1 + 
                128*Power(t1,2) + s1*(131 - 67*s2 + 97*t1 - 191*t2) - 
                24*t2 - 241*t1*t2 + 147*Power(t2,2) + 
                s2*(84 - 150*t1 + 92*t2)) + 
             s1*(38 + 18*Power(s2,2) - 62*Power(t1,2) + 
                s2*(-13 + 39*t1 - 111*t2) + 168*t2 - 185*Power(t2,2) + 
                t1*(-79 + 225*t2))) - 
          4*(74 + 30*Power(s,3) - 8*Power(s1,3) + 78*s2 - 
             40*Power(s2,2) - 18*Power(s2,3) - 108*t1 + 91*s2*t1 + 
             56*Power(s2,2)*t1 - 49*Power(t1,2) - 59*s2*Power(t1,2) + 
             21*Power(t1,3) - 
             2*Power(s,2)*(8 + 30*s1 + 16*s2 - 7*t1 - 22*t2) - 28*t2 - 
             44*s2*t2 - 34*Power(s2,2)*t2 + 146*t1*t2 + 16*s2*t1*t2 + 
             17*Power(t1,2)*t2 - 70*Power(t2,2) + 24*s2*Power(t2,2) - 
             98*t1*Power(t2,2) + 60*Power(t2,3) + 
             Power(s1,2)*(-64 + 30*s2 - 41*t1 + 66*t2) + 
             s*(-125 + 38*Power(s1,2) + 20*Power(s2,2) + 37*t1 + 
                86*Power(t1,2) + s1*(95 - 34*s2 + 59*t1 - 138*t2) - 
                24*t2 - 160*t1*t2 + 102*Power(t2,2) + 
                s2*(58 - 106*t1 + 56*t2)) + 
             s1*(12*Power(s2,2) - 33*Power(t1,2) + 
                2*s2*(7 + 9*t1 - 29*t2) + t1*(-74 + 131*t2) + 
                2*(7 + 70*t2 - 59*Power(t2,2)))) - 
          (2*(-2 + 6*s - 4*s1 - 6*s2 + 4*t1)*
             (6*Power(s,4)*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) + 
               (-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
               2*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - 
                  s2*t2 + s*(-2 - 2*s1 - 2*s2 + t1 + t2),2) - 
               2*Power(s,3)*(6 - 3*t1 - 5*Power(t1,2) + 
                  3*Power(t1,3) + 6*Power(s2,2)*(-1 + t1 - t2) - 
                  3*t2 + 10*t1*t2 - 3*Power(t1,2)*t2 - 
                  5*Power(t2,2) - 3*t1*Power(t2,2) + 3*Power(t2,3) + 
                  6*Power(s1,2)*(-1 + s2 - t1 + t2) - 
                  3*s2*(-3*t1 + 3*Power(t1,2) + t2 - 4*t1*t2 + 
                     Power(t2,2)) + 
                  3*s1*(2*Power(s2,2) - Power(t1,2) - 
                     3*(-1 + t2)*t2 - s2*(2 + t1 + t2) + 
                     t1*(-1 + 4*t2))) + 
               Power(s,2)*(-4*t1 + 5*Power(t1,2) - Power(t1,4) + 
                  6*Power(s2,3)*(-1 + t1 - t2) - 4*t2 - 4*t1*t2 - 
                  2*Power(t1,3)*t2 + 5*Power(t2,2) + 
                  6*Power(t1,2)*Power(t2,2) - 2*t1*Power(t2,3) - 
                  Power(t2,4) + 6*Power(s1,3)*(-1 + s2 - t1 + t2) + 
                  6*Power(s2,2)*
                   (-1 - 2*Power(t1,2) + Power(t2,2) + t1*(3 + t2)) + 
                  6*Power(s1,2)*
                   (-1 + 3*Power(s2,2) + Power(t1,2) + 3*t2 + 
                     t1*t2 - 2*Power(t2,2) + s2*(-2 - 4*t1 + t2)) + 
                  s2*(12 + 7*Power(t1,3) - 4*t2 - 9*Power(t2,2) + 
                     11*Power(t2,3) + Power(t1,2)*(-13 + 3*t2) + 
                     t1*(-6 + 16*t2 - 21*Power(t2,2))) + 
                  s1*(12 + 6*Power(s2,3) + 11*Power(t1,3) + 
                     6*Power(s2,2)*(-2 + t1 - 4*t2) - 6*t2 - 
                     13*Power(t2,2) + 7*Power(t2,3) - 
                     3*Power(t1,2)*(3 + 7*t2) + 
                     t1*(-4 + 16*t2 + 3*Power(t2,2)) + 
                     s2*(-6 - 23*Power(t1,2) + 14*t2 - 
                        23*Power(t2,2) + 2*t1*(7 + 26*t2)))) - 
               2*s*(-2 + 3*t1 + 2*Power(t1,2) - 5*Power(t1,3) + 
                  2*Power(t1,4) + 3*t2 - 6*t1*t2 + 5*Power(t1,2)*t2 - 
                  2*Power(t1,3)*t2 + 2*Power(t2,2) + 
                  5*t1*Power(t2,2) - 5*Power(t2,3) - 
                  2*t1*Power(t2,3) + 2*Power(t2,4) + 
                  3*Power(s2,3)*t2*(1 - t1 + t2) + 
                  3*Power(s1,3)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                  Power(s1,2)*
                   (3*Power(s2,3) + 2*Power(t1,3) + 
                     Power(t1,2)*(1 - 7*t2) + 3*Power(-1 + t2,2) - 
                     2*Power(s2,2)*(1 + 2*t1 + 2*t2) + 
                     s2*(-3 + t1 - Power(t1,2) + 10*t2 + 11*t1*t2 - 
                        7*Power(t2,2)) + t1*(1 - 6*t2 + 5*Power(t2,2))\
) + Power(s2,2)*(3 + t2 + Power(t2,2) + 2*Power(t2,3) + 
                     Power(t1,2)*(3 + 5*t2) - 
                     t1*(6 + 6*t2 + 7*Power(t2,2))) - 
                  s2*(Power(t1,3)*(5 + 2*t2) - 
                     Power(t1,2)*(10 + 7*t2 + 3*Power(t2,2)) + 
                     t1*(5 + 4*t2 + 5*Power(t2,2)) + 
                     t2*(1 + t2 - 3*Power(t2,2) + Power(t2,3))) + 
                  s1*(3*Power(t1,3) - Power(t1,4) + 
                     3*Power(s2,3)*(-1 + t1 - 2*t2) - 
                     5*Power(-1 + t2,2)*t2 + 
                     Power(t1,2)*(-1 - 5*t2 + 3*Power(t2,2)) - 
                     t1*(1 + 4*t2 - 7*Power(t2,2) + 2*Power(t2,3)) + 
                     Power(s2,2)*
                      (-3 - 7*Power(t1,2) + t2 - Power(t2,2) + 
                        t1*(10 + 11*t2)) + 
                     s2*(4 + 5*Power(t1,3) + t2 - 10*Power(t2,2) + 
                        5*Power(t2,3) - 5*Power(t1,2)*(2 + t2) + 
                        t1*(1 + 10*t2 - 5*Power(t2,2)))))))/
           (Power(s,2)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)) + 
          ((-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
               s*(-2 - 2*s1 - 2*s2 + t1 + t2))*
             (-4*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - 
                  s2*t2 + s*(-2 - 2*s1 - 2*s2 + t1 + t2),2) + 
               12*s*(-Power(s2,2) + 2*s2*t1 + 2*Power(s2,2)*t1 - 
                  Power(t1,2) - 4*s2*Power(t1,2) - 
                  Power(s2,2)*Power(t1,2) + 2*Power(t1,3) + 
                  2*s2*Power(t1,3) - Power(t1,4) - Power(s2,2)*t2 - 
                  Power(s2,3)*t2 + 3*s2*t1*t2 + 3*Power(s2,2)*t1*t2 + 
                  Power(s2,3)*t1*t2 - 2*Power(t1,2)*t2 - 
                  4*s2*Power(t1,2)*t2 - 2*Power(s2,2)*Power(t1,2)*t2 + 
                  2*Power(t1,3)*t2 + s2*Power(t1,3)*t2 - Power(t2,2) + 
                  2*s2*Power(t2,2) - 2*Power(s2,2)*Power(t2,2) - 
                  Power(s2,3)*Power(t2,2) - 2*t1*Power(t2,2) + 
                  4*s2*t1*Power(t2,2) + 3*Power(s2,2)*t1*Power(t2,2) - 
                  2*Power(t1,2)*Power(t2,2) - 
                  2*s2*Power(t1,2)*Power(t2,2) + 2*Power(t2,3) - 
                  2*s2*Power(t2,3) - Power(s2,2)*Power(t2,3) + 
                  2*t1*Power(t2,3) + s2*t1*Power(t2,3) - Power(t2,4) - 
                  Power(s1,3)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                  Power(s,3)*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) - 
                  Power(s1,2)*
                   (Power(s2,3) + Power(t1,3) + 
                     Power(t1,2)*(2 - 3*t2) + Power(-1 + t2,2) - 
                     Power(s2,2)*(t1 + t2) + 
                     t1*(1 - 3*t2 + 2*Power(t2,2)) - 
                     s2*(1 + Power(t1,2) + t1*(2 - 4*t2) - 3*t2 + 
                        2*Power(t2,2))) + 
                  Power(s,2)*
                   (-2 + t1 + Power(t1,2) - Power(t1,3) + t2 - 
                     2*t1*t2 + Power(t1,2)*t2 + Power(t2,2) + 
                     t1*Power(t2,2) - Power(t2,3) - 
                     2*Power(s1,2)*(-1 + s2 - t1 + t2) + 
                     Power(s2,2)*(2 - 2*t1 + 2*t2) + 
                     s2*(-3*t1 + 3*Power(t1,2) + t2 - 4*t1*t2 + 
                        Power(t2,2)) + 
                     s1*(-2*Power(s2,2) + t1 + Power(t1,2) - 4*t1*t2 + 
                        3*(-1 + t2)*t2 + s2*(2 + t1 + t2))) + 
                  s1*(Power(t1,3)*(-2 + t2) + 2*Power(-1 + t2,2)*t2 + 
                     Power(s2,3)*(1 - t1 + 2*t2) + 
                     Power(t1,2)*(2 + 4*t2 - 2*Power(t2,2)) + 
                     t1*t2*(3 - 4*t2 + Power(t2,2)) + 
                     Power(s2,2)*
                      (2*Power(t1,2) + Power(1 + t2,2) - 
                        t1*(3 + 4*t2)) - 
                     s2*(Power(t1,3) - Power(t1,2)*(4 + t2) + 
                        t1*(3 + 6*t2 - Power(t2,2)) + 
                        t2*(3 - 4*t2 + Power(t2,2)))) + 
                  s*(-2*t1 + 2*Power(t1,2) + 
                     Power(s2,3)*(-1 + t1 - t2) - 2*t2 - 3*t1*t2 - 
                     Power(t1,3)*t2 + 2*Power(t2,2) + 
                     2*Power(t1,2)*Power(t2,2) - t1*Power(t2,3) + 
                     Power(s1,3)*(-1 + s2 - t1 + t2) + 
                     Power(s2,2)*
                      (-1 - 2*Power(t1,2) + Power(t2,2) + t1*(3 + t2)) \
+ Power(s1,2)*(-1 + 3*Power(s2,2) + Power(t1,2) + 3*t2 + t1*t2 - 
                        2*Power(t2,2) + s2*(-2 - 4*t1 + t2)) + 
                     s1*(2 + Power(s2,3) + 2*Power(t1,3) + 
                        Power(s2,2)*(-2 + t1 - 4*t2) - t2 - 
                        4*Power(t1,2)*t2 - 2*Power(t2,2) + 
                        Power(t2,3) + 
                        s2*(-1 + t1 - 4*Power(t1,2) + t2 + 9*t1*t2 - 
                        4*Power(t2,2)) + t1*(2 + t2 + Power(t2,2))) + 
                     s2*(Power(t1,3) + Power(t1,2)*(-2 + t2) + 
                        t1*(-1 + t2 - 4*Power(t2,2)) + 
                        2*(1 + t2 + Power(t2,3)))))))/
           (Power(s,3)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)))/
        ((-1 + s1)*(-s + s1 - t2)) + 
       (64*s*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) + 
          (2*(-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
               s*(-2 - 2*s1 - 2*s2 + t1 + t2))*
             (-16 + 12*Power(s,2) + 9*Power(s1,2) - 3*s2 + 
               9*Power(s2,2) + 2*t1 - 5*s2*t1 + 
               s1*(-4 + 19*s2 - 6*t1 - 5*t2) + 5*t2 - 7*s2*t2 + 
               3*t1*t2 - Power(t2,2) + 
               s*(5 - 21*s1 - 21*s2 + 6*t1 + 7*t2)))/s + 
          16*(-25 + Power(s,3) - Power(s1,3) + 22*s2 + 5*Power(s2,2) - 
             Power(s2,3) - 2*Power(s1,2)*(6*s2 - 5*t1) + 18*t1 - 
             34*s2*t1 - 6*Power(s2,2)*t1 + 18*Power(t1,2) + 
             16*s2*Power(t1,2) - 9*Power(t1,3) + 15*t2 - 11*s2*t2 + 
             11*Power(s2,2)*t2 - 21*t1*t2 - 12*s2*t1*t2 + 
             11*Power(t1,2)*t2 + 11*Power(t2,2) - 6*s2*Power(t2,2) - 
             t1*Power(t2,2) - Power(t2,3) + 
             Power(s,2)*(17 - 19*s1 - 19*s2 + t1 + 15*t2) + 
             s1*(13 - 12*Power(s2,2) - 2*Power(t1,2) - 12*t2 + 
                2*Power(t2,2) - 2*t1*(9 + 4*t2) + 
                s2*(9 + 13*t1 + 18*t2)) + 
             s*(59 + 19*Power(s1,2) + 19*Power(s2,2) + 5*t1 - 
                66*Power(t1,2) + s2*(-78 + 56*t1 - 85*t2) + 3*t2 + 
                135*t1*t2 - 54*Power(t2,2) + 
                s1*(-73 + 98*s2 - 70*t1 + 36*t2))) - 
          8*(-56 + 6*Power(s,3) - 5*Power(s1,3) + 57*s2 - 
             5*Power(s2,3) + 49*t1 - 79*s2*t1 - 13*Power(s2,2)*t1 + 
             47*Power(t1,2) + 41*s2*Power(t1,2) - 24*Power(t1,3) + 
             33*t2 - 31*s2*t2 + 40*Power(s2,2)*t2 - 63*t1*t2 - 
             42*s2*t1*t2 + 31*Power(t1,2)*t2 + 36*Power(t2,2) - 
             12*s2*Power(t2,2) - 7*Power(t2,3) - 
             Power(s1,2)*(5 + 45*s2 - 34*t1 + 2*t2) + 
             Power(s,2)*(38 - 61*s1 - 61*s2 + 6*t1 + 40*t2) + 
             s*(99 + 60*Power(s1,2) + 60*Power(s2,2) + t1 - 
                130*Power(t1,2) + s2*(-136 + 98*t1 - 193*t2) + 11*t2 + 
                278*t1*t2 - 108*Power(t2,2) + 
                s1*(-131 + 241*s2 - 153*t1 + 53*t2)) + 
             s1*(34 - 45*Power(s2,2) - 7*Power(t1,2) - 36*t2 + 
                14*Power(t2,2) - t1*(41 + 30*t2) + 
                s2*(14 + 45*t1 + 57*t2))) - 
          4*(49 + Power(s1,3) - 10*s2 + 14*Power(s2,2) + Power(s2,3) - 
             53*t1 + 48*s2*t1 + 14*Power(s2,2)*t1 - 38*Power(t1,2) - 
             34*s2*Power(t1,2) + 20*Power(t1,3) - 40*t2 + 19*s2*t2 - 
             36*Power(s2,2)*t2 + 65*t1*t2 + 38*s2*t1*t2 - 
             26*Power(t1,2)*t2 - 36*Power(t2,2) + 8*s2*Power(t2,2) - 
             3*t1*Power(t2,2) + 9*Power(t2,3) + 
             Power(s1,2)*(14 + 33*s2 - 30*t1 + 8*t2) + 
             Power(s,2)*(41*s1 + 41*s2 - 2*(10 + t1 + 14*t2)) + 
             s1*(5 + 33*Power(s2,2) + 6*Power(t1,2) + 
                s2*(8 - 36*t1 - 45*t2) + 27*t2 - 18*Power(t2,2) + 
                t1*(20 + 29*t2)) + 
             s*(-95 - 42*Power(s1,2) - 42*Power(s2,2) + 15*t1 + 
                92*Power(t1,2) + s1*(65 - 172*s2 + 110*t1 - 38*t2) - 
                2*t2 - 201*t1*t2 + 79*Power(t2,2) + 
                s2*(65 - 70*t1 + 142*t2))) + 
          32*(3 - 2*s2 - Power(s2,2) + Power(s1,2)*(s2 - t1) - 2*t1 + 
             4*s2*t1 + Power(s2,2)*t1 - 2*Power(t1,2) - 
             2*s2*Power(t1,2) + Power(t1,3) + 
             2*Power(s,2)*(-1 + s1 + s2 - t2) - 2*t2 + s2*t2 - 
             Power(s2,2)*t2 + 2*t1*t2 + s2*t1*t2 - Power(t1,2)*t2 - 
             Power(t2,2) + s2*Power(t2,2) + 
             s1*(-1 + Power(s2,2) + t2 + t1*(2 + t2) - 
                s2*(1 + t1 + 2*t2)) - 
             s*(14 + 2*Power(s1,2) + 2*Power(s2,2) + t1 - 
                14*Power(t1,2) + 28*t1*t2 - 12*Power(t2,2) + 
                s1*(-15 + 17*s2 - 14*t1 + 10*t2) + 
                s2*(13*t1 - 16*(1 + t2)))) - 
          (2*(2 + 8*s - 7*s1 - 7*s2 + 2*t1 + 2*t2)*
             (6*Power(s,4)*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) + 
               (-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
               2*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - 
                  s2*t2 + s*(-2 - 2*s1 - 2*s2 + t1 + t2),2) - 
               2*Power(s,3)*(6 - 3*t1 - 5*Power(t1,2) + 
                  3*Power(t1,3) + 6*Power(s2,2)*(-1 + t1 - t2) - 
                  3*t2 + 10*t1*t2 - 3*Power(t1,2)*t2 - 
                  5*Power(t2,2) - 3*t1*Power(t2,2) + 3*Power(t2,3) + 
                  6*Power(s1,2)*(-1 + s2 - t1 + t2) - 
                  3*s2*(-3*t1 + 3*Power(t1,2) + t2 - 4*t1*t2 + 
                     Power(t2,2)) + 
                  3*s1*(2*Power(s2,2) - Power(t1,2) - 
                     3*(-1 + t2)*t2 - s2*(2 + t1 + t2) + 
                     t1*(-1 + 4*t2))) + 
               Power(s,2)*(-4*t1 + 5*Power(t1,2) - Power(t1,4) + 
                  6*Power(s2,3)*(-1 + t1 - t2) - 4*t2 - 4*t1*t2 - 
                  2*Power(t1,3)*t2 + 5*Power(t2,2) + 
                  6*Power(t1,2)*Power(t2,2) - 2*t1*Power(t2,3) - 
                  Power(t2,4) + 6*Power(s1,3)*(-1 + s2 - t1 + t2) + 
                  6*Power(s2,2)*
                   (-1 - 2*Power(t1,2) + Power(t2,2) + t1*(3 + t2)) + 
                  6*Power(s1,2)*
                   (-1 + 3*Power(s2,2) + Power(t1,2) + 3*t2 + 
                     t1*t2 - 2*Power(t2,2) + s2*(-2 - 4*t1 + t2)) + 
                  s2*(12 + 7*Power(t1,3) - 4*t2 - 9*Power(t2,2) + 
                     11*Power(t2,3) + Power(t1,2)*(-13 + 3*t2) + 
                     t1*(-6 + 16*t2 - 21*Power(t2,2))) + 
                  s1*(12 + 6*Power(s2,3) + 11*Power(t1,3) + 
                     6*Power(s2,2)*(-2 + t1 - 4*t2) - 6*t2 - 
                     13*Power(t2,2) + 7*Power(t2,3) - 
                     3*Power(t1,2)*(3 + 7*t2) + 
                     t1*(-4 + 16*t2 + 3*Power(t2,2)) + 
                     s2*(-6 - 23*Power(t1,2) + 14*t2 - 
                        23*Power(t2,2) + 2*t1*(7 + 26*t2)))) - 
               2*s*(-2 + 3*t1 + 2*Power(t1,2) - 5*Power(t1,3) + 
                  2*Power(t1,4) + 3*t2 - 6*t1*t2 + 5*Power(t1,2)*t2 - 
                  2*Power(t1,3)*t2 + 2*Power(t2,2) + 
                  5*t1*Power(t2,2) - 5*Power(t2,3) - 
                  2*t1*Power(t2,3) + 2*Power(t2,4) + 
                  3*Power(s2,3)*t2*(1 - t1 + t2) + 
                  3*Power(s1,3)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                  Power(s1,2)*
                   (3*Power(s2,3) + 2*Power(t1,3) + 
                     Power(t1,2)*(1 - 7*t2) + 3*Power(-1 + t2,2) - 
                     2*Power(s2,2)*(1 + 2*t1 + 2*t2) + 
                     s2*(-3 + t1 - Power(t1,2) + 10*t2 + 11*t1*t2 - 
                        7*Power(t2,2)) + t1*(1 - 6*t2 + 5*Power(t2,2))\
) + Power(s2,2)*(3 + t2 + Power(t2,2) + 2*Power(t2,3) + 
                     Power(t1,2)*(3 + 5*t2) - 
                     t1*(6 + 6*t2 + 7*Power(t2,2))) - 
                  s2*(Power(t1,3)*(5 + 2*t2) - 
                     Power(t1,2)*(10 + 7*t2 + 3*Power(t2,2)) + 
                     t1*(5 + 4*t2 + 5*Power(t2,2)) + 
                     t2*(1 + t2 - 3*Power(t2,2) + Power(t2,3))) + 
                  s1*(3*Power(t1,3) - Power(t1,4) + 
                     3*Power(s2,3)*(-1 + t1 - 2*t2) - 
                     5*Power(-1 + t2,2)*t2 + 
                     Power(t1,2)*(-1 - 5*t2 + 3*Power(t2,2)) - 
                     t1*(1 + 4*t2 - 7*Power(t2,2) + 2*Power(t2,3)) + 
                     Power(s2,2)*
                      (-3 - 7*Power(t1,2) + t2 - Power(t2,2) + 
                        t1*(10 + 11*t2)) + 
                     s2*(4 + 5*Power(t1,3) + t2 - 10*Power(t2,2) + 
                        5*Power(t2,3) - 5*Power(t1,2)*(2 + t2) + 
                        t1*(1 + 10*t2 - 5*Power(t2,2)))))))/
           (Power(s,2)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)) + 
          ((-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
               s*(-2 - 2*s1 - 2*s2 + t1 + t2))*
             (-4*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - 
                  s2*t2 + s*(-2 - 2*s1 - 2*s2 + t1 + t2),2) + 
               12*s*(-Power(s2,2) + 2*s2*t1 + 2*Power(s2,2)*t1 - 
                  Power(t1,2) - 4*s2*Power(t1,2) - 
                  Power(s2,2)*Power(t1,2) + 2*Power(t1,3) + 
                  2*s2*Power(t1,3) - Power(t1,4) - Power(s2,2)*t2 - 
                  Power(s2,3)*t2 + 3*s2*t1*t2 + 3*Power(s2,2)*t1*t2 + 
                  Power(s2,3)*t1*t2 - 2*Power(t1,2)*t2 - 
                  4*s2*Power(t1,2)*t2 - 2*Power(s2,2)*Power(t1,2)*t2 + 
                  2*Power(t1,3)*t2 + s2*Power(t1,3)*t2 - Power(t2,2) + 
                  2*s2*Power(t2,2) - 2*Power(s2,2)*Power(t2,2) - 
                  Power(s2,3)*Power(t2,2) - 2*t1*Power(t2,2) + 
                  4*s2*t1*Power(t2,2) + 3*Power(s2,2)*t1*Power(t2,2) - 
                  2*Power(t1,2)*Power(t2,2) - 
                  2*s2*Power(t1,2)*Power(t2,2) + 2*Power(t2,3) - 
                  2*s2*Power(t2,3) - Power(s2,2)*Power(t2,3) + 
                  2*t1*Power(t2,3) + s2*t1*Power(t2,3) - Power(t2,4) - 
                  Power(s1,3)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                  Power(s,3)*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) - 
                  Power(s1,2)*
                   (Power(s2,3) + Power(t1,3) + 
                     Power(t1,2)*(2 - 3*t2) + Power(-1 + t2,2) - 
                     Power(s2,2)*(t1 + t2) + 
                     t1*(1 - 3*t2 + 2*Power(t2,2)) - 
                     s2*(1 + Power(t1,2) + t1*(2 - 4*t2) - 3*t2 + 
                        2*Power(t2,2))) + 
                  Power(s,2)*
                   (-2 + t1 + Power(t1,2) - Power(t1,3) + t2 - 
                     2*t1*t2 + Power(t1,2)*t2 + Power(t2,2) + 
                     t1*Power(t2,2) - Power(t2,3) - 
                     2*Power(s1,2)*(-1 + s2 - t1 + t2) + 
                     Power(s2,2)*(2 - 2*t1 + 2*t2) + 
                     s2*(-3*t1 + 3*Power(t1,2) + t2 - 4*t1*t2 + 
                        Power(t2,2)) + 
                     s1*(-2*Power(s2,2) + t1 + Power(t1,2) - 4*t1*t2 + 
                        3*(-1 + t2)*t2 + s2*(2 + t1 + t2))) + 
                  s1*(Power(t1,3)*(-2 + t2) + 2*Power(-1 + t2,2)*t2 + 
                     Power(s2,3)*(1 - t1 + 2*t2) + 
                     Power(t1,2)*(2 + 4*t2 - 2*Power(t2,2)) + 
                     t1*t2*(3 - 4*t2 + Power(t2,2)) + 
                     Power(s2,2)*
                      (2*Power(t1,2) + Power(1 + t2,2) - 
                        t1*(3 + 4*t2)) - 
                     s2*(Power(t1,3) - Power(t1,2)*(4 + t2) + 
                        t1*(3 + 6*t2 - Power(t2,2)) + 
                        t2*(3 - 4*t2 + Power(t2,2)))) + 
                  s*(-2*t1 + 2*Power(t1,2) + 
                     Power(s2,3)*(-1 + t1 - t2) - 2*t2 - 3*t1*t2 - 
                     Power(t1,3)*t2 + 2*Power(t2,2) + 
                     2*Power(t1,2)*Power(t2,2) - t1*Power(t2,3) + 
                     Power(s1,3)*(-1 + s2 - t1 + t2) + 
                     Power(s2,2)*
                      (-1 - 2*Power(t1,2) + Power(t2,2) + t1*(3 + t2)) \
+ Power(s1,2)*(-1 + 3*Power(s2,2) + Power(t1,2) + 3*t2 + t1*t2 - 
                        2*Power(t2,2) + s2*(-2 - 4*t1 + t2)) + 
                     s1*(2 + Power(s2,3) + 2*Power(t1,3) + 
                        Power(s2,2)*(-2 + t1 - 4*t2) - t2 - 
                        4*Power(t1,2)*t2 - 2*Power(t2,2) + 
                        Power(t2,3) + 
                        s2*(-1 + t1 - 4*Power(t1,2) + t2 + 9*t1*t2 - 
                        4*Power(t2,2)) + t1*(2 + t2 + Power(t2,2))) + 
                     s2*(Power(t1,3) + Power(t1,2)*(-2 + t2) + 
                        t1*(-1 + t2 - 4*Power(t2,2)) + 
                        2*(1 + t2 + Power(t2,3)))))))/
           (Power(s,3)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)))/
        ((-1 + t1)*(-1 + t2)) + 
       (64*(s - s2 + t1)*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) + 
          (2*(-23 + 9*Power(s,2) + 6*Power(s1,2) + 4*s2 + 
               3*Power(s2,2) + 5*t1 - 2*Power(t1,2) + 
               s1*(7 + 11*s2 - 4*t1 - 6*t2) - 
               2*s*(1 + 8*s1 + 6*s2 - t1 - 4*t2) + 2*t2 - 7*s2*t2 + 
               3*t1*t2 + Power(t2,2))*
             (-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
               s*(-2 - 2*s1 - 2*s2 + t1 + t2)))/s + 
          32*(2 + 10*s2 - 12*Power(s2,2) + Power(s2,3) - 12*t1 + 
             12*s2*t1 + 10*Power(s2,2)*t1 - 23*s2*Power(t1,2) + 
             12*Power(t1,3) + Power(s,2)*(-1 + 2*s1 - s2 + t1 - t2) - 
             11*Power(s2,2)*t2 + 2*t1*t2 + 35*s2*t1*t2 - 
             24*Power(t1,2)*t2 - 3*Power(t2,2) - 11*s2*Power(t2,2) + 
             11*t1*Power(t2,2) + Power(t2,3) + 
             Power(s1,2)*(-1 + 2*s2 - 2*t1 + t2) + 
             s1*(12*Power(s2,2) + 11*Power(t1,2) + t1*(12 - 9*t2) - 
                2*Power(-1 + t2,2) + s2*(-12 - 23*t1 + 9*t2)) + 
             s*(-13 - 2*Power(s1,2) + 14*Power(t1,2) + 
                s1*(15 - 15*s2 + 14*t1 - 12*t2) - 2*t2 - 28*t1*t2 + 
                14*Power(t2,2) + s2*(15 - 14*t1 + 14*t2))) + 
          4*(-53 + Power(s,3) - 4*Power(s1,3) - 10*s2 + 
             21*Power(s2,2) - 3*Power(s2,3) + 40*t1 - 40*s2*t1 - 
             55*Power(s2,2)*t1 + 14*Power(t1,2) + 121*s2*Power(t1,2) - 
             63*Power(t1,3) + 16*t2 + 9*s2*t2 + 63*Power(s2,2)*t2 - 
             65*t1*t2 - 173*s2*t1*t2 + 116*Power(t1,2)*t2 + 
             55*Power(t2,2) + 40*s2*Power(t2,2) - 39*t1*Power(t2,2) - 
             14*Power(t2,3) + 
             Power(s,2)*(4 - 44*s1 + 7*s2 - 23*t1 + 27*t2) - 
             Power(s1,2)*(38*s2 + 5*(2 - 7*t1 + t2)) + 
             s1*(18 - 60*Power(s2,2) - 36*Power(t1,2) - 50*t2 + 
                23*Power(t2,2) + s2*(19 + 97*t1 + 4*t2) + 
                t1*(-42 + 5*t2)) + 
             s*(70 + 47*Power(s1,2) - 5*Power(s2,2) - 101*Power(t1,2) + 
                s2*(-65 + 110*t1 - 116*t2) + 6*t2 + 184*t1*t2 - 
                80*Power(t2,2) + s1*(-47 + 124*s2 - 86*t1 + 33*t2))) + 
          16*(-13 + Power(s,3) - Power(s1,3) - 26*s2 + 42*Power(s2,2) - 
             6*Power(s2,3) + 44*t1 - 42*s2*t1 - 40*Power(s2,2)*t1 - 
             Power(t1,2) + 98*s2*Power(t1,2) - 52*Power(t1,3) + 
             Power(s1,2)*(4 - 17*s2 + 16*t1 - 5*t2) - 2*t2 - 5*s2*t2 + 
             48*Power(s2,2)*t2 - 15*t1*t2 - 150*s2*t1*t2 + 
             103*Power(t1,2)*t2 + 24*Power(t2,2) + 44*s2*Power(t2,2) - 
             44*t1*Power(t2,2) - 7*Power(t2,3) + 
             Power(s,2)*(8 - 19*s1 + 4*s2 - 8*t1 + 10*t2) + 
             s1*(20 - 53*Power(s2,2) - 43*Power(t1,2) + 
                s2*(48 + 96*t1 - 26*t2) - 29*t2 + 13*Power(t2,2) + 
                t1*(-51 + 28*t2)) + 
             s*(49 + 19*Power(s1,2) + Power(s2,2) - 68*Power(t1,2) + 
                s2*(-69 + 68*t1 - 72*t2) + 16*t2 + 134*t1*t2 - 
                65*Power(t2,2) + s1*(-68 + 80*s2 - 67*t1 + 46*t2))) - 
          8*(-29 + 5*Power(s,3) - 5*Power(s1,3) - 14*s2 + 
             50*Power(s2,2) - 9*Power(s2,3) + 63*t1 - 61*s2*t1 - 
             77*Power(s2,2)*t1 + 4*Power(t1,2) + 182*s2*Power(t1,2) - 
             96*Power(t1,3) - 5*t2 - 9*s2*t2 + 94*Power(s2,2)*t2 - 
             52*t1*t2 - 272*s2*t1*t2 + 184*Power(t1,2)*t2 + 
             65*Power(t2,2) + 71*s2*Power(t2,2) - 71*t1*Power(t2,2) - 
             17*Power(t2,3) - Power(s1,2)*(3 + 49*s2 - 44*t1 + 7*t2) + 
             Power(s,2)*(14 - 59*s1 + 3*s2 - 23*t1 + 33*t2) + 
             s1*(54 - 98*Power(s2,2) - 67*Power(t1,2) + 
                s2*(63 + 165*t1 - 17*t2) - 69*t2 + 29*Power(t2,2) + 
                t1*(-85 + 28*t2)) + 
             s*(71 + 59*Power(s1,2) + Power(s2,2) - 140*Power(t1,2) + 
                2*s2*(-57 + 72*t1 - 80*t2) + 32*t2 + 266*t1*t2 - 
                121*Power(t2,2) + s1*(-106 + 179*s2 - 132*t1 + 63*t2))) \
- (2*(-1 + 7*s - 6*s1 - 5*s2 + t1 + 3*t2)*
             (6*Power(s,4)*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) + 
               (-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
               2*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - 
                  s2*t2 + s*(-2 - 2*s1 - 2*s2 + t1 + t2),2) - 
               2*Power(s,3)*(6 - 3*t1 - 5*Power(t1,2) + 
                  3*Power(t1,3) + 6*Power(s2,2)*(-1 + t1 - t2) - 
                  3*t2 + 10*t1*t2 - 3*Power(t1,2)*t2 - 
                  5*Power(t2,2) - 3*t1*Power(t2,2) + 3*Power(t2,3) + 
                  6*Power(s1,2)*(-1 + s2 - t1 + t2) - 
                  3*s2*(-3*t1 + 3*Power(t1,2) + t2 - 4*t1*t2 + 
                     Power(t2,2)) + 
                  3*s1*(2*Power(s2,2) - Power(t1,2) - 
                     3*(-1 + t2)*t2 - s2*(2 + t1 + t2) + 
                     t1*(-1 + 4*t2))) + 
               Power(s,2)*(-4*t1 + 5*Power(t1,2) - Power(t1,4) + 
                  6*Power(s2,3)*(-1 + t1 - t2) - 4*t2 - 4*t1*t2 - 
                  2*Power(t1,3)*t2 + 5*Power(t2,2) + 
                  6*Power(t1,2)*Power(t2,2) - 2*t1*Power(t2,3) - 
                  Power(t2,4) + 6*Power(s1,3)*(-1 + s2 - t1 + t2) + 
                  6*Power(s2,2)*
                   (-1 - 2*Power(t1,2) + Power(t2,2) + t1*(3 + t2)) + 
                  6*Power(s1,2)*
                   (-1 + 3*Power(s2,2) + Power(t1,2) + 3*t2 + 
                     t1*t2 - 2*Power(t2,2) + s2*(-2 - 4*t1 + t2)) + 
                  s2*(12 + 7*Power(t1,3) - 4*t2 - 9*Power(t2,2) + 
                     11*Power(t2,3) + Power(t1,2)*(-13 + 3*t2) + 
                     t1*(-6 + 16*t2 - 21*Power(t2,2))) + 
                  s1*(12 + 6*Power(s2,3) + 11*Power(t1,3) + 
                     6*Power(s2,2)*(-2 + t1 - 4*t2) - 6*t2 - 
                     13*Power(t2,2) + 7*Power(t2,3) - 
                     3*Power(t1,2)*(3 + 7*t2) + 
                     t1*(-4 + 16*t2 + 3*Power(t2,2)) + 
                     s2*(-6 - 23*Power(t1,2) + 14*t2 - 
                        23*Power(t2,2) + 2*t1*(7 + 26*t2)))) - 
               2*s*(-2 + 3*t1 + 2*Power(t1,2) - 5*Power(t1,3) + 
                  2*Power(t1,4) + 3*t2 - 6*t1*t2 + 5*Power(t1,2)*t2 - 
                  2*Power(t1,3)*t2 + 2*Power(t2,2) + 
                  5*t1*Power(t2,2) - 5*Power(t2,3) - 
                  2*t1*Power(t2,3) + 2*Power(t2,4) + 
                  3*Power(s2,3)*t2*(1 - t1 + t2) + 
                  3*Power(s1,3)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                  Power(s1,2)*
                   (3*Power(s2,3) + 2*Power(t1,3) + 
                     Power(t1,2)*(1 - 7*t2) + 3*Power(-1 + t2,2) - 
                     2*Power(s2,2)*(1 + 2*t1 + 2*t2) + 
                     s2*(-3 + t1 - Power(t1,2) + 10*t2 + 11*t1*t2 - 
                        7*Power(t2,2)) + t1*(1 - 6*t2 + 5*Power(t2,2))\
) + Power(s2,2)*(3 + t2 + Power(t2,2) + 2*Power(t2,3) + 
                     Power(t1,2)*(3 + 5*t2) - 
                     t1*(6 + 6*t2 + 7*Power(t2,2))) - 
                  s2*(Power(t1,3)*(5 + 2*t2) - 
                     Power(t1,2)*(10 + 7*t2 + 3*Power(t2,2)) + 
                     t1*(5 + 4*t2 + 5*Power(t2,2)) + 
                     t2*(1 + t2 - 3*Power(t2,2) + Power(t2,3))) + 
                  s1*(3*Power(t1,3) - Power(t1,4) + 
                     3*Power(s2,3)*(-1 + t1 - 2*t2) - 
                     5*Power(-1 + t2,2)*t2 + 
                     Power(t1,2)*(-1 - 5*t2 + 3*Power(t2,2)) - 
                     t1*(1 + 4*t2 - 7*Power(t2,2) + 2*Power(t2,3)) + 
                     Power(s2,2)*
                      (-3 - 7*Power(t1,2) + t2 - Power(t2,2) + 
                        t1*(10 + 11*t2)) + 
                     s2*(4 + 5*Power(t1,3) + t2 - 10*Power(t2,2) + 
                        5*Power(t2,3) - 5*Power(t1,2)*(2 + t2) + 
                        t1*(1 + 10*t2 - 5*Power(t2,2)))))))/
           (Power(s,2)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)) + 
          ((-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
               s*(-2 - 2*s1 - 2*s2 + t1 + t2))*
             (-4*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - 
                  s2*t2 + s*(-2 - 2*s1 - 2*s2 + t1 + t2),2) + 
               12*s*(-Power(s2,2) + 2*s2*t1 + 2*Power(s2,2)*t1 - 
                  Power(t1,2) - 4*s2*Power(t1,2) - 
                  Power(s2,2)*Power(t1,2) + 2*Power(t1,3) + 
                  2*s2*Power(t1,3) - Power(t1,4) - Power(s2,2)*t2 - 
                  Power(s2,3)*t2 + 3*s2*t1*t2 + 3*Power(s2,2)*t1*t2 + 
                  Power(s2,3)*t1*t2 - 2*Power(t1,2)*t2 - 
                  4*s2*Power(t1,2)*t2 - 2*Power(s2,2)*Power(t1,2)*t2 + 
                  2*Power(t1,3)*t2 + s2*Power(t1,3)*t2 - Power(t2,2) + 
                  2*s2*Power(t2,2) - 2*Power(s2,2)*Power(t2,2) - 
                  Power(s2,3)*Power(t2,2) - 2*t1*Power(t2,2) + 
                  4*s2*t1*Power(t2,2) + 3*Power(s2,2)*t1*Power(t2,2) - 
                  2*Power(t1,2)*Power(t2,2) - 
                  2*s2*Power(t1,2)*Power(t2,2) + 2*Power(t2,3) - 
                  2*s2*Power(t2,3) - Power(s2,2)*Power(t2,3) + 
                  2*t1*Power(t2,3) + s2*t1*Power(t2,3) - Power(t2,4) - 
                  Power(s1,3)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                  Power(s,3)*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) - 
                  Power(s1,2)*
                   (Power(s2,3) + Power(t1,3) + 
                     Power(t1,2)*(2 - 3*t2) + Power(-1 + t2,2) - 
                     Power(s2,2)*(t1 + t2) + 
                     t1*(1 - 3*t2 + 2*Power(t2,2)) - 
                     s2*(1 + Power(t1,2) + t1*(2 - 4*t2) - 3*t2 + 
                        2*Power(t2,2))) + 
                  Power(s,2)*
                   (-2 + t1 + Power(t1,2) - Power(t1,3) + t2 - 
                     2*t1*t2 + Power(t1,2)*t2 + Power(t2,2) + 
                     t1*Power(t2,2) - Power(t2,3) - 
                     2*Power(s1,2)*(-1 + s2 - t1 + t2) + 
                     Power(s2,2)*(2 - 2*t1 + 2*t2) + 
                     s2*(-3*t1 + 3*Power(t1,2) + t2 - 4*t1*t2 + 
                        Power(t2,2)) + 
                     s1*(-2*Power(s2,2) + t1 + Power(t1,2) - 4*t1*t2 + 
                        3*(-1 + t2)*t2 + s2*(2 + t1 + t2))) + 
                  s1*(Power(t1,3)*(-2 + t2) + 2*Power(-1 + t2,2)*t2 + 
                     Power(s2,3)*(1 - t1 + 2*t2) + 
                     Power(t1,2)*(2 + 4*t2 - 2*Power(t2,2)) + 
                     t1*t2*(3 - 4*t2 + Power(t2,2)) + 
                     Power(s2,2)*
                      (2*Power(t1,2) + Power(1 + t2,2) - 
                        t1*(3 + 4*t2)) - 
                     s2*(Power(t1,3) - Power(t1,2)*(4 + t2) + 
                        t1*(3 + 6*t2 - Power(t2,2)) + 
                        t2*(3 - 4*t2 + Power(t2,2)))) + 
                  s*(-2*t1 + 2*Power(t1,2) + 
                     Power(s2,3)*(-1 + t1 - t2) - 2*t2 - 3*t1*t2 - 
                     Power(t1,3)*t2 + 2*Power(t2,2) + 
                     2*Power(t1,2)*Power(t2,2) - t1*Power(t2,3) + 
                     Power(s1,3)*(-1 + s2 - t1 + t2) + 
                     Power(s2,2)*
                      (-1 - 2*Power(t1,2) + Power(t2,2) + t1*(3 + t2)) \
+ Power(s1,2)*(-1 + 3*Power(s2,2) + Power(t1,2) + 3*t2 + t1*t2 - 
                        2*Power(t2,2) + s2*(-2 - 4*t1 + t2)) + 
                     s1*(2 + Power(s2,3) + 2*Power(t1,3) + 
                        Power(s2,2)*(-2 + t1 - 4*t2) - t2 - 
                        4*Power(t1,2)*t2 - 2*Power(t2,2) + 
                        Power(t2,3) + 
                        s2*(-1 + t1 - 4*Power(t1,2) + t2 + 9*t1*t2 - 
                        4*Power(t2,2)) + t1*(2 + t2 + Power(t2,2))) + 
                     s2*(Power(t1,3) + Power(t1,2)*(-2 + t2) + 
                        t1*(-1 + t2 - 4*Power(t2,2)) + 
                        2*(1 + t2 + Power(t2,3)))))))/
           (Power(s,3)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)))/
        ((-1 + s2)*(-s + s2 - t1)) + 
       (64*(-1 + s + t1)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2) + 
          (2*(-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
               s*(-2 - 2*s1 - 2*s2 + t1 + t2))*
             (-18 + 9*Power(s,2) + 12*Power(s1,2) + s2 + 
               3*Power(s2,2) + t1 + 3*s2*t1 - Power(t1,2) + 
               s1*(14*s2 + 3*t1 - 15*t2) + 4*t2 - 10*s2*t2 - t1*t2 + 
               4*Power(t2,2) + s*(3 - 21*s1 - 12*s2 - 3*t1 + 13*t2)))/s - 
          16*(-32 + 18*Power(s1,3) + 34*s2 + Power(s2,2) - Power(s2,3) + 
             43*t1 - 85*s2*t1 + 39*Power(t1,2) + 53*s2*Power(t1,2) - 
             52*Power(t1,3) + Power(s1,2)*(-15 + 2*s2 + 33*t1 - 44*t2) + 
             Power(s,2)*(22 + 2*s1 - 25*s2 + 13*t1 - 13*t2) + 2*t2 + 
             28*s2*t2 + 17*Power(s2,2)*t2 - 82*t1*t2 - 74*s2*t1*t2 + 
             105*Power(t1,2)*t2 + 48*Power(t2,2) + 13*s2*Power(t2,2) - 
             45*t1*Power(t2,2) - 8*Power(t2,3) - 
             s1*(-43 + 18*Power(s2,2) + 46*Power(t1,2) + 
                t1*(22 - 13*t2) + 34*t2 - 34*Power(t2,2) + 
                s2*(34 - 81*t1 + 15*t2)) + 
             s*(33 - 20*Power(s1,2) + 26*Power(s2,2) + 10*t1 - 
                65*Power(t1,2) + s2*(-63 + 37*t1 - 70*t2) + 17*t2 + 
                152*t1*t2 - 87*Power(t2,2) + 
                3*s1*(-21 + 30*s2 - 35*t1 + 36*t2))) + 
          4*(5 + Power(s,3) - 40*Power(s1,3) - 32*s2 + 3*Power(s2,3) - 
             66*t1 + 55*s2*t1 + 7*Power(s2,2)*t1 - 21*Power(t1,2) - 
             75*s2*Power(t1,2) + 64*Power(t1,3) - 
             Power(s1,2)*(17 + s2 + 75*t1 - 99*t2) + 48*t2 - s2*t2 - 
             41*Power(s2,2)*t2 + 105*t1*t2 + 123*s2*t1*t2 - 
             129*Power(t1,2)*t2 - 96*Power(t2,2) - 32*s2*Power(t2,2) + 
             46*t1*Power(t2,2) + 19*Power(t2,3) + 
             Power(s,2)*(-50 - 2*s1 + 47*s2 - 23*t1 + 23*t2) + 
             s1*(-87 + 42*Power(s2,2) + 45*Power(t1,2) + 111*t2 - 
                78*Power(t2,2) + t1*(12 + 29*t2) + 
                s2*(7 - 124*t1 + 31*t2)) + 
             s*(14 + 41*Power(s1,2) - 51*Power(s2,2) + 22*t1 + 
                97*Power(t1,2) - 81*t2 - 237*t1*t2 + 139*Power(t2,2) - 
                2*s1*(-63 + 72*s2 - 88*t1 + 90*t2) + 
                s2*(70 - 42*t1 + 106*t2))) - 
          8*(41 + 3*Power(s,3) - 52*Power(s1,3) - 30*s2 - 
             5*Power(s2,2) + 4*Power(s2,3) - 89*t1 + 119*s2*t1 + 
             2*Power(s2,2)*t1 - 50*Power(t1,2) - 103*s2*Power(t1,2) + 
             96*Power(t1,3) + 12*t2 - 30*s2*t2 - 44*Power(s2,2)*t2 + 
             146*t1*t2 + 163*s2*t1*t2 - 195*Power(t1,2)*t2 - 
             105*Power(t2,2) - 41*s2*Power(t2,2) + 77*t1*Power(t2,2) + 
             22*Power(t2,3) + 
             Power(s,2)*(-55 - 14*s1 + 57*s2 - 35*t1 + 39*t2) + 
             Power(s1,2)*(9 - 11*s2 - 89*t1 + 125*t2) + 
             s1*(-72 + 47*Power(s2,2) + 79*Power(t1,2) + 96*t2 - 
                95*Power(t2,2) + 2*t1*(17 + 5*t2) + 
                s2*(37 - 175*t1 + 51*t2)) + 
             s*(-49 + 63*Power(s1,2) - 64*Power(s2,2) + 6*t1 + 
                130*Power(t1,2) + 
                s1*(-181*s2 + 3*(48 + 79*t1 - 85*t2)) - 68*t2 - 
                320*t1*t2 + 191*Power(t2,2) + s2*(109 - 56*t1 + 133*t2))) \
- 32*(10 - 2*Power(s1,3) - 11*s2 - 11*t1 + 23*s2*t1 - 11*Power(t1,2) - 
             12*s2*Power(t1,2) + 12*Power(t1,3) - 9*s2*t2 - 
             2*Power(s2,2)*t2 + 21*t1*t2 + 14*s2*t1*t2 - 
             24*Power(t1,2)*t2 - 11*Power(t2,2) - s2*Power(t2,2) + 
             11*t1*Power(t2,2) + Power(t2,3) + 
             Power(s,2)*(-3 + 3*s2 - t1 + t2) + 
             Power(s1,2)*(3 - 4*t1 + 5*t2) + 
             s*(-10 + 2*Power(s1,2) - 3*Power(s2,2) - 3*t1 + 
                14*Power(t1,2) + s1*(13 - 17*s2 + 18*t1 - 18*t2) - t2 - 
                30*t1*t2 + 16*Power(t2,2) + s2*(15 - 11*t1 + 15*t2)) + 
             s1*(2*Power(s2,2) + 11*Power(t1,2) + t1*(4 - 7*t2) + 
                s2*(10 - 15*t1 + t2) - 4*(3 - 2*t2 + Power(t2,2)))) - 
          (2*(1 + 7*s - 8*s1 - 5*s2 - t1 + 5*t2)*
             (6*Power(s,4)*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) + 
               (-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + s1*(s2 - t1) + t1 + t2 - s2*t2,2) - 
               2*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - 
                  s2*t2 + s*(-2 - 2*s1 - 2*s2 + t1 + t2),2) - 
               2*Power(s,3)*(6 - 3*t1 - 5*Power(t1,2) + 
                  3*Power(t1,3) + 6*Power(s2,2)*(-1 + t1 - t2) - 
                  3*t2 + 10*t1*t2 - 3*Power(t1,2)*t2 - 5*Power(t2,2) - 
                  3*t1*Power(t2,2) + 3*Power(t2,3) + 
                  6*Power(s1,2)*(-1 + s2 - t1 + t2) - 
                  3*s2*(-3*t1 + 3*Power(t1,2) + t2 - 4*t1*t2 + 
                     Power(t2,2)) + 
                  3*s1*(2*Power(s2,2) - Power(t1,2) - 3*(-1 + t2)*t2 - 
                     s2*(2 + t1 + t2) + t1*(-1 + 4*t2))) + 
               Power(s,2)*(-4*t1 + 5*Power(t1,2) - Power(t1,4) + 
                  6*Power(s2,3)*(-1 + t1 - t2) - 4*t2 - 4*t1*t2 - 
                  2*Power(t1,3)*t2 + 5*Power(t2,2) + 
                  6*Power(t1,2)*Power(t2,2) - 2*t1*Power(t2,3) - 
                  Power(t2,4) + 6*Power(s1,3)*(-1 + s2 - t1 + t2) + 
                  6*Power(s2,2)*
                   (-1 - 2*Power(t1,2) + Power(t2,2) + t1*(3 + t2)) + 
                  6*Power(s1,2)*
                   (-1 + 3*Power(s2,2) + Power(t1,2) + 3*t2 + t1*t2 - 
                     2*Power(t2,2) + s2*(-2 - 4*t1 + t2)) + 
                  s2*(12 + 7*Power(t1,3) - 4*t2 - 9*Power(t2,2) + 
                     11*Power(t2,3) + Power(t1,2)*(-13 + 3*t2) + 
                     t1*(-6 + 16*t2 - 21*Power(t2,2))) + 
                  s1*(12 + 6*Power(s2,3) + 11*Power(t1,3) + 
                     6*Power(s2,2)*(-2 + t1 - 4*t2) - 6*t2 - 
                     13*Power(t2,2) + 7*Power(t2,3) - 
                     3*Power(t1,2)*(3 + 7*t2) + 
                     t1*(-4 + 16*t2 + 3*Power(t2,2)) + 
                     s2*(-6 - 23*Power(t1,2) + 14*t2 - 
                        23*Power(t2,2) + 2*t1*(7 + 26*t2)))) - 
               2*s*(-2 + 3*t1 + 2*Power(t1,2) - 5*Power(t1,3) + 
                  2*Power(t1,4) + 3*t2 - 6*t1*t2 + 5*Power(t1,2)*t2 - 
                  2*Power(t1,3)*t2 + 2*Power(t2,2) + 5*t1*Power(t2,2) - 
                  5*Power(t2,3) - 2*t1*Power(t2,3) + 2*Power(t2,4) + 
                  3*Power(s2,3)*t2*(1 - t1 + t2) + 
                  3*Power(s1,3)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                  Power(s1,2)*
                   (3*Power(s2,3) + 2*Power(t1,3) + 
                     Power(t1,2)*(1 - 7*t2) + 3*Power(-1 + t2,2) - 
                     2*Power(s2,2)*(1 + 2*t1 + 2*t2) + 
                     s2*(-3 + t1 - Power(t1,2) + 10*t2 + 11*t1*t2 - 
                        7*Power(t2,2)) + t1*(1 - 6*t2 + 5*Power(t2,2))) \
+ Power(s2,2)*(3 + t2 + Power(t2,2) + 2*Power(t2,3) + 
                     Power(t1,2)*(3 + 5*t2) - 
                     t1*(6 + 6*t2 + 7*Power(t2,2))) - 
                  s2*(Power(t1,3)*(5 + 2*t2) - 
                     Power(t1,2)*(10 + 7*t2 + 3*Power(t2,2)) + 
                     t1*(5 + 4*t2 + 5*Power(t2,2)) + 
                     t2*(1 + t2 - 3*Power(t2,2) + Power(t2,3))) + 
                  s1*(3*Power(t1,3) - Power(t1,4) + 
                     3*Power(s2,3)*(-1 + t1 - 2*t2) - 
                     5*Power(-1 + t2,2)*t2 + 
                     Power(t1,2)*(-1 - 5*t2 + 3*Power(t2,2)) - 
                     t1*(1 + 4*t2 - 7*Power(t2,2) + 2*Power(t2,3)) + 
                     Power(s2,2)*
                      (-3 - 7*Power(t1,2) + t2 - Power(t2,2) + 
                        t1*(10 + 11*t2)) + 
                     s2*(4 + 5*Power(t1,3) + t2 - 10*Power(t2,2) + 
                        5*Power(t2,3) - 5*Power(t1,2)*(2 + t2) + 
                        t1*(1 + 10*t2 - 5*Power(t2,2)))))))/
           (Power(s,2)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)) + 
          ((-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - s2*t2 + 
               s*(-2 - 2*s1 - 2*s2 + t1 + t2))*
             (-4*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2)*
                Power(-1 + 2*Power(s,2) + s1*(s2 - t1) + t1 + t2 - 
                  s2*t2 + s*(-2 - 2*s1 - 2*s2 + t1 + t2),2) + 
               12*s*(-Power(s2,2) + 2*s2*t1 + 2*Power(s2,2)*t1 - 
                  Power(t1,2) - 4*s2*Power(t1,2) - 
                  Power(s2,2)*Power(t1,2) + 2*Power(t1,3) + 
                  2*s2*Power(t1,3) - Power(t1,4) - Power(s2,2)*t2 - 
                  Power(s2,3)*t2 + 3*s2*t1*t2 + 3*Power(s2,2)*t1*t2 + 
                  Power(s2,3)*t1*t2 - 2*Power(t1,2)*t2 - 
                  4*s2*Power(t1,2)*t2 - 2*Power(s2,2)*Power(t1,2)*t2 + 
                  2*Power(t1,3)*t2 + s2*Power(t1,3)*t2 - Power(t2,2) + 
                  2*s2*Power(t2,2) - 2*Power(s2,2)*Power(t2,2) - 
                  Power(s2,3)*Power(t2,2) - 2*t1*Power(t2,2) + 
                  4*s2*t1*Power(t2,2) + 3*Power(s2,2)*t1*Power(t2,2) - 
                  2*Power(t1,2)*Power(t2,2) - 
                  2*s2*Power(t1,2)*Power(t2,2) + 2*Power(t2,3) - 
                  2*s2*Power(t2,3) - Power(s2,2)*Power(t2,3) + 
                  2*t1*Power(t2,3) + s2*t1*Power(t2,3) - Power(t2,4) - 
                  Power(s1,3)*(s2 - t1)*(-1 + s2 - t1 + t2) + 
                  Power(s,3)*(-1 + s1 + t1 - t2)*(-1 + s2 - t1 + t2) - 
                  Power(s1,2)*
                   (Power(s2,3) + Power(t1,3) + 
                     Power(t1,2)*(2 - 3*t2) + Power(-1 + t2,2) - 
                     Power(s2,2)*(t1 + t2) + 
                     t1*(1 - 3*t2 + 2*Power(t2,2)) - 
                     s2*(1 + Power(t1,2) + t1*(2 - 4*t2) - 3*t2 + 
                        2*Power(t2,2))) + 
                  Power(s,2)*
                   (-2 + t1 + Power(t1,2) - Power(t1,3) + t2 - 
                     2*t1*t2 + Power(t1,2)*t2 + Power(t2,2) + 
                     t1*Power(t2,2) - Power(t2,3) - 
                     2*Power(s1,2)*(-1 + s2 - t1 + t2) + 
                     Power(s2,2)*(2 - 2*t1 + 2*t2) + 
                     s2*(-3*t1 + 3*Power(t1,2) + t2 - 4*t1*t2 + 
                        Power(t2,2)) + 
                     s1*(-2*Power(s2,2) + t1 + Power(t1,2) - 4*t1*t2 + 
                        3*(-1 + t2)*t2 + s2*(2 + t1 + t2))) + 
                  s1*(Power(t1,3)*(-2 + t2) + 2*Power(-1 + t2,2)*t2 + 
                     Power(s2,3)*(1 - t1 + 2*t2) + 
                     Power(t1,2)*(2 + 4*t2 - 2*Power(t2,2)) + 
                     t1*t2*(3 - 4*t2 + Power(t2,2)) + 
                     Power(s2,2)*
                      (2*Power(t1,2) + Power(1 + t2,2) - t1*(3 + 4*t2)) \
- s2*(Power(t1,3) - Power(t1,2)*(4 + t2) + 
                        t1*(3 + 6*t2 - Power(t2,2)) + 
                        t2*(3 - 4*t2 + Power(t2,2)))) + 
                  s*(-2*t1 + 2*Power(t1,2) + 
                     Power(s2,3)*(-1 + t1 - t2) - 2*t2 - 3*t1*t2 - 
                     Power(t1,3)*t2 + 2*Power(t2,2) + 
                     2*Power(t1,2)*Power(t2,2) - t1*Power(t2,3) + 
                     Power(s1,3)*(-1 + s2 - t1 + t2) + 
                     Power(s2,2)*
                      (-1 - 2*Power(t1,2) + Power(t2,2) + t1*(3 + t2)) + 
                     Power(s1,2)*
                      (-1 + 3*Power(s2,2) + Power(t1,2) + 3*t2 + 
                        t1*t2 - 2*Power(t2,2) + s2*(-2 - 4*t1 + t2)) + 
                     s1*(2 + Power(s2,3) + 2*Power(t1,3) + 
                        Power(s2,2)*(-2 + t1 - 4*t2) - t2 - 
                        4*Power(t1,2)*t2 - 2*Power(t2,2) + 
                        Power(t2,3) + 
                        s2*(-1 + t1 - 4*Power(t1,2) + t2 + 9*t1*t2 - 
                        4*Power(t2,2)) + t1*(2 + t2 + Power(t2,2))) + 
                     s2*(Power(t1,3) + Power(t1,2)*(-2 + t2) + 
                        t1*(-1 + t2 - 4*Power(t2,2)) + 
                        2*(1 + t2 + Power(t2,3)))))))/
           (Power(s,3)*(-1 + s1 + t1 - t2)*(1 - s2 + t1 - t2)))/
        ((-1 + s2)*(-1 + t1)))*
     B3(1 - s1 - t1 + t2,1 - s2 + t1 - t2,2 + s - s1 - s2))/(4.*Power(Pi,2));
     return a;
};
