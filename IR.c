#include <gsl/gsl_math.h>
#include <gsl/gsl_pow_int.h>
#include "nlo_functions.h"

#define Pi M_PI
#define Sqrt sqrt
#define Power gsl_pow_int
#define ln gsl_sf_log_abs
#define Log gsl_sf_log_abs

double IK(double s,double s1,double s2,double t1,double t2,double EE,double np1,double np2){
    double a=(((2 - 4*s + 6*s1 - 2*s2 + 2*s1*s2 - 4*t1 + 4*s1*t1 + 8*t2 + 8*s*t2 + 
          4*s2*t2 - 4*s1*s2*t2 - 
          (2*(2 + s*(-3 + s1) - s2 - Power(s1,2)*s2 + 2*s1*(3 + s2))*
             (-1 + t2))/(s - s1 + t2) + 
          ((-1 + s1)*(-5 + 2*s2 + 5*t1 + 21*t2 - 3*s2*t2 - 2*t1*t2 + 
               2*Power(t2,2) - 2*s2*Power(t2,2) - 
               s*(-2 + t1 - t2)*(-1 + 2*t2) + 
               s1*(8 - s2 - 5*t1 - 8*t2 + 2*s2*t2 + 2*t1*t2)))/
           (-1 + t1) + ((-1 + s1)*(-1 + t2)*
             (Power(s1,2)*(s2 - t1) + 
               s*(-2*(6 + 3*t1 - 5*t2) + s1*(4 + t1 - t2)) + 
               s1*(19 + 9*t1 + s2*(-8 + t2) - 9*t2) - 
               2*(-6 + 2*t1*(-2 + t2) + (9 + s2)*t2 - 4*Power(t2,2))))/
           ((s - s2 + t1)*(s - s1 + t2)) + 
          ((-1 + s1)*(-1 + t2)*
             (12 + 6*s2 + 7*t1 - 2*s2*t1 - Power(t1,2) + 
               s1*(6 + s2*(-6 + t1) + 2*t1 - Power(t1,2)) - 6*t2 - 
               6*s2*t2 + 3*t1*t2 + s2*t1*t2 + 
               s*(Power(t1,2) + 6*(-1 + t2) - t1*(2 + t2))))/
           ((-1 + s2)*(-s + s2 - t1)) + 
          ((-1 + s1)*(-1 + t2)*
             (2 + 5*t1 - 8*s1*t1 + Power(t1,2) + s1*Power(t1,2) + 
               18*t2 - 3*t1*t2 + 
               s2*(6 + s1*(-2 + t1) - 2*t2 - t1*(2 + t2)) - 
               s*(2 + Power(t1,2) + 6*t2 - t1*(4 + t2))))/
           ((-1 + s2)*(-1 + t1)))/(Power(-1 + s1,2)*Power(-1 + t2,2)) + 
       (20 - 22*s + 6*s1 + 12*s2 + 4*s*s2 + 4*s1*s2 - 8*Power(s2,2) - 
          12*t1 - 6*s*t1 - 4*Power(s,2)*t1 - 2*s1*t1 + 4*s*s1*t1 + 
          16*s2*t1 + 4*s*s2*t1 - 4*s1*s2*t1 - 8*Power(t1,2) - 
          4*s*Power(t1,2) + 4*s1*Power(t1,2) - 10*t2 + 8*s*t2 - 
          12*s2*t2 + 10*t1*t2 - 4*s*t1*t2 + 4*s2*t1*t2 - 
          4*Power(t1,2)*t2 + 
          ((s - s2 + t1)*(s - s1 + t2)*
             (13 + 4*s2 - 5*t1 + Power(s1,2)*(-s2 + t1) - 14*t2 + 
               s2*t2 + t1*t2 + 3*Power(t2,2) - s2*Power(t2,2) + 
               2*t1*Power(t2,2) + Power(s,2)*(-2 + t1 + t2) + 
               s1*(13 - 5*s2 + 8*t1 - 3*t2 + 2*s2*t2 - 3*t1*t2) + 
               s*(-11 + 2*s2 + s1*(10 + s2 - 2*t1 - t2) + 2*t2 - 
                  s2*t2 + 3*t1*t2 + Power(t2,2))))/
           ((-1 + s1)*(-1 + t2)) + 
          ((s - s1 + t2)*(-21 + 23*t1 - 9*s1*t1 - 2*Power(t1,2) + 
               2*s1*Power(t1,2) + 3*t2 + 2*t1*t2 + 
               Power(s2,2)*(4 - 2*s1 + 2*t2) + 
               Power(s,2)*(4 - 2*t1 + 2*t2) + 
               s*(24 + 3*t1 - 2*Power(t1,2) + 2*s1*(-4 + s2 + t1) + 
                  2*s2*(-6 + t1 - 2*t2) + t2 + 2*t1*t2) + 
               s2*(11*s1 - 2*t1*(1 + t2) - 5*(4 + t2))))/(-1 + s2) - 
          ((s - s2 + t1)*(s - s1 + t2)*
             (-12 - 9*s2 + 2*Power(s2,2) + 6*t1 - s2*t1 + 
               s1*(Power(s2,2) + 6*(-1 + t1) - s2*(2 + t1)) + 6*t2 + 
               5*s2*t2 - Power(s2,2)*t2 - 12*t1*t2 + 2*s2*t1*t2 + 
               s*(s2*(-2 + t1 + t2) - 6*(-1 + t1 + t2))))/
           ((-1 + t1)*(-1 + t2)) + 
          (2*(s - s2 + t1)*(-6 - 2*s2 + 2*t1 + Power(s,2)*t1 + 
               Power(s1,2)*t1 + 5*t2 - s2*t2 + t1*t2 - Power(t2,2) + 
               t1*Power(t2,2) - 
               s*(-7 + s2 + t1*(-1 + 2*s1 - 2*t2) + t2) + 
               s1*(-7 + s2 + t2 - t1*(1 + 2*t2))))/(-1 + s1) - 
          ((s - s2 + t1)*(s - s1 + t2)*
             (Power(s2,2)*(-s1 + t2) + 
               s*(12 - 10*t1 + s2*(-4 + t1 - t2) + 6*t2) - 
               s2*(19 + s1*(-8 + t1) - 9*t1 + 9*t2) + 
               2*(-6 - 4*Power(t1,2) - 4*t2 + t1*(9 + s1 + 2*t2))))/
           ((-1 + s2)*(-1 + t1)))/
        (Power(s - s2 + t1,2)*Power(s - s1 + t2,2)) + 
       (2*(3 - 3*s2 + 6*Power(s2,2) - Power(s2,3) - 6*t1 + 6*s2*t1 - 
             2*Power(s2,2)*t1 + 2*Power(t1,2) - s2*Power(t1,2) - 
             s1*(4 - 3*s2 + Power(s2,2) - 3*t1 + 2*s2*t1 + 
                Power(t1,2)) + 
             s*(1 + Power(s2,2) + 2*s2*(-2 + t1) - 3*t1 + 
                Power(t1,2) - 2*t2) + 5*t2 + 2*s2*t2 - t1*t2) + 
          ((-1 + t1)*(-5 + 4*Power(s2,3) - 3*t1 + 3*s1*t1 + 
               4*Power(t1,2) - 2*s1*Power(t1,2) + 
               2*Power(s2,2)*(-8 + 2*s1 + t1) + 
               2*Power(s,2)*(-2 + s2 + t1) + 
               s*(14*s2 - 6*Power(s2,2) + t1 - 4*s2*t1 + 
                  2*Power(t1,2) - 2*s1*(-2 + s2 + t1) - t2) + t2 - 
               2*t1*t2 + s2*(4 - 9*s1 - 2*t1 + 2*s1*t1 - 
                  2*Power(t1,2) + t2)))/(s - s2 + t1) + 
          (2*(-1 + s2)*(-1 + t1)*
             (5 + 2*Power(s1,2) + 2*s2 + Power(s2,2) - 6*t1 + 
               Power(t1,2) + 6*t2 + 2*s2*t2 - 2*t1*t2 + Power(t2,2) - 
               2*s1*(1 + t2) - Power(s,2)*(-2 + s2 - t1 + t2) + 
               s*(-4 + Power(s2,2) + s2*(-1 - t1 + t2) + 
                  s1*(-3 + s2 - t1 + t2))))/((-1 + s1)*(-s + s1 - t2)) \
+ (2*(-1 + s2)*(-1 + t1)*(3 + Power(s1,2) - s2 + Power(s2,2) + 2*t1 + 
               s2*t1 - Power(t1,2) + 2*t2 - s2*t2 + 2*t1*t2 - 
               Power(t2,2) + s1*(-1 + s2 - t1 + t2) - 
               Power(s,2)*(-1 + t1 + t2) + 
               s*(-t1 - t2 + s1*(-2 + t1 + t2) + s2*(-2 + t1 + t2))))/
           ((-1 + s1)*(-1 + t2)) - 
          ((-1 + s2)*(-5 + 4*s2 - 4*Power(s2,2) + 5*t1 - 8*s2*t1 + 
               2*Power(s2,2)*t1 - 6*Power(t1,2) + 2*s2*Power(t1,2) + 
               t2 - 7*s2*t2 + 2*Power(s2,2)*t2 - 2*t1*t2 + 
               2*s2*t1*t2 - s*
                (4 - 7*t1 + 2*Power(t1,2) - 3*t2 + 2*t1*t2 + 
                  2*s2*(-2 + t1 + t2)) + 
               s1*(4 - 7*t1 + 2*Power(t1,2) - 4*t2 + 2*t1*t2 + 
                  s2*(-3 + 2*t1 + 2*t2))))/(-1 + t2) - 
          (2*(-1 + s2)*(-1 + t1)*
             (-3 - Power(s1,2) - 4*s2 - 4*Power(s2,2) + 
               2*Power(s2,3) + 4*t1 + 3*s2*t1 - 3*Power(s2,2)*t1 - 
               Power(t1,2) + s2*Power(t1,2) - t2 - s2*t2 + 
               2*Power(s2,2)*t2 - t1*t2 - s2*t1*t2 + 
               Power(s,2)*(-2 + s2 - t1 + t2) + 
               s1*(-1 + 2*Power(s2,2) + Power(t1,2) - t1*(-2 + t2) + 
                  t2 + s2*(-3 - 3*t1 + 2*t2)) - 
               s*(-1 + 3*Power(s2,2) + 2*t1 + Power(t1,2) + t2 - 
                  t1*t2 + s1*(-4 + s2 - t1 + t2) + 
                  s2*(-5 - 4*t1 + 3*t2))))/((s - s2 + t1)*(s - s1 + t2))\
)/(Power(-1 + s2,2)*Power(-1 + t1,2)) + 
       (-2*(-2 + Power(s,2)*(-1 + s2) - 7*s2 - 6*Power(s2,2) + 
             Power(s2,3) + 7*s2*t1 - Power(s2,2)*t1 - Power(t1,2) + 
             s1*(-2 + 2*s2 + Power(s2,2) - t1 - s2*t1) - 
             s*(-4 + s1 - 6*s2 + s1*s2 + 2*Power(s2,2) + 2*t1 - 
                s2*t1) + t2 + s2*t2) + 
          ((-1 + s2)*(-7 - 2*Power(s1,2) - 21*s2 - 2*Power(s2,2) + 
               2*Power(s2,3) + 6*t1 + 9*s2*t1 - 4*Power(s2,2)*t1 - 
               5*Power(t1,2) + 2*s2*Power(t1,2) + 7*t2 + 8*s2*t2 + 
               Power(s2,2)*t2 - 3*t1*t2 - s2*t1*t2 - 2*Power(t2,2) + 
               Power(s,2)*(-8 + 2*s2 - t1 + t2) + 
               s1*(-6 + Power(s2,2) - t1 + Power(t1,2) - 
                  s2*(1 + 2*t1) + 4*t2) + 
               s*(17 - 4*Power(s2,2) - 8*t1 - Power(t1,2) + 
                  s1*(6 - s2 + t1) + s2*(12 + 5*t1 - 2*t2) - 10*t2 + 
                  t1*t2)))/(s - s1 + t2) + 
          ((s - s2 + t1)*(-2*Power(s,2) - 7*s2 - 8*Power(s2,2) - 
               2*t1 + s2*t1 + 2*Power(s2,2)*t1 + 2*Power(t1,2) + 
               s1*(-2 + 2*s2 + Power(s2,2) + 2*t1 + s2*t1) - 6*t2 - 
               s2*t2 - Power(s2,2)*t2 + 
               s*(4 - 4*t1 + 2*t2 + s2*(6 - t1 + t2))))/(-1 + t1) + 
          ((-1 + s2)*(-s + s2 - t1)*
             (10 + 4*s2 + 2*Power(s2,2) + 2*Power(s1,2)*(1 + s2) - 
               t1 - 4*s2*t1 + Power(t1,2) - t2 - 2*s2*t2 + 2*t1*t2 - 
               s2*t1*t2 + Power(t2,2) + s2*Power(t2,2) - 
               Power(s,2)*(-6 + t1 + t2) - 
               s*(15 + s1*(4 + 3*s2 - t1) - 5*t1 + Power(t1,2) - 
                  s2*(-4 + t2) - 5*t2 + Power(t2,2)) + 
               s1*(4 + 2*Power(s2,2) + Power(t1,2) - 4*t2 - 
                  t1*(2 + t2) - s2*(-2 + t1 + t2))))/
           ((-1 + s1)*(-s + s1 - t2)) + 
          ((-1 + s2)*(-s + s2 - t1)*
             (9 - 3*s2 - 2*Power(s2,2) + 
               2*Power(s1,2)*(1 + s2 - t1) + 10*t1 + 3*s2*t1 - 
               Power(t1,2) - Power(s,2)*(6 + t1 - 3*t2) - 6*t2 - 
               2*s2*t2 + Power(s2,2)*t2 - 2*t1*t2 + s2*t1*t2 - 
               Power(t2,2) - s2*Power(t2,2) + 
               s1*(4 + Power(s2,2) - Power(t1,2) + s2*(-1 + t2) - 
                  2*t2 + t1*(5 + t2)) + 
               s*(-3 - 8*t1 + Power(t1,2) + s2*(8 + t1 - 4*t2) + 
                  2*t2 + Power(t2,2) - s1*(-2 + s2 - 3*t1 + 2*t2))))/
           ((-1 + s1)*(-1 + t2)) - 
          ((-1 + s2)*(-s + s2 - t1)*
             (-9*s2 + 6*Power(s2,2) - 6*t1 - 9*s2*t1 - 2*Power(t1,2) + 
               2*s2*Power(t1,2) + Power(s1,2)*(-2 + s2 + t1) + 
               9*s2*t2 - 3*Power(s2,2)*t2 + 2*Power(t2,2) - 
               Power(s,2)*(-2 + t1 + t2) + 
               s1*(-3 + Power(s2,2) - 3*t1 + 2*Power(t1,2) + 
                  s2*(-2 + t1 - 3*t2) + t2) + 
               s*(3 + 3*t1 - 2*Power(t1,2) - 5*t2 + s1*(2 - s2 + t2) + 
                  s2*(-8 + t1 + 4*t2))))/((-1 + t1)*(-1 + t2)))/
        (Power(-1 + s2,2)*Power(s - s2 + t1,2)) + 
       (2*(1 + 3*s2 - 2*t1 - s2*t1 + 2*Power(t1,2) - 2*t2 + 5*t1*t2 - 
             s2*t1*t2 + 2*Power(t2,2) - s1*(-3 + t2 + t1*t2) + 
             s*(-2 + t1 + t2 + t1*t2)) + 
          ((-1 + t1)*(2*Power(s1,2) + 6*s2 + 2*t1 - 2*Power(t1,2) - 
               t2 - 4*s2*t2 + 7*t1*t2 + 3*Power(t2,2) - 
               s2*Power(t2,2) + s1*(-((-2 + t1)*t2) + s2*(2 + t2)) + 
               s*(-2 + 2*t2 + Power(t2,2) + t1*(2 + t2))))/(-1 + s1) + 
          ((-1 + t2)*(4 + 10*s2 - 2*Power(s2,2) - 7*t1 + 2*s2*t1 + 
               5*Power(t1,2) - 2*s2*Power(t1,2) - 
               s1*(-2 + 2*t1 + Power(t1,2) + s2*(2 + t1)) + 4*s2*t2 + 
               3*t1*t2 + s2*t1*t2 + 2*Power(t2,2) + 
               s*(Power(t1,2) - t1*(-4 + t2) - 2*(3 + t2))))/(-1 + s2) \
+ ((-1 + t1)*(-1 + t2)*(2 + 7*s2 + 2*Power(s2,2) - 4*t1 - 5*s2*t1 + 
               2*Power(t1,2) + 2*Power(s1,2)*(1 + t1) + 5*t2 + 
               7*s2*t2 - Power(s2,2)*t2 - 5*t1*t2 + 3*Power(t2,2) - 
               s2*Power(t2,2) + Power(s,2)*(-2 + t1 + t2) + 
               s*(1 + s2*(2 + s1 - t1) - 3*t2 + Power(t2,2) + 
                  t1*(-3 - 3*s1 + t2)) - 
               s1*(Power(s2,2) + t1*(-2 + t2) + s2*(2 - 3*t1 + t2))))/
           ((-1 + s1)*(-s + s1 - t2)) + 
          ((-1 + t1)*(-1 + t2)*
             (9*s2 - 6*Power(s2,2) + 6*t1 + 9*s2*t1 + 2*Power(t1,2) - 
               2*s2*Power(t1,2) - Power(s1,2)*(-2 + s2 + t1) - 
               9*s2*t2 + 3*Power(s2,2)*t2 - 2*Power(t2,2) + 
               Power(s,2)*(-2 + t1 + t2) - 
               s1*(-3 + Power(s2,2) - 3*t1 + 2*Power(t1,2) + 
                  s2*(-2 + t1 - 3*t2) + t2) + 
               s*(-3 - 3*t1 + 2*Power(t1,2) + s1*(-2 + s2 - t2) + 
                  5*t2 - s2*(-8 + t1 + 4*t2))))/
           ((-1 + s2)*(-s + s2 - t1)) - 
          ((-1 + t1)*(-1 + t2)*
             (-3 - 7*s2 + 6*Power(s2,2) + 
               2*Power(s1,2)*(-1 + s2 - t1) + 2*t1 - 11*s2*t1 + 
               Power(t1,2) + 8*s2*t2 - 3*Power(s2,2)*t2 - 4*t1*t2 + 
               5*s2*t1*t2 + Power(t2,2) - s2*Power(t2,2) - 
               Power(s,2)*(-2 + t1 + t2) + 
               s*(1 + Power(t1,2) + t1*(8 + 3*s1 - 4*t2) - 8*t2 + 
                  2*s1*t2 + Power(t2,2) + s2*(-8 - s1 + t1 + 4*t2)) + 
               s1*(-6 + Power(s2,2) - Power(t1,2) + 4*t2 - 
                  s2*(1 + 3*t2) + t1*(-7 + 5*t2))))/
           ((s - s2 + t1)*(s - s1 + t2)))/
        (Power(-1 + t1,2)*Power(-1 + t2,2)) - 
       ((-2*(-1 + s1)*(-s + s1 - t2)*
             (14 + Power(s1,2) + 2*s2 + s*(-2 + s1 - t1) - 11*t1 - 
               s2*t1 + 2*Power(t1,2) + s1*(-11 + s2 + 5*t1 - 3*t2) + 
               12*t2 - 5*t1*t2 + 2*Power(t2,2)))/((-1 + s2)*(-1 + t1)) + 
          (2*(-1 + s1)*(-s + s1 - t2)*
             (-16 - Power(s,2) + Power(s1,2) + 3*s2 + Power(s2,2) + 
               2*t1 - 3*s2*t1 + 2*Power(t1,2) + 
               s1*(3 - 5*s2 + 4*t1 - 3*t2) + 
               s*(3 + s1 + s2 - t1 - t2) + 2*t2 + 4*s2*t2 - 5*t1*t2 + 
               2*Power(t2,2)))/((-1 + s2)*(-s + s2 - t1)) + 
          (2*(-1 + s1)*(-s + s1 - t2)*
             (-9 + Power(s1,2) + s2 + 2*t1 + s*(-5 + s1 + 2*t1 - t2) - 
               4*t2 + t1*t2 + Power(t2,2) - s1*(-3 + t1 + 3*t2)))/
           ((-1 + t1)*(-1 + t2)) - 
          2*(4 - Power(s1,3) + 2*s2 - 2*t1 + 3*t2 - t1*t2 + 
             3*Power(t2,2) + 
             s*(1 + Power(s1,2) - t1 + s1*(-8 + t1 - t2) + 3*t2) + 
             s1*(1 - 2*s2 - 13*t2 - Power(t2,2) + t1*(3 + t2)) + 
             Power(s1,2)*(-t1 + 2*(5 + t2))) + 
          ((s - s1 + t2)*(9 + 2*s2 - t1 + 
               Power(s1,2)*(8 - s2 + t1 - 2*t2) + 5*t2 - s2*t2 + 
               2*t1*t2 - 6*Power(t2,2) + 
               s1*(-5 + s2*(-1 + t2) + 5*t2 - 2*t1*t2 + 
                  2*Power(t2,2)) + s*(8 + t1 - 5*t2 + s1*(4 - t1 + t2))\
))/(-1 + t2) + ((-1 + s1)*(24 - 12*s2 + 12*t1 + 
               Power(s1,2)*(6 - s2 + t1) - 27*t2 + 4*s2*t2 - t1*t2 + 
               11*Power(t2,2) - s2*Power(t2,2) + 
               Power(s,2)*(8 - 2*s1 - t1 + t2) - 
               s1*(-19 + 17*t2 - 2*s2*t2 + t1*(3 + t2)) + 
               s*(-15 + 2*Power(s1,2) + 4*s2 - t1 + 19*t2 - s2*t2 - 
                  t1*t2 + Power(t2,2) + s1*(s2 - 3*(6 + t2)))))/
           (s - s2 + t1))/(Power(-1 + s1,2)*Power(s - s1 + t2,2)))*F(np1))/
   (2.*Power(Pi,2)) + (((2 - 4*s + 6*s1 - 2*s2 + 2*s1*s2 - 4*t1 + 
          4*s1*t1 + 8*t2 + 8*s*t2 + 4*s2*t2 - 4*s1*s2*t2 - 
          (2*(2 + s*(-3 + s1) - s2 - Power(s1,2)*s2 + 2*s1*(3 + s2))*
             (-1 + t2))/(s - s1 + t2) + 
          ((-1 + s1)*(-5 + 2*s2 + 5*t1 + 21*t2 - 3*s2*t2 - 2*t1*t2 + 
               2*Power(t2,2) - 2*s2*Power(t2,2) - 
               s*(-2 + t1 - t2)*(-1 + 2*t2) + 
               s1*(8 - s2 - 5*t1 - 8*t2 + 2*s2*t2 + 2*t1*t2)))/
           (-1 + t1) + ((-1 + s1)*(-1 + t2)*
             (Power(s1,2)*(s2 - t1) + 
               s*(-2*(6 + 3*t1 - 5*t2) + s1*(4 + t1 - t2)) + 
               s1*(19 + 9*t1 + s2*(-8 + t2) - 9*t2) - 
               2*(-6 + 2*t1*(-2 + t2) + (9 + s2)*t2 - 4*Power(t2,2))))/
           ((s - s2 + t1)*(s - s1 + t2)) + 
          ((-1 + s1)*(-1 + t2)*
             (12 + 6*s2 + 7*t1 - 2*s2*t1 - Power(t1,2) + 
               s1*(6 + s2*(-6 + t1) + 2*t1 - Power(t1,2)) - 6*t2 - 
               6*s2*t2 + 3*t1*t2 + s2*t1*t2 + 
               s*(Power(t1,2) + 6*(-1 + t2) - t1*(2 + t2))))/
           ((-1 + s2)*(-s + s2 - t1)) + 
          ((-1 + s1)*(-1 + t2)*
             (2 + 5*t1 - 8*s1*t1 + Power(t1,2) + s1*Power(t1,2) + 
               18*t2 - 3*t1*t2 + 
               s2*(6 + s1*(-2 + t1) - 2*t2 - t1*(2 + t2)) - 
               s*(2 + Power(t1,2) + 6*t2 - t1*(4 + t2))))/
           ((-1 + s2)*(-1 + t1)))/(Power(-1 + s1,2)*Power(-1 + t2,2)) + 
       (20 - 22*s + 6*s1 + 12*s2 + 4*s*s2 + 4*s1*s2 - 8*Power(s2,2) - 
          12*t1 - 6*s*t1 - 4*Power(s,2)*t1 - 2*s1*t1 + 4*s*s1*t1 + 
          16*s2*t1 + 4*s*s2*t1 - 4*s1*s2*t1 - 8*Power(t1,2) - 
          4*s*Power(t1,2) + 4*s1*Power(t1,2) - 10*t2 + 8*s*t2 - 
          12*s2*t2 + 10*t1*t2 - 4*s*t1*t2 + 4*s2*t1*t2 - 
          4*Power(t1,2)*t2 + 
          ((s - s2 + t1)*(s - s1 + t2)*
             (13 + 4*s2 - 5*t1 + Power(s1,2)*(-s2 + t1) - 14*t2 + 
               s2*t2 + t1*t2 + 3*Power(t2,2) - s2*Power(t2,2) + 
               2*t1*Power(t2,2) + Power(s,2)*(-2 + t1 + t2) + 
               s1*(13 - 5*s2 + 8*t1 - 3*t2 + 2*s2*t2 - 3*t1*t2) + 
               s*(-11 + 2*s2 + s1*(10 + s2 - 2*t1 - t2) + 2*t2 - 
                  s2*t2 + 3*t1*t2 + Power(t2,2))))/
           ((-1 + s1)*(-1 + t2)) + 
          ((s - s1 + t2)*(-21 + 23*t1 - 9*s1*t1 - 2*Power(t1,2) + 
               2*s1*Power(t1,2) + 3*t2 + 2*t1*t2 + 
               Power(s2,2)*(4 - 2*s1 + 2*t2) + 
               Power(s,2)*(4 - 2*t1 + 2*t2) + 
               s*(24 + 3*t1 - 2*Power(t1,2) + 2*s1*(-4 + s2 + t1) + 
                  2*s2*(-6 + t1 - 2*t2) + t2 + 2*t1*t2) + 
               s2*(11*s1 - 2*t1*(1 + t2) - 5*(4 + t2))))/(-1 + s2) - 
          ((s - s2 + t1)*(s - s1 + t2)*
             (-12 - 9*s2 + 2*Power(s2,2) + 6*t1 - s2*t1 + 
               s1*(Power(s2,2) + 6*(-1 + t1) - s2*(2 + t1)) + 6*t2 + 
               5*s2*t2 - Power(s2,2)*t2 - 12*t1*t2 + 2*s2*t1*t2 + 
               s*(s2*(-2 + t1 + t2) - 6*(-1 + t1 + t2))))/
           ((-1 + t1)*(-1 + t2)) + 
          (2*(s - s2 + t1)*(-6 - 2*s2 + 2*t1 + Power(s,2)*t1 + 
               Power(s1,2)*t1 + 5*t2 - s2*t2 + t1*t2 - Power(t2,2) + 
               t1*Power(t2,2) - 
               s*(-7 + s2 + t1*(-1 + 2*s1 - 2*t2) + t2) + 
               s1*(-7 + s2 + t2 - t1*(1 + 2*t2))))/(-1 + s1) - 
          ((s - s2 + t1)*(s - s1 + t2)*
             (Power(s2,2)*(-s1 + t2) + 
               s*(12 - 10*t1 + s2*(-4 + t1 - t2) + 6*t2) - 
               s2*(19 + s1*(-8 + t1) - 9*t1 + 9*t2) + 
               2*(-6 - 4*Power(t1,2) - 4*t2 + t1*(9 + s1 + 2*t2))))/
           ((-1 + s2)*(-1 + t1)))/
        (Power(s - s2 + t1,2)*Power(s - s1 + t2,2)) + 
       (2*(3 - 3*s2 + 6*Power(s2,2) - Power(s2,3) - 6*t1 + 6*s2*t1 - 
             2*Power(s2,2)*t1 + 2*Power(t1,2) - s2*Power(t1,2) - 
             s1*(4 - 3*s2 + Power(s2,2) - 3*t1 + 2*s2*t1 + 
                Power(t1,2)) + 
             s*(1 + Power(s2,2) + 2*s2*(-2 + t1) - 3*t1 + 
                Power(t1,2) - 2*t2) + 5*t2 + 2*s2*t2 - t1*t2) + 
          ((-1 + t1)*(-5 + 4*Power(s2,3) - 3*t1 + 3*s1*t1 + 
               4*Power(t1,2) - 2*s1*Power(t1,2) + 
               2*Power(s2,2)*(-8 + 2*s1 + t1) + 
               2*Power(s,2)*(-2 + s2 + t1) + 
               s*(14*s2 - 6*Power(s2,2) + t1 - 4*s2*t1 + 
                  2*Power(t1,2) - 2*s1*(-2 + s2 + t1) - t2) + t2 - 
               2*t1*t2 + s2*(4 - 9*s1 - 2*t1 + 2*s1*t1 - 
                  2*Power(t1,2) + t2)))/(s - s2 + t1) + 
          (2*(-1 + s2)*(-1 + t1)*
             (5 + 2*Power(s1,2) + 2*s2 + Power(s2,2) - 6*t1 + 
               Power(t1,2) + 6*t2 + 2*s2*t2 - 2*t1*t2 + Power(t2,2) - 
               2*s1*(1 + t2) - Power(s,2)*(-2 + s2 - t1 + t2) + 
               s*(-4 + Power(s2,2) + s2*(-1 - t1 + t2) + 
                  s1*(-3 + s2 - t1 + t2))))/((-1 + s1)*(-s + s1 - t2)) \
+ (2*(-1 + s2)*(-1 + t1)*(3 + Power(s1,2) - s2 + Power(s2,2) + 2*t1 + 
               s2*t1 - Power(t1,2) + 2*t2 - s2*t2 + 2*t1*t2 - 
               Power(t2,2) + s1*(-1 + s2 - t1 + t2) - 
               Power(s,2)*(-1 + t1 + t2) + 
               s*(-t1 - t2 + s1*(-2 + t1 + t2) + s2*(-2 + t1 + t2))))/
           ((-1 + s1)*(-1 + t2)) - 
          ((-1 + s2)*(-5 + 4*s2 - 4*Power(s2,2) + 5*t1 - 8*s2*t1 + 
               2*Power(s2,2)*t1 - 6*Power(t1,2) + 2*s2*Power(t1,2) + 
               t2 - 7*s2*t2 + 2*Power(s2,2)*t2 - 2*t1*t2 + 
               2*s2*t1*t2 - s*
                (4 - 7*t1 + 2*Power(t1,2) - 3*t2 + 2*t1*t2 + 
                  2*s2*(-2 + t1 + t2)) + 
               s1*(4 - 7*t1 + 2*Power(t1,2) - 4*t2 + 2*t1*t2 + 
                  s2*(-3 + 2*t1 + 2*t2))))/(-1 + t2) - 
          (2*(-1 + s2)*(-1 + t1)*
             (-3 - Power(s1,2) - 4*s2 - 4*Power(s2,2) + 
               2*Power(s2,3) + 4*t1 + 3*s2*t1 - 3*Power(s2,2)*t1 - 
               Power(t1,2) + s2*Power(t1,2) - t2 - s2*t2 + 
               2*Power(s2,2)*t2 - t1*t2 - s2*t1*t2 + 
               Power(s,2)*(-2 + s2 - t1 + t2) + 
               s1*(-1 + 2*Power(s2,2) + Power(t1,2) - t1*(-2 + t2) + 
                  t2 + s2*(-3 - 3*t1 + 2*t2)) - 
               s*(-1 + 3*Power(s2,2) + 2*t1 + Power(t1,2) + t2 - 
                  t1*t2 + s1*(-4 + s2 - t1 + t2) + 
                  s2*(-5 - 4*t1 + 3*t2))))/((s - s2 + t1)*(s - s1 + t2))\
)/(Power(-1 + s2,2)*Power(-1 + t1,2)) + 
       (-2*(-2 + Power(s,2)*(-1 + s2) - 7*s2 - 6*Power(s2,2) + 
             Power(s2,3) + 7*s2*t1 - Power(s2,2)*t1 - Power(t1,2) + 
             s1*(-2 + 2*s2 + Power(s2,2) - t1 - s2*t1) - 
             s*(-4 + s1 - 6*s2 + s1*s2 + 2*Power(s2,2) + 2*t1 - 
                s2*t1) + t2 + s2*t2) + 
          ((-1 + s2)*(-7 - 2*Power(s1,2) - 21*s2 - 2*Power(s2,2) + 
               2*Power(s2,3) + 6*t1 + 9*s2*t1 - 4*Power(s2,2)*t1 - 
               5*Power(t1,2) + 2*s2*Power(t1,2) + 7*t2 + 8*s2*t2 + 
               Power(s2,2)*t2 - 3*t1*t2 - s2*t1*t2 - 2*Power(t2,2) + 
               Power(s,2)*(-8 + 2*s2 - t1 + t2) + 
               s1*(-6 + Power(s2,2) - t1 + Power(t1,2) - 
                  s2*(1 + 2*t1) + 4*t2) + 
               s*(17 - 4*Power(s2,2) - 8*t1 - Power(t1,2) + 
                  s1*(6 - s2 + t1) + s2*(12 + 5*t1 - 2*t2) - 10*t2 + 
                  t1*t2)))/(s - s1 + t2) + 
          ((s - s2 + t1)*(-2*Power(s,2) - 7*s2 - 8*Power(s2,2) - 
               2*t1 + s2*t1 + 2*Power(s2,2)*t1 + 2*Power(t1,2) + 
               s1*(-2 + 2*s2 + Power(s2,2) + 2*t1 + s2*t1) - 6*t2 - 
               s2*t2 - Power(s2,2)*t2 + 
               s*(4 - 4*t1 + 2*t2 + s2*(6 - t1 + t2))))/(-1 + t1) + 
          ((-1 + s2)*(-s + s2 - t1)*
             (10 + 4*s2 + 2*Power(s2,2) + 2*Power(s1,2)*(1 + s2) - 
               t1 - 4*s2*t1 + Power(t1,2) - t2 - 2*s2*t2 + 2*t1*t2 - 
               s2*t1*t2 + Power(t2,2) + s2*Power(t2,2) - 
               Power(s,2)*(-6 + t1 + t2) - 
               s*(15 + s1*(4 + 3*s2 - t1) - 5*t1 + Power(t1,2) - 
                  s2*(-4 + t2) - 5*t2 + Power(t2,2)) + 
               s1*(4 + 2*Power(s2,2) + Power(t1,2) - 4*t2 - 
                  t1*(2 + t2) - s2*(-2 + t1 + t2))))/
           ((-1 + s1)*(-s + s1 - t2)) + 
          ((-1 + s2)*(-s + s2 - t1)*
             (9 - 3*s2 - 2*Power(s2,2) + 
               2*Power(s1,2)*(1 + s2 - t1) + 10*t1 + 3*s2*t1 - 
               Power(t1,2) - Power(s,2)*(6 + t1 - 3*t2) - 6*t2 - 
               2*s2*t2 + Power(s2,2)*t2 - 2*t1*t2 + s2*t1*t2 - 
               Power(t2,2) - s2*Power(t2,2) + 
               s1*(4 + Power(s2,2) - Power(t1,2) + s2*(-1 + t2) - 
                  2*t2 + t1*(5 + t2)) + 
               s*(-3 - 8*t1 + Power(t1,2) + s2*(8 + t1 - 4*t2) + 
                  2*t2 + Power(t2,2) - s1*(-2 + s2 - 3*t1 + 2*t2))))/
           ((-1 + s1)*(-1 + t2)) - 
          ((-1 + s2)*(-s + s2 - t1)*
             (-9*s2 + 6*Power(s2,2) - 6*t1 - 9*s2*t1 - 2*Power(t1,2) + 
               2*s2*Power(t1,2) + Power(s1,2)*(-2 + s2 + t1) + 
               9*s2*t2 - 3*Power(s2,2)*t2 + 2*Power(t2,2) - 
               Power(s,2)*(-2 + t1 + t2) + 
               s1*(-3 + Power(s2,2) - 3*t1 + 2*Power(t1,2) + 
                  s2*(-2 + t1 - 3*t2) + t2) + 
               s*(3 + 3*t1 - 2*Power(t1,2) - 5*t2 + s1*(2 - s2 + t2) + 
                  s2*(-8 + t1 + 4*t2))))/((-1 + t1)*(-1 + t2)))/
        (Power(-1 + s2,2)*Power(s - s2 + t1,2)) + 
       (2*(1 + 3*s2 - 2*t1 - s2*t1 + 2*Power(t1,2) - 2*t2 + 5*t1*t2 - 
             s2*t1*t2 + 2*Power(t2,2) - s1*(-3 + t2 + t1*t2) + 
             s*(-2 + t1 + t2 + t1*t2)) + 
          ((-1 + t1)*(2*Power(s1,2) + 6*s2 + 2*t1 - 2*Power(t1,2) - 
               t2 - 4*s2*t2 + 7*t1*t2 + 3*Power(t2,2) - 
               s2*Power(t2,2) + s1*(-((-2 + t1)*t2) + s2*(2 + t2)) + 
               s*(-2 + 2*t2 + Power(t2,2) + t1*(2 + t2))))/(-1 + s1) + 
          ((-1 + t2)*(4 + 10*s2 - 2*Power(s2,2) - 7*t1 + 2*s2*t1 + 
               5*Power(t1,2) - 2*s2*Power(t1,2) - 
               s1*(-2 + 2*t1 + Power(t1,2) + s2*(2 + t1)) + 4*s2*t2 + 
               3*t1*t2 + s2*t1*t2 + 2*Power(t2,2) + 
               s*(Power(t1,2) - t1*(-4 + t2) - 2*(3 + t2))))/(-1 + s2) \
+ ((-1 + t1)*(-1 + t2)*(2 + 7*s2 + 2*Power(s2,2) - 4*t1 - 5*s2*t1 + 
               2*Power(t1,2) + 2*Power(s1,2)*(1 + t1) + 5*t2 + 
               7*s2*t2 - Power(s2,2)*t2 - 5*t1*t2 + 3*Power(t2,2) - 
               s2*Power(t2,2) + Power(s,2)*(-2 + t1 + t2) + 
               s*(1 + s2*(2 + s1 - t1) - 3*t2 + Power(t2,2) + 
                  t1*(-3 - 3*s1 + t2)) - 
               s1*(Power(s2,2) + t1*(-2 + t2) + s2*(2 - 3*t1 + t2))))/
           ((-1 + s1)*(-s + s1 - t2)) + 
          ((-1 + t1)*(-1 + t2)*
             (9*s2 - 6*Power(s2,2) + 6*t1 + 9*s2*t1 + 2*Power(t1,2) - 
               2*s2*Power(t1,2) - Power(s1,2)*(-2 + s2 + t1) - 
               9*s2*t2 + 3*Power(s2,2)*t2 - 2*Power(t2,2) + 
               Power(s,2)*(-2 + t1 + t2) - 
               s1*(-3 + Power(s2,2) - 3*t1 + 2*Power(t1,2) + 
                  s2*(-2 + t1 - 3*t2) + t2) + 
               s*(-3 - 3*t1 + 2*Power(t1,2) + s1*(-2 + s2 - t2) + 
                  5*t2 - s2*(-8 + t1 + 4*t2))))/
           ((-1 + s2)*(-s + s2 - t1)) - 
          ((-1 + t1)*(-1 + t2)*
             (-3 - 7*s2 + 6*Power(s2,2) + 
               2*Power(s1,2)*(-1 + s2 - t1) + 2*t1 - 11*s2*t1 + 
               Power(t1,2) + 8*s2*t2 - 3*Power(s2,2)*t2 - 4*t1*t2 + 
               5*s2*t1*t2 + Power(t2,2) - s2*Power(t2,2) - 
               Power(s,2)*(-2 + t1 + t2) + 
               s*(1 + Power(t1,2) + t1*(8 + 3*s1 - 4*t2) - 8*t2 + 
                  2*s1*t2 + Power(t2,2) + s2*(-8 - s1 + t1 + 4*t2)) + 
               s1*(-6 + Power(s2,2) - Power(t1,2) + 4*t2 - 
                  s2*(1 + 3*t2) + t1*(-7 + 5*t2))))/
           ((s - s2 + t1)*(s - s1 + t2)))/
        (Power(-1 + t1,2)*Power(-1 + t2,2)) - 
       ((-2*(-1 + s1)*(-s + s1 - t2)*
             (14 + Power(s1,2) + 2*s2 + s*(-2 + s1 - t1) - 11*t1 - 
               s2*t1 + 2*Power(t1,2) + s1*(-11 + s2 + 5*t1 - 3*t2) + 
               12*t2 - 5*t1*t2 + 2*Power(t2,2)))/((-1 + s2)*(-1 + t1)) + 
          (2*(-1 + s1)*(-s + s1 - t2)*
             (-16 - Power(s,2) + Power(s1,2) + 3*s2 + Power(s2,2) + 
               2*t1 - 3*s2*t1 + 2*Power(t1,2) + 
               s1*(3 - 5*s2 + 4*t1 - 3*t2) + 
               s*(3 + s1 + s2 - t1 - t2) + 2*t2 + 4*s2*t2 - 5*t1*t2 + 
               2*Power(t2,2)))/((-1 + s2)*(-s + s2 - t1)) + 
          (2*(-1 + s1)*(-s + s1 - t2)*
             (-9 + Power(s1,2) + s2 + 2*t1 + s*(-5 + s1 + 2*t1 - t2) - 
               4*t2 + t1*t2 + Power(t2,2) - s1*(-3 + t1 + 3*t2)))/
           ((-1 + t1)*(-1 + t2)) - 
          2*(4 - Power(s1,3) + 2*s2 - 2*t1 + 3*t2 - t1*t2 + 
             3*Power(t2,2) + 
             s*(1 + Power(s1,2) - t1 + s1*(-8 + t1 - t2) + 3*t2) + 
             s1*(1 - 2*s2 - 13*t2 - Power(t2,2) + t1*(3 + t2)) + 
             Power(s1,2)*(-t1 + 2*(5 + t2))) + 
          ((s - s1 + t2)*(9 + 2*s2 - t1 + 
               Power(s1,2)*(8 - s2 + t1 - 2*t2) + 5*t2 - s2*t2 + 
               2*t1*t2 - 6*Power(t2,2) + 
               s1*(-5 + s2*(-1 + t2) + 5*t2 - 2*t1*t2 + 
                  2*Power(t2,2)) + s*(8 + t1 - 5*t2 + s1*(4 - t1 + t2))\
))/(-1 + t2) + ((-1 + s1)*(24 - 12*s2 + 12*t1 + 
               Power(s1,2)*(6 - s2 + t1) - 27*t2 + 4*s2*t2 - t1*t2 + 
               11*Power(t2,2) - s2*Power(t2,2) + 
               Power(s,2)*(8 - 2*s1 - t1 + t2) - 
               s1*(-19 + 17*t2 - 2*s2*t2 + t1*(3 + t2)) + 
               s*(-15 + 2*Power(s1,2) + 4*s2 - t1 + 19*t2 - s2*t2 - 
                  t1*t2 + Power(t2,2) + s1*(s2 - 3*(6 + t2)))))/
           (s - s2 + t1))/(Power(-1 + s1,2)*Power(s - s1 + t2,2)))*F(np2))/
   (2.*Power(Pi,2)) - (((2 - 4*s + 6*s1 - 2*s2 + 2*s1*s2 - 4*t1 + 
          4*s1*t1 + 8*t2 + 8*s*t2 + 4*s2*t2 - 4*s1*s2*t2 - 
          (2*(2 + s*(-3 + s1) - s2 - Power(s1,2)*s2 + 2*s1*(3 + s2))*
             (-1 + t2))/(s - s1 + t2) + 
          ((-1 + s1)*(-5 + 2*s2 + 5*t1 + 21*t2 - 3*s2*t2 - 2*t1*t2 + 
               2*Power(t2,2) - 2*s2*Power(t2,2) - 
               s*(-2 + t1 - t2)*(-1 + 2*t2) + 
               s1*(8 - s2 - 5*t1 - 8*t2 + 2*s2*t2 + 2*t1*t2)))/(-1 + t1) \
+ ((-1 + s1)*(-1 + t2)*(Power(s1,2)*(s2 - t1) + 
               s*(-2*(6 + 3*t1 - 5*t2) + s1*(4 + t1 - t2)) + 
               s1*(19 + 9*t1 + s2*(-8 + t2) - 9*t2) - 
               2*(-6 + 2*t1*(-2 + t2) + (9 + s2)*t2 - 4*Power(t2,2))))/
           ((s - s2 + t1)*(s - s1 + t2)) + 
          ((-1 + s1)*(-1 + t2)*
             (12 + 6*s2 + 7*t1 - 2*s2*t1 - Power(t1,2) + 
               s1*(6 + s2*(-6 + t1) + 2*t1 - Power(t1,2)) - 6*t2 - 
               6*s2*t2 + 3*t1*t2 + s2*t1*t2 + 
               s*(Power(t1,2) + 6*(-1 + t2) - t1*(2 + t2))))/
           ((-1 + s2)*(-s + s2 - t1)) + 
          ((-1 + s1)*(-1 + t2)*
             (2 + 5*t1 - 8*s1*t1 + Power(t1,2) + s1*Power(t1,2) + 
               18*t2 - 3*t1*t2 + 
               s2*(6 + s1*(-2 + t1) - 2*t2 - t1*(2 + t2)) - 
               s*(2 + Power(t1,2) + 6*t2 - t1*(4 + t2))))/
           ((-1 + s2)*(-1 + t1)))/(Power(-1 + s1,2)*Power(-1 + t2,2)) + 
       (20 - 22*s + 6*s1 + 12*s2 + 4*s*s2 + 4*s1*s2 - 8*Power(s2,2) - 
          12*t1 - 6*s*t1 - 4*Power(s,2)*t1 - 2*s1*t1 + 4*s*s1*t1 + 
          16*s2*t1 + 4*s*s2*t1 - 4*s1*s2*t1 - 8*Power(t1,2) - 
          4*s*Power(t1,2) + 4*s1*Power(t1,2) - 10*t2 + 8*s*t2 - 
          12*s2*t2 + 10*t1*t2 - 4*s*t1*t2 + 4*s2*t1*t2 - 
          4*Power(t1,2)*t2 + ((s - s2 + t1)*(s - s1 + t2)*
             (13 + 4*s2 - 5*t1 + Power(s1,2)*(-s2 + t1) - 14*t2 + 
               s2*t2 + t1*t2 + 3*Power(t2,2) - s2*Power(t2,2) + 
               2*t1*Power(t2,2) + Power(s,2)*(-2 + t1 + t2) + 
               s1*(13 - 5*s2 + 8*t1 - 3*t2 + 2*s2*t2 - 3*t1*t2) + 
               s*(-11 + 2*s2 + s1*(10 + s2 - 2*t1 - t2) + 2*t2 - 
                  s2*t2 + 3*t1*t2 + Power(t2,2))))/((-1 + s1)*(-1 + t2)) \
+ ((s - s1 + t2)*(-21 + 23*t1 - 9*s1*t1 - 2*Power(t1,2) + 
               2*s1*Power(t1,2) + 3*t2 + 2*t1*t2 + 
               Power(s2,2)*(4 - 2*s1 + 2*t2) + 
               Power(s,2)*(4 - 2*t1 + 2*t2) + 
               s*(24 + 3*t1 - 2*Power(t1,2) + 2*s1*(-4 + s2 + t1) + 
                  2*s2*(-6 + t1 - 2*t2) + t2 + 2*t1*t2) + 
               s2*(11*s1 - 2*t1*(1 + t2) - 5*(4 + t2))))/(-1 + s2) - 
          ((s - s2 + t1)*(s - s1 + t2)*
             (-12 - 9*s2 + 2*Power(s2,2) + 6*t1 - s2*t1 + 
               s1*(Power(s2,2) + 6*(-1 + t1) - s2*(2 + t1)) + 6*t2 + 
               5*s2*t2 - Power(s2,2)*t2 - 12*t1*t2 + 2*s2*t1*t2 + 
               s*(s2*(-2 + t1 + t2) - 6*(-1 + t1 + t2))))/
           ((-1 + t1)*(-1 + t2)) + 
          (2*(s - s2 + t1)*(-6 - 2*s2 + 2*t1 + Power(s,2)*t1 + 
               Power(s1,2)*t1 + 5*t2 - s2*t2 + t1*t2 - Power(t2,2) + 
               t1*Power(t2,2) - 
               s*(-7 + s2 + t1*(-1 + 2*s1 - 2*t2) + t2) + 
               s1*(-7 + s2 + t2 - t1*(1 + 2*t2))))/(-1 + s1) - 
          ((s - s2 + t1)*(s - s1 + t2)*
             (Power(s2,2)*(-s1 + t2) + 
               s*(12 - 10*t1 + s2*(-4 + t1 - t2) + 6*t2) - 
               s2*(19 + s1*(-8 + t1) - 9*t1 + 9*t2) + 
               2*(-6 - 4*Power(t1,2) - 4*t2 + t1*(9 + s1 + 2*t2))))/
           ((-1 + s2)*(-1 + t1)))/
        (Power(s - s2 + t1,2)*Power(s - s1 + t2,2)) + 
       (2*(3 - 3*s2 + 6*Power(s2,2) - Power(s2,3) - 6*t1 + 6*s2*t1 - 
             2*Power(s2,2)*t1 + 2*Power(t1,2) - s2*Power(t1,2) - 
             s1*(4 - 3*s2 + Power(s2,2) - 3*t1 + 2*s2*t1 + 
                Power(t1,2)) + 
             s*(1 + Power(s2,2) + 2*s2*(-2 + t1) - 3*t1 + 
                Power(t1,2) - 2*t2) + 5*t2 + 2*s2*t2 - t1*t2) + 
          ((-1 + t1)*(-5 + 4*Power(s2,3) - 3*t1 + 3*s1*t1 + 
               4*Power(t1,2) - 2*s1*Power(t1,2) + 
               2*Power(s2,2)*(-8 + 2*s1 + t1) + 
               2*Power(s,2)*(-2 + s2 + t1) + 
               s*(14*s2 - 6*Power(s2,2) + t1 - 4*s2*t1 + 
                  2*Power(t1,2) - 2*s1*(-2 + s2 + t1) - t2) + t2 - 
               2*t1*t2 + s2*(4 - 9*s1 - 2*t1 + 2*s1*t1 - 
                  2*Power(t1,2) + t2)))/(s - s2 + t1) + 
          (2*(-1 + s2)*(-1 + t1)*
             (5 + 2*Power(s1,2) + 2*s2 + Power(s2,2) - 6*t1 + 
               Power(t1,2) + 6*t2 + 2*s2*t2 - 2*t1*t2 + Power(t2,2) - 
               2*s1*(1 + t2) - Power(s,2)*(-2 + s2 - t1 + t2) + 
               s*(-4 + Power(s2,2) + s2*(-1 - t1 + t2) + 
                  s1*(-3 + s2 - t1 + t2))))/((-1 + s1)*(-s + s1 - t2)) + 
          (2*(-1 + s2)*(-1 + t1)*
             (3 + Power(s1,2) - s2 + Power(s2,2) + 2*t1 + s2*t1 - 
               Power(t1,2) + 2*t2 - s2*t2 + 2*t1*t2 - Power(t2,2) + 
               s1*(-1 + s2 - t1 + t2) - Power(s,2)*(-1 + t1 + t2) + 
               s*(-t1 - t2 + s1*(-2 + t1 + t2) + s2*(-2 + t1 + t2))))/
           ((-1 + s1)*(-1 + t2)) - 
          ((-1 + s2)*(-5 + 4*s2 - 4*Power(s2,2) + 5*t1 - 8*s2*t1 + 
               2*Power(s2,2)*t1 - 6*Power(t1,2) + 2*s2*Power(t1,2) + 
               t2 - 7*s2*t2 + 2*Power(s2,2)*t2 - 2*t1*t2 + 
               2*s2*t1*t2 - s*
                (4 - 7*t1 + 2*Power(t1,2) - 3*t2 + 2*t1*t2 + 
                  2*s2*(-2 + t1 + t2)) + 
               s1*(4 - 7*t1 + 2*Power(t1,2) - 4*t2 + 2*t1*t2 + 
                  s2*(-3 + 2*t1 + 2*t2))))/(-1 + t2) - 
          (2*(-1 + s2)*(-1 + t1)*
             (-3 - Power(s1,2) - 4*s2 - 4*Power(s2,2) + 2*Power(s2,3) + 
               4*t1 + 3*s2*t1 - 3*Power(s2,2)*t1 - Power(t1,2) + 
               s2*Power(t1,2) - t2 - s2*t2 + 2*Power(s2,2)*t2 - t1*t2 - 
               s2*t1*t2 + Power(s,2)*(-2 + s2 - t1 + t2) + 
               s1*(-1 + 2*Power(s2,2) + Power(t1,2) - t1*(-2 + t2) + 
                  t2 + s2*(-3 - 3*t1 + 2*t2)) - 
               s*(-1 + 3*Power(s2,2) + 2*t1 + Power(t1,2) + t2 - 
                  t1*t2 + s1*(-4 + s2 - t1 + t2) + s2*(-5 - 4*t1 + 3*t2)\
)))/((s - s2 + t1)*(s - s1 + t2)))/(Power(-1 + s2,2)*Power(-1 + t1,2)) + 
       (-2*(-2 + Power(s,2)*(-1 + s2) - 7*s2 - 6*Power(s2,2) + 
             Power(s2,3) + 7*s2*t1 - Power(s2,2)*t1 - Power(t1,2) + 
             s1*(-2 + 2*s2 + Power(s2,2) - t1 - s2*t1) - 
             s*(-4 + s1 - 6*s2 + s1*s2 + 2*Power(s2,2) + 2*t1 - 
                s2*t1) + t2 + s2*t2) + 
          ((-1 + s2)*(-7 - 2*Power(s1,2) - 21*s2 - 2*Power(s2,2) + 
               2*Power(s2,3) + 6*t1 + 9*s2*t1 - 4*Power(s2,2)*t1 - 
               5*Power(t1,2) + 2*s2*Power(t1,2) + 7*t2 + 8*s2*t2 + 
               Power(s2,2)*t2 - 3*t1*t2 - s2*t1*t2 - 2*Power(t2,2) + 
               Power(s,2)*(-8 + 2*s2 - t1 + t2) + 
               s1*(-6 + Power(s2,2) - t1 + Power(t1,2) - 
                  s2*(1 + 2*t1) + 4*t2) + 
               s*(17 - 4*Power(s2,2) - 8*t1 - Power(t1,2) + 
                  s1*(6 - s2 + t1) + s2*(12 + 5*t1 - 2*t2) - 10*t2 + 
                  t1*t2)))/(s - s1 + t2) + 
          ((s - s2 + t1)*(-2*Power(s,2) - 7*s2 - 8*Power(s2,2) - 
               2*t1 + s2*t1 + 2*Power(s2,2)*t1 + 2*Power(t1,2) + 
               s1*(-2 + 2*s2 + Power(s2,2) + 2*t1 + s2*t1) - 6*t2 - 
               s2*t2 - Power(s2,2)*t2 + 
               s*(4 - 4*t1 + 2*t2 + s2*(6 - t1 + t2))))/(-1 + t1) + 
          ((-1 + s2)*(-s + s2 - t1)*
             (10 + 4*s2 + 2*Power(s2,2) + 2*Power(s1,2)*(1 + s2) - 
               t1 - 4*s2*t1 + Power(t1,2) - t2 - 2*s2*t2 + 2*t1*t2 - 
               s2*t1*t2 + Power(t2,2) + s2*Power(t2,2) - 
               Power(s,2)*(-6 + t1 + t2) - 
               s*(15 + s1*(4 + 3*s2 - t1) - 5*t1 + Power(t1,2) - 
                  s2*(-4 + t2) - 5*t2 + Power(t2,2)) + 
               s1*(4 + 2*Power(s2,2) + Power(t1,2) - 4*t2 - 
                  t1*(2 + t2) - s2*(-2 + t1 + t2))))/
           ((-1 + s1)*(-s + s1 - t2)) + 
          ((-1 + s2)*(-s + s2 - t1)*
             (9 - 3*s2 - 2*Power(s2,2) + 2*Power(s1,2)*(1 + s2 - t1) + 
               10*t1 + 3*s2*t1 - Power(t1,2) - 
               Power(s,2)*(6 + t1 - 3*t2) - 6*t2 - 2*s2*t2 + 
               Power(s2,2)*t2 - 2*t1*t2 + s2*t1*t2 - Power(t2,2) - 
               s2*Power(t2,2) + 
               s1*(4 + Power(s2,2) - Power(t1,2) + s2*(-1 + t2) - 
                  2*t2 + t1*(5 + t2)) + 
               s*(-3 - 8*t1 + Power(t1,2) + s2*(8 + t1 - 4*t2) + 
                  2*t2 + Power(t2,2) - s1*(-2 + s2 - 3*t1 + 2*t2))))/
           ((-1 + s1)*(-1 + t2)) - 
          ((-1 + s2)*(-s + s2 - t1)*
             (-9*s2 + 6*Power(s2,2) - 6*t1 - 9*s2*t1 - 2*Power(t1,2) + 
               2*s2*Power(t1,2) + Power(s1,2)*(-2 + s2 + t1) + 
               9*s2*t2 - 3*Power(s2,2)*t2 + 2*Power(t2,2) - 
               Power(s,2)*(-2 + t1 + t2) + 
               s1*(-3 + Power(s2,2) - 3*t1 + 2*Power(t1,2) + 
                  s2*(-2 + t1 - 3*t2) + t2) + 
               s*(3 + 3*t1 - 2*Power(t1,2) - 5*t2 + s1*(2 - s2 + t2) + 
                  s2*(-8 + t1 + 4*t2))))/((-1 + t1)*(-1 + t2)))/
        (Power(-1 + s2,2)*Power(s - s2 + t1,2)) + 
       (2*(1 + 3*s2 - 2*t1 - s2*t1 + 2*Power(t1,2) - 2*t2 + 5*t1*t2 - 
             s2*t1*t2 + 2*Power(t2,2) - s1*(-3 + t2 + t1*t2) + 
             s*(-2 + t1 + t2 + t1*t2)) + 
          ((-1 + t1)*(2*Power(s1,2) + 6*s2 + 2*t1 - 2*Power(t1,2) - 
               t2 - 4*s2*t2 + 7*t1*t2 + 3*Power(t2,2) - 
               s2*Power(t2,2) + s1*(-((-2 + t1)*t2) + s2*(2 + t2)) + 
               s*(-2 + 2*t2 + Power(t2,2) + t1*(2 + t2))))/(-1 + s1) + 
          ((-1 + t2)*(4 + 10*s2 - 2*Power(s2,2) - 7*t1 + 2*s2*t1 + 
               5*Power(t1,2) - 2*s2*Power(t1,2) - 
               s1*(-2 + 2*t1 + Power(t1,2) + s2*(2 + t1)) + 4*s2*t2 + 
               3*t1*t2 + s2*t1*t2 + 2*Power(t2,2) + 
               s*(Power(t1,2) - t1*(-4 + t2) - 2*(3 + t2))))/(-1 + s2) + 
          ((-1 + t1)*(-1 + t2)*
             (2 + 7*s2 + 2*Power(s2,2) - 4*t1 - 5*s2*t1 + 
               2*Power(t1,2) + 2*Power(s1,2)*(1 + t1) + 5*t2 + 
               7*s2*t2 - Power(s2,2)*t2 - 5*t1*t2 + 3*Power(t2,2) - 
               s2*Power(t2,2) + Power(s,2)*(-2 + t1 + t2) + 
               s*(1 + s2*(2 + s1 - t1) - 3*t2 + Power(t2,2) + 
                  t1*(-3 - 3*s1 + t2)) - 
               s1*(Power(s2,2) + t1*(-2 + t2) + s2*(2 - 3*t1 + t2))))/
           ((-1 + s1)*(-s + s1 - t2)) + 
          ((-1 + t1)*(-1 + t2)*
             (9*s2 - 6*Power(s2,2) + 6*t1 + 9*s2*t1 + 2*Power(t1,2) - 
               2*s2*Power(t1,2) - Power(s1,2)*(-2 + s2 + t1) - 
               9*s2*t2 + 3*Power(s2,2)*t2 - 2*Power(t2,2) + 
               Power(s,2)*(-2 + t1 + t2) - 
               s1*(-3 + Power(s2,2) - 3*t1 + 2*Power(t1,2) + 
                  s2*(-2 + t1 - 3*t2) + t2) + 
               s*(-3 - 3*t1 + 2*Power(t1,2) + s1*(-2 + s2 - t2) + 
                  5*t2 - s2*(-8 + t1 + 4*t2))))/
           ((-1 + s2)*(-s + s2 - t1)) - 
          ((-1 + t1)*(-1 + t2)*
             (-3 - 7*s2 + 6*Power(s2,2) + 
               2*Power(s1,2)*(-1 + s2 - t1) + 2*t1 - 11*s2*t1 + 
               Power(t1,2) + 8*s2*t2 - 3*Power(s2,2)*t2 - 4*t1*t2 + 
               5*s2*t1*t2 + Power(t2,2) - s2*Power(t2,2) - 
               Power(s,2)*(-2 + t1 + t2) + 
               s*(1 + Power(t1,2) + t1*(8 + 3*s1 - 4*t2) - 8*t2 + 
                  2*s1*t2 + Power(t2,2) + s2*(-8 - s1 + t1 + 4*t2)) + 
               s1*(-6 + Power(s2,2) - Power(t1,2) + 4*t2 - 
                  s2*(1 + 3*t2) + t1*(-7 + 5*t2))))/
           ((s - s2 + t1)*(s - s1 + t2)))/
        (Power(-1 + t1,2)*Power(-1 + t2,2)) - 
       ((-2*(-1 + s1)*(-s + s1 - t2)*
             (14 + Power(s1,2) + 2*s2 + s*(-2 + s1 - t1) - 11*t1 - 
               s2*t1 + 2*Power(t1,2) + s1*(-11 + s2 + 5*t1 - 3*t2) + 
               12*t2 - 5*t1*t2 + 2*Power(t2,2)))/((-1 + s2)*(-1 + t1)) + 
          (2*(-1 + s1)*(-s + s1 - t2)*
             (-16 - Power(s,2) + Power(s1,2) + 3*s2 + Power(s2,2) + 
               2*t1 - 3*s2*t1 + 2*Power(t1,2) + 
               s1*(3 - 5*s2 + 4*t1 - 3*t2) + 
               s*(3 + s1 + s2 - t1 - t2) + 2*t2 + 4*s2*t2 - 5*t1*t2 + 
               2*Power(t2,2)))/((-1 + s2)*(-s + s2 - t1)) + 
          (2*(-1 + s1)*(-s + s1 - t2)*
             (-9 + Power(s1,2) + s2 + 2*t1 + s*(-5 + s1 + 2*t1 - t2) - 
               4*t2 + t1*t2 + Power(t2,2) - s1*(-3 + t1 + 3*t2)))/
           ((-1 + t1)*(-1 + t2)) - 
          2*(4 - Power(s1,3) + 2*s2 - 2*t1 + 3*t2 - t1*t2 + 
             3*Power(t2,2) + s*
              (1 + Power(s1,2) - t1 + s1*(-8 + t1 - t2) + 3*t2) + 
             s1*(1 - 2*s2 - 13*t2 - Power(t2,2) + t1*(3 + t2)) + 
             Power(s1,2)*(-t1 + 2*(5 + t2))) + 
          ((s - s1 + t2)*(9 + 2*s2 - t1 + 
               Power(s1,2)*(8 - s2 + t1 - 2*t2) + 5*t2 - s2*t2 + 
               2*t1*t2 - 6*Power(t2,2) + 
               s1*(-5 + s2*(-1 + t2) + 5*t2 - 2*t1*t2 + 
                  2*Power(t2,2)) + s*(8 + t1 - 5*t2 + s1*(4 - t1 + t2)))\
)/(-1 + t2) + ((-1 + s1)*(24 - 12*s2 + 12*t1 + 
               Power(s1,2)*(6 - s2 + t1) - 27*t2 + 4*s2*t2 - t1*t2 + 
               11*Power(t2,2) - s2*Power(t2,2) + 
               Power(s,2)*(8 - 2*s1 - t1 + t2) - 
               s1*(-19 + 17*t2 - 2*s2*t2 + t1*(3 + t2)) + 
               s*(-15 + 2*Power(s1,2) + 4*s2 - t1 + 19*t2 - s2*t2 - 
                  t1*t2 + Power(t2,2) + s1*(s2 - 3*(6 + t2)))))/
           (s - s2 + t1))/(Power(-1 + s1,2)*Power(s - s1 + t2,2)))*( G(Sqrt(1 - 4/Power(-s + s1 + s2,2)),np1,(np1 - (2*np2)/(s - s1 - s2))/Sqrt(1 - 4/Power(-s + s1 + s2,2)))+Log(4.) ) )/
   Power(Pi,2);
   return a;
};
