#include <stdio.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_pow_int.h>
#include "nlo_functions.h"

#define Pi M_PI





double ir(double ss,double ss1,double ss2,double tt1,double tt2,double iEE,double inp1,double inp2)
{
    long double s=ss, s1=ss1, s2=ss2, t1=tt1, t2=tt2;
    long double EE=iEE, np1=inp1, np2=inp2;
    long double aG=1/Power(4.*Pi,2);
    long double a=2*aG*((4*((-2*(2 + s*(-3 + s1) - s2 - Power(s1,2)*s2 + 2*s1*(3 + s2))*
             (-1 + t2))/(s - s1 + t2) + 
          (2*(-1 + s1)*(-1 + t2)*
             (5 + Power(s1,2) - 2*s2 + 2*Power(s2,2) + 6*t1 - 
               2*s2*t1 + Power(t1,2) + 2*s1*(1 + t1) + 
               s*(-4 + Power(s1,2) + s2*(-3 + t1 - t2) + 
                  s1*(-1 + s2 + t1 - t2)) - 
               Power(s,2)*(-2 + s1 + t1 - t2) - 6*t2 - 2*t1*t2 + 
               Power(t2,2)))/((-1 + s2)*(-s + s2 - t1)) - 
          (2*(-1 + s1)*(2 + s2 - 4*t1 + s1*(-3 + t2) - 
               s*Power(-1 + t2,2) - 6*t2 - 2*s2*t2 + 3*t1*t2 + 
               s2*Power(t2,2) - t1*Power(t2,2)))/(-1 + t1) + 
          2*(7 - s2 + 3*t1 + 2*s*(-1 + t2) - 8*t2 + s2*t2 - t1*t2 + 
             3*Power(t2,2) + 
             s1*(-1 + s2 + t1*(-3 + t2) + 8*t2 - s2*t2 - Power(t2,2))) + 
          ((-1 + s1)*(-1 + t2)*
             (12 + s2 - 2*Power(s2,2) - 3*t1 + 3*s2*t1 + Power(t1,2) + 
               Power(s1,2)*(-2 + s2 + t1) - 3*t2 - s2*t2 + 
               Power(s2,2)*t2 + 2*t1*t2 - s2*t1*t2 + Power(t2,2) + 
               s2*Power(t2,2) - Power(s,2)*(-2 + t1 + t2) + 
               s1*(1 + Power(s2,2) + Power(t1,2) + 3*t2 - 
                  t1*(1 + t2)) + 
               s*(-13 + t1 - Power(t1,2) + s2*(2 + t1) + t2 - 
                  Power(t2,2) + s1*(2 - s2 + t2))))/
           ((-1 + s2)*(-1 + t1)) - 
          ((-1 + s1)*(-1 + t2)*
             (-3 - 3*s2 + 2*Power(s2,2) - 5*t1 - s2*t1 - 
               2*Power(t1,2) + Power(s1,2)*(-6 + s2 + t1) + 6*t2 - 
               Power(s2,2)*t2 + 5*t1*t2 - 3*Power(t2,2) + 
               s2*Power(t2,2) - 2*Power(s,2)*(1 + t2) - 
               s1*(7 + s2 + Power(s2,2) + 4*t1 - 3*s2*t1 - 5*t2 + 
                  4*s2*t2 - t1*t2) - 
               s*(-4 + s1*(-4 + t1 - 3*t2) + t2 - 3*s2*t2 + 
                  Power(t2,2) + t1*(1 + s2 + t2))))/
           ((s - s2 + t1)*(s - s1 + t2)))*
        (F(np1) + F(np2) - 
          2*(G(Sqrt(1 - 4/Power(-s + s1 + s2,2)),np1,
              (np1 - (2*np2)/(s - s1 - s2))/
               Sqrt(1 - 4/Power(-s + s1 + s2,2))) + Log(4))) + 
       8*(-18 + 8*s + 10*s1 + 6*s2 - 6*s1*s2 - 14*t1 + 14*s1*t1 + 40*t2 - 
          8*s*t2 - 24*s1*t2 - 6*s2*t2 + 6*s1*s2*t2 + 6*t1*t2 - 
          6*s1*t1*t2 - 14*Power(t2,2) + 6*s1*Power(t2,2) - 
          ((-1 + s1)*(-11 - 2*s2 + 7*t1 + 14*t2 + 5*s2*t2 - 7*t1*t2 - 
               7*Power(t2,2) - 3*s2*Power(t2,2) + 4*t1*Power(t2,2) + 
               s*(-1 + t2)*(-2 + t1 + 3*t2) + 
               s1*(4 + s2 + t1*(-5 + t2) - s2*t2)))/(-1 + t1) + 
          ((-1 + s1)*(-1 + t2)*
             (-1 + 2*s2 - 6*Power(s2,2) + 
               Power(s1,2)*(-2 + 3*s2 - t1) + 6*s2*t1 - 3*Power(t1,2) + 
               4*t2 - 3*s2*t2 + 2*Power(s2,2)*t2 + 2*t1*t2 - s2*t1*t2 - 
               3*Power(t2,2) + s2*Power(t2,2) + 
               Power(s,2)*(-4 + 2*s1 + t1 + t2) + 
               s1*(1 + 2*Power(s2,2) + Power(t1,2) - s2*(3 + t1) + 
                  t2 - t1*t2) - 
               s*(-3 + 2*Power(s1,2) + 4*t1 + Power(t1,2) + 
                  s1*(-2 + 3*s2 - t2) + 3*s2*(-2 + t2) + 2*t2 - 
                  4*t1*t2 + Power(t2,2))))/((-1 + s2)*(-s + s2 - t1)) - 
          ((-1 + t2)*(5 + 2*s2 - t1 + 
               Power(s1,2)*(4 + 3*s2 + t1 - 2*t2) + 5*t2 - s2*t2 + 
               2*t1*t2 - 6*Power(t2,2) + 
               s*(4 + t1 - s1*t1 + (-5 + s1)*t2) + 
               s1*(-13 + s2*(-5 + t2) + (5 - 2*t1)*t2 + 2*Power(t2,2))))/
           (s - s1 + t2) - ((-1 + s1)*(-1 + t2)*
             (12 - s2 - 2*Power(s2,2) - 11*t1 + s2*t1 + Power(t1,2) + 
               Power(s1,2)*(-2 + s2 + t1) - 11*t2 + 3*s2*t2 + 
               Power(s2,2)*t2 + 6*t1*t2 - s2*t1*t2 + Power(t2,2) + 
               s2*Power(t2,2) - 2*Power(s,2)*(-2 + t1 + t2) + 
               s1*(-1 + Power(s2,2) + Power(t1,2) - t1*(-3 + t2) + 
                  t2 - 2*s2*(t1 + t2)) + 
               s*(-14 + 6*t1 - Power(t1,2) + 6*t2 - 2*t1*t2 - 
                  Power(t2,2) + s2*(2 + t1 + t2) + 
                  s1*(2 - 2*s2 + t1 + t2))))/((-1 + s2)*(-1 + t1)) + 
          ((-1 + s1)*(-1 + t2)*
             (4 + s2 + 2*Power(s2,2) - 3*t1 - s2*t1 - 3*Power(t1,2) + 
               2*Power(s1,2)*(-2 + s2 + t1) - 6*t2 + s2*t2 - 
               Power(s2,2)*t2 + 3*t1*t2 - s2*t1*t2 + 2*Power(t2,2) + 
               2*t1*Power(t2,2) - 
               s*(6 + 2*s1*(-3 + t1) + 2*t1 + Power(t1,2) + 
                  s2*(2 + t1 - t2) - 4*t2 - 3*t1*t2) + 
               s1*(2 - Power(s2,2) + 2*t1 + Power(t1,2) + 2*t2 - 
                  4*t1*t2 + s2*(4*t1 - 2*(2 + t2)))))/
           ((s - s2 + t1)*(s - s1 + t2)) - 
          2*((-2*(2 + s*(-3 + s1) - s2 - Power(s1,2)*s2 + 
                  2*s1*(3 + s2))*(-1 + t2))/(s - s1 + t2) + 
             (2*(-1 + s1)*(-1 + t2)*
                (5 + Power(s1,2) - 2*s2 + 2*Power(s2,2) + 6*t1 - 
                  2*s2*t1 + Power(t1,2) + 2*s1*(1 + t1) + 
                  s*(-4 + Power(s1,2) + s2*(-3 + t1 - t2) + 
                     s1*(-1 + s2 + t1 - t2)) - 
                  Power(s,2)*(-2 + s1 + t1 - t2) - 6*t2 - 2*t1*t2 + 
                  Power(t2,2)))/((-1 + s2)*(-s + s2 - t1)) - 
             (2*(-1 + s1)*(2 + s2 - 4*t1 + s1*(-3 + t2) - 
                  s*Power(-1 + t2,2) - 6*t2 - 2*s2*t2 + 3*t1*t2 + 
                  s2*Power(t2,2) - t1*Power(t2,2)))/(-1 + t1) + 
             2*(7 - s2 + 3*t1 + 2*s*(-1 + t2) - 8*t2 + s2*t2 - t1*t2 + 
                3*Power(t2,2) + 
                s1*(-1 + s2 + t1*(-3 + t2) + 8*t2 - s2*t2 - Power(t2,2))\
) + ((-1 + s1)*(-1 + t2)*(12 + s2 - 2*Power(s2,2) - 3*t1 + 3*s2*t1 + 
                  Power(t1,2) + Power(s1,2)*(-2 + s2 + t1) - 3*t2 - 
                  s2*t2 + Power(s2,2)*t2 + 2*t1*t2 - s2*t1*t2 + 
                  Power(t2,2) + s2*Power(t2,2) - 
                  Power(s,2)*(-2 + t1 + t2) + 
                  s1*(1 + Power(s2,2) + Power(t1,2) + 3*t2 - 
                     t1*(1 + t2)) + 
                  s*(-13 + t1 - Power(t1,2) + s2*(2 + t1) + t2 - 
                     Power(t2,2) + s1*(2 - s2 + t2))))/
              ((-1 + s2)*(-1 + t1)) - 
             ((-1 + s1)*(-1 + t2)*
                (-3 - 3*s2 + 2*Power(s2,2) - 5*t1 - s2*t1 - 
                  2*Power(t1,2) + Power(s1,2)*(-6 + s2 + t1) + 6*t2 - 
                  Power(s2,2)*t2 + 5*t1*t2 - 3*Power(t2,2) + 
                  s2*Power(t2,2) - 2*Power(s,2)*(1 + t2) - 
                  s1*(7 + s2 + Power(s2,2) + 4*t1 - 3*s2*t1 - 5*t2 + 
                     4*s2*t2 - t1*t2) - 
                  s*(-4 + s1*(-4 + t1 - 3*t2) + t2 - 3*s2*t2 + 
                     Power(t2,2) + t1*(1 + s2 + t2))))/
              ((s - s2 + t1)*(s - s1 + t2)))*Log(EE))*
        (1 + (1 - 2/(2 - s + s1 + s2))*R2(2 + s - s1 - s2)))/
     (Power(-1 + s1,2)*Power(-1 + t2,2)) + 
    (4*((2*(-1 + t2)*(-2 + 3*s2 + s*Power(-1 + t1,2) - 
               s1*Power(-1 + t1,2) + 6*t1 - s2*t1 + 4*t2 - 3*t1*t2 + 
               Power(t1,2)*t2))/(-1 + s2) + 
          (2*(-1 + t1)*(-2 - s2 + 4*t1 - s1*(-3 + t2) + 
               s*Power(-1 + t2,2) + 6*t2 + 2*s2*t2 - 3*t1*t2 - 
               s2*Power(t2,2) + t1*Power(t2,2)))/(-1 + s1) + 
          ((-1 + t1)*(-1 + t2)*
             (2 + 2*Power(s2,2) + 5*t1 + 3*Power(t1,2) - 
               Power(s1,2)*(-2 + s2 + t1) - 4*t2 + 2*s2*t2 + 
               2*Power(s2,2)*t2 - 5*t1*t2 - s2*t1*t2 + 2*Power(t2,2) + 
               Power(s,2)*(-2 + t1 + t2) - 
               s1*(-7 - 7*t1 + Power(t1,2) + s2*(2 + t1 - 3*t2) + 
                  5*t2) + s*(1 + Power(t1,2) + s1*(2 + s2 - t2) + 
                  t1*(-3 + t2) - 3*t2 - 3*s2*t2)))/
           ((-1 + s2)*(-s + s2 - t1)) + 
          2*(s*(-1 + t1)*(-1 + t2) - 
             2*(-2 + s1*(-1 + t1) + t1 + s2*(-1 + t2) + t2 - 4*t1*t2)) + 
          ((-1 + t1)*(-1 + t2)*
             (2 + 7*s2 + 2*Power(s2,2) - 4*t1 - 5*s2*t1 + 
               2*Power(t1,2) + 2*Power(s1,2)*(1 + t1) + 5*t2 + 
               7*s2*t2 - Power(s2,2)*t2 - 5*t1*t2 + 3*Power(t2,2) - 
               s2*Power(t2,2) + Power(s,2)*(-2 + t1 + t2) + 
               s*(1 + s2*(2 + s1 - t1) - 3*t2 + Power(t2,2) + 
                  t1*(-3 - 3*s1 + t2)) - 
               s1*(Power(s2,2) + t1*(-2 + t2) + s2*(2 - 3*t1 + t2))))/
           ((-1 + s1)*(-s + s1 - t2)) - 
          (2*(-1 + t1)*(-1 + t2)*
             (-2*Power(s,2) - 3*s2 - 2*Power(s2,2) + 2*t1 + 2*s2*t1 + 
               Power(s2,2)*t1 - 2*Power(t1,2) - s2*Power(t1,2) + 2*t2 - 
               Power(s2,2)*t2 + 2*s2*t1*t2 - 2*Power(t2,2) - 
               s2*Power(t2,2) + Power(s1,2)*(-2 + s2 - t1 + t2) + 
               s1*(-3 - 4*s2 + Power(s2,2) - Power(t1,2) + 2*t2 + 
                  2*t1*t2 - Power(t2,2)) + 
               s*(-1 - 2*t1 + Power(t1,2) - 2*t2 - 2*t1*t2 + 
                  Power(t2,2) + s2*(3 - t1 + t2) - 
                  s1*(-3 + s2 - t1 + t2))))/((s - s2 + t1)*(s - s1 + t2))\
)*(F(np1) + F(np2) - 
          2*(G(Sqrt(1 - 4/Power(-s + s1 + s2,2)),np1,
              (np1 - (2*np2)/(s - s1 - s2))/
               Sqrt(1 - 4/Power(-s + s1 + s2,2))) + Log(4))) + 
       8*(-4 - 6*s - 8*s1 - 8*s2 + 12*t1 + 6*s*t1 + 8*s1*t1 + 12*t2 + 
          6*s*t2 + 8*s2*t2 - 20*t1*t2 - 6*s*t1*t2 - 
          ((-1 + t2)*(-11 + 4*s2 + 14*t1 - 7*Power(t1,2) - 
               s1*(-1 + t1)*(-2 + s2 + 3*t1) + 7*t2 - 5*s2*t2 - 
               7*t1*t2 + s2*t1*t2 + 4*Power(t1,2)*t2 + 
               s*(-1 + t1)*(-2 + 3*t1 + t2)))/(-1 + s2) - 
          ((-1 + t1)*(-11 - 2*s2 + 7*t1 + 14*t2 + 5*s2*t2 - 7*t1*t2 - 
               7*Power(t2,2) - 3*s2*Power(t2,2) + 4*t1*Power(t2,2) + 
               s*(-1 + t2)*(-2 + t1 + 3*t2) + 
               s1*(4 + s2 + t1*(-5 + t2) - s2*t2)))/(-1 + s1) - 
          ((-1 + t1)*(-1 + t2)*
             (1 + 4*s2 + 4*Power(s2,2) - 4*t1 - 8*s2*t1 + 
               3*Power(t1,2) + 5*t2 + 7*s2*t2 - 2*Power(s2,2)*t2 - 
               7*t1*t2 + s2*t1*t2 - 4*Power(t2,2) - 2*s2*Power(t2,2) + 
               2*t1*Power(t2,2) + Power(s,2)*(-2 + t1 + t2) + 
               s*(3 + s1*(2 + s2 - t1) - 2*t1 + Power(t1,2) - 8*t2 + 
                  t1*t2 + 2*Power(t2,2) + s2*(2 - 2*t1 + t2)) - 
               s1*(8 + 2*Power(s2,2) - 7*t1 + Power(t1,2) - 6*t2 + 
                  s2*(5 - 5*t1 + 2*t2))))/((-1 + s1)*(-s + s1 - t2)) - 
          ((-1 + t1)*(-1 + t2)*
             (1 - 8*s2 + 5*t1 + 6*s2*t1 - 4*Power(t1,2) - 
               2*Power(s1,2)*(-2 + s2 + t1) - 4*t2 + 7*s2*t2 - 
               7*t1*t2 + 2*Power(t1,2)*t2 + 3*Power(t2,2) - 
               s2*Power(t2,2) + Power(s,2)*(-2 + t1 + t2) + 
               s*(3 - 8*t1 + 2*Power(t1,2) + s1*(2 + s2 + t1 - 2*t2) - 
                  s2*(-2 + t2) - 2*t2 + t1*t2 + Power(t2,2)) + 
               s1*(4 - 2*Power(t1,2) - 8*t2 + t1*(7 + t2) + 
                  s2*(-5 - 2*t1 + 5*t2))))/((-1 + s2)*(-s + s2 - t1)) - 
          ((-1 + t1)*(-1 + t2)*
             (s2 + 4*Power(s2,2) - 2*t1 - 5*s2*t1 - 2*Power(s2,2)*t1 + 
               2*Power(t1,2) + 2*s2*Power(t1,2) - 2*t2 - s2*t2 - 
               Power(s2,2)*t2 + 2*s2*t1*t2 + 2*Power(t2,2) - 
               2*Power(s,2)*(-3 + t1 + t2) - 
               Power(s1,2)*(-4 + s2 + t1 + 2*t2) - 
               s1*(-1 + Power(s2,2) + t1 + 5*t2 - 2*t1*t2 - 
                  2*Power(t2,2) + s2*(-4 + t1 + t2)) + 
               s*(3*s1*(-2 + t1 + t2) + 3*s2*(-2 + t1 + t2) + 
                  4*(-1 + t1 + t2 - t1*t2))))/
           ((s - s2 + t1)*(s - s1 + t2)) - 
          2*((2*(-1 + t2)*(-2 + 3*s2 + s*Power(-1 + t1,2) - 
                  s1*Power(-1 + t1,2) + 6*t1 - s2*t1 + 4*t2 - 
                  3*t1*t2 + Power(t1,2)*t2))/(-1 + s2) + 
             (2*(-1 + t1)*(-2 - s2 + 4*t1 - s1*(-3 + t2) + 
                  s*Power(-1 + t2,2) + 6*t2 + 2*s2*t2 - 3*t1*t2 - 
                  s2*Power(t2,2) + t1*Power(t2,2)))/(-1 + s1) + 
             ((-1 + t1)*(-1 + t2)*
                (2 + 2*Power(s2,2) + 5*t1 + 3*Power(t1,2) - 
                  Power(s1,2)*(-2 + s2 + t1) - 4*t2 + 2*s2*t2 + 
                  2*Power(s2,2)*t2 - 5*t1*t2 - s2*t1*t2 + 
                  2*Power(t2,2) + Power(s,2)*(-2 + t1 + t2) - 
                  s1*(-7 - 7*t1 + Power(t1,2) + s2*(2 + t1 - 3*t2) + 
                     5*t2) + 
                  s*(1 + Power(t1,2) + s1*(2 + s2 - t2) + 
                     t1*(-3 + t2) - 3*t2 - 3*s2*t2)))/
              ((-1 + s2)*(-s + s2 - t1)) + 
             2*(s*(-1 + t1)*(-1 + t2) - 
                2*(-2 + s1*(-1 + t1) + t1 + s2*(-1 + t2) + t2 - 4*t1*t2)\
) + ((-1 + t1)*(-1 + t2)*(2 + 7*s2 + 2*Power(s2,2) - 4*t1 - 5*s2*t1 + 
                  2*Power(t1,2) + 2*Power(s1,2)*(1 + t1) + 5*t2 + 
                  7*s2*t2 - Power(s2,2)*t2 - 5*t1*t2 + 3*Power(t2,2) - 
                  s2*Power(t2,2) + Power(s,2)*(-2 + t1 + t2) + 
                  s*(1 + s2*(2 + s1 - t1) - 3*t2 + Power(t2,2) + 
                     t1*(-3 - 3*s1 + t2)) - 
                  s1*(Power(s2,2) + t1*(-2 + t2) + s2*(2 - 3*t1 + t2)))\
)/((-1 + s1)*(-s + s1 - t2)) - 
             (2*(-1 + t1)*(-1 + t2)*
                (-2*Power(s,2) - 3*s2 - 2*Power(s2,2) + 2*t1 + 
                  2*s2*t1 + Power(s2,2)*t1 - 2*Power(t1,2) - 
                  s2*Power(t1,2) + 2*t2 - Power(s2,2)*t2 + 2*s2*t1*t2 - 
                  2*Power(t2,2) - s2*Power(t2,2) + 
                  Power(s1,2)*(-2 + s2 - t1 + t2) + 
                  s1*(-3 - 4*s2 + Power(s2,2) - Power(t1,2) + 2*t2 + 
                     2*t1*t2 - Power(t2,2)) + 
                  s*(-1 - 2*t1 + Power(t1,2) - 2*t2 - 2*t1*t2 + 
                     Power(t2,2) + s2*(3 - t1 + t2) - 
                     s1*(-3 + s2 - t1 + t2))))/
              ((s - s2 + t1)*(s - s1 + t2)))*Log(EE))*
        (1 + (1 - 2/(2 - s + s1 + s2))*R2(2 + s - s1 - s2)))/
     (Power(-1 + t1,2)*Power(-1 + t2,2)) + 
    (4*((-2*(2 + s*(-3 + s2) - s1*Power(-1 + s2,2) + 6*s2)*
             (s - s2 + t1))/(-1 + t1) + 
          (2*(-1 + s2)*(-s + s2 - t1)*
             (5 + Power(s1,2) - 2*s2 + 2*Power(s2,2) + 6*t1 - 
               2*s2*t1 + Power(t1,2) + 2*s1*(1 + t1) + 
               s*(-4 + Power(s1,2) + s2*(-3 + t1 - t2) + 
                  s1*(-1 + s2 + t1 - t2)) - 
               Power(s,2)*(-2 + s1 + t1 - t2) - 6*t2 - 2*t1*t2 + 
               Power(t2,2)))/((-1 + s1)*(-1 + t2)) + 
          ((-1 + s2)*(-s + s2 - t1)*
             (2 + 2*Power(s2,2) + 5*t1 + 3*Power(t1,2) - 
               Power(s1,2)*(-2 + s2 + t1) - 4*t2 + 2*s2*t2 + 
               2*Power(s2,2)*t2 - 5*t1*t2 - s2*t1*t2 + 2*Power(t2,2) + 
               Power(s,2)*(-2 + t1 + t2) - 
               s1*(-7 - 7*t1 + Power(t1,2) + s2*(2 + t1 - 3*t2) + 
                  5*t2) + s*(1 + Power(t1,2) + s1*(2 + s2 - t2) + 
                  t1*(-3 + t2) - 3*t2 - 3*s2*t2)))/((-1 + t1)*(-1 + t2)) \
+ (2*(-1 + s2)*(-6 - 7*s2 + s1*(-2 + s2 - t1) + 5*t1 + s2*t1 - 
               Power(t1,2) + 2*t2 + Power(s,2)*t2 - s2*t2 + 
               Power(s2,2)*t2 + t1*t2 - 2*s2*t1*t2 + Power(t1,2)*t2 - 
               s*(-7 + s1 + t1 - t2 + 2*s2*t2 - 2*t1*t2)))/(s - s1 + t2) \
+ 2*(4 - 2*s1*(-1 + s2) + s2 + 10*Power(s2,2) - Power(s2,3) + 3*t1 - 
             13*s2*t1 + 2*Power(s2,2)*t1 + 3*Power(t1,2) - 
             s2*Power(t1,2) - 2*t2 + 3*s2*t2 - Power(s2,2)*t2 - t1*t2 + 
             s2*t1*t2 + s*(1 + Power(s2,2) + 3*t1 - t2 + 
                s2*(-8 - t1 + t2))) + 
          ((-1 + s2)*(-s + s2 - t1)*
             (10 + 4*s2 + 2*Power(s2,2) + 2*Power(s1,2)*(1 + s2) - t1 - 
               4*s2*t1 + Power(t1,2) - t2 - 2*s2*t2 + 2*t1*t2 - 
               s2*t1*t2 + Power(t2,2) + s2*Power(t2,2) - 
               Power(s,2)*(-6 + t1 + t2) - 
               s*(15 + s1*(4 + 3*s2 - t1) - 5*t1 + Power(t1,2) - 
                  s2*(-4 + t2) - 5*t2 + Power(t2,2)) + 
               s1*(4 + 2*Power(s2,2) + Power(t1,2) - 4*t2 - 
                  t1*(2 + t2) - s2*(-2 + t1 + t2))))/
           ((-1 + s1)*(-s + s1 - t2)))*
        (F(np1) + F(np2) - 
          2*(G(Sqrt(1 - 4/Power(-s + s1 + s2,2)),np1,
              (np1 - (2*np2)/(s - s1 - s2))/
               Sqrt(1 - 4/Power(-s + s1 + s2,2))) + Log(4))) + 
       8*(-10*s - 8*s1 + 18*s2 + 24*s*s2 + 8*s1*s2 - 32*Power(s2,2) - 
          6*s*Power(s2,2) + 6*Power(s2,3) - 18*t1 - 14*s*t1 + 46*s2*t1 + 
          6*s*s2*t1 - 12*Power(s2,2)*t1 - 14*Power(t1,2) + 
          6*s2*Power(t1,2) + 8*t2 + 6*s*t2 - 14*s2*t2 - 6*s*s2*t2 + 
          6*Power(s2,2)*t2 + 6*t1*t2 - 6*s2*t1*t2 - 
          ((s - s2 + t1)*(5 - 13*s2 + 4*Power(s2,2) + 5*t1 + 5*s2*t1 - 
               2*Power(s2,2)*t1 - 6*Power(t1,2) + 2*s2*Power(t1,2) + 
               s1*(-1 + s2)*(-2 + 3*s2 + t1) - t2 + Power(s2,2)*t2 + 
               2*t1*t2 - 2*s2*t1*t2 + s*(4 + (-5 + s2)*t1 + t2 - s2*t2))\
)/(-1 + t1) + ((-1 + s2)*(-s + s2 - t1)*
             (-1 + 2*s2 - 6*Power(s2,2) + 
               Power(s1,2)*(-2 + 3*s2 - t1) + 6*s2*t1 - 3*Power(t1,2) + 
               4*t2 - 3*s2*t2 + 2*Power(s2,2)*t2 + 2*t1*t2 - s2*t1*t2 - 
               3*Power(t2,2) + s2*Power(t2,2) + 
               Power(s,2)*(-4 + 2*s1 + t1 + t2) + 
               s1*(1 + 2*Power(s2,2) + Power(t1,2) - s2*(3 + t1) + 
                  t2 - t1*t2) - 
               s*(-3 + 2*Power(s1,2) + 4*t1 + Power(t1,2) + 
                  s1*(-2 + 3*s2 - t2) + 3*s2*(-2 + t2) + 2*t2 - 
                  4*t1*t2 + Power(t2,2))))/((-1 + s1)*(-1 + t2)) - 
          ((-1 + s2)*(-s + s2 - t1)*
             (-2 + 3*t1 - 2*s2*t1 + Power(t1,2) + 3*t2 - 10*s2*t2 + 
               6*t1*t2 - s2*t1*t2 + Power(t2,2) + s2*Power(t2,2) - 
               2*Power(s,2)*(-2 + t1 + t2) + 
               s*(2*s1*(-3 + t1) + 6*t1 - Power(t1,2) + 
                  2*s2*(-3 + t2) + 6*t2 - 2*t1*t2 - Power(t2,2)) + 
               s1*(Power(t1,2) - 2*t2 - t1*(10 + t2) + 
                  s2*(10 + t1 + t2))))/((-1 + s1)*(-s + s1 - t2)) - 
          ((-1 + s2)*(-s + s2 - t1)*
             (1 - 8*s2 + 5*t1 + 6*s2*t1 - 4*Power(t1,2) - 
               2*Power(s1,2)*(-2 + s2 + t1) - 4*t2 + 7*s2*t2 - 
               7*t1*t2 + 2*Power(t1,2)*t2 + 3*Power(t2,2) - 
               s2*Power(t2,2) + Power(s,2)*(-2 + t1 + t2) + 
               s*(3 - 8*t1 + 2*Power(t1,2) + s1*(2 + s2 + t1 - 2*t2) - 
                  s2*(-2 + t2) - 2*t2 + t1*t2 + Power(t2,2)) + 
               s1*(4 - 2*Power(t1,2) - 8*t2 + t1*(7 + t2) + 
                  s2*(-5 - 2*t1 + 5*t2))))/((-1 + t1)*(-1 + t2)) - 
          ((-1 + s2)*(-s2 + 2*Power(s2,2) + t1 - 5*s2*t1 + 
               3*Power(t1,2) - 
               s1*(4 + Power(s2,2) + Power(t1,2) - 2*s2*(2 + t1)) + 
               4*t2 - 3*s2*t2 + 5*Power(s2,2)*t2 - t1*t2 - 9*s2*t1*t2 + 
               4*Power(t1,2)*t2 + Power(s,2)*(4 - 2*s2 + t1 + 3*t2) + 
               s*(5 + 2*Power(s2,2) + Power(t1,2) + 
                  s2*(-10 + s1 - 3*t1 - 8*t2) - t2 + t1*(7 - s1 + 7*t2))\
))/(s - s1 + t2) - 2*((-2*(2 + s*(-3 + s2) - s1*Power(-1 + s2,2) + 
                  6*s2)*(s - s2 + t1))/(-1 + t1) + 
             (2*(-1 + s2)*(-s + s2 - t1)*
                (5 + Power(s1,2) - 2*s2 + 2*Power(s2,2) + 6*t1 - 
                  2*s2*t1 + Power(t1,2) + 2*s1*(1 + t1) + 
                  s*(-4 + Power(s1,2) + s2*(-3 + t1 - t2) + 
                     s1*(-1 + s2 + t1 - t2)) - 
                  Power(s,2)*(-2 + s1 + t1 - t2) - 6*t2 - 2*t1*t2 + 
                  Power(t2,2)))/((-1 + s1)*(-1 + t2)) + 
             ((-1 + s2)*(-s + s2 - t1)*
                (2 + 2*Power(s2,2) + 5*t1 + 3*Power(t1,2) - 
                  Power(s1,2)*(-2 + s2 + t1) - 4*t2 + 2*s2*t2 + 
                  2*Power(s2,2)*t2 - 5*t1*t2 - s2*t1*t2 + 
                  2*Power(t2,2) + Power(s,2)*(-2 + t1 + t2) - 
                  s1*(-7 - 7*t1 + Power(t1,2) + s2*(2 + t1 - 3*t2) + 
                     5*t2) + 
                  s*(1 + Power(t1,2) + s1*(2 + s2 - t2) + 
                     t1*(-3 + t2) - 3*t2 - 3*s2*t2)))/
              ((-1 + t1)*(-1 + t2)) + 
             (2*(-1 + s2)*(-6 - 7*s2 + s1*(-2 + s2 - t1) + 5*t1 + 
                  s2*t1 - Power(t1,2) + 2*t2 + Power(s,2)*t2 - s2*t2 + 
                  Power(s2,2)*t2 + t1*t2 - 2*s2*t1*t2 + 
                  Power(t1,2)*t2 - 
                  s*(-7 + s1 + t1 - t2 + 2*s2*t2 - 2*t1*t2)))/
              (s - s1 + t2) + 
             2*(4 - 2*s1*(-1 + s2) + s2 + 10*Power(s2,2) - 
                Power(s2,3) + 3*t1 - 13*s2*t1 + 2*Power(s2,2)*t1 + 
                3*Power(t1,2) - s2*Power(t1,2) - 2*t2 + 3*s2*t2 - 
                Power(s2,2)*t2 - t1*t2 + s2*t1*t2 + 
                s*(1 + Power(s2,2) + 3*t1 - t2 + s2*(-8 - t1 + t2))) + 
             ((-1 + s2)*(-s + s2 - t1)*
                (10 + 4*s2 + 2*Power(s2,2) + 2*Power(s1,2)*(1 + s2) - 
                  t1 - 4*s2*t1 + Power(t1,2) - t2 - 2*s2*t2 + 2*t1*t2 - 
                  s2*t1*t2 + Power(t2,2) + s2*Power(t2,2) - 
                  Power(s,2)*(-6 + t1 + t2) - 
                  s*(15 + s1*(4 + 3*s2 - t1) - 5*t1 + Power(t1,2) - 
                     s2*(-4 + t2) - 5*t2 + Power(t2,2)) + 
                  s1*(4 + 2*Power(s2,2) + Power(t1,2) - 4*t2 - 
                     t1*(2 + t2) - s2*(-2 + t1 + t2))))/
              ((-1 + s1)*(-s + s1 - t2)))*Log(EE))*
        (1 + (1 - 2/(2 - s + s1 + s2))*R2(2 + s - s1 - s2)))/
     (Power(-1 + s2,2)*Power(s - s2 + t1,2)) + 
    (4*((-2*(2 + s*(-3 + s1) - s2 - Power(s1,2)*s2 + 2*s1*(3 + s2))*
             (s - s1 + t2))/(-1 + t2) + 
          2*(4 - Power(s1,3) + 2*s2 - 2*t1 + 3*t2 - t1*t2 + 
             3*Power(t2,2) + 
             s*(1 + Power(s1,2) - t1 + s1*(-8 + t1 - t2) + 3*t2) + 
             s1*(1 - 2*s2 - 13*t2 - Power(t2,2) + t1*(3 + t2)) + 
             Power(s1,2)*(-t1 + 2*(5 + t2))) + 
          ((-1 + s1)*(-s + s1 - t2)*
             (2 + 7*s2 + 2*Power(s2,2) - 4*t1 - 5*s2*t1 + 
               2*Power(t1,2) + 2*Power(s1,2)*(1 + t1) + 5*t2 + 
               7*s2*t2 - Power(s2,2)*t2 - 5*t1*t2 + 3*Power(t2,2) - 
               s2*Power(t2,2) + Power(s,2)*(-2 + t1 + t2) + 
               s*(1 + s2*(2 + s1 - t1) - 3*t2 + Power(t2,2) + 
                  t1*(-3 - 3*s1 + t2)) - 
               s1*(Power(s2,2) + t1*(-2 + t2) + s2*(2 - 3*t1 + t2))))/
           ((-1 + t1)*(-1 + t2)) + 
          (2*(-1 + s1)*(-s + s1 - t2)*
             (5 + 2*Power(s1,2) + 2*s2 + Power(s2,2) - 6*t1 + 
               Power(t1,2) + 6*t2 + 2*s2*t2 - 2*t1*t2 + Power(t2,2) - 
               2*s1*(1 + t2) - Power(s,2)*(-2 + s2 - t1 + t2) + 
               s*(-4 + Power(s2,2) + s2*(-1 - t1 + t2) + 
                  s1*(-3 + s2 - t1 + t2))))/((-1 + s2)*(-1 + t1)) + 
          ((-1 + s1)*(-s + s1 - t2)*
             (10 + 4*s2 + 2*Power(s2,2) + 2*Power(s1,2)*(1 + s2) - 
               t1 - 4*s2*t1 + Power(t1,2) - t2 - 2*s2*t2 + 2*t1*t2 - 
               s2*t1*t2 + Power(t2,2) + s2*Power(t2,2) - 
               Power(s,2)*(-6 + t1 + t2) - 
               s*(15 + s1*(4 + 3*s2 - t1) - 5*t1 + Power(t1,2) - 
                  s2*(-4 + t2) - 5*t2 + Power(t2,2)) + 
               s1*(4 + 2*Power(s2,2) + Power(t1,2) - 4*t2 - 
                  t1*(2 + t2) - s2*(-2 + t1 + t2))))/
           ((-1 + s2)*(-s + s2 - t1)) + 
          (2*(-1 + s1)*(-6 - 2*s2 + 2*t1 + Power(s,2)*t1 + 
               Power(s1,2)*t1 + 5*t2 - s2*t2 + t1*t2 - Power(t2,2) + 
               t1*Power(t2,2) - 
               s*(-7 + s2 + t1*(-1 + 2*s1 - 2*t2) + t2) + 
               s1*(-7 + s2 + t2 - t1*(1 + 2*t2))))/(s - s2 + t1))*
        (F(np1) + F(np2) - 
          2*(G(Sqrt(1 - 4/Power(-s + s1 + s2,2)),np1,
              (np1 - (2*np2)/(s - s1 - s2))/
               Sqrt(1 - 4/Power(-s + s1 + s2,2))) + Log(4))) + 
       8*(-10*s + 18*s1 + 24*s*s1 - 32*Power(s1,2) - 6*s*Power(s1,2) + 
          6*Power(s1,3) - 8*s2 + 8*s1*s2 + 8*t1 + 6*s*t1 - 14*s1*t1 - 
          6*s*s1*t1 + 6*Power(s1,2)*t1 - 18*t2 - 14*s*t2 + 46*s1*t2 + 
          6*s*s1*t2 - 12*Power(s1,2)*t2 + 6*t1*t2 - 6*s1*t1*t2 - 
          14*Power(t2,2) + 6*s1*Power(t2,2) - 
          ((s - s1 + t2)*(5 + 2*s2 - t1 + 
               Power(s1,2)*(4 + 3*s2 + t1 - 2*t2) + 5*t2 - s2*t2 + 
               2*t1*t2 - 6*Power(t2,2) + 
               s*(4 + t1 - s1*t1 + (-5 + s1)*t2) + 
               s1*(-13 + s2*(-5 + t2) + (5 - 2*t1)*t2 + 2*Power(t2,2))))/
           (-1 + t2) - ((-1 + s1)*
             (-4*s2 + 4*t1 + Power(s1,2)*(2 - s2 + 5*t1) + t2 - t1*t2 + 
               3*Power(t2,2) - s2*Power(t2,2) + 4*t1*Power(t2,2) + 
               Power(s,2)*(4 - 2*s1 + 3*t1 + t2) + 
               s*(5 + 2*Power(s1,2) - t1 + 
                  s1*(-10 + s2 - 8*t1 - 3*t2) + 7*t2 - s2*t2 + 
                  7*t1*t2 + Power(t2,2)) + 
               s1*(-1 - 3*t1 - 5*t2 - 9*t1*t2 + 2*s2*(2 + t2))))/
           (s - s2 + t1) + ((-1 + s1)*(-s + s1 - t2)*
             (-1 + s2 - 2*Power(s2,2) + 4*t1 + s2*t1 - 3*Power(t1,2) + 
               2*Power(s1,2)*(-3 + s2 + t1) - Power(s2,2)*t2 + 
               2*t1*t2 - s2*t1*t2 - 3*Power(t2,2) + s2*Power(t2,2) + 
               Power(s,2)*(-4 + 2*s2 + t1 + t2) - 
               s*(-3 + 2*Power(s2,2) + 2*t1 + Power(t1,2) - 
                  s2*(2 + t1) + 3*s1*(-2 + s2 + t1) + 4*t2 - 4*t1*t2 + 
                  Power(t2,2)) + 
               s1*(2 + 3*Power(s2,2) + Power(t1,2) + 6*t2 - 
                  s2*(3 + t2) - t1*(3 + t2))))/((-1 + s2)*(-1 + t1)) - 
          ((-1 + s1)*(-s + s1 - t2)*
             (-2 + 3*t1 - 2*s2*t1 + Power(t1,2) + 3*t2 - 10*s2*t2 + 
               6*t1*t2 - s2*t1*t2 + Power(t2,2) + s2*Power(t2,2) - 
               2*Power(s,2)*(-2 + t1 + t2) + 
               s*(2*s1*(-3 + t1) + 6*t1 - Power(t1,2) + 
                  2*s2*(-3 + t2) + 6*t2 - 2*t1*t2 - Power(t2,2)) + 
               s1*(Power(t1,2) - 2*t2 - t1*(10 + t2) + 
                  s2*(10 + t1 + t2))))/((-1 + s2)*(-s + s2 - t1)) - 
          ((-1 + s1)*(-s + s1 - t2)*
             (1 + 4*s2 + 4*Power(s2,2) - 4*t1 - 8*s2*t1 + 
               3*Power(t1,2) + 5*t2 + 7*s2*t2 - 2*Power(s2,2)*t2 - 
               7*t1*t2 + s2*t1*t2 - 4*Power(t2,2) - 2*s2*Power(t2,2) + 
               2*t1*Power(t2,2) + Power(s,2)*(-2 + t1 + t2) + 
               s*(3 + s1*(2 + s2 - t1) - 2*t1 + Power(t1,2) - 8*t2 + 
                  t1*t2 + 2*Power(t2,2) + s2*(2 - 2*t1 + t2)) - 
               s1*(8 + 2*Power(s2,2) - 7*t1 + Power(t1,2) - 6*t2 + 
                  s2*(5 - 5*t1 + 2*t2))))/((-1 + t1)*(-1 + t2)) - 
          2*((-2*(2 + s*(-3 + s1) - s2 - Power(s1,2)*s2 + 
                  2*s1*(3 + s2))*(s - s1 + t2))/(-1 + t2) + 
             2*(4 - Power(s1,3) + 2*s2 - 2*t1 + 3*t2 - t1*t2 + 
                3*Power(t2,2) + 
                s*(1 + Power(s1,2) - t1 + s1*(-8 + t1 - t2) + 3*t2) + 
                s1*(1 - 2*s2 - 13*t2 - Power(t2,2) + t1*(3 + t2)) + 
                Power(s1,2)*(-t1 + 2*(5 + t2))) + 
             ((-1 + s1)*(-s + s1 - t2)*
                (2 + 7*s2 + 2*Power(s2,2) - 4*t1 - 5*s2*t1 + 
                  2*Power(t1,2) + 2*Power(s1,2)*(1 + t1) + 5*t2 + 
                  7*s2*t2 - Power(s2,2)*t2 - 5*t1*t2 + 3*Power(t2,2) - 
                  s2*Power(t2,2) + Power(s,2)*(-2 + t1 + t2) + 
                  s*(1 + s2*(2 + s1 - t1) - 3*t2 + Power(t2,2) + 
                     t1*(-3 - 3*s1 + t2)) - 
                  s1*(Power(s2,2) + t1*(-2 + t2) + s2*(2 - 3*t1 + t2)))\
)/((-1 + t1)*(-1 + t2)) + (2*(-1 + s1)*(-s + s1 - t2)*
                (5 + 2*Power(s1,2) + 2*s2 + Power(s2,2) - 6*t1 + 
                  Power(t1,2) + 6*t2 + 2*s2*t2 - 2*t1*t2 + 
                  Power(t2,2) - 2*s1*(1 + t2) - 
                  Power(s,2)*(-2 + s2 - t1 + t2) + 
                  s*(-4 + Power(s2,2) + s2*(-1 - t1 + t2) + 
                     s1*(-3 + s2 - t1 + t2))))/((-1 + s2)*(-1 + t1)) + 
             ((-1 + s1)*(-s + s1 - t2)*
                (10 + 4*s2 + 2*Power(s2,2) + 2*Power(s1,2)*(1 + s2) - 
                  t1 - 4*s2*t1 + Power(t1,2) - t2 - 2*s2*t2 + 
                  2*t1*t2 - s2*t1*t2 + Power(t2,2) + s2*Power(t2,2) - 
                  Power(s,2)*(-6 + t1 + t2) - 
                  s*(15 + s1*(4 + 3*s2 - t1) - 5*t1 + Power(t1,2) - 
                     s2*(-4 + t2) - 5*t2 + Power(t2,2)) + 
                  s1*(4 + 2*Power(s2,2) + Power(t1,2) - 4*t2 - 
                     t1*(2 + t2) - s2*(-2 + t1 + t2))))/
              ((-1 + s2)*(-s + s2 - t1)) + 
             (2*(-1 + s1)*(-6 - 2*s2 + 2*t1 + Power(s,2)*t1 + 
                  Power(s1,2)*t1 + 5*t2 - s2*t2 + t1*t2 - Power(t2,2) + 
                  t1*Power(t2,2) - 
                  s*(-7 + s2 + t1*(-1 + 2*s1 - 2*t2) + t2) + 
                  s1*(-7 + s2 + t2 - t1*(1 + 2*t2))))/(s - s2 + t1))*
           Log(EE))*(1 + (1 - 2/(2 - s + s1 + s2))*R2(2 + s - s1 - s2)))/
     (Power(-1 + s1,2)*Power(s - s1 + t2,2)) + 
    (4*((2*(s - s1 + t2)*(-6 - 7*s2 + s1*(-2 + s2 - t1) + 5*t1 + 
               s2*t1 - Power(t1,2) + 2*t2 + Power(s,2)*t2 - s2*t2 + 
               Power(s2,2)*t2 + t1*t2 - 2*s2*t1*t2 + Power(t1,2)*t2 - 
               s*(-7 + s1 + t1 - t2 + 2*s2*t2 - 2*t1*t2)))/(-1 + s2) + 
          2*(8 + Power(s,3) + 6*s2 + s1*(6 + 4*s2 - 6*t1) - 6*t1 - 
             6*t2 - 6*s2*t2 + 8*t1*t2 + 
             Power(s,2)*(8 - s1 - s2 + t1 + t2) + 
             s*(-12 + s1*(-6 + s2 - t1) + 8*t1 + 8*t2 + t1*t2 - 
                s2*(6 + t2))) + 
          ((s - s2 + t1)*(s - s1 + t2)*
             (3 + 3*s2 - 2*Power(s2,2) + 5*t1 + s2*t1 + 
               2*Power(t1,2) - Power(s1,2)*(-6 + s2 + t1) - 6*t2 + 
               Power(s2,2)*t2 - 5*t1*t2 + 3*Power(t2,2) - 
               s2*Power(t2,2) + 2*Power(s,2)*(1 + t2) + 
               s1*(7 + s2 + Power(s2,2) + 4*t1 - 3*s2*t1 - 5*t2 + 
                  4*s2*t2 - t1*t2) + 
               s*(-4 + s1*(-4 + t1 - 3*t2) + t2 - 3*s2*t2 + 
                  Power(t2,2) + t1*(1 + s2 + t2))))/
           ((-1 + s1)*(-1 + t2)) + 
          ((s - s2 + t1)*(s - s1 + t2)*
             (3 + 7*s2 + 6*Power(s2,2) - 6*t1 - 5*s2*t1 + 
               3*Power(t1,2) + 2*Power(s,2)*(1 + t1) + 
               Power(s1,2)*(-2 + s2 + t1) + 5*t2 + 4*s2*t2 - 
               Power(s2,2)*t2 - 5*t1*t2 - s2*t1*t2 + 2*Power(t2,2) + 
               s1*(3 + s2 - Power(s2,2) + 4*s2*t1 - Power(t1,2) + 
                  t2 - 3*s2*t2) + 
               s*(-4 + t1 - 3*s1*t1 + Power(t1,2) + t2 + s1*t2 + 
                  t1*t2 + s2*(-4 - 3*t1 + t2))))/((-1 + s2)*(-1 + t1)) + 
          (2*(s - s2 + t1)*(s - s1 + t2)*
             (2*Power(s,2) + 3*s2 + 2*Power(s2,2) - 2*t1 - 2*s2*t1 - 
               Power(s2,2)*t1 + 2*Power(t1,2) + s2*Power(t1,2) - 
               2*t2 + Power(s2,2)*t2 - 2*s2*t1*t2 + 2*Power(t2,2) + 
               s2*Power(t2,2) - Power(s1,2)*(-2 + s2 - t1 + t2) + 
               s1*(3 + 4*s2 - Power(s2,2) + Power(t1,2) - 2*t2 - 
                  2*t1*t2 + Power(t2,2)) + 
               s*(1 + 2*t1 - Power(t1,2) + s2*(-3 + t1 - t2) + 2*t2 + 
                  2*t1*t2 - Power(t2,2) + s1*(-3 + s2 - t1 + t2))))/
           ((-1 + t1)*(-1 + t2)) + 
          (2*(s - s2 + t1)*(-6 - 2*s2 + 2*t1 + Power(s,2)*t1 + 
               Power(s1,2)*t1 + 5*t2 - s2*t2 + t1*t2 - Power(t2,2) + 
               t1*Power(t2,2) - 
               s*(-7 + s2 + t1*(-1 + 2*s1 - 2*t2) + t2) + 
               s1*(-7 + s2 + t2 - t1*(1 + 2*t2))))/(-1 + s1))*
        (F(np1) + F(np2) - 
          2*(G(Sqrt(1 - 4/Power(-s + s1 + s2,2)),np1,
              (np1 - (2*np2)/(s - s1 - s2))/
               Sqrt(1 - 4/Power(-s + s1 + s2,2))) + Log(4))) + 
       8*(16*s - 20*Power(s,2) - 6*Power(s,3) - 8*s1 + 12*s*s1 + 
          6*Power(s,2)*s1 - 8*s2 + 12*s*s2 + 6*Power(s,2)*s2 - 4*s1*s2 - 
          6*s*s1*s2 + 8*t1 - 20*s*t1 - 6*Power(s,2)*t1 + 12*s1*t1 + 
          6*s*s1*t1 + 8*t2 - 20*s*t2 - 6*Power(s,2)*t2 + 12*s2*t2 + 
          6*s*s2*t2 - 20*t1*t2 - 6*s*t1*t2 - 
          ((s - s2 + t1)*(-4*s2 + 4*t1 + Power(s1,2)*(2 - s2 + 5*t1) + 
               t2 - t1*t2 + 3*Power(t2,2) - s2*Power(t2,2) + 
               4*t1*Power(t2,2) + Power(s,2)*(4 - 2*s1 + 3*t1 + t2) + 
               s*(5 + 2*Power(s1,2) - t1 + 
                  s1*(-10 + s2 - 8*t1 - 3*t2) + 7*t2 - s2*t2 + 
                  7*t1*t2 + Power(t2,2)) + 
               s1*(-1 - 3*t1 - 5*t2 - 9*t1*t2 + 2*s2*(2 + t2))))/
           (-1 + s1) - ((s - s2 + t1)*(s - s1 + t2)*
             (-4 - 2*s2 + 4*Power(s2,2) + 6*t1 - 2*s2*t1 - 
               2*Power(t1,2) + Power(s1,2)*(-2 + s2 + t1) + 
               s1*(-2*Power(s2,2) + 2*s2*(2 + t1 - 2*t2) + 
                  (1 + t1)*(-1 + t2)) + 3*t2 - 2*s2*t2 - 
               2*Power(s2,2)*t2 - 3*t1*t2 + 4*s2*t1*t2 - 
               2*Power(t1,2)*t2 + 3*Power(t2,2) - s2*Power(t2,2) + 
               s*(6 - 4*t1 + 2*s2*(-3 + t2) + 2*t2 - 3*t1*t2 + 
                  Power(t2,2) + s1*(2 - t1 + t2))))/((-1 + s2)*(-1 + t1)) \
- ((s - s2 + t1)*(s - s1 + t2)*
             (-4 - s2 - 2*Power(s2,2) + 3*t1 + s2*t1 + 3*Power(t1,2) - 
               2*Power(s1,2)*(-2 + s2 + t1) + 6*t2 - s2*t2 + 
               Power(s2,2)*t2 - 3*t1*t2 + s2*t1*t2 - 2*Power(t2,2) - 
               2*t1*Power(t2,2) + 
               s*(6 + 2*s1*(-3 + t1) + 2*t1 + Power(t1,2) + 
                  s2*(2 + t1 - t2) - 4*t2 - 3*t1*t2) + 
               s1*(Power(s2,2) - Power(t1,2) - 2*(1 + t2) + 
                  s2*(4 - 4*t1 + 2*t2) + t1*(-2 + 4*t2))))/
           ((-1 + s1)*(-1 + t2)) - 
          ((s - s1 + t2)*(-s2 + 2*Power(s2,2) + t1 - 5*s2*t1 + 
               3*Power(t1,2) - 
               s1*(4 + Power(s2,2) + Power(t1,2) - 2*s2*(2 + t1)) + 
               4*t2 - 3*s2*t2 + 5*Power(s2,2)*t2 - t1*t2 - 9*s2*t1*t2 + 
               4*Power(t1,2)*t2 + Power(s,2)*(4 - 2*s2 + t1 + 3*t2) + 
               s*(5 + 2*Power(s2,2) + Power(t1,2) + 
                  s2*(-10 + s1 - 3*t1 - 8*t2) - t2 + t1*(7 - s1 + 7*t2))\
))/(-1 + s2) + ((s - s2 + t1)*(s - s1 + t2)*
             (-s2 - 4*Power(s2,2) + 2*t1 + 5*s2*t1 + 2*Power(s2,2)*t1 - 
               2*Power(t1,2) - 2*s2*Power(t1,2) + 2*t2 + s2*t2 + 
               Power(s2,2)*t2 - 2*s2*t1*t2 - 2*Power(t2,2) + 
               2*Power(s,2)*(-3 + t1 + t2) + 
               Power(s1,2)*(-4 + s2 + t1 + 2*t2) + 
               s1*(-1 + Power(s2,2) + t1 + 5*t2 - 2*t1*t2 - 
                  2*Power(t2,2) + s2*(-4 + t1 + t2)) - 
               s*(3*s1*(-2 + t1 + t2) + 3*s2*(-2 + t1 + t2) + 
                  4*(-1 + t1 + t2 - t1*t2))))/((-1 + t1)*(-1 + t2)) - 
          2*((2*(s - s1 + t2)*
                (-6 - 7*s2 + s1*(-2 + s2 - t1) + 5*t1 + s2*t1 - 
                  Power(t1,2) + 2*t2 + Power(s,2)*t2 - s2*t2 + 
                  Power(s2,2)*t2 + t1*t2 - 2*s2*t1*t2 + 
                  Power(t1,2)*t2 - 
                  s*(-7 + s1 + t1 - t2 + 2*s2*t2 - 2*t1*t2)))/(-1 + s2) \
+ 2*(8 + Power(s,3) + 6*s2 + s1*(6 + 4*s2 - 6*t1) - 6*t1 - 6*t2 - 
                6*s2*t2 + 8*t1*t2 + 
                Power(s,2)*(8 - s1 - s2 + t1 + t2) + 
                s*(-12 + s1*(-6 + s2 - t1) + 8*t1 + 8*t2 + t1*t2 - 
                   s2*(6 + t2))) + 
             ((s - s2 + t1)*(s - s1 + t2)*
                (3 + 3*s2 - 2*Power(s2,2) + 5*t1 + s2*t1 + 
                  2*Power(t1,2) - Power(s1,2)*(-6 + s2 + t1) - 6*t2 + 
                  Power(s2,2)*t2 - 5*t1*t2 + 3*Power(t2,2) - 
                  s2*Power(t2,2) + 2*Power(s,2)*(1 + t2) + 
                  s1*(7 + s2 + Power(s2,2) + 4*t1 - 3*s2*t1 - 5*t2 + 
                     4*s2*t2 - t1*t2) + 
                  s*(-4 + s1*(-4 + t1 - 3*t2) + t2 - 3*s2*t2 + 
                     Power(t2,2) + t1*(1 + s2 + t2))))/
              ((-1 + s1)*(-1 + t2)) + 
             ((s - s2 + t1)*(s - s1 + t2)*
                (3 + 7*s2 + 6*Power(s2,2) - 6*t1 - 5*s2*t1 + 
                  3*Power(t1,2) + 2*Power(s,2)*(1 + t1) + 
                  Power(s1,2)*(-2 + s2 + t1) + 5*t2 + 4*s2*t2 - 
                  Power(s2,2)*t2 - 5*t1*t2 - s2*t1*t2 + 
                  2*Power(t2,2) + 
                  s1*(3 + s2 - Power(s2,2) + 4*s2*t1 - Power(t1,2) + 
                     t2 - 3*s2*t2) + 
                  s*(-4 + t1 - 3*s1*t1 + Power(t1,2) + t2 + s1*t2 + 
                     t1*t2 + s2*(-4 - 3*t1 + t2))))/
              ((-1 + s2)*(-1 + t1)) + 
             (2*(s - s2 + t1)*(s - s1 + t2)*
                (2*Power(s,2) + 3*s2 + 2*Power(s2,2) - 2*t1 - 
                  2*s2*t1 - Power(s2,2)*t1 + 2*Power(t1,2) + 
                  s2*Power(t1,2) - 2*t2 + Power(s2,2)*t2 - 
                  2*s2*t1*t2 + 2*Power(t2,2) + s2*Power(t2,2) - 
                  Power(s1,2)*(-2 + s2 - t1 + t2) + 
                  s1*(3 + 4*s2 - Power(s2,2) + Power(t1,2) - 2*t2 - 
                     2*t1*t2 + Power(t2,2)) + 
                  s*(1 + 2*t1 - Power(t1,2) + s2*(-3 + t1 - t2) + 
                     2*t2 + 2*t1*t2 - Power(t2,2) + 
                     s1*(-3 + s2 - t1 + t2))))/((-1 + t1)*(-1 + t2)) + 
             (2*(s - s2 + t1)*
                (-6 - 2*s2 + 2*t1 + Power(s,2)*t1 + Power(s1,2)*t1 + 
                  5*t2 - s2*t2 + t1*t2 - Power(t2,2) + t1*Power(t2,2) - 
                  s*(-7 + s2 + t1*(-1 + 2*s1 - 2*t2) + t2) + 
                  s1*(-7 + s2 + t2 - t1*(1 + 2*t2))))/(-1 + s1))*Log(EE))*
        (1 + (1 - 2/(2 - s + s1 + s2))*R2(2 + s - s1 - s2)))/
     (Power(s - s2 + t1,2)*Power(s - s1 + t2,2)) + 
    (4*((-2*(2 + s*(-3 + s2) - s1*Power(-1 + s2,2) + 6*s2)*(-1 + t1))/
           (s - s2 + t1) + 2*(7 - s2 + 2*s*(-1 + t1) - 8*t1 + 8*s2*t1 + 
             3*Power(t1,2) - s2*Power(t1,2) + 
             s1*(-1 + s2 + t1 - s2*t1) + 3*t2 - 3*s2*t2 - t1*t2 + 
             s2*t1*t2) - (2*(-1 + s2)*
             (2 - 3*s2 - s*Power(-1 + t1,2) + s1*Power(-1 + t1,2) - 
               6*t1 + s2*t1 - 4*t2 + 3*t1*t2 - Power(t1,2)*t2))/(-1 + t2) \
+ ((-1 + s2)*(-1 + t1)*(12 + s2 - 2*Power(s2,2) - 3*t1 + 3*s2*t1 + 
               Power(t1,2) + Power(s1,2)*(-2 + s2 + t1) - 3*t2 - 
               s2*t2 + Power(s2,2)*t2 + 2*t1*t2 - s2*t1*t2 + 
               Power(t2,2) + s2*Power(t2,2) - 
               Power(s,2)*(-2 + t1 + t2) + 
               s1*(1 + Power(s2,2) + Power(t1,2) + 3*t2 - 
                  t1*(1 + t2)) + 
               s*(-13 + t1 - Power(t1,2) + s2*(2 + t1) + t2 - 
                  Power(t2,2) + s1*(2 - s2 + t2))))/((-1 + s1)*(-1 + t2)) \
+ (2*(-1 + s2)*(-1 + t1)*(5 + 2*Power(s1,2) + 2*s2 + Power(s2,2) - 
               6*t1 + Power(t1,2) + 6*t2 + 2*s2*t2 - 2*t1*t2 + 
               Power(t2,2) - 2*s1*(1 + t2) - 
               Power(s,2)*(-2 + s2 - t1 + t2) + 
               s*(-4 + Power(s2,2) + s2*(-1 - t1 + t2) + 
                  s1*(-3 + s2 - t1 + t2))))/((-1 + s1)*(-s + s1 - t2)) - 
          ((-1 + s2)*(-1 + t1)*
             (-3 - 7*s2 - 6*Power(s2,2) + 6*t1 + 5*s2*t1 - 
               3*Power(t1,2) - 2*Power(s,2)*(1 + t1) - 
               Power(s1,2)*(-2 + s2 + t1) - 5*t2 - 4*s2*t2 + 
               Power(s2,2)*t2 + 5*t1*t2 + s2*t1*t2 - 2*Power(t2,2) - 
               s*(-4 + t1 - 3*s1*t1 + Power(t1,2) + t2 + s1*t2 + 
                  t1*t2 + s2*(-4 - 3*t1 + t2)) + 
               s1*(-3 + Power(s2,2) + Power(t1,2) - t2 + 
                  s2*(-1 - 4*t1 + 3*t2))))/((s - s2 + t1)*(s - s1 + t2)))*
        (F(np1) + F(np2) - 
          2*(G(Sqrt(1 - 4/Power(-s + s1 + s2,2)),np1,
              (np1 - (2*np2)/(s - s1 - s2))/
               Sqrt(1 - 4/Power(-s + s1 + s2,2))) + Log(4))) + 
       8*(-18 + 8*s + 6*s1 + 10*s2 - 6*s1*s2 + 40*t1 - 8*s*t1 - 6*s1*t1 - 
          24*s2*t1 + 6*s1*s2*t1 - 14*Power(t1,2) + 6*s2*Power(t1,2) - 
          14*t2 + 14*s2*t2 + 6*t1*t2 - 6*s2*t1*t2 - 
          ((-1 + s2)*(-11 + 4*s2 + 14*t1 - 7*Power(t1,2) - 
               s1*(-1 + t1)*(-2 + s2 + 3*t1) + 7*t2 - 5*s2*t2 - 
               7*t1*t2 + s2*t1*t2 + 4*Power(t1,2)*t2 + 
               s*(-1 + t1)*(-2 + 3*t1 + t2)))/(-1 + t2) - 
          ((-1 + t1)*(5 - 13*s2 + 4*Power(s2,2) + 5*t1 + 5*s2*t1 - 
               2*Power(s2,2)*t1 - 6*Power(t1,2) + 2*s2*Power(t1,2) + 
               s1*(-1 + s2)*(-2 + 3*s2 + t1) - t2 + Power(s2,2)*t2 + 
               2*t1*t2 - 2*s2*t1*t2 + s*(4 + (-5 + s2)*t1 + t2 - s2*t2)))/
           (s - s2 + t1) + ((-1 + s2)*(-1 + t1)*
             (-1 + s2 - 2*Power(s2,2) + 4*t1 + s2*t1 - 3*Power(t1,2) + 
               2*Power(s1,2)*(-3 + s2 + t1) - Power(s2,2)*t2 + 2*t1*t2 - 
               s2*t1*t2 - 3*Power(t2,2) + s2*Power(t2,2) + 
               Power(s,2)*(-4 + 2*s2 + t1 + t2) - 
               s*(-3 + 2*Power(s2,2) + 2*t1 + Power(t1,2) - 
                  s2*(2 + t1) + 3*s1*(-2 + s2 + t1) + 4*t2 - 4*t1*t2 + 
                  Power(t2,2)) + 
               s1*(2 + 3*Power(s2,2) + Power(t1,2) + 6*t2 - 
                  s2*(3 + t2) - t1*(3 + t2))))/((-1 + s1)*(-s + s1 - t2)) \
+ ((-1 + s2)*(-1 + t1)*(4 + 2*s2 - 4*Power(s2,2) - 6*t1 + 2*s2*t1 + 
               2*Power(t1,2) - Power(s1,2)*(-2 + s2 + t1) + 
               s1*(2*Power(s2,2) - 2*s2*(2 + t1 - 2*t2) - 
                  (1 + t1)*(-1 + t2)) - 3*t2 + 2*s2*t2 + 
               2*Power(s2,2)*t2 + 3*t1*t2 - 4*s2*t1*t2 + 
               2*Power(t1,2)*t2 - 3*Power(t2,2) + s2*Power(t2,2) - 
               s*(6 - 4*t1 + 2*s2*(-3 + t2) + 2*t2 - 3*t1*t2 + 
                  Power(t2,2) + s1*(2 - t1 + t2))))/
           ((s - s2 + t1)*(s - s1 + t2)) - 
          ((-1 + s2)*(-1 + t1)*
             (12 - s2 - 2*Power(s2,2) - 11*t1 + s2*t1 + Power(t1,2) + 
               Power(s1,2)*(-2 + s2 + t1) - 11*t2 + 3*s2*t2 + 
               Power(s2,2)*t2 + 6*t1*t2 - s2*t1*t2 + Power(t2,2) + 
               s2*Power(t2,2) - 2*Power(s,2)*(-2 + t1 + t2) + 
               s1*(-1 + Power(s2,2) + Power(t1,2) - t1*(-3 + t2) + t2 - 
                  2*s2*(t1 + t2)) + 
               s*(-14 + 6*t1 - Power(t1,2) + 6*t2 - 2*t1*t2 - 
                  Power(t2,2) + s2*(2 + t1 + t2) + 
                  s1*(2 - 2*s2 + t1 + t2))))/((-1 + s1)*(-1 + t2)) - 
          2*((-2*(2 + s*(-3 + s2) - s1*Power(-1 + s2,2) + 6*s2)*
                (-1 + t1))/(s - s2 + t1) + 
             2*(7 - s2 + 2*s*(-1 + t1) - 8*t1 + 8*s2*t1 + 
                3*Power(t1,2) - s2*Power(t1,2) + 
                s1*(-1 + s2 + t1 - s2*t1) + 3*t2 - 3*s2*t2 - t1*t2 + 
                s2*t1*t2) - (2*(-1 + s2)*
                (2 - 3*s2 - s*Power(-1 + t1,2) + s1*Power(-1 + t1,2) - 
                  6*t1 + s2*t1 - 4*t2 + 3*t1*t2 - Power(t1,2)*t2))/
              (-1 + t2) + ((-1 + s2)*(-1 + t1)*
                (12 + s2 - 2*Power(s2,2) - 3*t1 + 3*s2*t1 + 
                  Power(t1,2) + Power(s1,2)*(-2 + s2 + t1) - 3*t2 - 
                  s2*t2 + Power(s2,2)*t2 + 2*t1*t2 - s2*t1*t2 + 
                  Power(t2,2) + s2*Power(t2,2) - 
                  Power(s,2)*(-2 + t1 + t2) + 
                  s1*(1 + Power(s2,2) + Power(t1,2) + 3*t2 - 
                     t1*(1 + t2)) + 
                  s*(-13 + t1 - Power(t1,2) + s2*(2 + t1) + t2 - 
                     Power(t2,2) + s1*(2 - s2 + t2))))/
              ((-1 + s1)*(-1 + t2)) + 
             (2*(-1 + s2)*(-1 + t1)*
                (5 + 2*Power(s1,2) + 2*s2 + Power(s2,2) - 6*t1 + 
                  Power(t1,2) + 6*t2 + 2*s2*t2 - 2*t1*t2 + 
                  Power(t2,2) - 2*s1*(1 + t2) - 
                  Power(s,2)*(-2 + s2 - t1 + t2) + 
                  s*(-4 + Power(s2,2) + s2*(-1 - t1 + t2) + 
                     s1*(-3 + s2 - t1 + t2))))/((-1 + s1)*(-s + s1 - t2)) \
- ((-1 + s2)*(-1 + t1)*(-3 - 7*s2 - 6*Power(s2,2) + 6*t1 + 5*s2*t1 - 
                  3*Power(t1,2) - 2*Power(s,2)*(1 + t1) - 
                  Power(s1,2)*(-2 + s2 + t1) - 5*t2 - 4*s2*t2 + 
                  Power(s2,2)*t2 + 5*t1*t2 + s2*t1*t2 - 2*Power(t2,2) - 
                  s*(-4 + t1 - 3*s1*t1 + Power(t1,2) + t2 + s1*t2 + 
                     t1*t2 + s2*(-4 - 3*t1 + t2)) + 
                  s1*(-3 + Power(s2,2) + Power(t1,2) - t2 + 
                     s2*(-1 - 4*t1 + 3*t2))))/
              ((s - s2 + t1)*(s - s1 + t2)))*Log(EE))*
        (1 + (1 - 2/(2 - s + s1 + s2))*R2(2 + s - s1 - s2)))/
     (Power(-1 + s2,2)*Power(-1 + t1,2)));
    return a;
};
