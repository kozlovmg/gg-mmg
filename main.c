#include <stdio.h>

double nlo2(double s,double s1,double s2,double t1,double t2,double EE,double np1,double np2);

int main(void)
{
    double s, s1, s2, t1, t2, E1, np1, np2;
    s=400.0;
    s1=2.880;
    s2=107.061;
    t1=-145.969;
    t2=-198.059;
    np1=7.348;
    np2=9.95298;
    E1=0.16;
    double result=nlo2(s,s1,s2,t1,t2,E1,np1,np2);
    printf("%f\n",result);
    };