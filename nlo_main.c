#include <stdio.h>
#include "nlo.h"


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
    E1=0.167;
    double result=nlo(s,s1,s2,t1,t2,E1,np1,np2);
    //result=G(Sqrt(1 - 4/Power(-s + s1 + s2,2)),np1,(np1 - (2*np2)/(s - s1 - s2))/sqrt(1 - 4/Power(-s + s1 + s2,2)));
    //result=B1(1 - s + s1 - t2,s,s1);
    printf("%f\n",result);
};