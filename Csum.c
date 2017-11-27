#include <stdio.h>

long double Cpart1(long double s,long double s1,long double s2,long double t1,long double t2);
long double Cpart2(long double s,long double s1,long double s2,long double t1,long double t2);
long double Cpart3(long double s,long double s1,long double s2,long double t1,long double t2);
long double Cpart4(long double s,long double s1,long double s2,long double t1,long double t2);
long double Cpart5(long double s,long double s1,long double s2,long double t1,long double t2);
long double Cpart6(long double s,long double s1,long double s2,long double t1,long double t2);
long double Cpart7(long double s,long double s1,long double s2,long double t1,long double t2);


double C(double s,double s1,double s2,double t1,double t2){
    long double sl=s, s1l=s1, s2l=s2, t1l=t1, t2l=t2, a ;    
    long double a1,a2,a3,a4,a5,a6,a7;
    a1=Cpart1(sl,s1l,s2l,t1l,t2l);
    a2=Cpart2(sl,s1l,s2l,t1l,t2l);
    a3=Cpart3(sl,s1l,s2l,t1l,t2l);
    a4=Cpart4(sl,s1l,s2l,t1l,t2l);
    a5=Cpart5(sl,s1l,s2l,t1l,t2l);
    a6=Cpart6(sl,s1l,s2l,t1l,t2l);
    a7=Cpart7(sl,s1l,s2l,t1l,t2l);
      
    a=a1+a2+a3+a4+a5+a6+a7;
    //if(a>10000000.) printf("C=%Le a1=%Le  a2=%Le  a3=%Le  a4=%Le  a5=%Le  a6=%Le a7=%Le\n",a,a1,a2,a3,a4,a5,a6,a7);
    return a;
};