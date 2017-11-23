#include <stdio.h>

double IK(double s,double s1,double s2,double t1,double t2,double EE,double np1,double np2);
double C(double s,double s1,double s2,double t1,double t2,double EE,double np1,double np2);
double Box1(double s,double s1,double s2,double t1,double t2);
double Box2(double s,double s1,double s2,double t1,double t2);
double Box3(double s,double s1,double s2,double t1,double t2);
double Tri(double s,double s1,double s2,double t1,double t2);
double Rr(double s,double s1,double s2,double t1,double t2, double EE);



double nlo2(double *vars){
    double s=vars[0], s1=vars[1], s2=vars[2], t1=vars[3], t2=vars[4], np1=vars[5], np2=vars[6], EE=vars[7];
    double a,b1,b2,b3,tr,rr,cc,ik;
    b1=Box1(s,s1,s2,t1,t2);
    b2=Box2(s,s1,s2,t1,t2);
    b3=Box3(s,s1,s2,t1,t2);
    tr=Tri(s,s1,s2,t1,t2);
    rr=Rr(s,s1,s2,t1,t2,EE);
    cc=C(s,s1,s2,t1,t2,EE,np1,np2);
    ik=IK(s,s1,s2,t1,t2,EE,np1,np2);
    a=b1+b2+b3+tr+rr+cc+ik;
    //if(a>10000000.) printf("nlo2 function: b1=%e  b2=%e  b3=%e  tr=%e  rr=%e  cc=%e ik=%e\n",b1,b2,b3,tr,rr,cc,ik);
    //a=Box2(s,s1,s2,t1,t2)+Box3(s,s1,s2,t1,t2)+Tri(s,s1,s2,t1,t2)+Rr(s,s1,s2,t1,t2,EE)+C(s,s1,s2,t1,t2,EE,np1,np2)+Box1(s,s1,s2,t1,t2)+IK(s,s1,s2,t1,t2,EE,np1,np2);
    return a;
};