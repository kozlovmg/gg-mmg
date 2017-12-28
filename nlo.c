#include <stdio.h>

double IK(double s,double s1,double s2,double t1,double t2,double EE,double np1,double np2);
double C(double s,double s1,double s2,double t1,double t2);
double Box1(double s,double s1,double s2,double t1,double t2);
double Box2(double s,double s1,double s2,double t1,double t2);
double Box3(double s,double s1,double s2,double t1,double t2);
double Tri(double s,double s1,double s2,double t1,double t2);
double Rr(double s,double s1,double s2,double t1,double t2, double EE);

double lbl11(double,double,double,double,double,double);
double lbl12(double,double,double,double,double,double);
double lbl13(double,double,double,double,double,double);

double lbl21(double,double,double,double,double,double);
double lbl22(double,double,double,double,double,double);
double lbl23(double,double,double,double,double,double);

double lbl31(double,double,double,double,double,double);
double lbl32(double,double,double,double,double,double);
double lbl33(double,double,double,double,double,double);

double lbl41(double,double,double,double,double,double);
double lbl42(double,double,double,double,double,double);
double lbl43(double,double,double,double,double,double);

double lbl51(double,double,double,double,double,double);
double lbl52(double,double,double,double,double,double);
double lbl53(double,double,double,double,double,double);

double lbl61(double,double,double,double,double,double);
double lbl62(double,double,double,double,double,double);
double lbl63(double,double,double,double,double,double);



double nlo1(double *vars){
    double s=vars[0], s1=vars[1], s2=vars[2], t1=vars[3], t2=vars[4], np1=vars[5], np2=vars[6], EE=vars[7];
    double a,b1,b2,b3,tr,rr,cc,ik;
    b1=Box1(s,s1,s2,t1,t2);
    b2=Box2(s,s1,s2,t1,t2);
    b3=Box3(s,s1,s2,t1,t2);
    tr=Tri(s,s1,s2,t1,t2);
    rr=Rr(s,s1,s2,t1,t2,EE);
    cc=C(s,s1,s2,t1,t2);
    ik=IK(s,s1,s2,t1,t2,EE,np1,np2);
    a=(b1+b2+b3+tr+rr+cc+ik)/4.;// "/4" because nlo2=2*Re(Mnlo*(Mborn)^*)/8
    /*if(a>10000000.)
    { 
        printf("nlo2 function: b1=%e  b2=%e  b3=%e  tr=%e  rr=%e  cc=%e ik=%e\n",b1,b2,b3,tr,rr,cc,ik);
        printf("u1=%e  u2=%e\n",-(s+t2-s1-1),-(s+t1-s2-1));
    };*/
    
    return a;
};


double nlo2(double *vars){
    double s=vars[0], s1=vars[1], s2=vars[2], t1=vars[3], t2=vars[4], m=vars[5];
    double r1,r2,r3,r4,r5,r6,a;
    r1=lbl11(s,s1,s2,t1,t2,m)+lbl12(s,s1,s2,t1,t2,m)+lbl13(s,s1,s2,t1,t2,m);
    r2=lbl21(s,s1,s2,t1,t2,m)+lbl22(s,s1,s2,t1,t2,m)+lbl23(s,s1,s2,t1,t2,m);
    r3=lbl31(s,s1,s2,t1,t2,m)+lbl32(s,s1,s2,t1,t2,m)+lbl33(s,s1,s2,t1,t2,m);
    r4=lbl41(s,s1,s2,t1,t2,m)+lbl42(s,s1,s2,t1,t2,m)+lbl43(s,s1,s2,t1,t2,m);
    r5=lbl51(s,s1,s2,t1,t2,m)+lbl52(s,s1,s2,t1,t2,m)+lbl53(s,s1,s2,t1,t2,m);
    r6=lbl61(s,s1,s2,t1,t2,m)+lbl62(s,s1,s2,t1,t2,m)+lbl63(s,s1,s2,t1,t2,m);
    a=(r1+r2+r3+r4+r5+r6)/4.;
    return a;
};
