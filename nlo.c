#include <stdio.h>

double IK(double s,double s1,double s2,double t1,double t2,double EE,double np1,double np2);
double C(double s,double s1,double s2,double t1,double t2);
double Box1(double s,double s1,double s2,double t1,double t2);
double Box2(double s,double s1,double s2,double t1,double t2);
double Box3(double s,double s1,double s2,double t1,double t2);
double Tri(double s,double s1,double s2,double t1,double t2);
double Rr(double s,double s1,double s2,double t1,double t2, double EE);

double lbl1(double,double,double,double,double,double);
double lbl2(double,double,double,double,double,double);
double lbl3(double,double,double,double,double,double);
double lbl4(double,double,double,double,double,double);
double lbl5(double,double,double,double,double,double);
double lbl6(double,double,double,double,double,double);



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
    r1=lbl1(s,s1,s2,t1,t2,m);
    r2=lbl2(s,s1,s2,t1,t2,m);
    r3=lbl3(s,s1,s2,t1,t2,m);
    r4=lbl4(s,s1,s2,t1,t2,m);
    r5=lbl5(s,s1,s2,t1,t2,m);
    r6=lbl6(s,s1,s2,t1,t2,m);
    a=(r1+r2+r3+r4+r5+r6)/4.;
    return a;
};
