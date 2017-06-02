

double IK(double s,double s1,double s2,double t1,double t2,double EE,double np1,double np2);
double C(double s,double s1,double s2,double t1,double t2,double EE,double np1,double np2);
double Box1(double s,double s1,double s2,double t1,double t2);
double Box2(double s,double s1,double s2,double t1,double t2);
double Box3(double s,double s1,double s2,double t1,double t2);
double Tri(double s,double s1,double s2,double t1,double t2);
double Rr(double s,double s1,double s2,double t1,double t2, double EE);



double nlo2(double s,double s1,double s2,double t1,double t2,double EE,double np1,double np2){
    double a=C(s,s1,s2,t1,t2,EE,np1,np2)+IK(s,s1,s2,t1,t2,EE,np1,np2)+Box1(s,s1,s2,t1,t2)+Box2(s,s1,s2,t1,t2)+Box3(s,s1,s2,t1,t2)+Tri(s,s1,s2,t1,t2)+Rr(s,s1,s2,t1,t2,EE);
    return a;
};