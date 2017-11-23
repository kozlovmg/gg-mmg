
long double Cpart1(long double s,long double s1,long double s2,long double t1,long double t2);
long double Cpart2(long double s,long double s1,long double s2,long double t1,long double t2);
long double Cpart3(long double s,long double s1,long double s2,long double t1,long double t2);
long double Cpart4(long double s,long double s1,long double s2,long double t1,long double t2);
long double Cpart5(long double s,long double s1,long double s2,long double t1,long double t2);
long double Cpart6(long double s,long double s1,long double s2,long double t1,long double t2);


double C(double s,double s1,double s2,double t1,double t2,double EE,double np1,double np2){
    
    long double sl=s, s1l=s1, s2l=s2, t1l=t1, t2l=t2;    
    double a=Cpart1(sl,s1l,s2l,t1l,t2l)+Cpart3(sl,s1l,s2l,t1l,t2l)+Cpart3(sl,s1l,s2l,t1l,t2l)+Cpart4(sl,s1l,s2l,t1l,t2l)+Cpart5(sl,s1l,s2l,t1l,t2l)+Cpart6(sl,s1l,s2l,t1l,t2l);
    return a;
};