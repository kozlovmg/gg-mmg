
double Cpart1(double s,double s1,double s2,double t1,double t2);
double Cpart2(double s,double s1,double s2,double t1,double t2);
double Cpart3(double s,double s1,double s2,double t1,double t2);
double Cpart4(double s,double s1,double s2,double t1,double t2);
double Cpart5(double s,double s1,double s2,double t1,double t2);
double Cpart6(double s,double s1,double s2,double t1,double t2);


double C(double s,double s1,double s2,double t1,double t2,double EE,double np1,double np2){
        double a=Cpart1(s,s1,s2,t1,t2)+Cpart2(s,s1,s2,t1,t2)+Cpart3(s,s1,s2,t1,t2)+Cpart4(s,s1,s2,t1,t2)+Cpart5(s,s1,s2,t1,t2)+Cpart6(s,s1,s2,t1,t2);
        return a;
};