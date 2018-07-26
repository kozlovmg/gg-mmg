/* nlo_functions.h
 *
 *
 *
 */



long double Li2(long double x);
long double PolyLog(int,long double);
long double Log(long double);
__float128 L128(__float128);
long double Sqrt(long double);
double Power(double,int);

long double ht(long double);
long double R(long double,long double);
long double S3(long double, long double, long  double);
long double ff(long double x, long double y);
long double bb(long double x);
long double X(long double a);

long double pq(long double,unsigned int);
__float128 p128(__float128,unsigned int);


//Tests
int cuts(double, double, double,double, double);
int Tphys(double,double,double,double,double);
int inv_test(double, double, double, double, double);


/*Masters*/
//Bubles
long double R1(long double);
long double R2(long double);
long double R1n(long double,long double);
long double R2n(long double);

//Triangles
long double T1(long double);
long double T2(long double,long double);
long double T3(long double,long double);
long double T4(long double);
long double T5(long double);
long double T5n(long double);
long double T6(long double,long double);
long double T7(long double,long double);
long double T3n(long double,long double);


//Boxes
double B1(double s,double t,double m);
long double B2(long double s,long double t,long double m);
double B3(double s,double t,double m);
double B3n(double s,double t,double m);



//IR functions
double F(double x);
double integrand(double x, void * params);
double G(double a, double b, double c);


double GG(double x,double y,double z,double u,double v,double w);


