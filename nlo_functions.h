/* nlo_functions.h
 *
 *
 *
 */


double Li2(double x);
double ht(double x);
double R(double x, double y);
double S3(double x, double x1, double x2);
double ff(double x, double y);
double bb(double x);
double X(double a);

long double pw_n(long double,unsigned int);

/*Masters*/
//Bubles
long double R1(long double x);
long double R2(long double x);

//Triangles
double T1(double x);
double T2(double t1,double t2);
double T3(double t1,double t2);
double T4(double x);
double T5(double x);
double T6(double s,double t);

//Boxes
double B1(double s,double t,double m);
double B2(double s,double t,double m);
double B3(double s,double t,double m);



//IR functions
double F(double x);
double integrand(double x, void * params);
double G(double a, double b, double c);