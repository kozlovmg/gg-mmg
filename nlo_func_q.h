/* nlo_func_q.h
 *
 *
 *
 */


__float128 Li2q(__float128 x);
__float128 Logq(__float128);
__float128 Sqrtq(__float128);


__float128 htq(__float128);
__float128 Rq(__float128,__float128);
__float128 S3q(__float128, __float128,__float128);
__float128 ffq(__float128, __float128);
__float128 bbq(__float128);
__float128 Xq(__float128);


__float128 p128(__float128,unsigned int);


/*Masters*/
//Bubles
__float128 R1q(__float128);
__float128 R2q(__float128);
__float128 R1nq(__float128,__float128);
__float128 R2nq(__float128);

//Triangles
__float128 T1q(__float128);
__float128 T2q(__float128,__float128);
__float128 T3q(__float128,__float128);
__float128 T4q(__float128);
__float128 T5q(__float128);
__float128 T5nq(__float128);
__float128 T6q(__float128,__float128);
__float128 T7q(__float128,__float128);
__float128 T3nq(__float128,__float128);


//Boxes

__float128 B2q(__float128 s,__float128 t,__float128 m);

double B1(double s,double t,double m);
double B3(double s,double t,double m);
double B3n(double s,double t,double m);

#define Pi M_PIq
#define Power p128
#define Log Logq

#define R1  R1q
#define R2  R2q
#define R1n  R1nq
#define R2n  R2nq
#define T1  T1q
#define T2  T2q
#define T3  T3q
#define T3n  T3nq
#define T4  T4q
#define T5  T5q
#define T5n  T5nq
#define T6  T6q
#define T7  T7q
#define B2  B2q
