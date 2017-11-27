#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_pow_int.h>
#include "nlo_functions.h"


#define Pi M_PI
#define Power gsl_pow_int



double nlo2(double *);
double born(double *);

double test(double *);
int MK(double *, double ( *)(double *),long,double *);
double PhaseVol(double, double, double, double);
double t2_func(double,double,double,double,double);
int HevTheta(double);
double Gk(double,double,double,double,double,double);
double f_lambda(double, double, double);

double B1(double,double,double);
double B1N(double,double,double);

double C(double s,double s1,double s2,double t1,double t2);



int main(int argc, char *argv[])
{
    
    double mu, EminG, Ecut, Ecm;
    mu=105.6583745;     // Muon mass (MeV)
    EminG=50./mu;	// Minimum energy of gamma
    Ecut=10./mu;		// Cut for invisible gamma
    Ecm=500./mu;	// Energy of center mass
    
    long N=100000;
    

    
    

    double s = 4*Power(Ecm,2);
    double params[3];
    double result_born[2],result_nlo[2];
    
    params[0] = s;
    params[1] = EminG;
    params[2] = Ecut;
    
    printf("s=%f\n",s);
    
    printf("calc born\n");
    MK(params,born,N,result_born);
    
    printf("%e  %e\n",result_born[0],result_born[1]);

    printf("calc next to leading order\n");
    MK(params,nlo2,N,result_nlo);
      
    
    
    printf("%e  %e\n",result_nlo[0],result_nlo[1]);
    
    /*
    const gsl_rng_type * T;
    gsl_rng * r;
    
    gsl_rng_env_setup();

    T = gsl_rng_ranlxd2;
    r = gsl_rng_alloc (T);
    double s,s1,s2,t1,t2;
    double s1_min, s1_max, s2_min, s2_max, t1_min, t1_max, t2_min, t2_max;
    double x1,x2,x3,x4;
    int i;
    
    s=4*Power(500/105.,2);
    
    s1_min = 1;
    s1_max = Power(sqrt(s)-1,2);
    
    s2_min = 1;
    s2_max = Power(sqrt(s)-1,2);
    
    t1_min = (2-s-sqrt(s*(s-2)))/2;
    t1_max = (2-s+sqrt(s*(s-2)))/2;
    
    t2_min = (2-s-sqrt(s*(s-2)))/2;
    t2_max = (2-s+sqrt(s*(s-2)))/2;
    
    
    
    for(i=0;i<10;i++)
    {
        x1= gsl_rng_uniform(r);
	x2= gsl_rng_uniform(r);
	x3= gsl_rng_uniform(r);
	x4= gsl_rng_uniform(r);
	
	s1 = x1*s1_max + (1-x1)*s1_min;
	s2 = x2*s2_max + (1-x2)*s2_min;
	t1 = x3*t1_max + (1-x3)*t1_min;
        t2 = x4*t2_max + (1-x4)*t2_min;
    
    
    
    
    
    
    printf("B1[%f,%f,%f]= %e\n",1-(s-s1+t2),s,s1,B1(1-(s-s1+t2),s,s1) );
    //printf("B1[%f,%f,%f]= %e\n",1-(s-s1+t2),s,s1,B1N(1-(s-s1+t2),s,s1) );
    printf("B1[%f,%f,%f]= %e\n\n",s1,-(s2-t1+t2-1),1-(s-s1+t2),B1(s1,-(s2-t1+t2-1),1-(s-s1+t2)) );
    //printf("B1[%f,%f,%f]= %e\n\n",s1,-(s2-t1+t2-1),1-(s-s1+t2),B1N(s1,-(s2-t1+t2-1),1-(s-s1+t2)) );
   
    
    };
   
*/
};





//params={s/mu^2,EminG/mu,Ecut/mu}
int MK(double *params, double ( *fun)(double *),long N, double *result)
{
    double s, Em, Ec, s1, s2, t1, t2, lam, np1, np2 ;
    double s1_min, s1_max, s2_min, s2_max, t1_min, t1_max;
    double x1, x2, x3, x4;
    double nfunc,pvol;
    double *points;
    double *vars;
    double volume;
    double Mtmp=0, Stmp=0;
    long i=0,k=0,j=0, errs=0;
    //double a;
    
    points = (double *) malloc(N*sizeof(double));
    vars = (double *) malloc(8*sizeof(double));
    
    s=params[0];
    Em=params[1];
    Ec=params[2];
    
    /*initialisation of random generator*/
    
    const gsl_rng_type * T;
    gsl_rng * r;
    
    gsl_rng_env_setup();

    T = gsl_rng_ranlxd2;
    r = gsl_rng_alloc (T);
    
    gsl_rng_set(r, 10);
    
    /*end initialisation*/
   
    
    s1_min = 1;
    s1_max = Power(sqrt(s)-1,2);
    s2_min = 1;
    s2_max = Power(sqrt(s)-1,2);
    
    t1_min = (2-s-sqrt(s*(s-2)))/2;
    t1_max = (2-s+sqrt(s*(s-2)))/2;
    
    volume = (s1_max-s1_min)*(s2_max-s2_min)*(t1_max-t1_min)*2*Pi;
    
    // start calculation 
    vars[0] = s;
    vars[7] = Ec;
    
    
    
    for(i=0;i<N;i++)
    {
        x1= gsl_rng_uniform(r);
	x2= gsl_rng_uniform(r);
	x3= gsl_rng_uniform(r);
	x4= gsl_rng_uniform(r);
	
	s1 = x1*s1_max + (1-x1)*s1_min;
	s2 = x2*s2_max + (1-x2)*s2_min;
	t1 = x3*t1_max + (1-x3)*t1_min;
	lam = x4*2*Pi;
        
	
        np1 = (s-s2+1)/sqrt(s)/2;
        np2 = (s-s1+1)/sqrt(s)/2;
        
                
        vars[1] = s1;
        vars[2] = s2;
        vars[3] = t1;
        vars[5] = np1;
        vars[6] = np2;
        
        
        
	if( (s1+s2-2)/sqrt(s)>=Em && Gk(s,t1,s2,0.,0.,1.)<0 && Gk(s1,s2,s,0.,1.,1.)<0) 
        {
            k=k+1;
            t2 = t2_func(s,s1,s2,t1,lam);
            vars[4] = t2;
            nfunc = fun(vars);
            pvol=PhaseVol(s,s1,s2,t1);
            points[i] = nfunc*pvol*volume;
            if(nfunc>10000000.) j=j+1;
            if(isnanl(points[i])) errs=errs+1;
            //if(nfunc>10000000.) printf("point=%e, s1=%e, s2=%e, t1=%e, t2=%e\n",nfunc,s1,s2,t1,t2);
        }
        else points[i]=0;
        
    };
    
    printf("%ld/%ld , >10^7:  %ld, errors: %ld\n",k,N,j,errs);
    // Calc mean value
    
    for(i=0;i<N;i++)
    {
	Mtmp = Mtmp + points[i];
    };
	
    Mtmp = Mtmp/N;
    //printf("Mtmp=%e\n",Mtmp);
    
    result[0] = Mtmp;
    //Calc sigma2
    

    for(i=0;i<N;i++)
    {
	Stmp = Stmp + Power(points[i]-Mtmp,2);
    };

    Stmp = sqrt( Stmp/(N*(N-1)) );
    //printf("Stmp=%e\n",Stmp);
    result[1] = Stmp;

    // result is Mval +/- sigma2
    
    
    
    //free memory
    free(points);
    free(vars);
    
    // free random
    gsl_rng_free (r);
    
    return 1;
};






double PhaseVol(double s,double s1, double s2, double t1)
{
    return Pi/( 8.*sqrt(f_lambda(s,0.,0.)*f_lambda(s,s2,1.)) );
};


double t2_func(double s,double s1,double s2,double t1,double lam)
{
	return 1.+(2*sqrt(Gk(s,t1,s2,0.,0.,1.)*Gk(s1,s2,s,0.,1.,1.))*cos(lam) + (-1 + s1)*(-1 + s2)*(s2 - t1) + Power(s,2)*(-1 + t1) - s*(-1 + (2 + s1)*t1 + s2*(-4 + s1 + t1)) )/f_lambda(s,s2,1.)  ;
};

int HevTheta(double a)
{
    if (a<0) return 0;
    else return 1;
};

double Gk(double x,double y,double z,double u,double v,double w)
{
    return Power(v,2)*w + Power(u,2)*z + u*((w - x)*(-v + y) - (v + w + x + y)*z + Power(z,2)) + x*(y*(x + y - z) + w*(-y + z)) + v*(Power(w,2) + y*(-x + z) - w*(x + y + z)) ;
};

double f_lambda(double x, double y, double z)
{
    return Power(x-y-z,2)-4*y*z;
};





