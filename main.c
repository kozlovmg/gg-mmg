#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_pow_int.h>
#include "nlo_functions.h"


#define Pi M_PI
#define Power gsl_pow_int



double nlo1(double *);
double nlo2(double *);
double born(double *);

double test(double *);
int MK(double *, double ( *)(double *),long,double *);
int MK2(double *, double ( *)(double *),long,double *);
int MK_t(double *, double ( *)(double *),long,double *);
double PhaseVol(double, double, double, double);
double t2_func(double,double,double,double,double);
int HevTheta(double);
double Gk(double,double,double,double,double,double);
double f_lambda(double, double, double);
double delta4(double,double,double,double,double);

double B1(double,double,double);
double B1N(double,double,double);

double C(double s,double s1,double s2,double t1,double t2);

void display_results(char *, double, double);
double born_mc(double *, size_t, void *);
double born_mc_v3(double *, size_t, void *);

double nlo1_mc(double *, size_t, void *);
double nlo2_mc(double *, size_t, void *);




int main(int argc, char *argv[])
{
    
    double mu, EminG, Ecut, Ecm;
    mu=105.6583745;     // Muon mass (MeV)
    EminG=50./mu;	// Minimum energy of gamma
    Ecut=10./mu;		// Cut for invisible gamma
    Ecm=500./mu;	// Energy of center mass
    
    double s = 4*Power(Ecm,2);
    
    //long N=100000;
    
    double si_min, si_max, t1_min, t1_max;
    
    si_min = 1.;
    si_max = Power(sqrt(s)-1.,2);
    t1_min = (2-s-sqrt(s*(s-2)))/2;
    t1_max = (2-s+sqrt(s*(s-2)))/2;
    
    
    double res, err;
    
    // integration limits
    double xl[4], xu[4];
    // s1
    xl[0] = si_min;
    xu[0] = si_max;
    // s2
    xl[1] = si_min;
    xu[1] = si_max;
    // t1
    xl[2] = t1_min;
    xu[2] = t1_max;
    // lambda
    xl[3] = 0.;
    xu[3] = 2.*Pi;
    
    //t2
    //xl[3] = t1_min;
    //xu[3] = t1_max;
    
    
    
    // calls of function
    size_t N=50000;
    
    // initialisation of random numbers generator
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    
    gsl_monte_function func;
    
    double params[3];
    params[0] = s;
    params[1] = EminG;
    params[2] = Ecut;
    
    func.f = &born_mc;
    //func.f = &nlo1_mc;
    //func.f = &nlo2_mc;
    func.dim = 4;
    func.params = &params;
    
      
    // monte carlo 
    
    
    /*
    gsl_monte_miser_state *state = gsl_monte_miser_alloc (4);
    gsl_monte_miser_integrate (&func, xl, xu, 4, N, r, state,&res, &err);
    gsl_monte_miser_free (state);

    display_results ("miser", res, err);
  */
    
    gsl_monte_vegas_state *state = gsl_monte_vegas_alloc (4);

    gsl_monte_vegas_integrate (&func, xl, xu, 4, 50000, r, state, &res, &err);
    display_results ("vegas warm-up", res, err);

    printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (&func, xl, xu, 4, N, r, state, &res, &err);
        printf ("result = % .6f sigma = % .6f\n chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq(state));
      }
    while ( fabs(gsl_monte_vegas_chisq (state) - 1.0) > 0.5 );

    display_results ("vegas final", res, err);

    gsl_monte_vegas_free (state);
   
    
    
    gsl_rng_free (r);

    return 0;
};


void display_results(char *title, double result, double error)
{
  printf ("%s ==================\n", title);
  printf ("result = % .6f\n", result);
  printf ("sigma  = % .6f\n", error);
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
        
        
        
	if( (s1+s2-2)/sqrt(s)>=Em && Gk(s,t1,s2,0.,0.,1.)<0 && Gk(s1,s2,s,0.,1.,1.)<0 ) 
        {
            t2 = t2_func(s,s1,s2,t1,lam);
            if( (s+t2-s1-1.)*(s+t1-s2-1.)>1. && (-t1)*(-t2)>1. )
            {
                k=k+1;
                vars[4] = t2;
                nfunc = fun(vars);
                pvol=PhaseVol(s,s1,s2,t1);
                points[i] = nfunc*pvol*volume;
                if(nfunc>10000000.) j=j+1;
                //if(isnanl(points[i])) errs=errs+1;
                if(nfunc>10000000.) printf("point=%e, s1=%e, s2=%e, t1=%e, t2=%e\n\n",nfunc,s1,s2,t1,t2);
            };
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

double PhaseVol2(double s,double s2, double t1)
{
    return Pi/( 8.*sqrt(f_lambda(s,0.,0.)*f_lambda(s2,0.,t1)) );
};

double PhaseVol3(double s,double s1, double s2, double t1, double t2)
{
    return Pi/(16.*sqrt(f_lambda(s,0.,0.))*sqrt(-delta4(s,s1,s2,t1,t2)));
    
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



double s1_func(double s,double s2,double t1,double t2,double ph)
{
	return s+1.-( 2*sqrt(Gk(s,t1,s2,0.,0.,1.)*Gk(s2,t2,1.,t1,0.,0.))*cos(ph) + t1 - s2*(1 + t1 + s*t1) + Power(s2,2)*(1 + s - t2) + s*t1*(-2 + t2) - t1*t2 + s2*(1 + s + t1)*t2 )/f_lambda(s2,t1,0)  ;
};

double delta4(double s,double s1,double s2,double t1,double t2)
{
    return (Power(s1,2)*Power(s2 - t1,2) + Power(1 + (-1 + s)*t1,2) + 
   2*(-1 + s + s2 - 2*s*s2 - (-1 + Power(s,2) + s2 - s*(2 + s2))*t1)*t2 + 
   (Power(s,2) + Power(-1 + s2,2) - 2*s*(1 + s2))*Power(t2,2) + 
   2*s1*(-(Power(s2,2)*t2) - t1*(-1 + t1 + s*(2 + t1 - t2) + t2) + 
      s2*(-1 + t1 + s*t1 + (1 + s + t1)*t2)))/16.;
};



double born_mc(double *inv, size_t dim, void *params)
{
    (void)(dim);
    double *pars = (double *)params;
    double vars[5];
    double s1,s2,t1,lam; // variables
    double t2; // temporary variables
    double s, Em;
    double result;
        
    s1=inv[0];
    s2=inv[1];
    t1=inv[2];
    lam=inv[3];
    
    s=pars[0];
    Em=pars[1];
           
    if( (s1+s2-2)/sqrt(s)>=Em && Gk(s,t1,s2,0.,0.,1.)<0 && Gk(s1,s2,s,0.,1.,1.)<0 ) 
    {
        t2 = t2_func(s,s1,s2,t1,lam);
        vars[0] = s;
        vars[1] = s1;
        vars[2] = s2;
        vars[3] = t1;
        vars[4] = t2;
        
        result = born(vars)*PhaseVol(s,s1,s2,t1)/( 2*s*Power(2*Pi,5) );        
    }
    else
    {
        result = 0.;
    };
    
    return result;

};

double born_mc_v3(double *inv, size_t dim, void *params)
{
    (void)(dim);
    double *pars = (double *)params;
    double vars[5];
    double s1,s2,t1,t2; // variables
    double s, Em;
    double result;
        
    s1=inv[0];
    s2=inv[1];
    t1=inv[2];
    t2=inv[3];
    
    s=pars[0];
    Em=pars[1];
           
    if( (s1+s2-2)/sqrt(s)>=Em && delta4(s,s1,s2,t1,t2)<-0.2 ) 
    {
        vars[0] = s;
        vars[1] = s1;
        vars[2] = s2;
        vars[3] = t1;
        vars[4] = t2;
        //result = PhaseVol3(s,s1,s2,t1,t2)/( 2*s*Power(2*Pi,5) ); 
        result = born(vars)*PhaseVol3(s,s1,s2,t1,t2)/( 2*s*Power(2*Pi,5) );        
    }
    else
    {
        result = 0.;
    };
    
    return result;

};





double nlo1_mc(double *inv, size_t dim, void *params)
{
    (void)(dim);
    double *pars = (double *)params;
    double vars[8];
    double s1,s2,t1,lam; // variables
    double t2, np1, np2; // temporary variables
    double s, Ec, Em;
    double result;
        
    s1=inv[0];
    s2=inv[1];
    t1=inv[2];
    lam=inv[3];
    
    s=pars[0];
    Em=pars[1];
    Ec=pars[2];
    np1 = (s-s2+1)/sqrt(s)/2;
    np2 = (s-s1+1)/sqrt(s)/2;
        
    if( (s1+s2-2)/sqrt(s)>=Em && Gk(s,t1,s2,0.,0.,1.)<0 && Gk(s1,s2,s,0.,1.,1.)<0  ) 
    {
        t2 = t2_func(s,s1,s2,t1,lam);
        if( abs((s+t2-s1-1)-(s+t1-s2-1))>1. && abs(t1-t2)>1. )
        {
            vars[0] = s;
            vars[1] = s1;
            vars[2] = s2;
            vars[3] = t1;
            vars[4] = t2;
            vars[5] = np1;
            vars[6] = np2;
            vars[7] = Ec;
            result = nlo1(vars)*PhaseVol(s,s1,s2,t1)/( 2*s*Power(2*Pi,5) );        
        }
        else result=0.;
        
        
    }
    else
    {
        result = 0.;
    };
    
    return result;

};

double nlo2_mc(double *inv, size_t dim, void *params)
{
    (void)(dim);
    double *pars = (double *)params;
    double vars[6];
    double s1,s2,t1,lam; // variables
    double m, t2, np1, np2; // temporary variables
    double s, Ec, Em;
    double result;
        
    s1=inv[0];
    s2=inv[1];
    t1=inv[2];
    lam=inv[3];
    
    s=pars[0];
    Em=pars[1];
    Ec=pars[2];
    np1 = (s-s2+1)/sqrt(s)/2;
    np2 = (s-s1+1)/sqrt(s)/2;
        
    if( (s1+s2-2)/sqrt(s)>=Em && Gk(s,t1,s2,0.,0.,1.)<0 && Gk(s1,s2,s,0.,1.,1.)<0  ) 
    {
        t2 = t2_func(s,s1,s2,t1,lam);
        vars[0] = s;
        vars[1] = s1;
        vars[2] = s2;
        vars[3] = t1;
        vars[4] = t2;
        vars[5] = m;
        result = nlo2(vars)*PhaseVol(s,s1,s2,t1)/( 2*s*Power(2*Pi,5) );        
    }
    else
    {
        result = 0.;
    };
    
    return result;

};