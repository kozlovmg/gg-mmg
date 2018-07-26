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
double nlo3(double *);
double nlo_qcd(double *);
double born(double *);


double PhaseVol(double, double, double, double);
double t2_func(double,double,double,double,double);
int HevTheta(double);
double Gk(double,double,double,double,double,double);
double f_lambda(double, double, double);
double delta4(double,double,double,double,double);



double C(double s,double s1,double s2,double t1,double t2);

void display_results(char *, double, double);

int save(int,double *);



double born_mc(double *, size_t, void *);
double nlo1_mc(double *, size_t, void *);
double nlo2_mc(double *, size_t, void *);
double nlo3_mc(double *, size_t, void *);
double nlo_qcd_mc(double *, size_t, void *);






int main(int argc, char *argv[])
{
    size_t N;
    double mu, Ecut, Ecm, me, Const, mtau, mQ;
    mu=105.6583745;     // Muon mass (MeV)
    me=0.511/mu;       // Electron mass 
    mtau=1776.84/mu;    // Tau mass
    //m_pi=134.9766/mu;
    //m_u=2./mu;
    //m_d=4.5/mu;
    //m_s=104./mu;
    //m_c=1270./mu;
    //m_b=4200./mu;
    mQ=300./mu;
    Ecut=30./mu;	// Cut for invisible gamma
    
    //Const=Power(4*Pi/137.,3)/Power(mu,2)*0.00257*Power(10.,12); //  pb
    Const=24./137.*Power(Pi*me,2)*Power(10.,9)*0.6652; //nb
    
    if(argc==1)
    {
        N=50000;
        Ecm=500.;// Energy of center mass
    };
        
    if(argc==2)
    {
        N=50000;
        Ecm=atof(argv[1]);
        //printf("Ecm=%e\n",Ecm);
    };
    
    if(argc>2)
    {
        Ecm=atof(argv[1]);
        N=atoi(argv[2]);
        //printf("Ecm=%e  N=%i\n",Ecm,N);
    };
    

    
    double s = 4*Power(Ecm/mu,2);
    double si_cut=2.;
    
    double si_min, si_max, t1_min, t1_max, si_min_2, si_max_2;
    
    si_min = si_cut;
    si_max = Power(sqrt(s)-1.,2);
    t1_min = (2-s-sqrt(s*(s-2)))/2;
    t1_max = (2-s+sqrt(s*(s-2)))/2;
    
    si_min_2 = 1.;
    si_max_2 = si_cut;
    
    

    
    // integration limits
    double xl[4], xu[4], xl_1[4], xu_1[4], xl_2[4], xu_2[4], xl_nlo2[4], xu_nlo2[4];
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
    
    // s1<s_cut
    //=======================
    // s1
    xl_1[0] = si_min_2;
    xu_1[0] = si_max_2;
    // s2
    xl_1[1] = si_min;
    xu_1[1] = si_max;
    // t1
    xl_1[2] = t1_min;
    xu_1[2] = t1_max;
    // lambda
    xl_1[3] = 0.;
    xu_1[3] = 2.*Pi;
    //======================
    
    // s2<s_cut
    //=======================
    // s1
    xl_2[0] = si_min;
    xu_2[0] = si_max;
    // s2
    xl_2[1] = si_min_2;
    xu_2[1] = si_max_2;
    // t1
    xl_2[2] = t1_min;
    xu_2[2] = t1_max;
    // lambda
    xl_2[3] = 0.;
    xu_2[3] = 2.*Pi;
    //======================
    
    //===== for nlo2 ====================================
    
    
    // s1
    xl_nlo2[0] = 1.;
    xu_nlo2[0] = si_max;
    // s2
    xl_nlo2[1] = 1.;
    xu_nlo2[1] = si_max;
    // t1
    xl_nlo2[2] = t1_min;
    xu_nlo2[2] = t1_max;
    // lambda
    xl_nlo2[3] = 0.;
    xu_nlo2[3] = 2.*Pi;
    
    
    
    
    // initialisation of random numbers generator
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,10);
    
    
    double res=0., err=0., res1=0., err1=0., res2=0., err2=0.;
    double results[30];
    double params[4];
    gsl_monte_function func;
    gsl_monte_vegas_state *state;
    
    func.dim = 4;
    func.params = &params;
    
    params[0] = s; 
    params[1] = Ecut;
    
    //################## BORN #########################
    
    func.f = &born_mc;
    
    state = gsl_monte_vegas_alloc(4);
    gsl_monte_vegas_integrate (&func, xl_nlo2, xu_nlo2, 4, 10000, r, state, &res, &err);
    gsl_monte_vegas_integrate (&func, xl_nlo2, xu_nlo2, 4, N, r, state, &res, &err);
    
    gsl_monte_vegas_free(state);
    
    //printf("si>2mu: %.6f + pm %.6f nb\n",res*Const,err*Const);
    
       
    //state = gsl_monte_vegas_alloc(4);
    //gsl_monte_vegas_integrate (&func, xl_1, xu_1, 4, 10000, r, state, &res1, &err1);
    //gsl_monte_vegas_integrate (&func, xl_1, xu_1, 4, N, r, state, &res1, &err1);
    
    //gsl_monte_vegas_free(state);
    //printf("s1<2mu: %.6f + pm %.6f nb\n",res1*Const,err1*Const);
    
    //state = gsl_monte_vegas_alloc(4);
    //gsl_monte_vegas_integrate (&func, xl_2, xu_2, 4, 10000, r, state, &res2, &err2);
    //gsl_monte_vegas_integrate (&func, xl_2, xu_2, 4, N, r, state, &res2, &err2);
    
    //gsl_monte_vegas_free(state);
    //printf("s2<2mu: %.6f + pm %.6f nb\n",res2*Const,err2*Const);
    
    results[0]=Ecm;
    results[1]=res;
    results[2]=err;
    
    
    printf("==========Born============\n");
    printf("%.6f + pm %.6f nb\n",results[1]*Const,results[2]*Const);
    
    
    // =============== QCD =======================
    /*
    printf("======= NLO QCD ===========\n");
    params[1] = m_pi;
    func.f = &nlo_qcd_mc;
    
    state = gsl_monte_vegas_alloc(4);
    gsl_monte_vegas_integrate (&func, xl_nlo2, xu_nlo2, 4, 1000, r, state, &res, &err);
    gsl_monte_vegas_integrate (&func, xl_nlo2, xu_nlo2, 4, N, r, state, &res, &err);

    
    gsl_monte_vegas_free(state);
    
    results[11]=res/results[1]*4*Pi/137.;
    results[12]=err/results[1]*4*Pi/137.;
    
    
    printf("qcd %f    %f\n",results[11],results[12]);
   
    */
    
    //################## NLO_1 ###################################
   

    printf("==========NLO 1================\n");    
    
    func.f = &nlo1_mc;   
    
    
 
  
    state = gsl_monte_vegas_alloc(4);
    gsl_monte_vegas_integrate (&func, xl, xu, 4, 1000, r, state, &res, &err);
    gsl_monte_vegas_integrate (&func, xl, xu, 4, N/2, r, state, &res, &err);
    gsl_monte_vegas_free(state);
  
 
    
    state = gsl_monte_vegas_alloc(4);
    gsl_monte_vegas_integrate (&func, xl_1, xu_1, 4, 1000, r, state, &res1, &err1);
    gsl_monte_vegas_integrate (&func, xl_1, xu_1, 4, N/2, r, state, &res1, &err1);
    gsl_monte_vegas_free(state);
    
 
   
    state = gsl_monte_vegas_alloc(4);
    gsl_monte_vegas_integrate (&func, xl_2, xu_2, 4, 1000, r, state, &res2, &err2);
    gsl_monte_vegas_integrate (&func, xl_2, xu_2, 4, N/2, r, state, &res2, &err2);  
    gsl_monte_vegas_free(state);
   
    
    
    //printf("nlo1/born s1,s2>s_cut \n %.6f    %.6f\n", res/results[1]*4*Pi/137., err/results[1]*4*Pi/137.);
    printf("nlo1/born 1<s1<s_cut \n %.6f    %.6f\n", res1/results[1]*4*Pi/137., err1/results[1]*4*Pi/137.);
    printf("nlo1/born 1<s2<s_cut \n %.6f    %.6f\n", res2/results[1]*4*Pi/137., err2/results[1]*4*Pi/137.);
          
       
    results[3]=(res+res1+res2)/results[1]*4*Pi/137.;
    results[4]=sqrt(err*err+err1*err1+err2*err2)/results[1]*4*Pi/137.;
      
    
    printf("final result nlo1:\n");
    printf("%.6f + pm %.6f \n",results[3],results[4]);


    
    //################### NLO 3 #####################################
    /*
    params[1] = Ecut;
    func.f = &nlo3_mc;
    
    state = gsl_monte_vegas_alloc(4);
    gsl_monte_vegas_integrate (&func, xl_nlo2, xu_nlo2, 4, 10000, r, state, &res, &err);
    gsl_monte_vegas_integrate (&func, xl_nlo2, xu_nlo2, 4, N, r, state, &res, &err);
    
    gsl_monte_vegas_free(state);
    
    results[9]=(res)/results[1]*4*Pi/137.;
    results[10]=err/results[1]*4*Pi/137.;
    
    printf("==========NLO 3================\n");
    printf("%.6f + pm %.6f \n",results[9],results[10]);
*/
    
    
    
    //################### NLO_2###################################
 
    func.f = &nlo2_mc;


    params[1] = me;       
    state = gsl_monte_vegas_alloc(4);
    
    
    printf("========= NLO 2 ===========\n");
    
    
    gsl_monte_vegas_integrate (&func, xl_nlo2, xu_nlo2, 4, 1000, r, state, &res, &err);
    gsl_monte_vegas_integrate (&func, xl_nlo2, xu_nlo2, 4, N, r, state, &res, &err);
    
    results[5]=res/results[1]*4*Pi/137.;
    results[6]=err/results[1]*4*Pi/137.;

    
    printf("nlo2/born, e\n %.6f    %.6f\n",results[5],results[6]);
    gsl_monte_vegas_free(state);
    
    

 
    params[1] = 1.;
       
    
    state = gsl_monte_vegas_alloc(4);
    gsl_monte_vegas_integrate (&func, xl_nlo2, xu_nlo2, 4, 1000, r, state, &res, &err);
    gsl_monte_vegas_integrate (&func, xl_nlo2, xu_nlo2, 4, N, r, state, &res, &err);

    
    gsl_monte_vegas_free(state);
    
    results[7]=res/results[1]*4*Pi/137.;
    results[8]=err/results[1]*4*Pi/137.;
    
    printf("nlo2/born, mu\n %.6f    %.6f\n",results[7],results[8]);



    params[1] = mtau;
       
    
    state = gsl_monte_vegas_alloc(4);
    gsl_monte_vegas_integrate (&func, xl_nlo2, xu_nlo2, 4, 1000, r, state, &res, &err);
   
    gsl_monte_vegas_integrate (&func, xl_nlo2, xu_nlo2, 4, N, r, state, &res, &err);

    
    gsl_monte_vegas_free(state);
    
    results[9]=res/results[1]*4*Pi/137.;
    results[10]=err/results[1]*4*Pi/137.;
    
    printf("nlo2/born, tau\n %.6f    %.6f\n",results[9],results[10]);

    
    
 
     
   
    printf("======= NLO2 quarks ===========\n");
    params[1] = mQ;
       
    
    state = gsl_monte_vegas_alloc(4);
    gsl_monte_vegas_integrate (&func, xl_nlo2, xu_nlo2, 4, 1000, r, state, &res, &err);
    gsl_monte_vegas_integrate (&func, xl_nlo2, xu_nlo2, 4, N, r, state, &res, &err);

    
    gsl_monte_vegas_free(state);
    
    results[11]=res/results[1]*4*Pi/137.*(pow(2./3.,4)+pow(1./3.,4)+pow(1./3.,4));
    results[12]=err/results[1]*4*Pi/137.*(pow(2./3.,4)+pow(1./3.,4)+pow(1./3.,4));
    
    
    printf("3q with m=300 MeV %.6f    %.6f\n",results[11],results[12]);
 

    //================Test==================
   
    
    /*
 
   
    params[1] = m_c;
       
    
    state = gsl_monte_vegas_alloc(4);
    gsl_monte_vegas_integrate (&func, xl_nlo2, xu_nlo2, 4, 1000, r, state, &res, &err);
    gsl_monte_vegas_integrate (&func, xl_nlo2, xu_nlo2, 4, N, r, state, &res, &err);

    
    gsl_monte_vegas_free(state);
    
    results[17]=res/results[1]*4*Pi/137.*pow(2./3.,4);
    results[18]=err/results[1]*4*Pi/137.*pow(2./3.,4);
    
    
    printf("charm %.6f    %.6f\n",results[17],results[18]);
   
    params[1] = m_b;
       
    
    state = gsl_monte_vegas_alloc(4);
    gsl_monte_vegas_integrate (&func, xl_nlo2, xu_nlo2, 4, 1000, r, state, &res, &err);
    gsl_monte_vegas_integrate (&func, xl_nlo2, xu_nlo2, 4, N, r, state, &res, &err);

    
    gsl_monte_vegas_free(state);
    
    results[19]=res/results[1]*4*Pi/137.*pow(1./3.,4);
    results[20]=err/results[1]*4*Pi/137.*pow(1./3.,4);
    
    
    printf("beauty %.6f    %.6f\n",results[19],results[20]);
   
    
    
    printf("all quarks %.6f \n",results[11]+results[13]+results[15]+results[17]+results[19]);
    
    */
    //=============================
    
    save(13,results);
    
     
    
    
    
    

    
    gsl_rng_free (r);
    
    return 0;
};


void display_results(char *title, double result, double error)
{
  printf ("%s ==================\n", title);
  printf ("result = % .6f\n", result);
  printf ("sigma  = % .6f\n", error);
};



int save(int len,double *Mres)
{
    FILE *file;
    int i;
    file=fopen("num.txt","a");
    for(i=0;i<len;i++)
    {
        fprintf(file,"%.6f  ",Mres[i]);
    };
    fprintf(file,"\n");
    
    fclose(file);
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
    double s;
    double result;
    
        
    s1=inv[0];
    s2=inv[1];
    t1=inv[2];
    lam=inv[3];
    
    s=pars[0];
       
           
    if( Gk(s,t1,s2,0.,0.,1.)<0 && Gk(s1,s2,s,0.,1.,1.)<0 ) 
    {
        t2 = t2_func(s,s1,s2,t1,lam);
        if( cuts(s,s1,s2,t1,t2) )
        {
            vars[0] = s;
            vars[1] = s1;
            vars[2] = s2;
            vars[3] = t1;
            vars[4] = t2;
            result = born(vars)*PhaseVol(s,s1,s2,t1)/( 2*s*Power(2*Pi,5) );        
        }
        else result=0;
    }
    else  result = 0.;
    
    return result;
};





double nlo1_mc(double *inv, size_t dim, void *params)
{
    (void)(dim);
    double *pars = (double *)params;
    double vars[7];
    double s1,s2,t1,lam; // variables
    double t2; // temporary variables
    double s, omega;
    double result;
    
    
    s1=inv[0];
    s2=inv[1];
    t1=inv[2];
    lam=inv[3];
    
    s=pars[0];
    omega=pars[1];
    
    
    
        
    if( Gk(s,t1,s2,0.,0.,1.)<0 && Gk(s1,s2,s,0.,1.,1.)<0 ) 
    {
        t2 = t2_func(s,s1,s2,t1,lam);
        if( cuts(s,s1,s2,t1,t2) )
        {
            vars[0] = s;
            vars[1] = s1;
            vars[2] = s2;
            vars[3] = t1;
            vars[4] = t2;
            vars[5] = omega;
            result = nlo1(vars)*PhaseVol(s,s1,s2,t1)/( 2*s*Power(2*Pi,5) );   
            if(isnanl(result)!=0) printf("NAN nlo1 func: s1=%e  s2=%e  t1=%e  t2=%e  nlo1=%e\n",s1,s2,t1,t2, nlo1(vars) );
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
    double s, m, t2; // temporary variables
    double result;
    
        
    s1=inv[0];
    s2=inv[1];
    t1=inv[2];
    lam=inv[3];
    
    s=pars[0];
    m=pars[1];
    
        
    if( Gk(s,t1,s2,0.,0.,1.)<0 && Gk(s1,s2,s,0.,1.,1.)<0 ) 
    {
        t2 = t2_func(s,s1,s2,t1,lam);
        if( cuts(s,s1,s2,t1,t2) )
        {
            vars[0] = s;
            vars[1] = s1;
            vars[2] = s2;
            vars[3] = t1;
            vars[4] = t2;
            vars[5] = m;
            result = nlo2(vars)*PhaseVol(s,s1,s2,t1)/( 2*s*Power(2*Pi,5) );   
            //if(isnanl(result)!=0) printf("NAN nlo2 func: s1=%e  s2=%e  t1=%e  t2=%e  nlo2=%e\n",s1,s2,t1,t2, nlo2(vars) );
        }
        else result=0.;
    }
    else  result = 0.;
        
    return result;
};

double nlo3_mc(double *inv, size_t dim, void *params)
{
    (void)(dim);
    double *pars = (double *)params;
    double vars[8];
    double s1,s2,t1,lam; // variables
    double t2, np1, np2; // temporary variables
    double s, Ec;
    double result;
    
    
    s1=inv[0];
    s2=inv[1];
    t1=inv[2];
    lam=inv[3];
    
    s=pars[0];
    Ec=pars[1];
    
    
    np1 = (s-s2+1)/sqrt(s)/2;
    np2 = (s-s1+1)/sqrt(s)/2;
        
    if( Gk(s,t1,s2,0.,0.,1.)<0 && Gk(s1,s2,s,0.,1.,1.)<0 ) 
    {
        t2 = t2_func(s,s1,s2,t1,lam);
        if( cuts(s,s1,s2,t1,t2) )
        {
            vars[0] = s;
            vars[1] = s1;
            vars[2] = s2;
            vars[3] = t1;
            vars[4] = t2;
            vars[5] = np1;
            vars[6] = np2;
            vars[7] = Ec;
            result = nlo3(vars)*PhaseVol(s,s1,s2,t1)/( 2*s*Power(2*Pi,5) );   
            if(isnanl(result)!=0) printf("NAN nlo3 func: s1=%e  s2=%e  t1=%e  t2=%e  nlo1=%e\n",s1,s2,t1,t2, nlo1(vars) );
        }
        else result=0.;
        
        
    }
    else
    {
        result = 0.;
    };
    
    return result;

};

/*
double nlo_qcd_mc(double *inv, size_t dim, void *params)
{
    (void)(dim);
    double *pars = (double *)params;
    double vars[6];
    double s1,s2,t1,lam; // variables
    double s, m, t2; // temporary variables
    double result;
    
        
    s1=inv[0];
    s2=inv[1];
    t1=inv[2];
    lam=inv[3];
    
    s=pars[0];
    m=pars[1];
    
        
    if( Gk(s,t1,s2,0.,0.,1.)<0 && Gk(s1,s2,s,0.,1.,1.)<0 ) 
    {
        t2 = t2_func(s,s1,s2,t1,lam);
        if( cuts(s,s1,s2,t1,t2) )
        {
            vars[0] = s;
            vars[1] = s1;
            vars[2] = s2;
            vars[3] = t1;
            vars[4] = t2;
            vars[5] = m;
            result = nlo_qcd(vars)*PhaseVol(s,s1,s2,t1)/( 2*s*Power(2*Pi,5) );   
            //if(isnanl(result)!=0) printf("NAN nlo2 func: s1=%e  s2=%e  t1=%e  t2=%e  nlo2=%e\n",s1,s2,t1,t2, nlo2(vars) );
        }
        else result=0.;
    }
    else  result = 0.;
        
    return result;
};
*/