#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_pow_int.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_integration.h>

#define Pi M_PI
#define Sqrt sqrt
#define Power gsl_pow_int
#define ln gsl_sf_log_abs
#define Log gsl_sf_log_abs

/*some functions*/

double Li2(double x)
{
    if (x<=1) return gsl_sf_dilog(x);
    else return Power(Pi,2)/6.-ln(x)*ln(1-x)-gsl_sf_dilog(1-x);
   };

double ht(double x)
{
    if (x>0.) return 1.;
    else return 0.;
    };

double R(double x, double y)
{
    return Li2(x/(x-y))-Li2((x-1)/(x-y));
    };
    
double S3(double x, double x1, double x2)
{
    return R(x,x1)+R(x,x2);
    };

double ff(double x, double y)
{
    double a;
    a=Power(Pi,2)*(ht(-x) + ht(-y))*ht(-1 + x*y) + ht(-(x*y))*(Power(Pi,2)/6.-Li2(x*y))-(ln(x) + ln(y))*ln(1-x*y)+
    ht(x*y)*(Li2(1-x*y)+ln(x*y)*ln(1- x*y));
    return a;
    };

double bb(double x)
{
    return sqrt(1-4./x);
    };

double X(double a)
{
    return (bb(a)-1)/(bb(a)+1);
    };

/*Simple functions for master integrals*/    
    
double R1(double x)
{
    return (1/x-1)*ln(1-x);
    };

double R2(double x)
{
    return -(sqrt(1 - 4/x)*ln((1 + sqrt(1 - 4/x))/(1 - sqrt(1 - 4/x))));
    };

double T1(double x)
{
    return  (Power(M_PI,2)/6.-Li2(1+x))/x;
    };

double T3(double t1, double t2)
{
    return -(Li2(t1)-Li2(t2))/(t1-t2);
    };
    
double T4(double x)
{
    return (Li2(1-x)+ln(x)*ln(x-1))/(x-1);
    };
    
double T5(double x)
{
    return (Power(ln((1+sqrt(1-4/x))/(1-sqrt(1-4/x))),2)-Power(M_PI,2)*ht(x))/2/x;
    };

double T6(double s, double t)
{
    return (-(Power(Pi,2)*(ht(s)-ht(t)))+Power(ln((1+sqrt(1-4/s))/(1-sqrt(1-4/s))),2)-Power(ln((1+Sqrt(1-4/t))/(1-Sqrt(1-4/t))),2))/(2.*(s - t));
    };

double F(double x)
{
    return (x*ln((1 + sqrt(1-Power(x,-2)))/(1-sqrt(1-Power(x,-2)))))/sqrt(-1+Power(x,2));
    };



/*Complicated functions*/

double T2(double t1,double t2)
{
    double x1,x2,x3,x31,x32,D3;
    D3=sqrt(Power((t1+t2-1),2)-4.*t1*t2);
    x1=(1.-Power(t1,2)-t2+t1*t2+(1.+t1)*D3)/2./t1/D3;
    x2=-(-2.+2.*t1+3.*t2+t1*t2-Power(t2,2)+(t2-2.)*D3)/(D3*(1.+t1-t2+D3));
    x3=-(-2.+2.*t1+3.*t2+t1*t2-Power(t2,2)+(t2-2.)*D3)/(D3*(1.-t1-t2+D3));
    x31=(-sqrt((-4.+t2)*t2)+t2)/(2.*t2);
    x32=(sqrt((-4.+t2)*t2)+t2)/(2.*t2);
    return (S3(x1,1.,1./t1)-S3(x2,1.,1.)+S3(x3,x31,x32))/D3;        
    };

/* B1 ????? */
double B1(double s,double t,double m)
{
    double D1,D2,y1p,y1m,y2p,y2m,b,x1p,x1m,x2p,x2m,x3p,x3m,x4p,x4m,tmp;
    D1=4*t+Power(1-m+t,2);
    D2=-4*m+s*(4+4*m-2*t)+Power(s,2)*(-4+t)+9*t;
    b=-(1+t-m+sqrt(D1))/2;
    y1p=-((-2.*t+s*t+Power(t,2)+t*sqrt(D1)-sqrt(t*D2)-t*m)/(-1.+s-3.*t+s*t+(1.-s)*sqrt(D1)+m-s*m));
    y1m=-((-2.*t+s*t+Power(t,2)+t*sqrt(D1)+sqrt(t*D2)-t*m)/(-1.+s-3.*t+s*t+(1.-s)*sqrt(D1)+m-s*m));
    y2p=(t+s*t+sqrt(t*D2))/(2.*(Power(s,2)-2.*t+s*(-1.+t-m)+m));
    y2m=(t+s*t-sqrt(t*D2))/(2.*(Power(s,2)-2.*t+s*(-1.+t-m)+m));
    x1p=-1.+sqrt(2.);
    x1m=-1.-sqrt(2.);
    x2p=(1.+m+sqrt(8.+Power(m-1.,2)))/2./(m-2.);
    x2m=(1.+m-sqrt(8.+Power(m-1.,2)))/2./(m-2.);
    x3p=(1.+sqrt(1.-4./t))/2.;
    x3m=(1.-sqrt(1.-4./t))/2.;
    x4p=(1.+s+sqrt(8.+Power(s-1.,2)))/2./(s-2.);
    x4m=(1.+s-sqrt(8.+Power(s-1.,2)))/2./(s-2.);
    tmp=( -R(y1m/(1.-b),t/(-1.+s+t))+R(-(y1m/b),1.) + R(b + y1m,0.) + R(y1p/(1.-b),t/(-1.+s+t)) - R(-(y1p/b),1.) - R(b+y1p,0.) - R(y2m,0.) + R(y2m,t/(-1.+s+t)) + R(y2p,0.) - R(y2p,t/(-1.+s+t)) + S3(y1m/(1.-b),x2p,x2m) - S3(-(y1m/b),x3p,x3m) - S3(b+y1m,x1p,x1m) - S3(y1p/(1.-b),x2p,x2m) + S3(-(y1p/b),x3p,x3m) + S3(b +y1p,x1p,x1m) - S3(y2m,x2p,x2m) + S3(y2m,x4p,x4m) + S3(y2p,x2p,x2m) - S3(y2p,x4p,x4m))/sqrt(D2*t);
    return tmp;
    };

double B2(double s,double t,double m)
{
    return (Power(Pi,2)/6.-2.*(ff(X(s),1./X(m)) + ff(X(s),X(m))) + Power(Pi,2)*ht(-X(m)) - Li2(Power(X(s),2)) - Power(ln(X(m)),2) - 2.*(-(Power(Pi,2)*ht(-1.+t)*ht(-X(s))) + ln(1.-t)*ln(X(s))) - 2.*ln(X(s))*ln(1.-Power(X(s),2)))/(s*(-1.+t)*bb(s));    
    };

double B3(double s,double t,double m)
{
    double y1p,y1m,x1p,x1m,x2p,x2m,y2p,y2m,x3p,x3m,tmp;
    y1p=(-1.+sqrt(1.+4.*(m-s-t)/(s*t)))/2.;
    y1m=(-1.-sqrt(1.+4.*(m-s-t)/(s*t)))/2.;
    y2p=t/(s+t-m)/2.*(1+sqrt(1.+4.*(m-s-t)/(s*t)));
    y2m=t/(s+t-m)/2.*(1-sqrt(1.+4.*(m-s-t)/(s*t)));
    x1p=(1.+sqrt(1.-4./m))/2.;
    x1m=(1.-sqrt(1.-4./m))/2.;
    x2p=(1.+sqrt(1.-4./t))/2.;
    x2m=(1.-sqrt(1.-4./t))/2.;
    x3p=(1.+sqrt(1.-4./s))/2.;
    x3m=(1.-sqrt(1.-4./s))/2.;
    tmp=(R(-y1m,1.) + R(1.+y1m,0.) - R(-y1p,1.) - R(1.+y1p,0.) - R(y2m,0.) + R(y2m,-(t/(m-s-t))) + R(y2p,0.) - R(y2p,-(t/(m-s-t))) - S3(-y1m,x2p,x2m) - S3(1.+y1m,x1p,x1m) + S3(-y1p,x2p,x2m) + S3(1.+y1p,x1p,x1m) + S3(y2m,x3p,x3m) - S3(y2p,x3p,x3m))/(s*sqrt(1.+(4.*(m-s-t))/(s*t))*t);   
    return tmp;
    };


double integrand(double x, void * params)
{
    double a=*(double *) params;
    double b=*((double *) params+1);
    double c=*((double *) params+2);
    double f=ln((b-c*x+sqrt(1-Power(x,2)+Power(b*x-c,2)))/2./sqrt(1-Power(x,2)))/(1-a*x);
    return f;
};
    
    
    
double G(double a, double b, double c)
{
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
    
    double result, error;
    double pars[3];
    
    pars[0]=a;
    pars[1]=b;
    pars[2]=c;
    
    gsl_function func;
    func.function=&integrand;
    func.params=pars;
    
    gsl_integration_qags(&func,-1,1,0,1e-7,1000,w,&result,&error);
    gsl_integration_workspace_free(w);
    
    return result;
};



