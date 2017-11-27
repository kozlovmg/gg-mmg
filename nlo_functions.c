

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_pow_int.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#define Pi M_PI
#define Sqrt sqrt
#define Power gsl_pow_int
#define ln gsl_sf_log_abs
#define Log gsl_sf_log_abs
#define arctan atan

#define plus gsl_complex_add
#define div gsl_complex_div

double B1_inty(double , void *);
double B1_integrand(double , void *);

long double pw_n(long double x,unsigned int a)
{
    long double tmp=1.;
    do
    {
        if (a & 1) tmp=tmp*x;
        a>>=1;
        x=x*x;
    }
    while(a);
    return tmp;
};

/*some functions*/

double ln2b(double x)
{
    if(x<0) 
    {
        return Power( Log((1+sqrt(1-4/x))/(sqrt(1-4/x)-1)) ,2);
    };
    if(x>4)
    {
        return Power( Log((1+sqrt(1-4/x))/(sqrt(1-4/x)-1)) ,2)-Power(Pi,2);
    }
    else
    {
        return -Power(Pi-2*atan(sqrt(4/x-1)),2);
    };   
    if(x==0) return 0.;
    return 0.;
};




double Li2(double x)
{
    if (x<=1) return gsl_sf_dilog(x);
    else return Power(Pi,2)/6.-ln(x)*ln(1-x)-gsl_sf_dilog(1-x);
   };
   
double Li2im(double re,double im)
{
    gsl_sf_result result_re, result_im;
    double r,theta;
    double result;
    
    r=sqrt(re*re+im*im);
    theta=atan(im/re);
    
    gsl_sf_complex_dilog_e(r,theta,&result_re,&result_im);
    
    result=result_im.val;
    
    return result;
};

double Li2re(double re,double im)
{
    gsl_sf_result result_re, result_im;
    double r,theta;
    double result;
    
    r=sqrt(re*re+im*im);
    theta=atan(im/re);
    
    gsl_sf_complex_dilog_e(r,theta,&result_re,&result_im);
    
    result=result_re.val;
    
    return result;
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
    
    
double Rv2(double a,double b, double y)
{
    double re1,re2,im1,im2;
    if(b<0)
    {
        re1=(a*(a-y)-b)/(Power(a-y,2)-b);
        im1=-y*sqrt(-b)/(Power(a-y,2)-b);
        re2=((a-1)*(a-y)-b)/(Power(a-y,2)-b);
        im2=(1-y)*sqrt(-b)/(Power(a-y,2)-b);
        return 2*(Li2im(re1,im1)-Li2im(re2,im2));
    }
    else
    {
        return R(a+sqrt(b),y)-R(a-sqrt(b),y);
    };
    return 0;
};

double Rv3(double x,double a,double b)
{
    double re1,re2,im1,im2;
    if(b<0)
    {
        re1=(x*(x-a))/(Power(x-a,2)-b);
        im1=x*sqrt(-b)/(Power(x-a,2)-b);
        re2=(x-1)*(x-a)/(Power(x-a,2)-b);
        im2=(x-1)*sqrt(-b)/(Power(x-a,2)-b);
              
        return 2*(Li2re(re1,im1)-Li2re(re2,im2));
    }
    else
    {
        return R(x,a+sqrt(b))+R(x,a-sqrt(b));
    };
    return 0;
    
};

double reRc(double a,double b,double y)
{
    double re1,re2,im1,im2;
    
    re1=(a*a+b*b-y*a)/(Power(a-y,2)+b*b);
    im1=-b*y/(Power(a-y,2)+b*b);
    re2=((a-1)*(a-y)+b*b)/(Power(a-y,2)+b*b);
    im2=b*(1-y)/(Power(a-y,2)+b*b);
    
    return Li2re(re1,im1)-Li2re(re2,im2);    
};

double imRc(double a,double b,double y)
{
    double re1,re2,im1,im2;
    
    re1=(a*a+b*b-y*a)/(Power(a-y,2)+b*b);
    im1=-b*y/(Power(a-y,2)+b*b);
    re2=((a-1)*(a-y)+b*b)/(Power(a-y,2)+b*b);
    im2=b*(1-y)/(Power(a-y,2)+b*b);
    
    return Li2im(re1,im1)-Li2im(re2,im2);    
};



double RD(double a, double b, double y, double D)
{
    if(D>=0 && b<=0) return 0.;
    if(D>0 && b>0) return Rv2(a,b,y)/sqrt(D);
    if(D<0 && b<0) return Rv2(a,b,y)/sqrt(-D);
    return 0;
};
    

    
    
double S3(double x, double x1, double x2)
{
    return R(x,x1)+R(x,x2);
    };
    

double reS3c(double a,double b, double x1, double x2)
{
    return reRc(a,b,x1)+reRc(a,b,x2);
};

double imS3c(double a,double b, double x1, double x2)
{
    return imRc(a,b,x1)+imRc(a,b,x2);
};


    
double S3v2(double ya,double yb, double xa, double xb)
{
    double re1,re2,im1,im2,re3,re4,im3,im4;
    if(xb>=0)
    {
        return Rv2(ya,yb,xa+sqrt(xb))+Rv2(ya,yb,xa-sqrt(xb));
    };
    if(yb>=0)
    {
        return Rv3(ya+sqrt(yb),xa,xb)-Rv3(ya-sqrt(yb),xa,xb);
    };
    if(xb<0 && yb<0)
    {
        re1=(ya*(ya-xa)-yb-sqrt(xb*yb))/(Power(ya-xa,2)+Power(sqrt(-xb)-sqrt(-yb),2));
        im1=(ya*sqrt(-xb)-xa*sqrt(-yb))/(Power(ya-xa,2)+Power(sqrt(-xb)-sqrt(-yb),2));
        re2=((ya-1)*(ya-xa)-yb-sqrt(xb*yb))/(Power(ya-xa,2)+Power(sqrt(-xb)-sqrt(-yb),2));
        im2=((1-xa)*sqrt(-yb)+(ya-1)*sqrt(-xb))/(Power(ya-xa,2)+Power(sqrt(-xb)-sqrt(-yb),2));
        
        re3=(ya*(ya-xa)-yb+sqrt(xb*yb))/(Power(ya-xa,2)+Power(sqrt(-xb)+sqrt(-yb),2));
        im3=(-ya*sqrt(-xb)-xa*sqrt(-yb))/(Power(ya-xa,2)+Power(sqrt(-xb)+sqrt(-yb),2));
        re4=((ya-1)*(ya-xa)-yb+sqrt(xb*yb))/(Power(ya-xa,2)+Power(sqrt(-xb)+sqrt(-yb),2));
        im4=((1-xa)*sqrt(-yb)+(1-ya)*sqrt(-xb))/(Power(ya-xa,2)+Power(sqrt(-xb)+sqrt(-yb),2));
        
        return 2*(Li2im(re1,im1)-Li2im(re2,im2)+Li2im(re3,im3)-Li2im(re4,im4));
    };
    return 0;
};
    

double S3D(double ya,double yb, double xa, double xb,double D)
{
    if(yb>=0 && xb>=0 && D>=0)
    {
        return S3v2(ya,yb,xa,xb)/sqrt(D);
    };
    if(yb>0 && xb>0 && D<0) return 0.;
    if(yb<0 && xb>0 && D>0) return 0.;
    if(yb>0 && xb<0 && D>0)
    { 
        return S3v2(ya,yb,xa,xb)/sqrt(D);
    };
    if(yb<0 && xb<0 && D>0) return 0.;
    if(yb<0 && xb>0 && D<0) 
    {
        return S3v2(ya,yb,xa,xb)/sqrt(-D);
    };
    if(yb>0 && xb<0 && D<0) return 0.;
    if(yb<0 && xb<0 && D<0) 
    {
        return S3v2(ya,yb,xa,xb)/sqrt(-D);
    };
    return 0;
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
    
long double R1(long double x)
{
    return (1./x-1.)*logl(fabsl(1.-x));
    };

long double R2(long double x)
{
    if(x<0. || x>4.)
    { 
        return -sqrtl(1. - 4./x)*logl( fabsl( (1. + sqrtl(1. - 4./x))/(1. - sqrtl(1. - 4./x)) ) );
    }
    else
    {
        return sqrtl(4./x-1.)*(2.*atanl(sqrtl(4./x-1.))-Pi);
    };
};

double T1(double x)
{
    return  (Power(Pi,2)/6.-Li2(1+x))/x;
    };

double T3(double t1,double t2)
{
    return -(Li2(t1)-Li2(t2))/(t1-t2);
    };
    
double T4(double x)
{
    return (Li2(1-x)+ln(x)*ln(x-1))/(x-1);
     };
    
double T5(double x)
{
    return ln2b(x)/2/x;
    };

double T6(double s,double t)
{
    return (ln2b(s)-ln2b(t))/(2.*(s - t));
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
/*double B1(double s,double t,double m)
{
    double D1,D2,y1p,y1m,y2p,y2m,b,x1p,x1m,x2p,x2m,x3p,x3m,x4p,x4m,tmp;
    D1=4*t+Power(1-m+t,2);
    D2=(-4*m+s*(4+4*m-2*t)+Power(s,2)*(-4+t)+9*t)*t;
    b=-(1+t-m+sqrt(D1))/2;
    y1p=-((-2.*t+s*t+Power(t,2)+t*sqrt(D1)-sqrt(D2)-t*m)/(-1.+s-3.*t+s*t+(1.-s)*sqrt(D1)+m-s*m));
    y1m=-((-2.*t+s*t+Power(t,2)+t*sqrt(D1)+sqrt(D2)-t*m)/(-1.+s-3.*t+s*t+(1.-s)*sqrt(D1)+m-s*m));
    y2p=(t+s*t+sqrt(D2))/(2.*(Power(s,2)-2.*t+s*(-1.+t-m)+m));
    y2m=(t+s*t-sqrt(D2))/(2.*(Power(s,2)-2.*t+s*(-1.+t-m)+m));
    x1p=-1.+sqrt(2.);
    x1m=-1.-sqrt(2.);
    x2p=(1.+m+sqrt(8.+Power(m-1.,2)))/2./(m-2.);
    x2m=(1.+m-sqrt(8.+Power(m-1.,2)))/2./(m-2.);
    x3p=(1.+sqrt(1.-4./t))/2.;
    x3m=(1.-sqrt(1.-4./t))/2.;
    x4p=(1.+s+sqrt(8.+Power(s-1.,2)))/2./(s-2.);
    x4m=(1.+s-sqrt(8.+Power(s-1.,2)))/2./(s-2.);
    tmp=( -R(y1m/(1.-b),t/(-1.+s+t))+R(-(y1m/b),1.) + R(b + y1m,0.) + R(y1p/(1.-b),t/(-1.+s+t)) - R(-(y1p/b),1.) - R(b+y1p,0.) - R(y2m,0.) + R(y2m,t/(-1.+s+t)) + R(y2p,0.) - R(y2p,t/(-1.+s+t)) + S3(y1m/(1.-b),x2p,x2m) - S3(-(y1m/b),x3p,x3m) - S3(b+y1m,x1p,x1m) - S3(y1p/(1.-b),x2p,x2m) + S3(-(y1p/b),x3p,x3m) + S3(b +y1p,x1p,x1m) - S3(y2m,x2p,x2m) + S3(y2m,x4p,x4m) + S3(y2p,x2p,x2m) - S3(y2p,x4p,x4m))/sqrt(D2);
    return tmp;
    };*/
    
/* B1 ????? */
double B1(double s,double t,double m)
{
    double D1,D2,x1p,x1m,x2p,x2m,x3p,x3m, tmp;
    double a2,b2,a4,b4,y2a,y2b;
    D1=4*t+Power(1-m+t,2);
    D2=(-4*m+s*(4+4*m-2*t)+Power(s,2)*(-4+t)+9*t)*t;
    
    x2p=(1.+m+sqrt(8.+Power(m-1.,2)))/2./(m-2.);
    x2m=(1.+m-sqrt(8.+Power(m-1.,2)))/2./(m-2.);
            
    y2a=(t+s*t)/(2.*(Power(s,2)-2.*t+s*(-1.+t-m)+m));
    y2b=D2/Power(2.*(Power(s,2)-2.*t+s*(-1.+t-m)+m),2);
    
    x1p=-1.+sqrt(2.);
    x1m=-1.-sqrt(2.);
    x3p=(1.+sqrt(1.-4./t))/2.;
    x3m=(1.-sqrt(1.-4./t))/2.;
    
        
    a2=(1+m)/(2*(m-2));
    b2=(8+Power(m-1,2))/Power(2*(m-2),2);
    a4=(1+s)/(2*(s-2));
    b4=(8+Power(s-1,2))/Power(2*(s-2),2);
    
    if(D1>=0 && D2>=0)
    {
        double b,y1p,y1m;
        b=-(1+t-m+sqrt(D1))/2;
        y1p=-((-2.*t+s*t+Power(t,2)+t*sqrt(D1)-sqrt(D2)-t*m)/(-1.+s-3.*t+s*t+(1.-s)*sqrt(D1)+m-s*m));
        y1m=-((-2.*t+s*t+Power(t,2)+t*sqrt(D1)+sqrt(D2)-t*m)/(-1.+s-3.*t+s*t+(1.-s)*sqrt(D1)+m-s*m));
        tmp=(-R(y1m/(1.-b),t/(-1.+s+t))+R(-(y1m/b),1.) + R(b + y1m,0.) + R(y1p/(1.-b),t/(-1.+s+t)) - R(-(y1p/b),1.) - R(b+y1p,0.)   + S3(y1m/(1.-b),x2p,x2m) - S3(-(y1m/b),x3p,x3m) - S3(b+y1m,x1p,x1m) - S3(y1p/(1.-b),x2p,x2m) + S3(-(y1p/b),x3p,x3m) + S3(b +y1p,x1p,x1m))/sqrt(D2) - RD(y2a,y2b,t/(-1+s+t),D2) + RD(y2a,y2b,0,D2) + S3D(y2a,y2b,a2,b2,D2) - S3D(y2a,y2b,a4,b4,D2);
        return tmp;    
    };
    if(D1>=0 && D2<0)
    {
        double r1,i1,r2,i2,r3,i3,b,ry1,iy1;
        b=-(1+t-m+sqrt(D1))/2;
        ry1=-((-2.*t+s*t+Power(t,2)+t*sqrt(D1)-t*m)/(-1.+s-3.*t+s*t+(1.-s)*sqrt(D1)+m-s*m));
        iy1=sqrt(-D2)/(-1.+s-3.*t+s*t+(1.-s)*sqrt(D1)+m-s*m);
        r1=(ry1)/(1-b);
        i1=(iy1)/(1-b);
        r2=-(ry1)/b;
        i2=-(iy1)/b;
        r3=(ry1)+b;
        i3=(iy1)+b;           
        tmp=( -imRc(r1,i1,t/(-1.+s+t)) + imRc(r2,i2,1) + imRc(r3,i3,0) +imRc(r1,-i1,t/(-1.+s+t)) - imRc(r2,-i2,1) - imRc(r3,-i3,0) + imS3c(r1,i1,x2p,x2m) - imS3c(r2,i2,x3p,x3m) - imS3c(r3,i3,x1p,x1m) - imS3c(r1,-i1,x2p,x2m) + imS3c(r2,-i2,x3p,x3m) + imS3c(r3,-i3,x1p,x1m) )/sqrt(-D2)- RD(y2a,y2b,t/(-1+s+t),D2) + RD(y2a,y2b,0,D2) + S3D(y2a,y2b,a2,b2,D2) - S3D(y2a,y2b,a4,b4,D2);
        return tmp;
    };
    if(D1<0 && D2>=0)
    {
        double r1m,i1m,r1p,i1p,r2m,i2m,r2p,i2p,r3m,i3m,r3p,i3p;
        gsl_complex bb, bb1, bb2, np, nm,den, arg1m, arg2m, arg3m, arg1p, arg2p, arg3p, y1p, y1m;
        
        
        bb=gsl_complex_rect(-(1+t-m)/2,-(sqrt(-D1))/2);
        bb1=gsl_complex_rect(1+(1+t-m)/2,(sqrt(-D1))/2);
        bb2=gsl_complex_rect((1+t-m)/2,(sqrt(-D1))/2);
        
        den=gsl_complex_rect(-1.+s-3.*t+s*t+m-s*m,(1.-s)*sqrt(-D1));
        np=gsl_complex_rect(-(-2.*t+s*t+Power(t,2)-sqrt(D2)-t*m),-t*sqrt(-D1));
        nm=gsl_complex_rect(-(-2.*t+s*t+Power(t,2)+sqrt(D2)-t*m),-t*sqrt(-D1));
        
        
                
        y1p=div(np,den);
        y1m=div(nm,den);

        arg1m=div(y1m,bb1);
        arg1p=div(y1p,bb1);
        
        arg2m=div(y1m,bb2);
        arg2p=div(y1p,bb2);
        
        arg3m=plus(y1m,bb);
        arg3p=plus(y1p,bb);
                
        
        r1m= arg1m.dat[0];
        i1m= arg1m.dat[1];
        r1p= arg1p.dat[0];
        i1p= arg1p.dat[1];
        
        r2m= arg2m.dat[0];
        i2m= arg2m.dat[1];
        r2p= arg2p.dat[0];
        i2p= arg2p.dat[1];
        
        r3m= arg3m.dat[0];
        i3m= arg3m.dat[1];
        r3p= arg3p.dat[0];
        i3p= arg3p.dat[1];
              

        tmp=( -reRc(r1m,i1m,t/(-1.+s+t)) + reRc(r2m,i2m,1) + reRc(r3m,i3m,0) + reRc(r1p,i1p,t/(-1.+s+t)) - reRc(r2p,i2p,1) - reRc(r3p,i3p,0) + reS3c(r1m,i1m,x2p,x2m) - reS3c(r2m,i2m,x3p,x3m) - reS3c(r3m,i3m,x1p,x1m) - reS3c(r1p,i1p,x2p,x2m) + reS3c(r2p,i2p,x3p,x3m) + reS3c(r3p,i3p,x1p,x1m) )/sqrt(D2)- RD(y2a,y2b,t/(-1+s+t),D2) + RD(y2a,y2b,0,D2) + S3D(y2a,y2b,a2,b2,D2) - S3D(y2a,y2b,a4,b4,D2);
        return tmp;
    };
    if(D1<0 && D2<0)
    {
        double r1m,i1m,r1p,i1p,r2m,i2m,r2p,i2p,r3m,i3m,r3p,i3p;
        gsl_complex bb, bb1, bb2, np, nm,den, arg1m, arg2m, arg3m, arg1p, arg2p, arg3p, y1p, y1m;
        
        
        bb=gsl_complex_rect(-(1+t-m)/2,-(sqrt(-D1))/2);
        bb1=gsl_complex_rect(1+(1+t-m)/2,(sqrt(-D1))/2);
        bb2=gsl_complex_rect((1+t-m)/2,(sqrt(-D1))/2);
        
        den=gsl_complex_rect(-1.+s-3.*t+s*t+m-s*m,(1.-s)*sqrt(-D1));
        np=gsl_complex_rect(-(-2.*t+s*t+Power(t,2)-t*m),-t*sqrt(-D1)+sqrt(-D2));
        nm=gsl_complex_rect(-(-2.*t+s*t+Power(t,2)-t*m),-t*sqrt(-D1)-sqrt(-D2));
                
        y1p=div(np,den);
        y1m=div(nm,den);

        arg1m=div(y1m,bb1);
        arg1p=div(y1p,bb1);
        
        arg2m=div(y1m,bb2);
        arg2p=div(y1p,bb2);
        
        arg3m=plus(y1m,bb);
        arg3p=plus(y1p,bb);
                
        
        r1m= arg1m.dat[0];
        i1m= arg1m.dat[1];
        r1p= arg1p.dat[0];
        i1p= arg1p.dat[1];
        
        r2m= arg2m.dat[0];
        i2m= arg2m.dat[1];
        r2p= arg2p.dat[0];
        i2p= arg2p.dat[1];
        
        r3m= arg3m.dat[0];
        i3m= arg3m.dat[1];
        r3p= arg3p.dat[0];
        i3p= arg3p.dat[1];
              

        tmp=( -imRc(r1m,i1m,t/(-1.+s+t)) + imRc(r2m,i2m,1) + imRc(r3m,i3m,0) + imRc(r1p,i1p,t/(-1.+s+t)) - imRc(r2p,i2p,1) - imRc(r3p,i3p,0) + imS3c(r1m,i1m,x2p,x2m) - imS3c(r2m,i2m,x3p,x3m) - imS3c(r3m,i3m,x1p,x1m) - imS3c(r1p,i1p,x2p,x2m) + imS3c(r2p,i2p,x3p,x3m) + imS3c(r3p,i3p,x1p,x1m) )/sqrt(-D2)- RD(y2a,y2b,t/(-1+s+t),D2) + RD(y2a,y2b,0,D2) + S3D(y2a,y2b,a2,b2,D2) - S3D(y2a,y2b,a4,b4,D2);
        return tmp;
    };
    return 0;
};
    


double B1N(double s, double t, double m)
{
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
    
    double result, error;
    double pars[3];
    
    pars[0]=s;
    pars[1]=t;
    pars[2]=m;
    
    
    gsl_function func;
    func.function=&B1_inty;
    func.params=pars;
    
    gsl_integration_qags(&func,-1,1,0,1e-7,1000,w,&result,&error);
    gsl_integration_workspace_free(w);
    
    return result;
};

double B1_inty(double x, void * params)
{
    double s=*(double *) params;
    double t=*((double *) params+1);
    double m=*((double *) params+2);
    
    
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
    
    double result, error;
    double pars[4];
    
    pars[0]=x;
    pars[1]=s;
    pars[2]=t;
    pars[3]=m;
    
    gsl_function func;
    func.function=&B1_integrand;
    func.params=pars;
    
    gsl_integration_qags(&func,-1,1,0,1e-7,1000,w,&result,&error);
    gsl_integration_workspace_free(w);
    
    return result;
};

double B1_integrand(double y, void * params)
{
    double x=*(double *) params;
    double s=*((double *) params+1);
    double t=*((double *) params+2);
    double m=*((double *) params+3);
    return (x-y)/(-1+t*x-t*Power(x,2)+Power(y,2)+ x*y*(1+t-m)+y*(1-t+m))/(-1+(2-s)*Power(y,2)+x*y*(s-m)+y*(1+m));
};



double B2(double s,double t,double m)
{
    return (Power(Pi,2)/6.-2.*(ff(X(s),1./X(m)) + ff(X(s),X(m))) + Power(Pi,2)*ht(-X(m)) - Li2(Power(X(s),2)) - Power(ln(X(m)),2) - 2.*(-(Power(Pi,2)*ht(-1.+t)*ht(-X(s))) + ln(1.-t)*ln(X(s))) - 2.*ln(X(s))*ln(1.-Power(X(s),2)))/(s*(-1.+t)*bb(s));    
    };

/*double B3(double s,double t,double m)
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
    };*/

double B3(double s,double t,double m)
{
    
    double a1,a2,a3,b1,b2,b3;
    double y1a,y1b,y2a,y2b;
    double D, tmp;
    D=1.+4.*(m-s-t)/(s*t);
    
    y1a=-1/2;
    y1b=D/4;
    
    y2a=t/(s+t-m)/2.;
    y2b=Power(t/(s+t-m)/2.,2)*D;
    
    a1=1/2;
    b1=(1.-4./m)/4;
    a2=1/2;
    b2=(1.-4./t)/4;
    a3=1/2;
    b3=(1.-4./s)/4;
    
    
    tmp=( RD(-y1a,y1b,1,D)-RD(1+y1a,y1b,0,D)+RD(y2a,y2b,0,D)-RD(y2a,y2b,-t/(m-s-t),D)+S3D(1+y1a,y1b,a1,b1,D)-S3D(-y1a,y1b,a2,b2,D)-S3D(y2a,y2b,a3,b3,D) )/(s*t);
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



