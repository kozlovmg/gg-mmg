#include <stdio.h>
#include <math.h>
#include "lbl.h"



double ir(double,double,double,double,double,double,double,double);
double bubl(double *);
double triangle(double *);
double box(double *);
double pent(double,double,double,double,double);







double nlo1(double *vars){
    double s=vars[0], s1=vars[1], s2=vars[2], t1=vars[3], t2=vars[4], omega=vars[5];
    double a=0., pp=0., bb=0., tt=0., bu=0., ii=0.;
    double np1 = (s-s2+1)/sqrt(s)/2;
    double np2 = (s-s1+1)/sqrt(s)/2;
    
    bb=box(vars);
    tt = triangle(vars);
    bu = bubl(vars);
    pp = pent(s,s1,s2,t1,t2);
    ii=ir(s,s1,s2,t1,t2,omega,np1,np2);
    
    
    
  
    
    
    //a=(pp+bb+tt+bu)/4.;// "/4" because nlo2=2*Re(Mnlo*(Mborn)^*)/4
    a=(tt+bu+pp+bb+ii)/2.;
    //printf("nlo1=%e\n",a);
    if( isnanl(a)!=0 )
    { 
        printf("nlo1 function is NAN:\n");
        printf("pent=%e box=%e  tri=%e  buble=%e ir=%e \n",pp,bb,tt,bu,ii);
    };
    
 
    /*
    if(a>1000000.)
    { 
        
        printf("nlo1 function:t1-t2=%e  s1-s2+t1-t2=%e\n",t1-t2,s1-s2+t1-t2);
        printf("s=%0.6f s1=%0.6f  s2=%0.6f  t1=%0.6f  t2=%0.6f \n ",s,s1,s2,t1,t2);
        printf("box=%e  tri=%e  buble=%e \n",bb,tt,bu);
        printf("u1=%e  u2=%e\n",-(s+t2-s1-1),-(s+t1-s2-1));
    };
 */
  return a;
  
};






double nlo2(double *vars){
    double s=vars[0], s1=vars[1], s2=vars[2], t1=vars[3], t2=vars[4], m=vars[5];
    long double r1,r2,r3,r4,r5,r6,a;
    r1=lbl1_q(s,s1,s2,t1,t2,m);
    //if(r1>10000000) printf("r1\n");
    r2=lbl2_q(s,s1,s2,t1,t2,m);
    //if(r2>10000000) printf("r2\n");
    r3=lbl3_q(s,s1,s2,t1,t2,m);
    //if(r3>10000000) printf("r3\n");
    r4=lbl4_q(s,s1,s2,t1,t2,m);
    //if(r4>10000000) printf("r4\n");
    r5=lbl5_q(s,s1,s2,t1,t2,m);
    //if(r5>10000000) printf("r5\n");
    r6=lbl6_q(s,s1,s2,t1,t2,m);
    //if(r6>10000000) printf("r6\n");
    a=(r1+r2+r3+r4+r5+r6)/2.;
    if( isnanl(a)!=0  )
    { 
        printf("nlo2 function is NAN: \n");
        printf("s=%0.6f s1=%0.6f  s2=%0.6f  t1=%0.6f  t2=%0.6f  m=%0.6f\n ",s,s1,s2,t1,t2,m);
        printf("r1=%Le  r2=%Le  r3=%Le  r4=%Le  r5=%Le  r6=%Le\n",r1,r2,r3,r4,r5,r6);
        //printf("u1=%e  u2=%e\n",-(s+t2-s1-1),-(s+t1-s2-1));
    };
    /*
    if(a>1000000.)
    { 
       
        printf("nlo2 function:\n");
        printf("s=%0.6f s1=%0.6f  s2=%0.6f  t1=%0.6f  t2=%0.6f \n ",s,s1,s2,t1,t2);
        printf("r1=%Le  r2=%Le  r3=%Le  r4=%Le  r5=%Le r6=%Le\n",r1,r2,r3,r4,r5,r6);
        //printf("cos(gamma)=%0.6f\n",(s1-s2+2*(t1-t2))/(s1+s2-2));
        //printf("u1=%e  u2=%e\n",-(s+t2-s1-1),-(s+t1-s2-1));
    };
    */
    return a;
};


double nlo3(double *vars){
    double s=vars[0], s1=vars[1], s2=vars[2], t1=vars[3], t2=vars[4], np1=vars[5], np2=vars[6], EE=vars[7];
    double a, ii;
    

    
    ii=ir(s,s1,s2,t1,t2,EE,np1,np2);
    
  
    
    
    
    a=(ii)/2.;
    //printf("nlo1=%e\n",a);
    if( isnanl(a)!=0 )
    { 
        printf("nlo3 function is NAN:\n");
        //printf("pent=%e  box=%e  tri=%e  buble=%e  ir=%e\n",pp,bb,tt,bu,ii);
        //printf("u1=%e  u2=%e\n",-(s+t2-s1-1),-(s+t1-s2-1));
    };
    
  
    /*
    if(a>10000000.)
    { 
        
        printf("nlo1 function:t1-t2=%e  s1-s2+t1-t2=%e\n",t1-t2,s1-s2+t1-t2);
        printf("s=%0.6f s1=%0.6f  s2=%0.6f  t1=%0.6f  t2=%0.6f \n ",s,s1,s2,t1,t2);
        printf("pent=%e  box=%e  tri=%e  buble=%e  ir=%e\n",pp,bb,tt,bu,ii);
        //printf("u1=%e  u2=%e\n",-(s+t2-s1-1),-(s+t1-s2-1));
    };
    */
  return a;
  
};


/*

double nlo_qcd(double *vars){
    double s=vars[0], s1=vars[1], s2=vars[2], t1=vars[3], t2=vars[4], m=vars[5];
    double a;
    a=(qcd(s,s1,s2,t1,t2,m))/2.;
    if( isnanl(a) )
    { 
        printf("nlo_qcd function is NAN:\n");
        printf("s=%0.6f s1=%0.6f  s2=%0.6f  t1=%0.6f  t2=%0.6f \n ",s,s1,s2,t1,t2);
                //printf("u1=%e  u2=%e\n",-(s+t2-s1-1),-(s+t1-s2-1));
    };
       
    return a;
};

*/