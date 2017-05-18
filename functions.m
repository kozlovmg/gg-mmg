
(*Simple functions*)
R1[x_]:=(1/x-1) Log[Abs[1-x]];
R2[x_]:=-Sqrt[1-4/x] Log[Abs[(1+Sqrt[1-4/x])/(1-Sqrt[1-4/x])]];

T1[x_]:=(Pi^2/6-PolyLog[2,1+x])/x;
T3[t1_,t2_]:=-(PolyLog[2,t1]-PolyLog[2,t2])/(t1-t2);
T4[x_]:=(PolyLog[2,1-x]+Log[x] Log[Abs[x-1]])/(x-1);
T5[x_]:=(Log[(1+Sqrt[1-4/x])/(1-Sqrt[1-4/x])]^2-Pi^2 HeavisideTheta[x])/2/x;
T6[s_,t_]:=(Log[Abs[(1+Sqrt[1-4/s])/(1-Sqrt[1-4/s])]]^2-Log[Abs[(1+Sqrt[1-4/s])/(1-Sqrt[1-4/s])]]^2+Pi^2(HeavisideTheta[s]-HeavisideTheta[t]))/2/(s+t);

F[x_]:=x Log[Abs[(1+Sqrt[1-1/x^2])/(1-Sqrt[1-1/x^2])]]/Sqrt[x^2-1];


(*Complicated functions*)

(*special functions*)

R[x_,y_]:=PolyLog[2,x/(x-y)]-PolyLog[2,(x-1)/(x-y)];
S3[x_,x1_,x2_]:=R[x,x1]+R[x,x2];

ff[x_, y_] := (PolyLog[2,1-x y]+Log[Abs[1-x y]] Log[Abs[x y]]) HeavisideTheta[x y] + (Pi^2/6-PolyLog[2, x y]) HeavisideTheta[-x y]-Log[Abs[1-x y]] (Log[Abs[x]]+Log[Abs[y]])+Pi^2 (HeavisideTheta[-x]+HeavisideTheta[-y]) HeavisideTheta[x y-1];



T2[t1_,t2_]:=Module[{x1,x2,x3,x31,x32,D3,tmp},
    D3=Sqrt[(t1+t2-1)^2-4 t1 t2];
    x1=(1-t2^2-t2+t1 t2+(1+t1) D3)/2/t1/D3;
    x2=-(-2+2*t1+3*t2+t1*t2-t2^2+(t2-2) D3)/(D3*(1+t1-t2+D3));
    x3=-(-2+2*t1+3*t2+t1*t2-t2^2+(t2-2) D3)/(D3*(1-t1-t2+D3));
    x31=(-Sqrt[(-4 + t2) t2] + t2)/(2*t2);
    x32=(Sqrt[(-4 + t2) t2] + t2)/(2*t2);
    tmp=(S3[x1,1,1/t1]-S3[x2,1,1]+S3[x3,x31,x32])/D3;
    tmp 
        ];

(* B1 ????? *)
B1[s_,t_,m_]:=Module[{D1,D2,y1p,y1m,y2p,y2m,b,x1p,x1m,x2p,x2m,x3p,x3m,x4p,x4m,tmp},
    D1=4t+(1+t-m)^2;
    D2=s^2(-4+t)+9 t-4 m+s(4-2 t+4 m);
    b=-(1+t-m+Sqrt[D1])/2;
    y1p=-((-2 t+s t+t^2+t Sqrt[D1]-Sqrt[t D2]-t m)/(-1+s-3 t+s t+(1-s) Sqrt[D1]+m-s m));
    y1m=-((-2 t+s t+t^2+t Sqrt[D1]+Sqrt[t D2]-t m)/(-1+s-3 t+s t+(1-s) Sqrt[D1]+m-s m));
    y2p=(t+s*t+Sqrt[t D2])/(2*(s^2-2*t+s*(-1+t-m)+m));
    y2m=(t+s*t-Sqrt[t D2])/(2*(s^2-2*t+s*(-1+t-m)+m));
    x1p=-1+Sqrt[2];
    x1m=-1-Sqrt[2];
    x2p=(1+m+Sqrt[8+(m-1)^2])/2/(m-2);
    x2m=(1+m-Sqrt[8+(m-1)^2])/2/(m-2);
    x3p=(1+Sqrt[1-4/t])/2;
    x3m=(1-Sqrt[1-4/t])/2;
    x4p=(1+s+Sqrt[8+(s-1)^2])/2/(s-2);
    x4m=(1+s-Sqrt[8+(s-1)^2])/2/(s-2);
    tmp=1/Sqrt[t D2] (-(R[y1p+b, 0]-S3[y1p+b,x1p,x1m]-R[y1p/(1-b),t/(s+t-1)]+S3[y1p/(1-b),x2p,x2m]+ R[-y1p/b,1]-S3[-y1p/b,x3p,x3m])+(R[y2p,0]-S3[y2p,x4p,x4m]-R[y2p,t/(s+t-1)]+S3[y2p,x2p,x2m])+(R[y1m+b,0]-S3[y1m+b,x1p,x1m]-R[y1m/(1-b),t/(s+t-1)]+S3[y1m/(1-b),x2p,x2m]+R[-y1m/b,1]-S3[-y1m/b, x3p,x3m])-(R[y2m,0]-S3[y2m,x4p,x4m]-R[y2m, t/(s+t-1)] + S3[y2m,x2p,x2m]));
    tmp
];

B2[s_,t_,m_]:=Module[{b,X,tmp},
    b[a_]:=Sqrt[1-4/a];
    X[a_]:=(b[a]-1)/(b[a]+1);
    tmp=(2 ff[X[s],X[m]]+2 ff[X[s],1/X[m]]-PolyLog[2,X[s]^2]-2 Log[Abs[X[s]]] Log[1-X[s]^2]-2 Log[Abs[X[s]]] Log[1-t]-Log[X[m]]^2+Pi^2/6)/b[s]/s/(t-1);
    tmp
];

B3[s_,t_,m_]:=Module[{y1p,y1m,x1p,x1m,x2p,x2m,y2p,y2m,x3p,x3m,tmp},
    y1p=(-1+Sqrt[1+4 (m-s-t)/(s t)])/2;
    y1m=(-1-Sqrt[1+4 (m-s-t)/(s t)])/2;
    y2p=t/(s+t-m)/2 (1+Sqrt[1+4 (m-s-t)/(s t)]);
    y2m=t/(s+t-m)/2 (1-Sqrt[1+4 (m-s-t)/(s t)]);
    x1p=(1+Sqrt[1-4/m])/2;
    x1m=(1-Sqrt[1-4/m])/2;
    x2p=(1+Sqrt[1-4/t])/2;
    x2m=(1-Sqrt[1-4/t])/2;
    x3p=(1+Sqrt[1-4/s])/2;
    x3m=(1-Sqrt[1-4/s])/2;
    tmp=1/(s t Sqrt[1+4 (m-s-t)/(s t)]) (-(R[1+y1p,0]-S3[1+y1p,x1p,x1m]+R[-y1p,1]-S3[-y1p,x2p,x2m])+(R[y2p,0]- S3[y2p,x3p,x3m]-R[y2p,-t/(m-t-s)]) + (R[1+y1m,0]-S3[1+y1m,x1p,x1m]+R[-y1m,1]-S3[-y1m,x2p,x2m])-(R[y2m,0]- S3[y2m,x3p,x3m]-R[y2m,-t/(m-t-s)]));
    tmp
];





(*Integral functions*)

G[a_,b_,c_]:=NIntegrate[Log[(b-c x+Sqrt[1-x^2+(b x-c)^2])/2/Sqrt[1-x^2]]/(1-a x),{x,-1,1}];