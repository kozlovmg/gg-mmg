(*special functions*)



Li2 = Compile[{x}, 
  If[x <= 1, PolyLog[2, x], 
   Pi^2/6 - Log[x] Log[Abs[1 - x]] - PolyLog[2, 1 - x]]];

ln=Compile[{{x,_Real}},Log[Abs[x]]];
   
ht=Compile[{{x,_Real}},HeavisideTheta[x]];

R=Compile[{x,y}, Li2[x/(x-y)]-Li2[(x-1)/(x-y)] ];
S3=Compile[{x,x1,x2}, R[x,x1]+R[x,x2] ];

ff=Compile[{x,y}, (Li2[1-x y]+ln[1-x y] ln[x y]) ht[x y] + (Pi^2/6-Li2[x y]) ht[-x y]-ln[1-x y] (ln[x]+ln[y])+Pi^2 (ht[-x]+ht[-y]) ht[x y-1] ];

    
bb=Compile[{{a,_Real}}, Sqrt[1-4/a] ];
X=Compile[{a}, (bb[a]-1)/(bb[a]+1) ];



(***********************************************************************************************************)
(*Simple functions*)
R1=Compile[{{x,_Real}}, (1/x-1) ln[1-x] ];
R2=Compile[{{x,_Real}}, -Sqrt[1-4/x] ln[(1+Sqrt[1-4/x])/(1-Sqrt[1-4/x])] ];

T1=Compile[{{x,_Real}}, (Pi^2/6-Li2[1+x])/x ];
T3=Compile[{{t1,_Real},{t2,_Real}}, -(Li2[t1]-Li2[t2])/(t1-t2) ];
T4=Compile[{{x,_Real}}, (Li2[1-x]+ln[x] ln[x-1])/(x-1) ];
T5=Compile[{{x,_Real}}, (ln[(1+Sqrt[1-4/x])/(1-Sqrt[1-4/x])]^2-Pi^2 HeavisideTheta[x])/2/x ];
T6=Compile[{{s,_Real},{t,_Real}}, (ln[(1+Sqrt[1-4/s])/(1-Sqrt[1-4/s])]^2-ln[(1+Sqrt[1-4/t])/(1-Sqrt[1-4/t])]^2-Pi^2(HeavisideTheta[s]-HeavisideTheta[t]))/2/(s-t) ];

F=Compile[{{x,_Real}}, x ln[(1+Sqrt[1-1/x^2])/(1-Sqrt[1-1/x^2])]/Sqrt[x^2-1] ];


(*Complicated functions*)




T2=Compile[{{t1,_Real},{t2,_Real}}, 
    Module[{x1,x2,x3,x31,x32,D3,tmp},
        D3=Sqrt[(t1+t2-1)^2-4 t1 t2];
        x1=(1-t1^2-t2+t1 t2+(1+t1) D3)/2/t1/D3;
        x2=-(-2+2*t1+3*t2+t1*t2-t2^2+(t2-2) D3)/(D3*(1+t1-t2+D3));
        x3=-(-2+2*t1+3*t2+t1*t2-t2^2+(t2-2) D3)/(D3*(1-t1-t2+D3));
        x31=(-Sqrt[(-4 + t2) t2] + t2)/(2*t2);
        x32=(Sqrt[(-4 + t2) t2] + t2)/(2*t2);
        tmp=(S3[x1,1,1/t1]-S3[x2,1,1]+S3[x3,x31,x32])/D3;
        tmp 
        ] 
    ];

(* B1 ????? *)
B1=Compile[{{s,_Real},{t,_Real},{m,_Real}},         
    Module[{D1,D2,y1p,y1m,y2p,y2m,b,x1p,x1m,x2p,x2m,x3p,x3m,x4p,x4m,tmp},
        D1=4 t+(1+t-m)^2;
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
        ]
    ];

  
    
(*B2=Compile[{{s,_Real},{t,_Real},{m,_Real}},         
    Module[{tmp},
        tmp=(2 ff[X[s],X[m]]+2 ff[X[s],1/X[m]]-Li2[X[s]^2]-2 Log[Abs[X[s]]] Log[1-X[s]^2]-2 Log[Abs[X[s]]] Log[1-t]-Log[X[m]]^2+Pi^2/6)/bb[s]/s/(t-1);
        tmp
    ]
];*)

B2=Compile[{{s,_Real},{t,_Real},{m,_Real}},         
    Module[{tmp},
        tmp=(-2(ff[X[s],X[m]]+ff[X[s],1/X[m]])-Li2[X[s]^2]-2 ln[X[s]] ln[1-X[s]^2]-2 (ln[X[s]] ln[1-t]-Pi^2 ht[t-1] ht[-X[s]])-(ln[X[m]]^2-Pi^2 ht[-X[m]])+Pi^2/6)/bb[s]/s/(t-1);
        tmp
    ]
];



B3=Compile[{{s,_Real},{t,_Real},{m,_Real}},         
    Module[{y1p,y1m,x1p,x1m,x2p,x2m,y2p,y2m,x3p,x3m,tmp},
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
    ]
];





(*Integral functions*)

G=Compile[{{a,_Real},{b,_Real},{c,_Real}},         
        NIntegrate[Log[(b-c x+Sqrt[1-x^2+(b x-c)^2])/2/Sqrt[1-x^2]]/(1-a x),{x,-1,1}]   ];