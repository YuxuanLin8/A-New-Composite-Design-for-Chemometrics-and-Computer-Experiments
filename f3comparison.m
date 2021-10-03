L2=[1 	1 
1 	2 
1 	3 
1 	4 
2 	1 
2 	2 
2 	3 
2 	4 
3 	1 
3 	2 
3 	3 
3 	4 
4 	1 
4 	2 
4 	3 
4 	4 ];

%L
Lprac1f3=L2(:,1);
Lprac1f3(Lprac1f3==1)=-2;
Lprac1f3(Lprac1f3==2)=-2/3;
Lprac1f3(Lprac1f3==3)=2/3;
Lprac1f3(Lprac1f3==4)=2;
Lprac2f3=L2(:,2);
Lprac2f3(Lprac2f3==1)=-2;
Lprac2f3(Lprac2f3==2)=-2/3;
Lprac2f3(Lprac2f3==3)=2/3;
Lprac2f3(Lprac2f3==4)=2;

%S---> [-3/2,3/2] x [-3/2,3/2]
Sprac1f3=L2(:,1);
Sprac1f3(Sprac1f3==1)=-3/2;
Sprac1f3(Sprac1f3==2)=-3/2+1;
Sprac1f3(Sprac1f3==3)=-3/2+2;
Sprac1f3(Sprac1f3==4)=3/2;
Sprac2f3=L2(:,2);
Sprac2f3(Sprac2f3==1)=-3/2;
Sprac2f3(Sprac2f3==2)=-3/2+1;
Sprac2f3(Sprac2f3==3)=-3/2+2;
Sprac2f3(Sprac2f3==4)=3/2;

%D
Dprac1f3=L2(:,1);
Dprac1f3(Dprac1f3==1)=-2;
Dprac1f3(Dprac1f3==2)=-2+2*(1-1/sqrt(5));
Dprac1f3(Dprac1f3==3)=-2+2*(1+1/sqrt(5));
Dprac1f3(Dprac1f3==4)=2;
Dprac2f3=L2(:,2);
Dprac2f3(Dprac2f3==1)=-2;
Dprac2f3(Dprac2f3==2)=-2+2*(1-1/sqrt(5));
Dprac2f3(Dprac2f3==3)=-2+2*(1+1/sqrt(5));
Dprac2f3(Dprac2f3==4)=2;

U2=[1 	10 
2 	3 
3 	15 
4 	6 
5 	8 
6 	13 
7 	1 
8 	12 
9 	5 
10 	16 
11 	4 
12 	9 
13 	11 
14 	2 
15 	14 
16 	7 ];

Uprac1f3=-2+4/15*(U2(:,1)-1);
Uprac2f3=-2+4/15*(U2(:,2)-1);

Rprac1f3=-15/8+1/4*(U2(:,1)-1);
Rprac2f3=-15/8+1/4*(U2(:,2)-1);


y_L3=f3(Lprac1f3,Lprac2f3);
y_S3=f3(Sprac1f3,Sprac2f3);
y_D3=f3(Dprac1f3,Dprac2f3);
y_U3=f3(Uprac1f3,Uprac2f3);
y_R3=f3(Rprac1f3,Rprac2f3);

theta = [10 10]; lob = [1e-1 1e-1]; upb = [20 20];
x1newf3=-2+4*rand(1000,1);
x2newf3=-2+4*rand(1000,1);
newdataf3=[x1newf3,x2newf3];

ynewf3=f3(x1newf3,x2newf3);

%L poly0
[dmodel_L0f3, perf_L0f3] = dacefit([Lprac1f3,Lprac2f3],y_L3, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_L0f3] = predictor(newdataf3, dmodel_L0f3);
sMSE_L0f3=sqrt(1/1000*sum((ynewf3-Yhat_L0f3).^2))

%L poly1
[dmodel_L1f3, perf_L1f3] = dacefit([Lprac1f3,Lprac2f3],y_L3, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_L1f3] = predictor(newdataf3, dmodel_L1f3);
sMSE_L1f3=sqrt(1/1000*sum((ynewf3-Yhat_L1f3).^2))

%L poly2
[dmodel_L2f3, perf_L2f3] = dacefit([Lprac1f3,Lprac2f3],y_L3, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_L2f3] = predictor(newdataf3, dmodel_L2f3);
sMSE_L2f3=sqrt(1/1000*sum((ynewf3-Yhat_L2f3).^2))

%L quadratic
beta0=[1;1;1;1;1];
quafac2=@(b,x)b(1)+b(2)*x(:,1)+b(3)*x(:,1).^2+b(4)*x(:,2)+b(5)*x(:,2).^2;
[betaLf3]=nlinfit([Lprac1f3,Lprac2f3],y_L3,quafac2,beta0);
Yhat_Lquaf3=quafac2(betaLf3,newdataf3);
sMSE_Lquaf3=sqrt(1/1000*sum((ynewf3-Yhat_Lquaf3).^2))


%S poly0
[dmodel_S0f3, perf_S0f3] = dacefit([Sprac1f3,Sprac2f3],y_S3, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_S0f3] = predictor(newdataf3, dmodel_S0f3);
sMSE_S0f3=sqrt(1/1000*sum((ynewf3-Yhat_S0f3).^2))

%S poly1
[dmodel_S1f3, perf_S1f3] = dacefit([Sprac1f3,Sprac2f3],y_S3, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_S1f3] = predictor(newdataf3, dmodel_S1f3);
sMSE_S1f3=sqrt(1/1000*sum((ynewf3-Yhat_S1f3).^2))

%S poly2
[dmodel_S2f3, perf_S2f3] = dacefit([Sprac1f3,Sprac2f3],y_S3, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_S2f3] = predictor(newdataf3, dmodel_S2f3);
sMSE_S2f3=sqrt(1/1000*sum((ynewf3-Yhat_S2f3).^2))

%S quadratic
[betaSf3]=nlinfit([Sprac1f3,Sprac2f3],y_S3,quafac2,beta0);
Yhat_Squaf3=quafac2(betaSf3,newdataf3);
sMSE_Squaf3=sqrt(1/1000*sum((ynewf3-Yhat_Squaf3).^2))


%U poly0
[dmodel_U0f3, perf_U0f3] = dacefit([Uprac1f3,Uprac2f3],y_U3, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_U0f3] = predictor(newdataf3, dmodel_U0f3);
sMSE_U0f3=sqrt(1/1000*sum((ynewf3-Yhat_U0f3).^2))

%U poly1
[dmodel_U1f3, perf_U1f3] = dacefit([Uprac1f3,Uprac2f3],y_U3, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_U1f3] = predictor(newdataf3, dmodel_U1f3);
sMSE_U1f3=sqrt(1/1000*sum((ynewf3-Yhat_U1f3).^2))

%U poly2
[dmodel_U2f3, perf_U2f3] = dacefit([Uprac1f3,Uprac2f3],y_U3, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_U2f3] = predictor(newdataf3, dmodel_U2f3);
sMSE_U2f3=sqrt(1/1000*sum((ynewf3-Yhat_U2f3).^2))

%U quadratic
[betaUf3]=nlinfit([Uprac1f3,Uprac2f3],y_U3,quafac2,beta0);
Yhat_Uquaf3=quafac2(betaUf3,newdataf3);
sMSE_Uquaf3=sqrt(1/1000*sum((ynewf3-Yhat_Uquaf3).^2))

%R poly0
[dmodel_R0f3, perf_R0f3] = dacefit([Rprac1f3,Rprac2f3],y_R3, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_R0f3] = predictor(newdataf3, dmodel_R0f3);
sMSE_R0f3=sqrt(1/1000*sum((ynewf3-Yhat_R0f3).^2))

%R poly1
[dmodel_R1f3, perf_R1f3] = dacefit([Rprac1f3,Rprac2f3],y_R3, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_R1f3] = predictor(newdataf3, dmodel_R1f3);
sMSE_R1f3=sqrt(1/1000*sum((ynewf3-Yhat_R1f3).^2))

%R poly2
[dmodel_R2f3, perf_R2f3] = dacefit([Rprac1f3,Rprac2f3],y_R3, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_R2f3] = predictor(newdataf3, dmodel_R2f3);
sMSE_R2f3=sqrt(1/1000*sum((ynewf3-Yhat_R2f3).^2))

%R quadratic
[betaRf3]=nlinfit([Rprac1f3,Rprac2f3],y_R3,quafac2,beta0);
Yhat_Rquaf3=quafac2(betaRf3,newdataf3);
sMSE_Rquaf3=sqrt(1/1000*sum((ynewf3-Yhat_Rquaf3).^2))

%D poly0
[dmodel_D0f3, perf_D0f3] = dacefit([Dprac1f3,Dprac2f3],y_D3, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_D0f3] = predictor(newdataf3, dmodel_D0f3);
sMSE_D0f3=sqrt(1/1000*sum((ynewf3-Yhat_D0f3).^2))

%D poly1
[dmodel_D1f3, perf_D1f3] = dacefit([Dprac1f3,Dprac2f3],y_D3, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_D1f3] = predictor(newdataf3, dmodel_D1f3);
sMSE_D1f3=sqrt(1/1000*sum((ynewf3-Yhat_D1f3).^2))

%D poly2
[dmodel_D2f3, perf_D2f3] = dacefit([Dprac1f3,Dprac2f3],y_D3, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_D2f3] = predictor(newdataf3, dmodel_D2f3);
sMSE_D2f3=sqrt(1/1000*sum((ynewf3-Yhat_D2f3).^2))

%D quadratic
[betaDf3]=nlinfit([Dprac1f3,Dprac2f3],y_D3,quafac2,beta0);
Yhat_Dquaf3=quafac2(betaDf3,newdataf3);
sMSE_Dquaf3=sqrt(1/1000*sum((ynewf3-Yhat_Dquaf3).^2))

%L+D poly0/poly1/poly2
LDdataf3=[Lprac1f3,Lprac2f3;Dprac1f3,Dprac2f3];
Y_LDf3=[y_L3;y_D3];
[LDdataf3_rev,ILDf3]=unique(LDdataf3,'rows');
YLDf3_rev=Y_LDf3(ILDf3);

[dmodel_LD0f3, perf_LD0f3] = dacefit(LDdataf3_rev,YLDf3_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_LD0f3] = predictor(newdataf3, dmodel_LD0f3);
sMSE_LD0f3=sqrt(1/1000*sum((ynewf3-Yhat_LD0f3).^2))


[dmodel_LD1f3, perf_LD1f3] = dacefit(LDdataf3_rev,YLDf3_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_LD1f3] = predictor(newdataf3, dmodel_LD1f3);
sMSE_LD1f3=sqrt(1/1000*sum((ynewf3-Yhat_LD1f3).^2))


[dmodel_LD2f3, perf_LD2f3] = dacefit(LDdataf3_rev,YLDf3_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_LD2f3] = predictor(newdataf3, dmodel_LD2f3);
sMSE_LD2f3=sqrt(1/1000*sum((ynewf3-Yhat_LD2f3).^2))

%L+D quadratic
[betaLDf3]=nlinfit(LDdataf3_rev,YLDf3_rev,quafac2,beta0);
Yhat_LDquaf3=quafac2(betaLDf3,newdataf3);
sMSE_LDquaf3=sqrt(1/1000*sum((ynewf3-Yhat_LDquaf3).^2))

%S+D poly0/poly1/poly2
SDdataf3=[Sprac1f3,Sprac2f3;Dprac1f3,Dprac2f3];
Y_SDf3=[y_S3;y_D3];
[SDdataf3_rev,ISDf3]=unique(SDdataf3,'rows');
YSDf3_rev=Y_SDf3(ISDf3);

[dmodel_SD0f3, perf_SD0f3] = dacefit(SDdataf3_rev,YSDf3_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_SD0f3] = predictor(newdataf3, dmodel_SD0f3);
sMSE_SD0f3=sqrt(1/1000*sum((ynewf3-Yhat_SD0f3).^2))


[dmodel_SD1f3, perf_SD1f3] = dacefit(SDdataf3_rev,YSDf3_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_SD1f3] = predictor(newdataf3, dmodel_SD1f3);
sMSE_SD1f3=sqrt(1/1000*sum((ynewf3-Yhat_SD1f3).^2))


[dmodel_SD2f3, perf_SD2f3] = dacefit(SDdataf3_rev,YSDf3_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_SD2f3] = predictor(newdataf3, dmodel_SD2f3);
sMSE_SD2f3=sqrt(1/1000*sum((ynewf3-Yhat_SD2f3).^2))

%S+D quadratic
[betaSDf3]=nlinfit(SDdataf3_rev,YSDf3_rev,quafac2,beta0);
Yhat_SDquaf3=quafac2(betaSDf3,newdataf3);
sMSE_SDquaf3=sqrt(1/1000*sum((ynewf3-Yhat_SDquaf3).^2))

%U+D poly0/poly1/poly2
UDdataf3=[Uprac1f3,Uprac2f3;Dprac1f3,Dprac2f3];
Y_UDf3=[y_U3;y_D3];
[UDdataf3_rev,IUDf3]=unique(UDdataf3,'rows');
YUDf3_rev=Y_UDf3(IUDf3);

[dmodel_UD0f3, perf_UD0f3] = dacefit(UDdataf3_rev,YUDf3_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_UD0f3] = predictor(newdataf3, dmodel_UD0f3);
sMSE_UD0f3=sqrt(1/1000*sum((ynewf3-Yhat_UD0f3).^2))


[dmodel_UD1f3, perf_UD1f3] = dacefit(UDdataf3_rev,YUDf3_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_UD1f3] = predictor(newdataf3, dmodel_UD1f3);
sMSE_UD1f3=sqrt(1/1000*sum((ynewf3-Yhat_UD1f3).^2))


[dmodel_UD2f3, perf_UD2f3] = dacefit(UDdataf3_rev,YUDf3_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_UD2f3] = predictor(newdataf3, dmodel_UD2f3);
sMSE_UD2f3=sqrt(1/1000*sum((ynewf3-Yhat_UD2f3).^2))

%U+D quadratic
[betaUDf3]=nlinfit(UDdataf3_rev,YUDf3_rev,quafac2,beta0);
Yhat_UDquaf3=quafac2(betaUDf3,newdataf3);
sMSE_UDquaf3=sqrt(1/1000*sum((ynewf3-Yhat_UDquaf3).^2))


%R+D poly0/poly1/poly2
RDdataf3=[Rprac1f3,Rprac2f3;Dprac1f3,Dprac2f3];
Y_RDf3=[y_R3;y_D3];
[RDdataf3_rev,IRDf3]=unique(RDdataf3,'rows');
YRDf3_rev=Y_RDf3(IRDf3);

[dmodel_RD0f3, perf_RD0f3] = dacefit(RDdataf3_rev,YRDf3_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_RD0f3] = predictor(newdataf3, dmodel_RD0f3);
sMSE_RD0f3=sqrt(1/1000*sum((ynewf3-Yhat_RD0f3).^2))


[dmodel_RD1f3, perf_RD1f3] = dacefit(RDdataf3_rev,YRDf3_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_RD1f3] = predictor(newdataf3, dmodel_RD1f3);
sMSE_RD1f3=sqrt(1/1000*sum((ynewf3-Yhat_RD1f3).^2))


[dmodel_RD2f3, perf_RD2f3] = dacefit(RDdataf3_rev,YRDf3_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_RD2f3] = predictor(newdataf3, dmodel_RD2f3);
sMSE_RD2f3=sqrt(1/1000*sum((ynewf3-Yhat_RD2f3).^2))

%R+D quadratic
[betaRDf3]=nlinfit(RDdataf3_rev,YRDf3_rev,quafac2,beta0);
Yhat_RDquaf3=quafac2(betaRDf3,newdataf3);
sMSE_RDquaf3=sqrt(1/1000*sum((ynewf3-Yhat_RDquaf3).^2))

%R+E poly0/poly1/poly2
RLdataf3=[Rprac1f3,Rprac2f3;Lprac1f3,Lprac2f3];
Y_RLf3=[y_R3;y_L3];
[RLdataf3_rev,IRLf3]=unique(RLdataf3,'rows');
YRLf3_rev=Y_RLf3(IRLf3);

[dmodel_RL0f3, perf_RL0f3] = dacefit(RLdataf3_rev,YRLf3_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_RL0f3] = predictor(newdataf3, dmodel_RL0f3);
sMSE_RL0f3=sqrt(1/1000*sum((ynewf3-Yhat_RL0f3).^2))


[dmodel_RL1f3, perf_RL1f3] = dacefit(RLdataf3_rev,YRLf3_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_RL1f3] = predictor(newdataf3, dmodel_RL1f3);
sMSE_RL1f3=sqrt(1/1000*sum((ynewf3-Yhat_RL1f3).^2))


[dmodel_RL2f3, perf_RL2f3] = dacefit(RLdataf3_rev,YRLf3_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_RL2f3] = predictor(newdataf3, dmodel_RL2f3);
sMSE_RL2f3=sqrt(1/1000*sum((ynewf3-Yhat_RL2f3).^2))

%R+E quadratic
[betaRLf3]=nlinfit(RLdataf3_rev,YRLf3_rev,quafac2,beta0);
Yhat_RLquaf3=quafac2(betaRLf3,newdataf3);
sMSE_RLquaf3=sqrt(1/1000*sum((ynewf3-Yhat_RLquaf3).^2))

%R+S poly0/poly1/poly2
RSdataf3=[Rprac1f3,Rprac2f3;Sprac1f3,Sprac2f3];
Y_RSf3=[y_R3;y_S3];
[RSdataf3_rev,IRSf3]=unique(RSdataf3,'rows');
YRSf3_rev=Y_RSf3(IRSf3);

[dmodel_RS0f3, perf_RS0f3] = dacefit(RSdataf3_rev,YRSf3_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_RS0f3] = predictor(newdataf3, dmodel_RS0f3);
sMSE_RS0f3=sqrt(1/1000*sum((ynewf3-Yhat_RS0f3).^2))


[dmodel_RS1f3, perf_RS1f3] = dacefit(RSdataf3_rev,YRSf3_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_RS1f3] = predictor(newdataf3, dmodel_RS1f3);
sMSE_RS1f3=sqrt(1/1000*sum((ynewf3-Yhat_RS1f3).^2))


[dmodel_RS2f3, perf_RS2f3] = dacefit(RSdataf3_rev,YRSf3_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_RS2f3] = predictor(newdataf3, dmodel_RS2f3);
sMSE_RS2f3=sqrt(1/1000*sum((ynewf3-Yhat_RS2f3).^2))

%R+S quadratic
[betaRSf3]=nlinfit(RSdataf3_rev,YRSf3_rev,quafac2,beta0);
Yhat_RSquaf3=quafac2(betaRSf3,newdataf3);
sMSE_RSquaf3=sqrt(1/1000*sum((ynewf3-Yhat_RSquaf3).^2))

%U+S poly0/poly1/poly2
USdataf3=[Uprac1f3,Uprac2f3;Sprac1f3,Sprac2f3];
Y_USf3=[y_U3;y_S3];
[USdataf3_rev,IUSf3]=unique(USdataf3,'rows');
YUSf3_rev=Y_USf3(IUSf3);
[dmodel_US0f3, perf_US0f3] = dacefit(USdataf3_rev,YUSf3_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_US0f3] = predictor(newdataf3, dmodel_US0f3);
sMSE_US0f3=sqrt(1/1000*sum((ynewf3-Yhat_US0f3).^2))

[dmodel_US1f3, perf_US1f3] = dacefit(USdataf3_rev,YUSf3_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_US1f3] = predictor(newdataf3, dmodel_US1f3);
sMSE_US1f3=sqrt(1/1000*sum((ynewf3-Yhat_US1f3).^2))
[dmodel_US2f3, perf_US2f3] = dacefit(USdataf3_rev,YUSf3_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_US2f3] = predictor(newdataf3, dmodel_US2f3);
sMSE_US2f3=sqrt(1/1000*sum((ynewf3-Yhat_US2f3).^2))

%U+S quadratic
[betaUSf3]=nlinfit(USdataf3_rev,YUSf3_rev,quafac2,beta0);
Yhat_USquaf3=quafac2(betaUSf3,newdataf3);
sMSE_USquaf3=sqrt(1/1000*sum((ynewf3-Yhat_USquaf3).^2))

%U+E poly0/poly1/poly2
UEdataf3=[Uprac1f3,Uprac2f3;Lprac1f3,Lprac2f3];
Y_UEf3=[y_U3;y_L3];
[UEdataf3_rev,IUEf3]=unique(UEdataf3,'rows');
YUEf3_rev=Y_UEf3(IUEf3);

[dmodel_UE0f3, perf_UE0f3] = dacefit(UEdataf3_rev,YUEf3_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_UE0f3] = predictor(newdataf3, dmodel_UE0f3);
sMSE_UE0f3=sqrt(1/1000*sum((ynewf3-Yhat_UE0f3).^2))


[dmodel_UE1f3, perf_UE1f3] = dacefit(UEdataf3_rev,YUEf3_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_UE1f3] = predictor(newdataf3, dmodel_UE1f3);
sMSE_UE1f3=sqrt(1/1000*sum((ynewf3-Yhat_UE1f3).^2))


[dmodel_UE2f3, perf_UE2f3] = dacefit(UEdataf3_rev,YUEf3_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_UE2f3] = predictor(newdataf3, dmodel_UE2f3);
sMSE_UE2f3=sqrt(1/1000*sum((ynewf3-Yhat_UE2f3).^2))

%U+E quadratic
[betaUEf3]=nlinfit(UEdataf3_rev,YUEf3_rev,quafac2,beta0);
Yhat_UEquaf3=quafac2(betaUEf3,newdataf3);
sMSE_UEquaf3=sqrt(1/1000*sum((ynewf3-Yhat_UEquaf3).^2))

