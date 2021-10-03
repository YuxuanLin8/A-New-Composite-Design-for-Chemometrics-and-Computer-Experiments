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

%L[0,14]^2
Lprac1f4=L2(:,1);
Lprac1f4(Lprac1f4==1)=0;
Lprac1f4(Lprac1f4==2)=14/3;
Lprac1f4(Lprac1f4==3)=28/3;
Lprac1f4(Lprac1f4==4)=14;
Lprac2f4=L2(:,2);
Lprac2f4(Lprac2f4==1)=0;
Lprac2f4(Lprac2f4==2)=14/3;
Lprac2f4(Lprac2f4==3)=28/3;
Lprac2f4(Lprac2f4==4)=14;

%S---> [7/4,49/4]^2 
Sprac1f4=L2(:,1);
Sprac1f4(Sprac1f4==1)=7/4;
Sprac1f4(Sprac1f4==2)=7/4+7/2;
Sprac1f4(Sprac1f4==3)=7/4+7;
Sprac1f4(Sprac1f4==4)=49/4;
Sprac2f4=L2(:,2);
Sprac2f4(Sprac2f4==1)=7/4;
Sprac2f4(Sprac2f4==2)=7/4+7/2;
Sprac2f4(Sprac2f4==3)=7/4+7;
Sprac2f4(Sprac2f4==4)=49/4;

%D
Dprac1f4=L2(:,1);
Dprac1f4(Dprac1f4==1)=0;
Dprac1f4(Dprac1f4==2)=7-7/sqrt(5);
Dprac1f4(Dprac1f4==3)=7+7/sqrt(5);
Dprac1f4(Dprac1f4==4)=14;
Dprac2f4=L2(:,2);
Dprac2f4(Dprac2f4==1)=0;
Dprac2f4(Dprac2f4==2)=7-7/sqrt(5);
Dprac2f4(Dprac2f4==3)=7+7/sqrt(5);
Dprac2f4(Dprac2f4==4)=14;

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

Uprac1f4=14/15*(U2(:,1)-1);
Uprac2f4=14/15*(U2(:,2)-1);

Rprac1f4=7/16+7/8*(U2(:,1)-1);
Rprac2f4=7/16+7/8*(U2(:,2)-1);


y_L4=f4(Lprac1f4,Lprac2f4);
y_S4=f4(Sprac1f4,Sprac2f4);
y_D4=f4(Dprac1f4,Dprac2f4);
y_U4=f4(Uprac1f4,Uprac2f4);
y_R4=f4(Rprac1f4,Rprac2f4);

theta = [10 10]; lob = [1e-1 1e-1]; upb = [20 20];
x1newf4=-2+4*rand(1000,1);
x2newf4=-2+4*rand(1000,1);
newdataf4=[x1newf4,x2newf4];

ynewf4=f4(x1newf4,x2newf4);

%L poly0
[dmodel_L0f4, perf_L0f4] = dacefit([Lprac1f4,Lprac2f4],y_L4, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_L0f4] = predictor(newdataf4, dmodel_L0f4);
sMSE_L0f4=sqrt(1/1000*sum((ynewf4-Yhat_L0f4).^2))

%L poly1
[dmodel_L1f4, perf_L1f4] = dacefit([Lprac1f4,Lprac2f4],y_L4, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_L1f4] = predictor(newdataf4, dmodel_L1f4);
sMSE_L1f4=sqrt(1/1000*sum((ynewf4-Yhat_L1f4).^2))

%L poly2
[dmodel_L2f4, perf_L2f4] = dacefit([Lprac1f4,Lprac2f4],y_L4, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_L2f4] = predictor(newdataf4, dmodel_L2f4);
sMSE_L2f4=sqrt(1/1000*sum((ynewf4-Yhat_L2f4).^2))

%L quadratic
% beta0=[1;1;1;1;1];
% quafac2=@(b,x)b(1)+b(2)*x(:,1)+b(3)*x(:,1).^2+b(4)*x(:,2)+b(5)*x(:,2).^2;
% [betaLf4]=nlinfit([Lprac1f4,Lprac2f4],y_L4,quafac2,beta0);
% Yhat_Lquaf4=quafac2(betaLf4,newdataf4);
% sMSE_Lquaf4=sqrt(1/1000*sum((ynewf4-Yhat_Lquaf4).^2))


%S poly0
[dmodel_S0f4, perf_S0f4] = dacefit([Sprac1f4,Sprac2f4],y_S4, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_S0f4] = predictor(newdataf4, dmodel_S0f4);
sMSE_S0f4=sqrt(1/1000*sum((ynewf4-Yhat_S0f4).^2))

%S poly1
[dmodel_S1f4, perf_S1f4] = dacefit([Sprac1f4,Sprac2f4],y_S4, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_S1f4] = predictor(newdataf4, dmodel_S1f4);
sMSE_S1f4=sqrt(1/1000*sum((ynewf4-Yhat_S1f4).^2))

%S poly2
[dmodel_S2f4, perf_S2f4] = dacefit([Sprac1f4,Sprac2f4],y_S4, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_S2f4] = predictor(newdataf4, dmodel_S2f4);
sMSE_S2f4=sqrt(1/1000*sum((ynewf4-Yhat_S2f4).^2))

%S quadratic
% [betaSf4]=nlinfit([Sprac1f4,Sprac2f4],y_S4,quafac2,beta0);
% Yhat_Squaf4=quafac2(betaSf4,newdataf4);
% sMSE_Squaf4=sqrt(1/1000*sum((ynewf4-Yhat_Squaf4).^2))


%U poly0
[dmodel_U0f4, perf_U0f4] = dacefit([Uprac1f4,Uprac2f4],y_U4, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_U0f4] = predictor(newdataf4, dmodel_U0f4);
sMSE_U0f4=sqrt(1/1000*sum((ynewf4-Yhat_U0f4).^2))

%U poly1
[dmodel_U1f4, perf_U1f4] = dacefit([Uprac1f4,Uprac2f4],y_U4, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_U1f4] = predictor(newdataf4, dmodel_U1f4);
sMSE_U1f4=sqrt(1/1000*sum((ynewf4-Yhat_U1f4).^2))

%U poly2
[dmodel_U2f4, perf_U2f4] = dacefit([Uprac1f4,Uprac2f4],y_U4, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_U2f4] = predictor(newdataf4, dmodel_U2f4);
sMSE_U2f4=sqrt(1/1000*sum((ynewf4-Yhat_U2f4).^2))

%U quadratic
% [betaUf4]=nlinfit([Uprac1f4,Uprac2f4],y_U4,quafac2,beta0);
% Yhat_Uquaf4=quafac2(betaUf4,newdataf4);
% sMSE_Uquaf4=sqrt(1/1000*sum((ynewf4-Yhat_Uquaf4).^2))

%R poly0
[dmodel_R0f4, perf_R0f4] = dacefit([Rprac1f4,Rprac2f4],y_R4, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_R0f4] = predictor(newdataf4, dmodel_R0f4);
sMSE_R0f4=sqrt(1/1000*sum((ynewf4-Yhat_R0f4).^2))

%R poly1
[dmodel_R1f4, perf_R1f4] = dacefit([Rprac1f4,Rprac2f4],y_R4, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_R1f4] = predictor(newdataf4, dmodel_R1f4);
sMSE_R1f4=sqrt(1/1000*sum((ynewf4-Yhat_R1f4).^2))

%R poly2
[dmodel_R2f4, perf_R2f4] = dacefit([Rprac1f4,Rprac2f4],y_R4, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_R2f4] = predictor(newdataf4, dmodel_R2f4);
sMSE_R2f4=sqrt(1/1000*sum((ynewf4-Yhat_R2f4).^2))

%R quadratic
% [betaRf4]=nlinfit([Rprac1f4,Rprac2f4],y_R4,quafac2,beta0);
% Yhat_Rquaf4=quafac2(betaRf4,newdataf4);
% sMSE_Rquaf4=sqrt(1/1000*sum((ynewf4-Yhat_Rquaf4).^2))

%D poly0
[dmodel_D0f4, perf_D0f4] = dacefit([Dprac1f4,Dprac2f4],y_D4, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_D0f4] = predictor(newdataf4, dmodel_D0f4);
sMSE_D0f4=sqrt(1/1000*sum((ynewf4-Yhat_D0f4).^2))

%D poly1
[dmodel_D1f4, perf_D1f4] = dacefit([Dprac1f4,Dprac2f4],y_D4, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_D1f4] = predictor(newdataf4, dmodel_D1f4);
sMSE_D1f4=sqrt(1/1000*sum((ynewf4-Yhat_D1f4).^2))

%D poly2
[dmodel_D2f4, perf_D2f4] = dacefit([Dprac1f4,Dprac2f4],y_D4, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_D2f4] = predictor(newdataf4, dmodel_D2f4);
sMSE_D2f4=sqrt(1/1000*sum((ynewf4-Yhat_D2f4).^2))

%D quadratic
% [betaDf4]=nlinfit([Dprac1f4,Dprac2f4],y_D4,quafac2,beta0);
% Yhat_Dquaf4=quafac2(betaDf4,newdataf4);
% sMSE_Dquaf4=sqrt(1/1000*sum((ynewf4-Yhat_Dquaf4).^2))

%L+D poly0/poly1/poly2
LDdataf4=[Lprac1f4,Lprac2f4;Dprac1f4,Dprac2f4];
Y_LDf4=[y_L4;y_D4];
[LDdataf4_rev,ILDf4]=unique(LDdataf4,'rows');
YLDf4_rev=Y_LDf4(ILDf4);

[dmodel_LD0f4, perf_LD0f4] = dacefit(LDdataf4_rev,YLDf4_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_LD0f4] = predictor(newdataf4, dmodel_LD0f4);
sMSE_LD0f4=sqrt(1/1000*sum((ynewf4-Yhat_LD0f4).^2))


[dmodel_LD1f4, perf_LD1f4] = dacefit(LDdataf4_rev,YLDf4_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_LD1f4] = predictor(newdataf4, dmodel_LD1f4);
sMSE_LD1f4=sqrt(1/1000*sum((ynewf4-Yhat_LD1f4).^2))


[dmodel_LD2f4, perf_LD2f4] = dacefit(LDdataf4_rev,YLDf4_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_LD2f4] = predictor(newdataf4, dmodel_LD2f4);
sMSE_LD2f4=sqrt(1/1000*sum((ynewf4-Yhat_LD2f4).^2))

%L+D quadratic
% [betaLDf4]=nlinfit(LDdataf4_rev,YLDf4_rev,quafac2,beta0);
% Yhat_LDquaf4=quafac2(betaLDf4,newdataf4);
% sMSE_LDquaf4=sqrt(1/1000*sum((ynewf4-Yhat_LDquaf4).^2))

%S+D poly0/poly1/poly2
SDdataf4=[Sprac1f4,Sprac2f4;Dprac1f4,Dprac2f4];
Y_SDf4=[y_S4;y_D4];
[SDdataf4_rev,ISDf4]=unique(SDdataf4,'rows');
YSDf4_rev=Y_SDf4(ISDf4);

[dmodel_SD0f4, perf_SD0f4] = dacefit(SDdataf4_rev,YSDf4_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_SD0f4] = predictor(newdataf4, dmodel_SD0f4);
sMSE_SD0f4=sqrt(1/1000*sum((ynewf4-Yhat_SD0f4).^2))


[dmodel_SD1f4, perf_SD1f4] = dacefit(SDdataf4_rev,YSDf4_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_SD1f4] = predictor(newdataf4, dmodel_SD1f4);
sMSE_SD1f4=sqrt(1/1000*sum((ynewf4-Yhat_SD1f4).^2))


[dmodel_SD2f4, perf_SD2f4] = dacefit(SDdataf4_rev,YSDf4_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_SD2f4] = predictor(newdataf4, dmodel_SD2f4);
sMSE_SD2f4=sqrt(1/1000*sum((ynewf4-Yhat_SD2f4).^2))

%S+D quadratic
% [betaSDf4]=nlinfit(SDdataf4_rev,YSDf4_rev,quafac2,beta0);
% Yhat_SDquaf4=quafac2(betaSDf4,newdataf4);
% sMSE_SDquaf4=sqrt(1/1000*sum((ynewf4-Yhat_SDquaf4).^2))

%U+D poly0/poly1/poly2
UDdataf4=[Uprac1f4,Uprac2f4;Dprac1f4,Dprac2f4];
Y_UDf4=[y_U4;y_D4];
[UDdataf4_rev,IUDf4]=unique(UDdataf4,'rows');
YUDf4_rev=Y_UDf4(IUDf4);

[dmodel_UD0f4, perf_UD0f4] = dacefit(UDdataf4_rev,YUDf4_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_UD0f4] = predictor(newdataf4, dmodel_UD0f4);
sMSE_UD0f4=sqrt(1/1000*sum((ynewf4-Yhat_UD0f4).^2))


[dmodel_UD1f4, perf_UD1f4] = dacefit(UDdataf4_rev,YUDf4_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_UD1f4] = predictor(newdataf4, dmodel_UD1f4);
sMSE_UD1f4=sqrt(1/1000*sum((ynewf4-Yhat_UD1f4).^2))


[dmodel_UD2f4, perf_UD2f4] = dacefit(UDdataf4_rev,YUDf4_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_UD2f4] = predictor(newdataf4, dmodel_UD2f4);
sMSE_UD2f4=sqrt(1/1000*sum((ynewf4-Yhat_UD2f4).^2))

%U+D quadratic
% [betaUDf4]=nlinfit(UDdataf4_rev,YUDf4_rev,quafac2,beta0);
% Yhat_UDquaf4=quafac2(betaUDf4,newdataf4);
% sMSE_UDquaf4=sqrt(1/1000*sum((ynewf4-Yhat_UDquaf4).^2))


%R+D poly0/poly1/poly2
RDdataf4=[Rprac1f4,Rprac2f4;Dprac1f4,Dprac2f4];
Y_RDf4=[y_R4;y_D4];
[RDdataf4_rev,IRDf4]=unique(RDdataf4,'rows');
YRDf4_rev=Y_RDf4(IRDf4);

[dmodel_RD0f4, perf_RD0f4] = dacefit(RDdataf4_rev,YRDf4_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_RD0f4] = predictor(newdataf4, dmodel_RD0f4);
sMSE_RD0f4=sqrt(1/1000*sum((ynewf4-Yhat_RD0f4).^2))


[dmodel_RD1f4, perf_RD1f4] = dacefit(RDdataf4_rev,YRDf4_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_RD1f4] = predictor(newdataf4, dmodel_RD1f4);
sMSE_RD1f4=sqrt(1/1000*sum((ynewf4-Yhat_RD1f4).^2))


[dmodel_RD2f4, perf_RD2f4] = dacefit(RDdataf4_rev,YRDf4_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_RD2f4] = predictor(newdataf4, dmodel_RD2f4);
sMSE_RD2f4=sqrt(1/1000*sum((ynewf4-Yhat_RD2f4).^2))

%R+D quadratic
% [betaRDf4]=nlinfit(RDdataf4_rev,YRDf4_rev,quafac2,beta0);
% Yhat_RDquaf4=quafac2(betaRDf4,newdataf4);
% sMSE_RDquaf4=sqrt(1/1000*sum((ynewf4-Yhat_RDquaf4).^2))

%R+E poly0/poly1/poly2
RLdataf4=[Rprac1f4,Rprac2f4;Lprac1f4,Lprac2f4];
Y_RLf4=[y_R4;y_L4];
[RLdataf4_rev,IRLf4]=unique(RLdataf4,'rows');
YRLf4_rev=Y_RLf4(IRLf4);

[dmodel_RL0f4, perf_RL0f4] = dacefit(RLdataf4_rev,YRLf4_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_RL0f4] = predictor(newdataf4, dmodel_RL0f4);
sMSE_RL0f4=sqrt(1/1000*sum((ynewf4-Yhat_RL0f4).^2))


[dmodel_RL1f4, perf_RL1f4] = dacefit(RLdataf4_rev,YRLf4_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_RL1f4] = predictor(newdataf4, dmodel_RL1f4);
sMSE_RL1f4=sqrt(1/1000*sum((ynewf4-Yhat_RL1f4).^2))


[dmodel_RL2f4, perf_RL2f4] = dacefit(RLdataf4_rev,YRLf4_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_RL2f4] = predictor(newdataf4, dmodel_RL2f4);
sMSE_RL2f4=sqrt(1/1000*sum((ynewf4-Yhat_RL2f4).^2))

%R+E quadratic
% [betaRLf4]=nlinfit(RLdataf4_rev,YRLf4_rev,quafac2,beta0);
% Yhat_RLquaf4=quafac2(betaRLf4,newdataf4);
% sMSE_RLquaf4=sqrt(1/1000*sum((ynewf4-Yhat_RLquaf4).^2))

%R+S poly0/poly1/poly2
RSdataf4=[Rprac1f4,Rprac2f4;Sprac1f4,Sprac2f4];
Y_RSf4=[y_R4;y_S4];
[RSdataf4_rev,IRSf4]=unique(RSdataf4,'rows');
YRSf4_rev=Y_RSf4(IRSf4);

[dmodel_RS0f4, perf_RS0f4] = dacefit(RSdataf4_rev,YRSf4_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_RS0f4] = predictor(newdataf4, dmodel_RS0f4);
sMSE_RS0f4=sqrt(1/1000*sum((ynewf4-Yhat_RS0f4).^2))


[dmodel_RS1f4, perf_RS1f4] = dacefit(RSdataf4_rev,YRSf4_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_RS1f4] = predictor(newdataf4, dmodel_RS1f4);
sMSE_RS1f4=sqrt(1/1000*sum((ynewf4-Yhat_RS1f4).^2))


[dmodel_RS2f4, perf_RS2f4] = dacefit(RSdataf4_rev,YRSf4_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_RS2f4] = predictor(newdataf4, dmodel_RS2f4);
sMSE_RS2f4=sqrt(1/1000*sum((ynewf4-Yhat_RS2f4).^2))

%R+S quadratic
% [betaRSf4]=nlinfit(RSdataf4_rev,YRSf4_rev,quafac2,beta0);
% Yhat_RSquaf4=quafac2(betaRSf4,newdataf4);
% sMSE_RSquaf4=sqrt(1/1000*sum((ynewf4-Yhat_RSquaf4).^2))

%U+S poly0/poly1/poly2
USdataf4=[Uprac1f4,Uprac2f4;Sprac1f4,Sprac2f4];
Y_USf4=[y_U4;y_S4];
[USdataf4_rev,IUSf4]=unique(USdataf4,'rows');
YUSf4_rev=Y_USf4(IUSf4);

[dmodel_US0f4, perf_US0f4] = dacefit(USdataf4_rev,YUSf4_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_US0f4] = predictor(newdataf4, dmodel_US0f4);
sMSE_US0f4=sqrt(1/1000*sum((ynewf4-Yhat_US0f4).^2))

[dmodel_US1f4, perf_US1f4] = dacefit(USdataf4_rev,YUSf4_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_US1f4] = predictor(newdataf4, dmodel_US1f4);
sMSE_US1f4=sqrt(1/1000*sum((ynewf4-Yhat_US1f4).^2))

[dmodel_US2f4, perf_US2f4] = dacefit(USdataf4_rev,YUSf4_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_US2f4] = predictor(newdataf4, dmodel_US2f4);
sMSE_US2f4=sqrt(1/1000*sum((ynewf4-Yhat_US2f4).^2))

%U+S quadratic
% [betaUSf4]=nlinfit(USdataf4_rev,YUSf4_rev,quafac2,beta0);
% Yhat_USquaf4=quafac2(betaUSf4,newdataf4);
% sMSE_USquaf4=sqrt(1/1000*sum((ynewf4-Yhat_USquaf4).^2))

%U+E poly0/poly1/poly2
UEdataf4=[Uprac1f4,Uprac2f4;Lprac1f4,Lprac2f4];
Y_UEf4=[y_U4;y_L4];
[UEdataf4_rev,IUEf4]=unique(UEdataf4,'rows');
YUEf4_rev=Y_UEf4(IUEf4);

[dmodel_UE0f4, perf_UE0f4] = dacefit(UEdataf4_rev,YUEf4_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_UE0f4] = predictor(newdataf4, dmodel_UE0f4);
sMSE_UE0f4=sqrt(1/1000*sum((ynewf4-Yhat_UE0f4).^2))


[dmodel_UE1f4, perf_UE1f4] = dacefit(UEdataf4_rev,YUEf4_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_UE1f4] = predictor(newdataf4, dmodel_UE1f4);
sMSE_UE1f4=sqrt(1/1000*sum((ynewf4-Yhat_UE1f4).^2))


[dmodel_UE2f4, perf_UE2f4] = dacefit(UEdataf4_rev,YUEf4_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_UE2f4] = predictor(newdataf4, dmodel_UE2f4);
sMSE_UE2f4=sqrt(1/1000*sum((ynewf4-Yhat_UE2f4).^2))

%U+E quadratic
% [betaUEf4]=nlinfit(UEdataf4_rev,YUEf4_rev,quafac2,beta0);
% Yhat_UEquaf4=quafac2(betaUEf4,newdataf4);
% sMSE_UEquaf4=sqrt(1/1000*sum((ynewf4-Yhat_UEquaf4).^2))

