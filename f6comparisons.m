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
Lprac1f6=L2(:,1);
Lprac1f6(Lprac1f6==1)=-10;
Lprac1f6(Lprac1f6==2)=-13/3;
Lprac1f6(Lprac1f6==3)=4/3;
Lprac1f6(Lprac1f6==4)=7;
Lprac2f6=L2(:,2);
Lprac2f6(Lprac2f6==1)=-6;
Lprac2f6(Lprac2f6==2)=-5/3;
Lprac2f6(Lprac2f6==3)=8/3;
Lprac2f6(Lprac2f6==4)=7;

%S---> [-63/8,39/8] x [-35/8,43/8]
Sprac1f6=L2(:,1);
Sprac1f6(Sprac1f6==1)=-63/8;
Sprac1f6(Sprac1f6==2)=-29/8;
Sprac1f6(Sprac1f6==3)=5/8;
Sprac1f6(Sprac1f6==4)=39/8;
Sprac2f6=L2(:,2);
Sprac2f6(Sprac2f6==1)=-35/8;
Sprac2f6(Sprac2f6==2)=-9/8;
Sprac2f6(Sprac2f6==3)=17/8;
Sprac2f6(Sprac2f6==4)=43/8;

%D
Dprac1f6=L2(:,1);
Dprac1f6(Dprac1f6==1)=-10;
Dprac1f6(Dprac1f6==2)=-10+17/2*(1-1/sqrt(5));
Dprac1f6(Dprac1f6==3)=-10+17/2*(1+1/sqrt(5));
Dprac1f6(Dprac1f6==4)=7;
Dprac2f6=L2(:,2);
Dprac2f6(Dprac2f6==1)=-6;
Dprac2f6(Dprac2f6==2)=-6+13/2*(1-1/sqrt(5));
Dprac2f6(Dprac2f6==3)=-6+13/2*(1+1/sqrt(5));
Dprac2f6(Dprac2f6==4)=7;

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

Uprac1f6=-10+17/15*(U2(:,1)-1);
Uprac2f6=-6+13/15*(U2(:,2)-1);

Rprac1f6=-303/32+34/32*(U2(:,1)-1);
Rprac2f6=-179/32+26/32*(U2(:,2)-1);


y_L6=f6(Lprac1f6,Lprac2f6);
y_S6=f6(Sprac1f6,Sprac2f6);
y_D6=f6(Dprac1f6,Dprac2f6);
y_U6=f6(Uprac1f6,Uprac2f6);
y_R6=f6(Rprac1f6,Rprac2f6);

theta = [10 10]; lob = [1e-1 1e-1]; upb = [20 20];
x1newf6=-10+17*rand(1000,1);
x2newf6=-6+13*rand(1000,1);
newdataf6=[x1newf6,x2newf6];

ynewf6=f6(x1newf6,x2newf6);

%L poly0
[dmodel_L0f6, perf_L0f6] = dacefit([Lprac1f6,Lprac2f6],y_L6, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_L0f6] = predictor(newdataf6, dmodel_L0f6);
sMSE_L0f6=sqrt(1/1000*sum((ynewf6-Yhat_L0f6).^2))

%L poly1
[dmodel_L1f6, perf_L1f6] = dacefit([Lprac1f6,Lprac2f6],y_L6, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_L1f6] = predictor(newdataf6, dmodel_L1f6);
sMSE_L1f6=sqrt(1/1000*sum((ynewf6-Yhat_L1f6).^2))

%L poly2
[dmodel_L2f6, perf_L2f6] = dacefit([Lprac1f6,Lprac2f6],y_L6, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_L2f6] = predictor(newdataf6, dmodel_L2f6);
sMSE_L2f6=sqrt(1/1000*sum((ynewf6-Yhat_L2f6).^2))

%L quadratic
% % beta0=[1;1;1;1;1];
% % quafac2=@(b,x)b(1)+b(2)*x(:,1)+b(3)*x(:,1).^2+b(4)*x(:,2)+b(5)*x(:,2).^2;
% % [betaLf6]=nlinfit([Lprac1f6,Lprac2f6],y_L6,quafac2,beta0);
% % Yhat_Lquaf6=quafac2(betaLf6,newdataf6);
% % sMSE_Lquaf6=sqrt(1/1000*sum((ynewf6-Yhat_Lquaf6).^2))
% % 

%S poly0
[dmodel_S0f6, perf_S0f6] = dacefit([Sprac1f6,Sprac2f6],y_S6, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_S0f6] = predictor(newdataf6, dmodel_S0f6);
sMSE_S0f6=sqrt(1/1000*sum((ynewf6-Yhat_S0f6).^2))

%S poly1
[dmodel_S1f6, perf_S1f6] = dacefit([Sprac1f6,Sprac2f6],y_S6, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_S1f6] = predictor(newdataf6, dmodel_S1f6);
sMSE_S1f6=sqrt(1/1000*sum((ynewf6-Yhat_S1f6).^2))

%S poly2
[dmodel_S2f6, perf_S2f6] = dacefit([Sprac1f6,Sprac2f6],y_S6, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_S2f6] = predictor(newdataf6, dmodel_S2f6);
sMSE_S2f6=sqrt(1/1000*sum((ynewf6-Yhat_S2f6).^2))

%S quadratic
[betaSf6]=nlinfit([Sprac1f6,Sprac2f6],y_S6,quafac2,beta0);
Yhat_Squaf6=quafac2(betaSf6,newdataf6);
sMSE_Squaf6=sqrt(1/1000*sum((ynewf6-Yhat_Squaf6).^2))


%U poly0
[dmodel_U0f6, perf_U0f6] = dacefit([Uprac1f6,Uprac2f6],y_U6, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_U0f6] = predictor(newdataf6, dmodel_U0f6);
sMSE_U0f6=sqrt(1/1000*sum((ynewf6-Yhat_U0f6).^2))

%U poly1
[dmodel_U1f6, perf_U1f6] = dacefit([Uprac1f6,Uprac2f6],y_U6, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_U1f6] = predictor(newdataf6, dmodel_U1f6);
sMSE_U1f6=sqrt(1/1000*sum((ynewf6-Yhat_U1f6).^2))

%U poly2
[dmodel_U2f6, perf_U2f6] = dacefit([Uprac1f6,Uprac2f6],y_U6, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_U2f6] = predictor(newdataf6, dmodel_U2f6);
sMSE_U2f6=sqrt(1/1000*sum((ynewf6-Yhat_U2f6).^2))

%U quadratic
[betaUf6]=nlinfit([Uprac1f6,Uprac2f6],y_U6,quafac2,beta0);
Yhat_Uquaf6=quafac2(betaUf6,newdataf6);
sMSE_Uquaf6=sqrt(1/1000*sum((ynewf6-Yhat_Uquaf6).^2))

%R poly0
[dmodel_R0f6, perf_R0f6] = dacefit([Rprac1f6,Rprac2f6],y_R6, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_R0f6] = predictor(newdataf6, dmodel_R0f6);
sMSE_R0f6=sqrt(1/1000*sum((ynewf6-Yhat_R0f6).^2))

%R poly1
[dmodel_R1f6, perf_R1f6] = dacefit([Rprac1f6,Rprac2f6],y_R6, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_R1f6] = predictor(newdataf6, dmodel_R1f6);
sMSE_R1f6=sqrt(1/1000*sum((ynewf6-Yhat_R1f6).^2))

%R poly2
[dmodel_R2f6, perf_R2f6] = dacefit([Rprac1f6,Rprac2f6],y_R6, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_R2f6] = predictor(newdataf6, dmodel_R2f6);
sMSE_R2f6=sqrt(1/1000*sum((ynewf6-Yhat_R2f6).^2))

%R quadratic
[betaRf6]=nlinfit([Rprac1f6,Rprac2f6],y_R6,quafac2,beta0);
Yhat_Rquaf6=quafac2(betaRf6,newdataf6);
sMSE_Rquaf6=sqrt(1/1000*sum((ynewf6-Yhat_Rquaf6).^2))

%D poly0
[dmodel_D0f6, perf_D0f6] = dacefit([Dprac1f6,Dprac2f6],y_D6, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_D0f6] = predictor(newdataf6, dmodel_D0f6);
sMSE_D0f6=sqrt(1/1000*sum((ynewf6-Yhat_D0f6).^2))

%D poly1
[dmodel_D1f6, perf_D1f6] = dacefit([Dprac1f6,Dprac2f6],y_D6, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_D1f6] = predictor(newdataf6, dmodel_D1f6);
sMSE_D1f6=sqrt(1/1000*sum((ynewf6-Yhat_D1f6).^2))

%D poly2
[dmodel_D2f6, perf_D2f6] = dacefit([Dprac1f6,Dprac2f6],y_D6, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_D2f6] = predictor(newdataf6, dmodel_D2f6);
sMSE_D2f6=sqrt(1/1000*sum((ynewf6-Yhat_D2f6).^2))

%D quadratic
[betaDf6]=nlinfit([Dprac1f6,Dprac2f6],y_D6,quafac2,beta0);
Yhat_Dquaf6=quafac2(betaDf6,newdataf6);
sMSE_Dquaf6=sqrt(1/1000*sum((ynewf6-Yhat_Dquaf6).^2))

%L+D poly0/poly1/poly2
LDdataf6=[Lprac1f6,Lprac2f6;Dprac1f6,Dprac2f6];
Y_LDf6=[y_L6;y_D6];
[LDdataf6_rev,ILDf6]=unique(LDdataf6,'rows');
YLDf6_rev=Y_LDf6(ILDf6);

[dmodel_LD0f6, perf_LD0f6] = dacefit(LDdataf6_rev,YLDf6_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_LD0f6] = predictor(newdataf6, dmodel_LD0f6);
sMSE_LD0f6=sqrt(1/1000*sum((ynewf6-Yhat_LD0f6).^2))


[dmodel_LD1f6, perf_LD1f6] = dacefit(LDdataf6_rev,YLDf6_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_LD1f6] = predictor(newdataf6, dmodel_LD1f6);
sMSE_LD1f6=sqrt(1/1000*sum((ynewf6-Yhat_LD1f6).^2))


[dmodel_LD2f6, perf_LD2f6] = dacefit(LDdataf6_rev,YLDf6_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_LD2f6] = predictor(newdataf6, dmodel_LD2f6);
sMSE_LD2f6=sqrt(1/1000*sum((ynewf6-Yhat_LD2f6).^2))

%L+D quadratic
[betaLDf6]=nlinfit(LDdataf6_rev,YLDf6_rev,quafac2,beta0);
Yhat_LDquaf6=quafac2(betaLDf6,newdataf6);
sMSE_LDquaf6=sqrt(1/1000*sum((ynewf6-Yhat_LDquaf6).^2))

%S+D poly0/poly1/poly2
SDdataf6=[Sprac1f6,Sprac2f6;Dprac1f6,Dprac2f6];
Y_SDf6=[y_S6;y_D6];
[SDdataf6_rev,ISDf6]=unique(SDdataf6,'rows');
YSDf6_rev=Y_SDf6(ISDf6);

[dmodel_SD0f6, perf_SD0f6] = dacefit(SDdataf6_rev,YSDf6_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_SD0f6] = predictor(newdataf6, dmodel_SD0f6);
sMSE_SD0f6=sqrt(1/1000*sum((ynewf6-Yhat_SD0f6).^2))


[dmodel_SD1f6, perf_SD1f6] = dacefit(SDdataf6_rev,YSDf6_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_SD1f6] = predictor(newdataf6, dmodel_SD1f6);
sMSE_SD1f6=sqrt(1/1000*sum((ynewf6-Yhat_SD1f6).^2))


[dmodel_SD2f6, perf_SD2f6] = dacefit(SDdataf6_rev,YSDf6_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_SD2f6] = predictor(newdataf6, dmodel_SD2f6);
sMSE_SD2f6=sqrt(1/1000*sum((ynewf6-Yhat_SD2f6).^2))

%S+D quadratic
[betaSDf6]=nlinfit(SDdataf6_rev,YSDf6_rev,quafac2,beta0);
Yhat_SDquaf6=quafac2(betaSDf6,newdataf6);
sMSE_SDquaf6=sqrt(1/1000*sum((ynewf6-Yhat_SDquaf6).^2))

%U+D poly0/poly1/poly2
UDdataf6=[Uprac1f6,Uprac2f6;Dprac1f6,Dprac2f6];
Y_UDf6=[y_U6;y_D6];
[UDdataf6_rev,IUDf6]=unique(UDdataf6,'rows');
YUDf6_rev=Y_UDf6(IUDf6);

[dmodel_UD0f6, perf_UD0f6] = dacefit(UDdataf6_rev,YUDf6_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_UD0f6] = predictor(newdataf6, dmodel_UD0f6);
sMSE_UD0f6=sqrt(1/1000*sum((ynewf6-Yhat_UD0f6).^2))


[dmodel_UD1f6, perf_UD1f6] = dacefit(UDdataf6_rev,YUDf6_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_UD1f6] = predictor(newdataf6, dmodel_UD1f6);
sMSE_UD1f6=sqrt(1/1000*sum((ynewf6-Yhat_UD1f6).^2))


[dmodel_UD2f6, perf_UD2f6] = dacefit(UDdataf6_rev,YUDf6_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_UD2f6] = predictor(newdataf6, dmodel_UD2f6);
sMSE_UD2f6=sqrt(1/1000*sum((ynewf6-Yhat_UD2f6).^2))

%U+D quadratic
[betaUDf6]=nlinfit(UDdataf6_rev,YUDf6_rev,quafac2,beta0);
Yhat_UDquaf6=quafac2(betaUDf6,newdataf6);
sMSE_UDquaf6=sqrt(1/1000*sum((ynewf6-Yhat_UDquaf6).^2))


%R+D poly0/poly1/poly2
RDdataf6=[Rprac1f6,Rprac2f6;Dprac1f6,Dprac2f6];
Y_RDf6=[y_R6;y_D6];
[RDdataf6_rev,IRDf6]=unique(RDdataf6,'rows');
YRDf6_rev=Y_RDf6(IRDf6);

[dmodel_RD0f6, perf_RD0f6] = dacefit(RDdataf6_rev,YRDf6_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_RD0f6] = predictor(newdataf6, dmodel_RD0f6);
sMSE_RD0f6=sqrt(1/1000*sum((ynewf6-Yhat_RD0f6).^2))


[dmodel_RD1f6, perf_RD1f6] = dacefit(RDdataf6_rev,YRDf6_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_RD1f6] = predictor(newdataf6, dmodel_RD1f6);
sMSE_RD1f6=sqrt(1/1000*sum((ynewf6-Yhat_RD1f6).^2))


[dmodel_RD2f6, perf_RD2f6] = dacefit(RDdataf6_rev,YRDf6_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_RD2f6] = predictor(newdataf6, dmodel_RD2f6);
sMSE_RD2f6=sqrt(1/1000*sum((ynewf6-Yhat_RD2f6).^2))

%R+D quadratic
[betaRDf6]=nlinfit(RDdataf6_rev,YRDf6_rev,quafac2,beta0);
Yhat_RDquaf6=quafac2(betaRDf6,newdataf6);
sMSE_RDquaf6=sqrt(1/1000*sum((ynewf6-Yhat_RDquaf6).^2))

%R+E poly0/poly1/poly2
RLdataf6=[Rprac1f6,Rprac2f6;Lprac1f6,Lprac2f6];
Y_RLf6=[y_R6;y_L6];
[RLdataf6_rev,IRLf6]=unique(RLdataf6,'rows');
YRLf6_rev=Y_RLf6(IRLf6);

[dmodel_RL0f6, perf_RL0f6] = dacefit(RLdataf6_rev,YRLf6_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_RL0f6] = predictor(newdataf6, dmodel_RL0f6);
sMSE_RL0f6=sqrt(1/1000*sum((ynewf6-Yhat_RL0f6).^2))


[dmodel_RL1f6, perf_RL1f6] = dacefit(RLdataf6_rev,YRLf6_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_RL1f6] = predictor(newdataf6, dmodel_RL1f6);
sMSE_RL1f6=sqrt(1/1000*sum((ynewf6-Yhat_RL1f6).^2))


[dmodel_RL2f6, perf_RL2f6] = dacefit(RLdataf6_rev,YRLf6_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_RL2f6] = predictor(newdataf6, dmodel_RL2f6);
sMSE_RL2f6=sqrt(1/1000*sum((ynewf6-Yhat_RL2f6).^2))

%R+E quadratic
[betaRLf6]=nlinfit(RLdataf6_rev,YRLf6_rev,quafac2,beta0);
Yhat_RLquaf6=quafac2(betaRLf6,newdataf6);
sMSE_RLquaf6=sqrt(1/1000*sum((ynewf6-Yhat_RLquaf6).^2))

%R+S poly0/poly1/poly2
RSdataf6=[Rprac1f6,Rprac2f6;Sprac1f6,Sprac2f6];
Y_RSf6=[y_R6;y_S6];
[RSdataf6_rev,IRSf6]=unique(RSdataf6,'rows');
YRSf6_rev=Y_RSf6(IRSf6);

[dmodel_RS0f6, perf_RS0f6] = dacefit(RSdataf6_rev,YRSf6_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_RS0f6] = predictor(newdataf6, dmodel_RS0f6);
sMSE_RS0f6=sqrt(1/1000*sum((ynewf6-Yhat_RS0f6).^2))


[dmodel_RS1f6, perf_RS1f6] = dacefit(RSdataf6_rev,YRSf6_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_RS1f6] = predictor(newdataf6, dmodel_RS1f6);
sMSE_RS1f6=sqrt(1/1000*sum((ynewf6-Yhat_RS1f6).^2))


[dmodel_RS2f6, perf_RS2f6] = dacefit(RSdataf6_rev,YRSf6_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_RS2f6] = predictor(newdataf6, dmodel_RS2f6);
sMSE_RS2f6=sqrt(1/1000*sum((ynewf6-Yhat_RS2f6).^2))

%R+S quadratic
[betaRSf6]=nlinfit(RSdataf6_rev,YRSf6_rev,quafac2,beta0);
Yhat_RSquaf6=quafac2(betaRSf6,newdataf6);
sMSE_RSquaf6=sqrt(1/1000*sum((ynewf6-Yhat_RSquaf6).^2))

%U+S poly0/poly1/poly2
USdataf6=[Uprac1f6,Uprac2f6;Sprac1f6,Sprac2f6];
Y_USf6=[y_U6;y_S6];
[USdataf6_rev,IUSf6]=unique(USdataf6,'rows');
YUSf6_rev=Y_USf6(IUSf6);

[dmodel_US0f6, perf_US0f6] = dacefit(USdataf6_rev,YUSf6_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_US0f6] = predictor(newdataf6, dmodel_US0f6);
sMSE_US0f6=sqrt(1/1000*sum((ynewf6-Yhat_US0f6).^2))

[dmodel_US1f6, perf_US1f6] = dacefit(USdataf6_rev,YUSf6_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_US1f6] = predictor(newdataf6, dmodel_US1f6);
sMSE_US1f6=sqrt(1/1000*sum((ynewf6-Yhat_US1f6).^2))

[dmodel_US2f6, perf_US2f6] = dacefit(USdataf6_rev,YUSf6_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_US2f6] = predictor(newdataf6, dmodel_US2f6);
sMSE_US2f6=sqrt(1/1000*sum((ynewf6-Yhat_US2f6).^2))

%U+S quadratic
[betaUSf6]=nlinfit(USdataf6_rev,YUSf6_rev,quafac2,beta0);
Yhat_USquaf6=quafac2(betaUSf6,newdataf6);
sMSE_USquaf6=sqrt(1/1000*sum((ynewf6-Yhat_USquaf6).^2))

%U+E poly0/poly1/poly2
UEdataf6=[Uprac1f6,Uprac2f6;Lprac1f6,Lprac2f6];
Y_UEf6=[y_U6;y_L6];
[UEdataf6_rev,IUEf]=unique(UEdataf6,'rows');
YUEf6_rev=Y_UEf6(IUEf4);

[dmodel_UE0f6, perf_UE0f6] = dacefit(UEdataf6_rev,YUEf6_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_UE0f6] = predictor(newdataf6, dmodel_UE0f6);
sMSE_UE0f6=sqrt(1/1000*sum((ynewf6-Yhat_UE0f6).^2))


[dmodel_UE1f6, perf_UE1f6] = dacefit(UEdataf6_rev,YUEf6_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_UE1f6] = predictor(newdataf6, dmodel_UE1f6);
sMSE_UE1f6=sqrt(1/1000*sum((ynewf6-Yhat_UE1f6).^2))


[dmodel_UE2f6, perf_UE2f6] = dacefit(UEdataf6_rev,YUEf6_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_UE2f6] = predictor(newdataf6, dmodel_UE2f6);
sMSE_UE2f6=sqrt(1/1000*sum((ynewf6-Yhat_UE2f6).^2))

%U+E quadratic
[betaUEf6]=nlinfit(UEdataf6_rev,YUEf6_rev,quafac2,beta0);
Yhat_UEquaf6=quafac2(betaUEf6,newdataf6);
sMSE_UEquaf6=sqrt(1/1000*sum((ynewf6-Yhat_UEquaf6).^2))

