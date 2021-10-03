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

UL2=[
    2     1
     2     0
     2     2
     2     3
     3     1
     3     0
     3     2
     3     3
     0     1
     0     0
     0     2
     0     3
     1     1
     1     0
     1     2
     1     3];

%L32
L32_2=[0,0;0,0;0,1;0,1;0,2;0,2;0,3;0,3;1,0;1,0;1,1;1,1;1,2;1,2;1,3;1,3;2,0;2,0;2,1;2,1;2,2;2,2;2,3;2,3;3,0;3,0;3,1;3,1;3,2;3,2;3,3;3,3];

%L
Lprac1f2=L2(:,1);
Lprac1f2(Lprac1f2==0)=-3;
Lprac1f2(Lprac1f2==1)=-1;
Lprac1f2(Lprac1f2==2)=1;
Lprac1f2(Lprac1f2==3)=3;
Lprac2f2=L2(:,2);
Lprac2f2(Lprac2f2==0)=-2;
Lprac2f2(Lprac2f2==1)=-2/3;
Lprac2f2(Lprac2f2==2)=2/3;
Lprac2f2(Lprac2f2==3)=2;

%S---> [-9/4,9/4] x [-3/2,3/2]
Sprac1f2=L2(:,1);
Sprac1f2(Sprac1f2==1)=-9/4;
Sprac1f2(Sprac1f2==2)=-9/4+3/2;
Sprac1f2(Sprac1f2==3)=-9/4+6/2;
Sprac1f2(Sprac1f2==4)=9/4;
Sprac2f2=L2(:,2);
Sprac2f2(Sprac2f2==1)=-3/2;
Sprac2f2(Sprac2f2==2)=-3/2+1;
Sprac2f2(Sprac2f2==3)=-3/2+2;
Sprac2f2(Sprac2f2==4)=3/2;

%D32
Dprac1f2=L32_2(:,1);
Dprac1f2(Dprac1f2==0)=-3;
Dprac1f2(Dprac1f2==1)=-3+3*(1-1/sqrt(5));
Dprac1f2(Dprac1f2==2)=-3+3*(1+1/sqrt(5));
Dprac1f2(Dprac1f2==3)=3;
Dprac2f2=L32_2(:,2);
Dprac2f2(Dprac2f2==0)=-2;
Dprac2f2(Dprac2f2==1)=-2+2*(1-1/sqrt(5));
Dprac2f2(Dprac2f2==2)=-2+2*(1+1/sqrt(5));
Dprac2f2(Dprac2f2==3)=2;

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

Uprac1f2=-3+6/15*(U2(:,1)-1);
Uprac2f2=-2+4/15*(U2(:,2)-1);

Rprac1f2=-45/16+3/8*(U2(:,1)-1);
Rprac2f2=-15/8+1/4*(U2(:,2)-1);


y_L2=f2(Lprac1f2,Lprac2f2);
y_S2=f2(Sprac1f2,Sprac2f2);
y_D2=f2(Dprac1f2,Dprac2f2);
y_U2=f2(Uprac1f2,Uprac2f2);
y_R2=f2(Rprac1f2,Rprac2f2);

theta = [10 10]; lob = [1e-1 1e-1]; upb = [20 20];
x1newf2=-3+6*rand(1000,1);
x2newf2=-2+4*rand(1000,1);
newdataf2=[x1newf2,x2newf2];

ynewf2=f2(x1newf2,x2newf2);

%L poly0
[dmodel_L0f2, perf_L0f2] = dacefit([Lprac1f2,Lprac2f2],y_L2, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_L0f2] = predictor(newdataf2, dmodel_L0f2);
MSE_L0f2=1/1000*sum((ynewf2-Yhat_L0f2).^2)

%L poly1
[dmodel_L1f2, perf_L1f2] = dacefit([Lprac1f2,Lprac2f2],y_L2, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_L1f2] = predictor(newdataf2, dmodel_L1f2);
MSE_L1f2=1/1000*sum((ynewf2-Yhat_L1f2).^2)

%L poly2
[dmodel_L2f2, perf_L2f2] = dacefit([Lprac1f2,Lprac2f2],y_L2, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_L2f2] = predictor(newdataf2, dmodel_L2f2);
MSE_L2f2=1/1000*sum((ynewf2-Yhat_L2f2).^2)

%L quadratic
beta0=[1;1;1;1;1];
quafac2=@(b,x)b(1)+b(2)*x(:,1)+b(3)*x(:,1).^2+b(4)*x(:,2)+b(5)*x(:,2).^2;
[betaLf2]=nlinfit([Lprac1f2,Lprac2f2],y_L2,quafac2,beta0);
Yhat_Lquaf2=quafac2(betaLf2,newdataf2);
sMSE_Lquaf2=sqrt(1/1000*sum((ynewf2-Yhat_Lquaf2).^2))


%S poly0
[dmodel_S0f2, perf_S0f2] = dacefit([Sprac1f2,Sprac2f2],y_S2, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_S0f2] = predictor(newdataf2, dmodel_S0f2);
sMSE_S0f2=sqrt(1/1000*sum((ynewf2-Yhat_S0f2).^2))

%S poly1
[dmodel_S1f2, perf_S1f2] = dacefit([Sprac1f2,Sprac2f2],y_S2, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_S1f2] = predictor(newdataf2, dmodel_S1f2);
sMSE_S1f2=sqrt(1/1000*sum((ynewf2-Yhat_S1f2).^2))

%S poly2
[dmodel_S2f2, perf_S2f2] = dacefit([Sprac1f2,Sprac2f2],y_S2, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_S2f2] = predictor(newdataf2, dmodel_S2f2);
sMSE_S2f2=sqrt(1/1000*sum((ynewf2-Yhat_S2f2).^2))

%S quadratic
[betaSf2]=nlinfit([Sprac1f2,Sprac2f2],y_S2,quafac2,beta0);
Yhat_Squaf2=quafac2(betaSf2,newdataf2);
sMSE_Squaf2=sqrt(1/1000*sum((ynewf2-Yhat_Squaf2).^2))


%U poly0
[dmodel_U0f2, perf_U0f2] = dacefit([Uprac1f2,Uprac2f2],y_U2, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_U0f2] = predictor(newdataf2, dmodel_U0f2);
sMSE_U0f2=sqrt(1/1000*sum((ynewf2-Yhat_U0f2).^2))

%U poly1
[dmodel_U1f2, perf_U1f2] = dacefit([Uprac1f2,Uprac2f2],y_U2, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_U1f2] = predictor(newdataf2, dmodel_U1f2);
sMSE_U1f2=sqrt(1/1000*sum((ynewf2-Yhat_U1f2).^2))

%U poly2
[dmodel_U2f2, perf_U2f2] = dacefit([Uprac1f2,Uprac2f2],y_U2, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_U2f2] = predictor(newdataf2, dmodel_U2f2);
sMSE_U2f2=sqrt(1/1000*sum((ynewf2-Yhat_U2f2).^2))

%U quadratic
[betaUf2]=nlinfit([Uprac1f2,Uprac2f2],y_U2,quafac2,beta0);
Yhat_Uquaf2=quafac2(betaUf2,newdataf2);
sMSE_Uquaf2=sqrt(1/1000*sum((ynewf2-Yhat_Uquaf2).^2))

%R poly0
[dmodel_R0f2, perf_R0f2] = dacefit([Rprac1f2,Rprac2f2],y_R2, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_R0f2] = predictor(newdataf2, dmodel_R0f2);
sMSE_R0f2=sqrt(1/1000*sum((ynewf2-Yhat_R0f2).^2))

%R poly1
[dmodel_R1f2, perf_R1f2] = dacefit([Rprac1f2,Rprac2f2],y_R2, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_R1f2] = predictor(newdataf2, dmodel_R1f2);
sMSE_R1f2=sqrt(1/1000*sum((ynewf2-Yhat_R1f2).^2))

%R poly2
[dmodel_R2f2, perf_R2f2] = dacefit([Rprac1f2,Rprac2f2],y_R2, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_R2f2] = predictor(newdataf2, dmodel_R2f2);
sMSE_R2f2=sqrt(1/1000*sum((ynewf2-Yhat_R2f2).^2))

%R quadratic
[betaRf2]=nlinfit([Rprac1f2,Rprac2f2],y_R2,quafac2,beta0);
Yhat_Rquaf2=quafac2(betaRf2,newdataf2);
sMSE_Rquaf2=sqrt(1/1000*sum((ynewf2-Yhat_Rquaf2).^2))

%D poly0
[dmodel_D0f2, perf_D0f2] = dacefit([Dprac1f2,Dprac2f2],y_D2, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_D0f2] = predictor(newdataf2, dmodel_D0f2);
sMSE_D0f2=sqrt(1/1000*sum((ynewf2-Yhat_D0f2).^2))

%D poly1
[dmodel_D1f2, perf_D1f2] = dacefit([Dprac1f2,Dprac2f2],y_D2, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_D1f2] = predictor(newdataf2, dmodel_D1f2);
sMSE_D1f2=sqrt(1/1000*sum((ynewf2-Yhat_D1f2).^2))

%D poly2
[dmodel_D2f2, perf_D2f2] = dacefit([Dprac1f2,Dprac2f2],y_D2, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_D2f2] = predictor(newdataf2, dmodel_D2f2);
sMSE_D2f2=sqrt(1/1000*sum((ynewf2-Yhat_D2f2).^2))

%D quadratic
[betaDf2]=nlinfit([Dprac1f2,Dprac2f2],y_D2,quafac2,beta0);
Yhat_Dquaf2=quafac2(betaDf2,newdataf2);
sMSE_Dquaf2=sqrt(1/1000*sum((ynewf2-Yhat_Dquaf2).^2))

%L+D poly0/poly1/poly2
LDdataf2=[Lprac1f2,Lprac2f2;Dprac1f2,Dprac2f2];
Y_LDf2=[y_L2;y_D2];
[LDdataf2_rev,ILDf2]=unique(LDdataf2,'rows');
YLDf2_rev=Y_LDf2(ILDf2);

[dmodel_LD0f2, perf_LD0f2] = dacefit(LDdataf2_rev,YLDf2_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_LD0f2] = predictor(newdataf2, dmodel_LD0f2);
sMSE_LD0f2=sqrt(1/1000*sum((ynewf2-Yhat_LD0f2).^2))


[dmodel_LD1f2, perf_LD1f2] = dacefit(LDdataf2_rev,YLDf2_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_LD1f2] = predictor(newdataf2, dmodel_LD1f2);
sMSE_LD1f2=sqrt(1/1000*sum((ynewf2-Yhat_LD1f2).^2))


[dmodel_LD2f2, perf_LD2f2] = dacefit(LDdataf2_rev,YLDf2_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_LD2f2] = predictor(newdataf2, dmodel_LD2f2);
sMSE_LD2f2=sqrt(1/1000*sum((ynewf2-Yhat_LD2f2).^2))

%L+D quadratic
[betaLDf2]=nlinfit(LDdataf2_rev,YLDf2_rev,quafac2,beta0);
Yhat_LDquaf2=quafac2(betaLDf2,newdataf2);
sMSE_LDquaf2=sqrt(1/1000*sum((ynewf2-Yhat_LDquaf2).^2))

%S+D poly0/poly1/poly2
SDdataf2=[Sprac1f2,Sprac2f2;Dprac1f2,Dprac2f2];
Y_SDf2=[y_S2;y_D2];
[SDdataf2_rev,ISDf2]=unique(SDdataf2,'rows');
YSDf2_rev=Y_SDf2(ISDf2);

[dmodel_SD0f2, perf_SD0f2] = dacefit(SDdataf2_rev,YSDf2_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_SD0f2] = predictor(newdataf2, dmodel_SD0f2);
sMSE_SD0f2=sqrt(1/1000*sum((ynewf2-Yhat_SD0f2).^2))


[dmodel_SD1f2, perf_SD1f2] = dacefit(SDdataf2_rev,YSDf2_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_SD1f2] = predictor(newdataf2, dmodel_SD1f2);
sMSE_SD1f2=sqrt(1/1000*sum((ynewf2-Yhat_SD1f2).^2))


[dmodel_SD2f2, perf_SD2f2] = dacefit(SDdataf2_rev,YSDf2_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_SD2f2] = predictor(newdataf2, dmodel_SD2f2);
sMSE_SD2f2=sqrt(1/1000*sum((ynewf2-Yhat_SD2f2).^2))

%S+D quadratic
[betaSDf2]=nlinfit(SDdataf2_rev,YSDf2_rev,quafac2,beta0);
Yhat_SDquaf2=quafac2(betaSDf2,newdataf2);
sMSE_SDquaf2=sqrt(1/1000*sum((ynewf2-Yhat_SDquaf2).^2))

%U+D poly0/poly1/poly2
UDdataf2=[Uprac1f2,Uprac2f2;Dprac1f2,Dprac2f2];
Y_UDf2=[y_U2;y_D2];
[UDdataf2_rev,IUDf2]=unique(UDdataf2,'rows');
YUDf2_rev=Y_UDf2(IUDf2);

[dmodel_UD0f2, perf_UD0f2] = dacefit(UDdataf2_rev,YUDf2_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_UD0f2] = predictor(newdataf2, dmodel_UD0f2);
sMSE_UD0f2=sqrt(1/1000*sum((ynewf2-Yhat_UD0f2).^2))


[dmodel_UD1f2, perf_UD1f2] = dacefit(UDdataf2_rev,YUDf2_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_UD1f2] = predictor(newdataf2, dmodel_UD1f2);
sMSE_UD1f2=sqrt(1/1000*sum((ynewf2-Yhat_UD1f2).^2))


[dmodel_UD2f2, perf_UD2f2] = dacefit(UDdataf2_rev,YUDf2_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_UD2f2] = predictor(newdataf2, dmodel_UD2f2);
sMSE_UD2f2=sqrt(1/1000*sum((ynewf2-Yhat_UD2f2).^2))

%U+D quadratic
[betaUDf2]=nlinfit(UDdataf2_rev,YUDf2_rev,quafac2,beta0);
Yhat_UDquaf2=quafac2(betaUDf2,newdataf2);
sMSE_UDquaf2=sqrt(1/1000*sum((ynewf2-Yhat_UDquaf2).^2))


%R+D poly0/poly1/poly2
RDdataf2=[Rprac1f2,Rprac2f2;Dprac1f2,Dprac2f2];
Y_RDf2=[y_R2;y_D2];
[RDdataf2_rev,IRDf2]=unique(RDdataf2,'rows');
YRDf2_rev=Y_RDf2(IRDf2);

[dmodel_RD0f2, perf_RD0f2] = dacefit(RDdataf2_rev,YRDf2_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_RD0f2] = predictor(newdataf2, dmodel_RD0f2);
sMSE_RD0f2=sqrt(1/1000*sum((ynewf2-Yhat_RD0f2).^2))


[dmodel_RD1f2, perf_RD1f2] = dacefit(RDdataf2_rev,YRDf2_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_RD1f2] = predictor(newdataf2, dmodel_RD1f2);
sMSE_RD1f2=sqrt(1/1000*sum((ynewf2-Yhat_RD1f2).^2))


[dmodel_RD2f2, perf_RD2f2] = dacefit(RDdataf2_rev,YRDf2_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_RD2f2] = predictor(newdataf2, dmodel_RD2f2);
sMSE_RD2f2=sqrt(1/1000*sum((ynewf2-Yhat_RD2f2).^2))

%R+D quadratic
[betaRDf2]=nlinfit(RDdataf2_rev,YRDf2_rev,quafac2,beta0);
Yhat_RDquaf2=quafac2(betaRDf2,newdataf2);
sMSE_RDquaf2=sqrt(1/1000*sum((ynewf2-Yhat_RDquaf2).^2))

%R+E poly0/poly1/poly2
RLdataf2=[Rprac1f2,Rprac2f2;Lprac1f2,Lprac2f2];
Y_RLf2=[y_R2;y_L2];
[RLdataf2_rev,IRLf2]=unique(RLdataf2,'rows');
YRLf2_rev=Y_RLf2(IRLf2);

[dmodel_RL0f2, perf_RL0f2] = dacefit(RLdataf2_rev,YRLf2_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_RL0f2] = predictor(newdataf2, dmodel_RL0f2);
sMSE_RL0f2=sqrt(1/1000*sum((ynewf2-Yhat_RL0f2).^2))


[dmodel_RL1f2, perf_RL1f2] = dacefit(RLdataf2_rev,YRLf2_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_RL1f2] = predictor(newdataf2, dmodel_RL1f2);
sMSE_RL1f2=sqrt(1/1000*sum((ynewf2-Yhat_RL1f2).^2))


[dmodel_RL2f2, perf_RL2f2] = dacefit(RLdataf2_rev,YRLf2_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_RL2f2] = predictor(newdataf2, dmodel_RL2f2);
sMSE_RL2f2=sqrt(1/1000*sum((ynewf2-Yhat_RL2f2).^2))

%R+E quadratic
[betaRLf2]=nlinfit(RLdataf2_rev,YRLf2_rev,quafac2,beta0);
Yhat_RLquaf2=quafac2(betaRLf2,newdataf2);
sMSE_RLquaf2=sqrt(1/1000*sum((ynewf2-Yhat_RLquaf2).^2))

%R+S poly0/poly1/poly2
RSdataf2=[Rprac1f2,Rprac2f2;Sprac1f2,Sprac2f2];
Y_RSf2=[y_R2;y_S2];
[RSdataf2_rev,IRSf2]=unique(RSdataf2,'rows');
YRSf2_rev=Y_RSf2(IRSf2);

[dmodel_RS0f2, perf_RS0f2] = dacefit(RSdataf2_rev,YRSf2_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_RS0f2] = predictor(newdataf2, dmodel_RS0f2);
sMSE_RS0f2=sqrt(1/1000*sum((ynewf2-Yhat_RS0f2).^2))


[dmodel_RS1f2, perf_RS1f2] = dacefit(RSdataf2_rev,YRSf2_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_RS1f2] = predictor(newdataf2, dmodel_RS1f2);
sMSE_RS1f2=sqrt(1/1000*sum((ynewf2-Yhat_RS1f2).^2))


[dmodel_RS2f2, perf_RS2f2] = dacefit(RSdataf2_rev,YRSf2_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_RS2f2] = predictor(newdataf2, dmodel_RS2f2);
sMSE_RS2f2=sqrt(1/1000*sum((ynewf2-Yhat_RS2f2).^2))

%R+S quadratic
[betaRSf2]=nlinfit(RSdataf2_rev,YRSf2_rev,quafac2,beta0);
Yhat_RSquaf2=quafac2(betaRSf2,newdataf2);
sMSE_RSquaf2=sqrt(1/1000*sum((ynewf2-Yhat_RSquaf2).^2))

%U+S poly0/poly1/poly2
USdataf2=[Uprac1f2,Uprac2f2;Sprac1f2,Sprac2f2];
Y_USf2=[y_U2;y_S2];
[USdataf2_rev,IUSf2]=unique(USdataf2,'rows');
YUSf2_rev=Y_USf2(IUSf2);
[dmodel_US0f2, perf_US0f2] = dacefit(USdataf2_rev,YUSf2_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_US0f2] = predictor(newdataf2, dmodel_US0f2);
sMSE_US0f2=sqrt(1/1000*sum((ynewf2-Yhat_US0f2).^2))
[dmodel_US1f2, perf_US1f2] = dacefit(USdataf2_rev,YUSf2_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_US1f2] = predictor(newdataf2, dmodel_US1f2);
sMSE_US1f2=sqrt(1/1000*sum((ynewf2-Yhat_US1f2).^2))
[dmodel_US2f2, perf_US2f2] = dacefit(USdataf2_rev,YUSf2_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_US2f2] = predictor(newdataf2, dmodel_US2f2);
sMSE_US2f2=sqrt(1/1000*sum((ynewf2-Yhat_US2f2).^2))
%U+S quadratic
[betaUSf2]=nlinfit(USdataf2_rev,YUSf2_rev,quafac2,beta0);
Yhat_USquaf2=quafac2(betaUSf2,newdataf2);
sMSE_USquaf2=sqrt(1/1000*sum((ynewf2-Yhat_USquaf2).^2))

%U+E poly0/poly1/poly2
UEdataf2=[Uprac1f2,Uprac2f2;Lprac1f2,Lprac2f2];
Y_UEf2=[y_U2;y_L2];
[UEdataf2_rev,IUEf2]=unique(UEdataf2,'rows');
YUEf2_rev=Y_UEf2(IUEf2);

[dmodel_UE0f2, perf_UE0f2] = dacefit(UEdataf2_rev,YUEf2_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_UE0f2] = predictor(newdataf2, dmodel_UE0f2);
sMSE_UE0f2=sqrt(1/1000*sum((ynewf2-Yhat_UE0f2).^2))


[dmodel_UE1f2, perf_UE1f2] = dacefit(UEdataf2_rev,YUEf2_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_UE1f2] = predictor(newdataf2, dmodel_UE1f2);
sMSE_UE1f2=sqrt(1/1000*sum((ynewf2-Yhat_UE1f2).^2))


[dmodel_UE2f2, perf_UE2f2] = dacefit(UEdataf2_rev,YUEf2_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_UE2f2] = predictor(newdataf2, dmodel_UE2f2);
sMSE_UE2f2=sqrt(1/1000*sum((ynewf2-Yhat_UE2f2).^2))

%U+E quadratic
[betaUEf2]=nlinfit(UEdataf2_rev,YUEf2_rev,quafac2,beta0);
Yhat_UEquaf2=quafac2(betaUEf2,newdataf2);
sMSE_UEquaf2=sqrt(1/1000*sum((ynewf2-Yhat_UEquaf2).^2))

