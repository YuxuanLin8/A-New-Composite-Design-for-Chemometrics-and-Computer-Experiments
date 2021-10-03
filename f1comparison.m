L3=[0 0 0
0 1 1
0 2 2
0 3 3
1 0 1
1 1 0
1 2 3
1 3 2
2 0 2
2 1 3
2 2 0
2 3 1
3 0 3
3 1 2
3 2 1
3 3 0];

%Uniform L3
UL3=[
     0     1     2
     0     2     0
     0     3     3
     0     0     1
     1     1     0
     1     2     2
     1     3     1
     1     0     3
     2     1     3
     2     2     1
     2     3     2
     2     0     0
     3     1     1
     3     2     3
     3     3     0
     3     0     2];

%L32_3
L32_3=[
    0 0 0
0 0 1
0 1 0
0 1 2
0 2 1
0 2 3
0 3 2
0 3 3
1 0 0
1 0 2
1 1 2
1 1 3
1 2 0
1 2 1
1 3 1
1 3 3
2 0 1
2 0 3
2 1 0
2 1 1
2 2 2
2 2 3
2 3 0
2 3 2
3 0 2
3 0 3
3 1 1
3 1 3
3 2 0
3 2 2
3 3 0
3 3 1];

%L
Lprac1=L3;
Lprac1(Lprac1==0)=-2;
Lprac1(Lprac1==1)=-2/3;
Lprac1(Lprac1==2)=2/3;
Lprac1(Lprac1==3)=2;
Lprac2=L3;
Lprac2(Lprac2==0)=-1;
Lprac2(Lprac2==1)=2/3;
Lprac2(Lprac2==2)=7/3;
Lprac2(Lprac2==3)=4;

%L32
Lprac1=L32_3;
Lprac1(Lprac1==0)=-2;
Lprac1(Lprac1==1)=-2/3;
Lprac1(Lprac1==2)=2/3;
Lprac1(Lprac1==3)=2;
Lprac2=L32_3;
Lprac2(Lprac2==0)=-1;
Lprac2(Lprac2==1)=2/3;
Lprac2(Lprac2==2)=7/3;
Lprac2(Lprac2==3)=4;

%S--->[-3/2,3/2]^2 x [-3/8,27/8]
Sprac1=L3;
Sprac1(Sprac1==0)=-3/2;
Sprac1(Sprac1==1)=-3/2+1;
Sprac1(Sprac1==2)=-3/2+2;
Sprac1(Sprac1==3)=3/2;
Sprac2=L3;
Sprac2(Sprac2==0)=-3/8;
Sprac2(Sprac2==1)=-3/8+5/4;
Sprac2(Sprac2==2)=-3/8+5/2;
Sprac2(Sprac2==3)=27/8;


%D
Dprac1=L3;
Dprac1(Dprac1==0)=-2;
Dprac1(Dprac1==1)=-2+2*(1-1/sqrt(5));
Dprac1(Dprac1==2)=-2+2*(1+1/sqrt(5));
Dprac1(Dprac1==3)=2;
Dprac2=L3;
Dprac2(Dprac2==0)=-1;
Dprac2(Dprac2==1)=-1+2.5*(1-1/sqrt(5));
Dprac2(Dprac2==2)=-1+2.5*(1+1/sqrt(5));
Dprac2(Dprac2==3)=4;

%D32
Dprac1=L32_3;
Dprac1(Dprac1==0)=-2;
Dprac1(Dprac1==1)=-2+2*(1-1/sqrt(5));
Dprac1(Dprac1==2)=-2+2*(1+1/sqrt(5));
Dprac1(Dprac1==3)=2;
Dprac2=L32_3;
Dprac2(Dprac2==0)=-1;
Dprac2(Dprac2==1)=-1+2.5*(1-1/sqrt(5));
Dprac2(Dprac2==2)=-1+2.5*(1+1/sqrt(5));
Dprac2(Dprac2==3)=4;

U3=[
1	12	7
2	6	12
3	4	4
4	14	15
5	9	1
6	1	10
7	16	5
8	8	16
9	5	8
10	11	13
11	2	2
12	15	11
13	7	6
14	13	3
15	3	14
16	10	9
];

%U32_3
U32_3=[
    ];

Uprac1=-2+4/15*(U3-1);
Uprac2=-1+1/3*(U3-1);

Rprac1=-15/8+1/4*(U3-1);
Rprac2=-27/32+(5/16)*(U3-1);



y_L1=f1(Lprac1(:,1),Lprac1(:,2),Lprac2(:,3));
y_S1=f1(Sprac1(:,1),Sprac1(:,2),Sprac2(:,3));
y_D1=f1(Dprac1(:,1),Dprac1(:,2),Dprac2(:,3));
y_U1=f1(Uprac1(:,1),Uprac1(:,2),Uprac2(:,3));
y_R1=f1(Rprac1(:,1),Rprac1(:,2),Rprac2(:,3));

theta = [10 10 10]; lob = [1e-1 1e-1 1e-1]; upb = [20 20 20];
x12newf1=-2+4*rand(1000,2);
x3newf1=-1+rand(1000,1)*5;
newdataf1=[x12newf1,x3newf1];

ynewf1=f1(x12newf1(:,1),x12newf1(:,2),x3newf1);
%L poly0
[dmodel_L0f1, perf_L0f1] = dacefit([Lprac1(:,1:2),Lprac2(:,3)],y_L1, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_L0f1] = predictor(newdataf1, dmodel_L0f1);
MSE_L0f1=1/1000*sum((ynewf1-Yhat_L0f1).^2)
%sMSE_L0f1=sqrt(MSE_L0f1)

%L poly1
[dmodel_L1f1, perf_L1f1] = dacefit([Lprac1(:,1:2),Lprac2(:,3)],y_L1, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_L1f1] = predictor(newdataf1, dmodel_L1f1);
MSE_L1f1=1/1000*sum((ynewf1-Yhat_L1f1).^2)
%sMSE_L1f1=sqrt(MSE_L1f1)


%L poly2
[dmodel_L2f1, perf_L2f1] = dacefit([Lprac1(:,1:2),Lprac2(:,3)],y_L1, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_L2f1] = predictor(newdataf1, dmodel_L2f1);
MSE_L2f1=1/1000*sum((ynewf1-Yhat_L2f1).^2)
%sMSE_L2f1=sqrt(MSE_L2f1)

%L quadratic
% beta0=[1;1;1;1;1;1;1];
% quafac3=@(b,x)b(1)+b(2)*x(:,1)+b(3)*x(:,1).^2+b(4)*x(:,2)+b(5)*x(:,2).^2+b(6)*x(:,3)+b(7)*x(:,3).^2;
% [betaL1]=nlinfit([Lprac1(:,1:2),Lprac2(:,3)],y_L1,quafac3,beta0);
% Yhat_Lquaf1=quafac3(betaL1,newdataf1);
% sMSE_Lquaf1=sqrt(1/1000*sum((ynewf1-Yhat_Lquaf1).^2))


%S poly0
[dmodel_S0f1, perf_S0f1] = dacefit([Sprac1(:,1:2),Sprac2(:,3)],y_S1, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_S0f1] = predictor(newdataf1, dmodel_S0f1);
sMSE_S0f1=sqrt(1/1000*sum((ynewf1-Yhat_S0f1).^2))

%S poly1
[dmodel_S1f1, perf_S1f1] = dacefit([Sprac1(:,1:2),Sprac2(:,3)],y_S1, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_S1f1] = predictor(newdataf1, dmodel_S1f1);
sMSE_S1f1=sqrt(1/1000*sum((ynewf1-Yhat_S1f1).^2))

%S poly2
[dmodel_S2f1, perf_S2f1] = dacefit([Sprac1(:,1:2),Sprac2(:,3)],y_S1, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_S2f1] = predictor(newdataf1, dmodel_S2f1);
sMSE_S2f1=sqrt(1/1000*sum((ynewf1-Yhat_S2f1).^2))

%S quadratic
% [betaS1]=nlinfit([Sprac1(:,1:2),Sprac2(:,3)],y_S1,quafac3,beta0);
% Yhat_Squaf1=quafac3(betaS1,newdataf1);
% sMSE_Squaf1=sqrt(1/1000*sum((ynewf1-Yhat_Squaf1).^2))


%U poly0
[dmodel_U0f1, perf_U0f1] = dacefit([Uprac1(:,1:2),Uprac2(:,3)],y_U1, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_U0f1] = predictor(newdataf1, dmodel_U0f1);
MSE_U0f1=1/1000*sum((ynewf1-Yhat_U0f1).^2);
sMSE_U0f1=sqrt(MSE_U0f1)


%U poly1
[dmodel_U1f1, perf_U1f1] = dacefit([Uprac1(:,1:2),Uprac2(:,3)],y_U1, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_U1f1] = predictor(newdataf1, dmodel_U1f1);
MSE_U1f1=1/1000*sum((ynewf1-Yhat_U1f1).^2);
sMSE_U1f1=sqrt(MSE_U1f1)

%U poly2
[dmodel_U2f1, perf_U2f1] = dacefit([Uprac1(:,1:2),Uprac2(:,3)],y_U1, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_U2f1] = predictor(newdataf1, dmodel_U2f1);
MSE_U2f1=1/1000*sum((ynewf1-Yhat_U2f1).^2);
sMSE_U2f1=sqrt(MSE_U2f1)

%U quadratic
% [betaU1]=nlinfit([Uprac1(:,1:2),Uprac2(:,3)],y_U1,quafac3,beta0);
% Yhat_Uquaf1=quafac3(betaU1,newdataf1);
% sMSE_Uquaf1=sqrt(1/1000*sum((ynewf1-Yhat_Uquaf1).^2))

%R poly0
[dmodel_R0f1, perf_R0f1] = dacefit([Rprac1(:,1:2),Rprac2(:,3)],y_R1, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_R0f1] = predictor(newdataf1, dmodel_R0f1);
sMSE_R0f1=sqrt(1/1000*sum((ynewf1-Yhat_R0f1).^2))

%R poly1
[dmodel_R1f1, perf_R1f1] = dacefit([Rprac1(:,1:2),Rprac2(:,3)],y_R1, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_R1f1] = predictor(newdataf1, dmodel_R1f1);
sMSE_R1f1=sqrt(1/1000*sum((ynewf1-Yhat_R1f1).^2))

%R poly2
[dmodel_R2f1, perf_R2f1] = dacefit([Rprac1(:,1:2),Rprac2(:,3)],y_R1, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_R2f1] = predictor(newdataf1, dmodel_R2f1);
sMSE_R2f1=sqrt(1/1000*sum((ynewf1-Yhat_R2f1).^2))

%R quadratic
% [betaR1]=nlinfit([Rprac1(:,1:2),Rprac2(:,3)],y_R1,quafac3,beta0);
% Yhat_Rquaf1=quafac3(betaR1,newdataf1);
% sMSE_Rquaf1=sqrt(1/1000*sum((ynewf1-Yhat_Rquaf1).^2))

%D poly0
[dmodel_D0f1, perf_D0f1] = dacefit([Dprac1(:,1:2),Dprac2(:,3)],y_D1, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_D0f1] = predictor(newdataf1, dmodel_D0f1);
MSE_D0f1=1/1000*sum((ynewf1-Yhat_D0f1).^2)
%sMSE_D0f1=sqrt(MSE_D0f1)

%D poly1
[dmodel_D1f1, perf_D1f1] = dacefit([Dprac1(:,1:2),Dprac2(:,3)],y_D1, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_D1f1] = predictor(newdataf1, dmodel_D1f1);
MSE_D1f1=1/1000*sum((ynewf1-Yhat_D1f1).^2)
%sMSE_D1f1=sqrt(MSE_D1f1)

%D poly2
[dmodel_D2f1, perf_D2f1] = dacefit([Dprac1(:,1:2),Dprac2(:,3)],y_D1, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_D2f1] = predictor(newdataf1, dmodel_D2f1);
MSE_D2f1=1/1000*sum((ynewf1-Yhat_D2f1).^2)
%sMSE_D2f1=sqrt(MSE_D2f1)

%D quadratic
% [betaD1]=nlinfit([Dprac1(:,1:2),Dprac2(:,3)],y_D1,quafac3,beta0);
% Yhat_Dquaf1=quafac3(betaD1,newdataf1);
% sMSE_Dquaf1=sqrt(1/1000*sum((ynewf1-Yhat_Dquaf1).^2))

%L+D poly0/poly1/poly2
LDdataf1=[Lprac1(:,1:2),Lprac2(:,3);Dprac1(:,1:2),Dprac2(:,3)];
Y_L_D1=[y_L1;y_D1];
[LDdataf1_rev,ILDf1]=unique(LDdataf1,'rows');
YLDf1_rev=Y_L_D1(ILDf1);
[dmodel_LD0f1, perf_LD0f1] = dacefit(LDdataf1_rev,YLDf1_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_LD0f1] = predictor(newdataf1, dmodel_LD0f1);
MSE_LD0f1=1/1000*sum((ynewf1-Yhat_LD0f1).^2);
sMSE_LD0f1=sqrt(MSE_LD0f1)


[dmodel_LD1f1, perf_LD1f1] = dacefit(LDdataf1_rev,YLDf1_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_LD1f1] = predictor(newdataf1, dmodel_LD1f1);
MSE_LD1f1=1/1000*sum((ynewf1-Yhat_LD1f1).^2);
sMSE_LD1f1=sqrt(MSE_LD1f1)


[dmodel_LD2f1, perf_LD2f1] = dacefit(LDdataf1_rev,YLDf1_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_LD2f1] = predictor(newdataf1, dmodel_LD2f1);
MSE_LD2f1=1/1000*sum((ynewf1-Yhat_LD2f1).^2);
sMSE_LD2f1=sqrt(MSE_LD2f1)

%L+D quadratic
% [betaLDf1]=nlinfit(LDdataf1_rev,YLDf1_rev,quafac3,beta0);
% Yhat_LDquaf1=quafac3(betaLDf1,newdataf1);
% sMSE_LDquaf1=sqrt(1/1000*sum((ynewf1-Yhat_LDquaf1).^2))

%S+D poly0/poly1/poly2
SDdataf1=[Sprac1(:,1:2),Sprac2(:,3);Dprac1(:,1:2),Dprac2(:,3)];
Y_S_D1=[y_S1;y_D1];
[SDdataf1_rev,ISDf1]=unique(SDdataf1,'rows');
YSDf1_rev=Y_S_D1(ISDf1);
[dmodel_SD0f1, perf_SD0f1] = dacefit(SDdataf1_rev,YSDf1_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_SD0f1] = predictor(newdataf1, dmodel_SD0f1);
MSE_SD0f1=1/1000*sum((ynewf1-Yhat_SD0f1).^2);
sMSE_SD0f1=sqrt(MSE_SD0f1)

[dmodel_SD1f1, perf_SD1f1] = dacefit(SDdataf1_rev,YSDf1_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_SD1f1] = predictor(newdataf1, dmodel_SD1f1);
MSE_SD1f1=1/1000*sum((ynewf1-Yhat_SD1f1).^2);
sMSE_SD1f1=sqrt(MSE_SD1f1)

[dmodel_SD2f1, perf_SD2f1] = dacefit(SDdataf1_rev,YSDf1_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_SD2f1] = predictor(newdataf1, dmodel_SD2f1);
MSE_SD2f1=1/1000*sum((ynewf1-Yhat_SD2f1).^2);
sMSE_SD2f1=sqrt(MSE_SD2f1)

%S+D quadratic
% [betaSDf1]=nlinfit(SDdataf1_rev,YSDf1_rev,quafac3,beta0);
% Yhat_SDquaf1=quafac3(betaSDf1,newdataf1);
% sMSE_SDquaf1=sqrt(1/1000*sum((ynewf1-Yhat_SDquaf1).^2))

%U+D poly0/poly1/poly2
UDdataf1=[Uprac1(:,1:2),Uprac2(:,3);Dprac1(:,1:2),Dprac2(:,3)];
Y_U_D1=[y_U1;y_D1];
[UDdataf1_rev,IUDf1]=unique(UDdataf1,'rows');
YUDf1_rev=Y_U_D1(IUDf1);
[dmodel_UD0f1, perf_UD0f1] = dacefit(UDdataf1_rev,YUDf1_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_UD0f1] = predictor(newdataf1, dmodel_UD0f1);
MSE_UD0f1=1/1000*sum((ynewf1-Yhat_UD0f1).^2);
sMSE_UD0f1=sqrt(MSE_UD0f1)


[dmodel_UD1f1, perf_UD1f1] = dacefit(UDdataf1_rev,YUDf1_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_UD1f1] = predictor(newdataf1, dmodel_UD1f1);
MSE_UD1f1=1/1000*sum((ynewf1-Yhat_UD1f1).^2);
sMSE_UD1f1=sqrt(MSE_UD1f1)

[dmodel_UD2f1, perf_UD2f1] = dacefit(UDdataf1_rev,YUDf1_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_UD2f1] = predictor(newdataf1, dmodel_UD2f1);
MSE_UD2f1=1/1000*sum((ynewf1-Yhat_UD2f1).^2);
sMSE_UD2f1=sqrt(MSE_UD2f1)

%U+D quadratic
% [betaUDf1]=nlinfit(UDdataf1_rev,YUDf1_rev,quafac3,beta0);
% Yhat_UDquaf1=quafac3(betaUDf1,newdataf1);
% sMSE_UDquaf1=sqrt(1/1000*sum((ynewf1-Yhat_UDquaf1).^2))

%R+D poly0/poly1/poly2
RDdataf1=[Rprac1(:,1:2),Rprac2(:,3);Dprac1(:,1:2),Dprac2(:,3)];
Y_R_D1=[y_R1;y_D1];
[RDdataf1_rev,IRDf1]=unique(RDdataf1,'rows');
YRDf1_rev=Y_R_D1(IRDf1);
[dmodel_RD0f1, perf_RD0f1] = dacefit(RDdataf1_rev,YRDf1_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_RD0f1] = predictor(newdataf1, dmodel_RD0f1);
MSE_RD0f1=1/1000*sum((ynewf1-Yhat_RD0f1).^2);
sMSE_RD0f1=sqrt(MSE_RD0f1)

[dmodel_RD1f1, perf_RD1f1] = dacefit(RDdataf1_rev,YRDf1_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_RD1f1] = predictor(newdataf1, dmodel_RD1f1);
sMSE_RD1f1=sqrt(1/1000*sum((ynewf1-Yhat_RD1f1).^2))

[dmodel_RD2f1, perf_RD2f1] = dacefit(RDdataf1_rev,YRDf1_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_RD2f1] = predictor(newdataf1, dmodel_RD2f1);
sMSE_RD2f1=sqrt(1/1000*sum((ynewf1-Yhat_RD2f1).^2))

%R+D quadratic
% [betaRDf1]=nlinfit(RDdataf1_rev,YRDf1_rev,quafac3,beta0);
% Yhat_RDquaf1=quafac3(betaRDf1,newdataf1);
% sMSE_RDquaf1=sqrt(1/1000*sum((ynewf1-Yhat_RDquaf1).^2))


%L+U poly0/poly1/poly2
LUdataf1=[Lprac1(:,1:2),Lprac2(:,3);Uprac1(:,1:2),Uprac2(:,3)];
Y_L_U1=[y_L1;y_U1];
[LUdataf1_rev,ILUf1]=unique(LUdataf1,'rows');
YLUf1_rev=Y_L_U1(ILUf1);
[dmodel_LU0f1, perf_LU0f1] = dacefit(LUdataf1_rev,YLUf1_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_LU0f1] = predictor(newdataf1, dmodel_LU0f1);
MSE_LU0f1=1/1000*sum((ynewf1-Yhat_LU0f1).^2);
sMSE_LU0f1=sqrt(MSE_LU0f1)

[dmodel_LU1f1, perf_LU1f1] = dacefit(LUdataf1_rev,YLUf1_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_LU1f1] = predictor(newdataf1, dmodel_LU1f1);
MSE_LU1f1=1/1000*sum((ynewf1-Yhat_LU1f1).^2);
sMSE_LU1f1=sqrt(MSE_LU1f1)

[dmodel_LU2f1, perf_LU2f1] = dacefit(LUdataf1_rev,YLUf1_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_LU2f1] = predictor(newdataf1, dmodel_LU2f1);
MSE_LU2f1=1/1000*sum((ynewf1-Yhat_LU2f1).^2);
sMSE_LU2f1=sqrt(MSE_LU2f1)

%L+U quadratic
% [betaLUf1]=nlinfit(LUdataf1_rev,YLUf1_rev,quafac3,beta0);
% Yhat_LUquaf1=quafac3(betaLUf1,newdataf1);
% sMSE_LUquaf1=sqrt(1/1000*sum((ynewf1-Yhat_LUquaf1).^2))


%L+R poly0/poly1/poly2
LRdataf1=[Lprac1(:,1:2),Lprac2(:,3);Rprac1(:,1:2),Rprac2(:,3)];
Y_L_R1=[y_L1;y_R1];
[LRdataf1_rev,ILRf1]=unique(LRdataf1,'rows');
YLRf1_rev=Y_L_R1(ILRf1);
[dmodel_LR0f1, perf_LR0f1] = dacefit(LRdataf1_rev,YLRf1_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_LR0f1] = predictor(newdataf1, dmodel_LR0f1);
sMSE_LR0f1=sqrt(1/1000*sum((ynewf1-Yhat_LR0f1).^2))

[dmodel_LR1f1, perf_LR1f1] = dacefit(LRdataf1_rev,YLRf1_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_LR1f1] = predictor(newdataf1, dmodel_LR1f1);
sMSE_LR1f1=sqrt(1/1000*sum((ynewf1-Yhat_LR1f1).^2))

[dmodel_LR2f1, perf_LR2f1] = dacefit(LRdataf1_rev,YLRf1_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_LR2f1] = predictor(newdataf1, dmodel_LR2f1);
sMSE_LR2f1=sqrt(1/1000*sum((ynewf1-Yhat_LR2f1).^2))

%L+R quadratic
% [betaLRf1]=nlinfit(LRdataf1_rev,YLRf1_rev,quafac3,beta0);
% Yhat_LRquaf1=quafac3(betaLRf1,newdataf1);
% sMSE_LRquaf1=sqrt(1/1000*sum((ynewf1-Yhat_LRquaf1).^2))

%S+R poly0/poly1/poly2
SRdataf1=[Sprac1(:,1:2),Sprac2(:,3);Rprac1(:,1:2),Rprac2(:,3)];
Y_S_R1=[y_S1;y_R1];
[SRdataf1_rev,ISRf1]=unique(SRdataf1,'rows');
YSRf1_rev=Y_S_R1(ISRf1);
[dmodel_SR0f1, perf_SR0f1] = dacefit(SRdataf1_rev,YSRf1_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_SR0f1] = predictor(newdataf1, dmodel_SR0f1);
sMSE_SR0f1=sqrt(1/1000*sum((ynewf1-Yhat_SR0f1).^2))

[dmodel_SR1f1, perf_SR1f1] = dacefit(SRdataf1_rev,YSRf1_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_SR1f1] = predictor(newdataf1, dmodel_SR1f1);
sMSE_SR1f1=sqrt(1/1000*sum((ynewf1-Yhat_SR1f1).^2))

[dmodel_SR2f1, perf_SR2f1] = dacefit(SRdataf1_rev,YSRf1_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_SR2f1] = predictor(newdataf1, dmodel_SR2f1);
sMSE_SR2f1=sqrt(1/1000*sum((ynewf1-Yhat_SR2f1).^2))

%S+R quadratic
% [betaSRf1]=nlinfit(SRdataf1_rev,YSRf1_rev,quafac3,beta0);
% Yhat_SRquaf1=quafac3(betaSRf1,newdataf1);
% sMSE_SRquaf1=sqrt(1/1000*sum((ynewf1-Yhat_SRquaf1).^2))

%S+U poly0/poly1/poly2
SUdataf1=[Sprac1(:,1:2),Sprac2(:,3);Uprac1(:,1:2),Uprac2(:,3)];
Y_S_U1=[y_S1;y_U1];
[SUdataf1_rev,ISUf1]=unique(SUdataf1,'rows');
YSUf1_rev=Y_S_U1(ISUf1);
[dmodel_SU0f1, perf_SU0f1] = dacefit(SUdataf1_rev,YSUf1_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_SU0f1] = predictor(newdataf1, dmodel_SU0f1);
sMSE_SU0f1=sqrt(1/1000*sum((ynewf1-Yhat_SU0f1).^2))

[dmodel_SU1f1, perf_SU1f1] = dacefit(SUdataf1_rev,YSUf1_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_SU1f1] = predictor(newdataf1, dmodel_SU1f1);
sMSE_SU1f1=sqrt(1/1000*sum((ynewf1-Yhat_SU1f1).^2))

[dmodel_SU2f1, perf_SU2f1] = dacefit(SUdataf1_rev,YSUf1_rev, @regpoly2, @corrgauss, theta, lob, upb);
[Yhat_SU2f1] = predictor(newdataf1, dmodel_SU2f1);
sMSE_SU2f1=sqrt(1/1000*sum((ynewf1-Yhat_SU2f1).^2))


%S+U quadratic
% [betaSUf1]=nlinfit(SUdataf1_rev,YSUf1_rev,quafac3,beta0);
% Yhat_SUquaf1=quafac3(betaSUf1,newdataf1);
% sMSE_SUquaf1=sqrt(1/1000*sum((ynewf1-Yhat_SUquaf1).^2))

%optima estimation
%funLf1=@(x)quafac3(betaL1,x);
%[xLf1,Lf1min]=fminsearch(funLf1,[-1,1,-1])


%Empirical estimation:
%template
%[yminf1,imin]=min(Yhat_R2f1)
%xminf1=[x12newf1(imin,:),x3newf1(imin)]
%e=-7-yminf1
