U4=[
1	10	4	6
2	4	13	15
3	13	10	10
4	8	7	1
5	6	1	12
6	15	15	4
7	1	11	7
8	16	8	14
9	3	3	3
10	7	16	9
11	11	5	16
12	12	12	2
13	2	6	11
14	14	2	8
15	9	14	13
16	5	9	5
];

    U32=[
    30     7     8    20
    17     8    32    19
    31    17    30    12
     5    13     6     2
    12     5    11    13
    13     1    22    23
     2    12    27    10
    28    28     3    21
    15    29     7     9
     6    23    31    24
     8    16    16     7
    21    25    26    29
    22     6    15    16
    19    19     1    11
    18     3    25     3
    11    27    29     6
    29    30    13     4
    26    21    24    14
    14    32    19    30
     3     2    14    27
    16    18    12    22
     1    26     9    15
     7    22     5    25
    25    20    10    31
    20    24    18     1
    27    11    21     8
     9    31    23    18
    32    14    17    28
    10     9     2    32
    23     4     4     5
     4    15    20    17
    24    10    28    26];
Upracf5=-2+4/31*(U32-1);

L4=[0 0 0 0
0 1 1 1
0 2 2 2
0 3 3 3
1 0 1 2
1 1 0 3
1 2 3 0
1 3 2 1
2 0 2 3
2 1 3 2
2 2 0 1
2 3 1 0
3 0 3 1
3 1 2 0
3 2 1 3
3 3 0 2];

%L32
L32_7=[0 0 0 0 0 0 0
0 0 1 1 1 1 1
0 1 0 1 2 2 2
0 1 1 0 3 3 3
0 2 2 2 0 1 2
0 2 3 3 1 0 3
0 3 2 3 2 3 0
0 3 3 2 3 2 1
1 0 0 2 1 2 3
1 0 1 3 0 3 2
1 1 0 3 3 0 1
1 1 1 2 2 1 0
1 2 2 0 1 3 1
1 2 3 1 0 2 0
1 3 2 1 3 1 3
1 3 3 0 2 0 2
2 0 2 0 3 2 2
2 0 3 1 2 3 3
2 1 2 1 1 0 0
2 1 3 0 0 1 1
2 2 0 2 3 3 0
2 2 1 3 2 2 1
2 3 0 3 1 1 2
2 3 1 2 0 0 3
3 0 2 2 2 0 1
3 0 3 3 3 1 0
3 1 2 3 0 2 3
3 1 3 2 1 3 2
3 2 0 0 2 1 3
3 2 1 1 3 0 2
3 3 0 1 0 3 1
3 3 1 0 1 2 0];
L4=L32_7(:,1:4);
%f5
Lpracf5=L4;
Lpracf5(Lpracf5==0)=-2;
Lpracf5(Lpracf5==1)=-2/3;
Lpracf5(Lpracf5==2)=2/3;
Lpracf5(Lpracf5==3)=2;

%S---> [-3/2,3/2] x [-3/2,3/2]
Spracf5=L4;
Spracf5(Spracf5==1)=-3/2;
Spracf5(Spracf5==2)=-3/2+1;
Spracf5(Spracf5==3)=-3/2+2;
Spracf5(Spracf5==4)=3/2;



%D
Dpracf5=L4;
Dpracf5(Dpracf5==0)=-2;
Dpracf5(Dpracf5==1)=-2+2*(1-1/sqrt(5));
Dpracf5(Dpracf5==2)=-2+2*(1+1/sqrt(5));
Dpracf5(Dpracf5==3)=2;

Upracf5=-2+4/15*(U4-1);

Rpracf5=-15/8+1/4*(U4-1);

y_L5=f5(Lpracf5(:,1),Lpracf5(:,2),Lpracf5(:,3),Lpracf5(:,4));
y_S5=f5(Spracf5(:,1),Spracf5(:,2),Spracf5(:,3),Spracf5(:,4));
y_D5=f5(Dpracf5(:,1),Dpracf5(:,2),Dpracf5(:,3),Dpracf5(:,4));
y_U5=f5(Upracf5(:,1),Upracf5(:,2),Upracf5(:,3),Upracf5(:,4));
y_R5=f5(Rpracf5(:,1),Rpracf5(:,2),Rpracf5(:,3),Rpracf5(:,4));

theta = [10 10 10 10]; lob = [1e-1 1e-1 1e-1 1e-1]; upb = [20 20 20 20];
xnewf5=-2+4*rand(1000,4);
ynewf5=f5(xnewf5(:,1),xnewf5(:,2),xnewf5(:,3),xnewf5(:,4));

%L poly0
[dmodel_L0f5, perf_L0f5] = dacefit(Lpracf5,y_L5, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_L0f5] = predictor(xnewf5, dmodel_L0f5);
sMSE_L0f5=sqrt(1/1000*sum((ynewf5-Yhat_L0f5).^2))

%L poly1
[dmodel_L1f5, perf_L1f5] = dacefit(Lpracf5,y_L5, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_L1f5] = predictor(xnewf5, dmodel_L1f5);
sMSE_L1f5=sqrt(1/1000*sum((ynewf5-Yhat_L1f5).^2))

%L poly2
%[dmodel_L2f5, perf_L2f5] = dacefit(Lpracf5,y_L5, @regpoly2, @corrgauss, theta, lob, upb);
%[Yhat_L2f5] = predictor(xnewf5, dmodel_L2f5);
%sMSE_L2f5=sqrt(1/1000*sum((ynewf5-Yhat_L2f5).^2))

%L quadratic
beta0fac4=[1;1;1;1;1;1;1;1;1];
quafac4=@(b,x)b(1)+b(2)*x(:,1)+b(3)*x(:,1).^2+b(4)*x(:,2)+b(5)*x(:,2).^2+b(6)*x(:,3)+b(7)*x(:,3).^2+b(8)*x(:,4)+b(9)*x(:,4).^2;
[betaLf5]=nlinfit(Lpracf5,y_L5,quafac4,beta0fac4);
Yhat_Lquaf5=quafac4(betaLf5,xnewf5);
sMSE_Lquaf5=sqrt(1/1000*sum((ynewf5-Yhat_Lquaf5).^2))



%S poly0
[dmodel_S0f5, perf_S0f5] = dacefit(Spracf5,y_S5, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_S0f5] = predictor(xnewf5, dmodel_S0f5);
sMSE_S0f5=sqrt(1/1000*sum((ynewf5-Yhat_S0f5).^2))

%S poly1
[dmodel_S1f5, perf_S1f5] = dacefit(Spracf5,y_S5, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_S1f5] = predictor(xnewf5, dmodel_S1f5);
sMSE_S1f5=sqrt(1/1000*sum((ynewf5-Yhat_S1f5).^2))

%S poly2
%[dmodel_S2f5, perf_S2f5] = dacefit(Spracf5,y_S5, @regpoly2, @corrgauss, theta, lob, upb);
%[Yhat_S2f5] = predictor(xnewf5, dmodel_S2f5);
%sMSE_S2f5=sqrt(1/1000*sum((ynewf5-Yhat_S2f5).^2))

%S quadratic
[betaSf5]=nlinfit(Spracf5,y_S5,quafac4,beta0fac4);
Yhat_Squaf5=quafac4(betaSf5,xnewf5);
sMSE_Squaf5=sqrt(1/1000*sum((ynewf5-Yhat_Squaf5).^2))


%U poly0
[dmodel_U0f5, perf_U0f5] = dacefit(Upracf5,y_U5, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_U0f5] = predictor(xnewf5, dmodel_U0f5);
sMSE_U0f5=sqrt(1/1000*sum((ynewf5-Yhat_U0f5).^2))

%U poly1
[dmodel_U1f5, perf_U1f5] = dacefit(Upracf5,y_U5, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_U1f5] = predictor(xnewf5, dmodel_U1f5);
sMSE_U1f5=sqrt(1/1000*sum((ynewf5-Yhat_U1f5).^2))

%U poly2
%[dmodel_U2f5, perf_U2f5] = dacefit(Upracf5,y_U5, @regpoly2, @corrgauss, theta, lob, upb);
%[Yhat_U2f5] = predictor(xnewf5, dmodel_U2f5);
%sMSE_U2f5=sqrt(1/1000*sum((ynewf5-Yhat_U2f5).^2))

%U quadratic
[betaUf5]=nlinfit(Upracf5,y_U5,quafac4,beta0fac4);
Yhat_Uquaf5=quafac4(betaUf5,xnewf5);
sMSE_Uquaf5=sqrt(1/1000*sum((ynewf5-Yhat_Uquaf5).^2))


%R poly0
[dmodel_R0f5, perf_R0f5] = dacefit(Rpracf5,y_R5, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_R0f5] = predictor(xnewf5, dmodel_R0f5);
sMSE_R0f5=sqrt(1/1000*sum((ynewf5-Yhat_R0f5).^2))

%R poly1
[dmodel_R1f5, perf_R1f5] = dacefit(Rpracf5,y_R5, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_R1f5] = predictor(xnewf5, dmodel_R1f5);
sMSE_R1f5=sqrt(1/1000*sum((ynewf5-Yhat_R1f5).^2))


%R quadratic
[betaRf5]=nlinfit(Rpracf5,y_R5,quafac4,beta0fac4);
Yhat_Rquaf5=quafac4(betaRf5,xnewf5);
sMSE_Rquaf5=sqrt(1/1000*sum((ynewf5-Yhat_Rquaf5).^2))



%D poly0
[dmodel_D0f5, perf_D0f5] = dacefit(Dpracf5,y_D5, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_D0f5] = predictor(xnewf5, dmodel_D0f5);
sMSE_D0f5=sqrt(1/1000*sum((ynewf5-Yhat_D0f5).^2))

%D poly1
[dmodel_D1f5, perf_D1f5] = dacefit(Dpracf5,y_D5, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_D1f5] = predictor(xnewf5, dmodel_D1f5);
sMSE_D1f5=sqrt(1/1000*sum((ynewf5-Yhat_D1f5).^2))

%D poly2
%[dmodel_D2f5, perf_D2f5] = dacefit(Dpracf5,y_D5, @regpoly2, @corrgauss, theta, lob, upb);
%[Yhat_D2f5] = predictor(xnewf5, dmodel_D2f5);
%sMSE_D2f5=sqrt(1/1000*sum((ynewf5-Yhat_D2f5).^2))

%D quadratic
[betaDf5]=nlinfit(Dpracf5,y_D5,quafac4,beta0fac4);
Yhat_Dquaf5=quafac4(betaDf5,xnewf5);
sMSE_Dquaf5=sqrt(1/1000*sum((ynewf5-Yhat_Dquaf5).^2))


%L+D poly0/poly1/poly2/
LDdataf5=[Lpracf5;Dpracf5];
Y_L_D5=[y_L5;y_D5];
[LDdataf5_rev,ILDf5]=unique(LDdataf5,'rows');
YLDf5_rev=Y_L_D5(ILDf5);
[dmodel_LD0f5, perf_LD0f5] = dacefit(LDdataf5_rev,YLDf5_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_LD0f5] = predictor(xnewf5, dmodel_LD0f5);
sMSE_LD0f5=sqrt(1/1000*sum((ynewf5-Yhat_LD0f5).^2))

[dmodel_LD1f5, perf_LD1f5] = dacefit(LDdataf5_rev,YLDf5_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_LD1f5] = predictor(xnewf5, dmodel_LD1f5);
sMSE_LD1f5=sqrt(1/1000*sum((ynewf5-Yhat_LD1f5).^2))

%[dmodel_LD2f5, perf_LD2f5] = dacefit(LDdataf5_rev,YLDf5_rev, @regpoly2, @corrgauss, theta, lob, upb);
%[Yhat_LD2f5] = predictor(xnewf5, dmodel_LD2f5);
%sMSE_LD2f5=sqrt(1/1000*sum((ynewf5-Yhat_LD2f5).^2))

%L+D quadratic
[betaLDf5]=nlinfit(LDdataf5_rev,YLDf5_rev,quafac4,beta0fac4);
Yhat_LDquaf5=quafac4(betaLDf5,xnewf5);
sMSE_LDquaf5=sqrt(1/1000*sum((ynewf5-Yhat_LDquaf5).^2))


%S+D poly0/poly1/poly2
SDdataf5=[Spracf5;Dpracf5];
Y_SDf5=[y_S5;y_D5];
[SDdataf5_rev,ISDf5]=unique(SDdataf5,'rows');
YSDf5_rev=Y_SDf5(ISDf5);

[dmodel_SD0f5, perf_SD0f5] = dacefit(SDdataf5_rev,YSDf5_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_SD0f5] = predictor(xnewf5, dmodel_SD0f5);
sMSE_SD0f5=sqrt(1/1000*sum((ynewf5-Yhat_SD0f5).^2))


[dmodel_SD1f5, perf_SD1f5] = dacefit(SDdataf5_rev,YSDf5_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_SD1f5] = predictor(xnewf5, dmodel_SD1f5);
sMSE_SD1f5=sqrt(1/1000*sum((ynewf5-Yhat_SD1f5).^2))


%[dmodel_SD2f5, perf_SD2f5] = dacefit(SDdataf5_rev,YSDf5_rev, @regpoly2, @corrgauss, theta, lob, upb);
%[Yhat_SD2f5] = predictor(xnewf5, dmodel_SD2f5);
%sMSE_SD2f5=sqrt(1/1000*sum((ynewf5-Yhat_SD2f5).^2))

%S+D quadratic
[betaSDf5]=nlinfit(SDdataf5_rev,YSDf5_rev,quafac4,beta0fac4);
Yhat_SDquaf5=quafac4(betaSDf5,xnewf5);
sMSE_SDquaf5=sqrt(1/1000*sum((ynewf5-Yhat_SDquaf5).^2))


%U+D poly0/poly1/poly2
UDdataf5=[Upracf5;Dpracf5];
Y_U_D5=[y_U5;y_D5];
[UDdataf5_rev,IUDf5]=unique(UDdataf5,'rows');
YUDf5_rev=Y_U_D5(IUDf5);
[dmodel_UD0f5, perf_UD0f5] = dacefit(UDdataf5_rev,YUDf5_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_UD0f5] = predictor(xnewf5, dmodel_UD0f5);
sMSE_UD0f5=sqrt(1/1000*sum((ynewf5-Yhat_UD0f5).^2))

[dmodel_UD1f5, perf_UD1f5] = dacefit(UDdataf5_rev,YUDf5_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_UD1f5] = predictor(xnewf5, dmodel_UD1f5);
sMSE_UD1f5=sqrt(1/1000*sum((ynewf5-Yhat_UD1f5).^2))

%[dmodel_UD2f5, perf_UD2f5] = dacefit(UDdataf5_rev,YUDf5_rev, @regpoly2, @corrgauss, theta, lob, upb);
%[Yhat_UD2f5] = predictor(xnewf5, dmodel_UD2f5);
%sMSE_UD2f5=sqrt(1/1000*sum((ynewf5-Yhat_UD2f5).^2))

%U+D quadratic
[betaUDf5]=nlinfit(UDdataf5_rev,YUDf5_rev,quafac4,beta0fac4);
Yhat_UDquaf5=quafac4(betaUDf5,xnewf5);
sMSE_UDquaf5=sqrt(1/1000*sum((ynewf5-Yhat_UDquaf5).^2))


%R+D poly0/poly1/poly2
RDdataf5=[Rpracf5;Dpracf5];
Y_RDf5=[y_R5;y_D5];
[RDdataf5_rev,IRDf5]=unique(RDdataf5,'rows');
YRDf5_rev=Y_RDf5(IRDf5);

[dmodel_RD0f5, perf_RD0f5] = dacefit(RDdataf5_rev,YRDf5_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_RD0f5] = predictor(xnewf5, dmodel_RD0f5);
sMSE_RD0f5=sqrt(1/1000*sum((ynewf5-Yhat_RD0f5).^2))


[dmodel_RD1f5, perf_RD1f5] = dacefit(RDdataf5_rev,YRDf5_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_RD1f5] = predictor(xnewf5, dmodel_RD1f5);
sMSE_RD1f5=sqrt(1/1000*sum((ynewf5-Yhat_RD1f5).^2))


%[dmodel_RD2f5, perf_RD2f5] = dacefit(RDdataf5_rev,YRDf5_rev, @regpoly2, @corrgauss, theta, lob, upb);
%[Yhat_RD2f5] = predictor(xnewf5, dmodel_RD2f5);
%sMSE_RD2f5=sqrt(1/1000*sum((ynewf5-Yhat_RD2f5).^2))

%R+D quadratic
[betaRDf5]=nlinfit(RDdataf5_rev,YRDf5_rev,quafac4,beta0fac4);
Yhat_RDquaf5=quafac4(betaRDf5,xnewf5);
sMSE_RDquaf5=sqrt(1/1000*sum((ynewf5-Yhat_RDquaf5).^2))


%R+E poly0/poly1/poly2
RLdataf5=[Rpracf5;Lpracf5];
Y_RLf5=[y_R5;y_L5];
[RLdataf5_rev,IRLf5]=unique(RLdataf5,'rows');
YRLf5_rev=Y_RLf5(IRLf5);

[dmodel_RL0f5, perf_RL0f5] = dacefit(RLdataf5_rev,YRLf5_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_RL0f5] = predictor(xnewf5, dmodel_RL0f5);
sMSE_RL0f5=sqrt(1/1000*sum((ynewf5-Yhat_RL0f5).^2))


[dmodel_RL1f5, perf_RL1f5] = dacefit(RLdataf5_rev,YRLf5_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_RL1f5] = predictor(xnewf5, dmodel_RL1f5);
sMSE_RL1f5=sqrt(1/1000*sum((ynewf5-Yhat_RL1f5).^2))


%[dmodel_RL2f5, perf_RL2f5] = dacefit(RLdataf5_rev,YRLf5_rev, @regpoly2, @corrgauss, theta, lob, upb);
%[Yhat_RL2f5] = predictor(xnewf5, dmodel_RL2f5);
%sMSE_RL2f5=sqrt(1/1000*sum((ynewf5-Yhat_RL2f5).^2))

%R+E quadratic
[betaRLf5]=nlinfit(RLdataf5_rev,YRLf5_rev,quafac4,beta0fac4);
Yhat_RLquaf5=quafac4(betaRLf5,xnewf5);
sMSE_RLquaf5=sqrt(1/1000*sum((ynewf5-Yhat_RLquaf5).^2))


%R+S poly0/poly1/poly2
RSdataf5=[Rpracf5;Spracf5];
Y_RSf5=[y_R5;y_S5];
[RSdataf5_rev,IRSf5]=unique(RSdataf5,'rows');
YRSf5_rev=Y_RSf5(IRSf5);

[dmodel_RS0f5, perf_RS0f5] = dacefit(RSdataf5_rev,YRSf5_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_RS0f5] = predictor(xnewf5, dmodel_RS0f5);
sMSE_RS0f5=sqrt(1/1000*sum((ynewf5-Yhat_RS0f5).^2))


[dmodel_RS1f5, perf_RS1f5] = dacefit(RSdataf5_rev,YRSf5_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_RS1f5] = predictor(xnewf5, dmodel_RS1f5);
sMSE_RS1f5=sqrt(1/1000*sum((ynewf5-Yhat_RS1f5).^2))


%[dmodel_RS2f5, perf_RS2f5] = dacefit(RSdataf5_rev,YRSf5_rev, @regpoly2, @corrgauss, theta, lob, upb);
%[Yhat_RS2f5] = predictor(xnewf5, dmodel_RS2f5);
%sMSE_RS2f5=sqrt(1/1000*sum((ynewf5-Yhat_RS2f5).^2))

%R+S quadratic
[betaRSf5]=nlinfit(RSdataf5_rev,YRSf5_rev,quafac4,beta0fac4);
Yhat_RSquaf5=quafac4(betaRSf5,xnewf5);
sMSE_RSquaf5=sqrt(1/1000*sum((ynewf5-Yhat_RSquaf5).^2))


%U+S poly0/poly1/poly2
USdataf5=[Upracf5;Spracf5];
Y_USf5=[y_U5;y_S5];
[USdataf5_rev,IUSf5]=unique(USdataf5,'rows');
YUSf5_rev=Y_USf5(IUSf5);

[dmodel_US0f5, perf_US0f5] = dacefit(USdataf5_rev,YUSf5_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_US0f5] = predictor(xnewf5, dmodel_US0f5);
sMSE_US0f5=sqrt(1/1000*sum((ynewf5-Yhat_US0f5).^2))


[dmodel_US1f5, perf_US1f5] = dacefit(USdataf5_rev,YUSf5_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_US1f5] = predictor(xnewf5, dmodel_US1f5);
sMSE_US1f5=sqrt(1/1000*sum((ynewf5-Yhat_US1f5).^2))


%[dmodel_US2f5, perf_US2f5] = dacefit(USdataf5_rev,YUSf5_rev, @regpoly2, @corrgauss, theta, lob, upb);
%[Yhat_US2f5] = predictor(xnewf5, dmodel_US2f5);
%sMSE_US2f5=sqrt(1/1000*sum((ynewf5-Yhat_US2f5).^2))

%U+S quadratic
[betaUSf5]=nlinfit(USdataf5_rev,YUSf5_rev,quafac4,beta0fac4);
Yhat_USquaf5=quafac4(betaUSf5,xnewf5);
sMSE_USquaf5=sqrt(1/1000*sum((ynewf5-Yhat_USquaf5).^2))


%U+E poly0/poly1/poly2
UEdataf5=[Upracf5;Lpracf5];
Y_UEf5=[y_U5;y_L5];
[UEdataf5_rev,IUEf5]=unique(UEdataf5,'rows');
YUEf5_rev=Y_UEf5(IUEf5);

[dmodel_UE0f5, perf_UE0f5] = dacefit(UEdataf5_rev,YUEf5_rev, @regpoly0, @corrgauss, theta, lob, upb);
[Yhat_UE0f5] = predictor(xnewf5, dmodel_UE0f5);
sMSE_UE0f5=sqrt(1/1000*sum((ynewf5-Yhat_UE0f5).^2))


[dmodel_UE1f5, perf_UE1f5] = dacefit(UEdataf5_rev,YUEf5_rev, @regpoly1, @corrgauss, theta, lob, upb);
[Yhat_UE1f5] = predictor(xnewf5, dmodel_UE1f5);
sMSE_UE1f5=sqrt(1/1000*sum((ynewf5-Yhat_UE1f5).^2))


%[dmodel_UE2f5, perf_UE2f5] = dacefit(UEdataf5_rev,YUEf5_rev, @regpoly2, @corrgauss, theta, lob, upb);
%[Yhat_UE2f5] = predictor(xnewf5, dmodel_UE2f5);
%sMSE_UE2f5=sqrt(1/1000*sum((ynewf5-Yhat_UE2f5).^2))

%U+E quadratic
[betaUEf5]=nlinfit(UEdataf5_rev,YUEf5_rev,quafac4,beta0fac4);
Yhat_UEquaf5=quafac4(betaUEf5,xnewf5);
sMSE_UEquaf5=sqrt(1/1000*sum((ynewf5-Yhat_UEquaf5).^2))

