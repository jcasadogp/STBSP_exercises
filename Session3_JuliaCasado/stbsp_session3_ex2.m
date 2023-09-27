%high anxiety pregnant women (ID 518, 542, 547, 562, 576)
load('RR518.mat')
RR1 = RR;
load('RR542.mat')
RR2 = RR;
load('RR547.mat')
RR3 = RR;
load('RR562.mat')
RR4 = RR;
load('RR576.mat')
RR5 = RR;

%low anxiety pregnant women (ID 507, 514, 571, 619, 621),
load('RR507.mat')
RR6 = RR;
load('RR514.mat')
RR7 = RR;
load('RR571.mat')
RR8 = RR;
load('RR619.mat')
RR9 = RR;
load('RR621.mat')
RR10 = RR;
%% Detrended Fluctuation Analysis
%DFA for every signal in order to detect the differences in a1 and a2
[a1_1 , a2_1, Fall1] = DFA(RR1);
[a1_2 , a2_2 ,Fall2] = DFA(RR2);
[a1_3 , a2_3, Fall3] = DFA(RR3);
[a1_4 , a2_4, Fall4] = DFA(RR4);
[a1_5 , a2_5, Fall5] = DFA(RR5);
[a1_6 , a2_6, Fall6] = DFA(RR6);
[a1_7 , a2_7, Fall7] = DFA(RR7);
[a1_8 , a2_8, Fall8] = DFA(RR8);
[a1_9 , a2_9, Fall9] = DFA(RR9);
[a1_10 ,a2_10, Fall10] = DFA(RR10);

%Matrices which contain the a1 and a values for the two populations
a1_high = [a1_1 a1_2 a1_3 a1_4 a1_5];
a1_low = [a1_6 a1_7 a1_8 a1_9 a1_10];
a2_high = [a2_1 a2_2 a2_3 a2_4 a2_5];
a2_low = [a2_6 a2_7 a2_8 a2_9 a2_10];

%% 1/f slope

%resampled signals
RR1_ = resample(RR1,20,1);  %length(resampled) = p / q * length(RR)
RR2_ = resample(RR2,20,1);
RR3_ = resample(RR3,20,1);
RR4_ = resample(RR4,20,1);
RR5_ = resample(RR5,20,1);
RR6_ = resample(RR6,20,1);
RR7_ = resample(RR7,20,1);
RR8_ = resample(RR8,20,1);
RR9_ = resample(RR9,20,1);
RR10_ = resample(RR10,20,1);

subplot(3,2,1); [sl1] = slope(RR1_);
subplot(3,2,2); sl2 = slope(RR2_);
subplot(3,2,3); sl3 =slope(RR3_);
subplot(3,2,4); sl4 = slope(RR4_);
subplot(3,2,5); sl5 = slope(RR5_);
figure;
subplot(3,2,1); sl6 = slope(RR6_);
subplot(3,2,2); sl7 = slope(RR7_);
subplot(3,2,3); sl8 = slope(RR8_);
subplot(3,2,4); sl9 = slope(RR9_);
subplot(3,2,5); sl10 = slope(RR10_);

sl_high = [sl1 sl2 sl3 sl4 sl6];
sl_low = [sl6 sl7 sl8 sl9 sl10];

%% Poincaré plot
%high anxiety pregnant women (ID 518, 542, 547, 562, 576)
SD1all_high = [];
SD2all_high = [];
SDRRall_high = [];

subplot(3,2,1);
[SD1,SD2,SDRR] = poincare_plot(RR1);
SD1all_high = [SD1all_high,SD1];
SD2all_high = [SD2all_high,SD2];
SDRRall_high = [SDRRall_high,SDRR];

subplot(3,2,2);
[SD1,SD2,SDRR] = poincare_plot(RR2);
SD1all_high = [SD1all_high,SD1];
SD2all_high = [SD2all_high,SD2];
SDRRall_high = [SDRRall_high,SDRR];

subplot(3,2,3);
[SD1,SD2,SDRR] = poincare_plot(RR3);
SD1all_high = [SD1all_high,SD1];
SD2all_high = [SD2all_high,SD2];
SDRRall_high = [SDRRall_high,SDRR];

subplot(3,2,4);
[SD1,SD2,SDRR] = poincare_plot(RR4);
SD1all_high = [SD1all_high,SD1];
SD2all_high = [SD2all_high,SD2];
SDRRall_high = [SDRRall_high,SDRR];

subplot(3,2,5);
[SD1,SD2,SDRR] = poincare_plot(RR5);
SD1all_high = [SD1all_high,SD1];
SD2all_high = [SD2all_high,SD2];
SDRRall_high = [SDRRall_high,SDRR];

SD1high_mean = mean(SD1all_high);
SD2high_mean = mean(SD2all_high);
SDRRhigh_mean = mean(SDRRall_high);

%low anxiety pregnant women (ID 507, 514, 571, 619, 621),
figure
SD1all_low = [];
SD2all_low = [];
SDRRall_low = [];

subplot(3,2,1);
[SD1,SD2,SDRR] = poincare_plot(RR6);
SD1all_low = [SD1all_low,SD1];
SD2all_low = [SD2all_low,SD2];
SDRRall_low = [SDRRall_low,SDRR];

subplot(3,2,2);
[SD1,SD2,SDRR] = poincare_plot(RR7);
SD1all_low = [SD1all_low,SD1];
SD2all_low = [SD2all_low,SD2];
SDRRall_low = [SDRRall_low,SDRR];

subplot(3,2,3);
[SD1,SD2,SDRR] = poincare_plot(RR8);
SD1all_low = [SD1all_low,SD1];
SD2all_low = [SD2all_low,SD2];
SDRRall_low = [SDRRall_low,SDRR];

subplot(3,2,4);
[SD1,SD2,SDRR] = poincare_plot(RR9);
SD1all_low = [SD1all_low,SD1];
SD2all_low = [SD2all_low,SD2];
SDRRall_low = [SDRRall_low,SDRR];

subplot(3,2,5);
[SD1,SD2,SDRR] = poincare_plot(RR10);
SD1all_low = [SD1all_low,SD1];
SD2all_low = [SD2all_low,SD2];
SDRRall_low = [SDRRall_low,SDRR];

SD1low_mean = mean(SD1all_low);
SD2low_mean = mean(SD2all_low);
SDRRlow_mean = mean(SDRRall_low);

%%  Nonlinear HRV measures
%FD = fractal dimension
%SE = Sample Entropy
%LE = Lyapunov Exponent
%CD = Correlation Dimension
FD_high = [];
SE_high = [];
LE_high = [];
CD_high = [];
FD_low = [];
SE_low = [];
LE_low = [];
CD_low = [];

[FD,SE,LE,CD] = nonlin_measures(RR1);
FD_high = [FD_high,FD];
SE_high = [SE_high,SE];
LE_high = [LE_high,LE];
CD_high = [CD_high,CD];

[FD,SE,LE,CD] = nonlin_measures(RR2);
FD_high = [FD_high,FD];
SE_high = [SE_high,SE];
LE_high = [LE_high,LE];
CD_high = [CD_high,CD];

[FD,SE,LE,CD] = nonlin_measures(RR3);
FD_high = [FD_high,FD];
SE_high = [SE_high,SE];
LE_high = [LE_high,LE];
CD_high = [CD_high,CD];

[FD,SE,LE,CD] = nonlin_measures(RR4);
FD_high = [FD_high,FD];
SE_high = [SE_high,SE];
LE_high = [LE_high,LE];
CD_high = [CD_high,CD];

[FD,SE,LE,CD] = nonlin_measures(RR5);
FD_high = [FD_high,FD];
SE_high = [SE_high,SE];
LE_high = [LE_high,LE];
CD_high = [CD_high,CD];

[FD,SE,LE,CD] = nonlin_measures(RR6);
FD_low = [FD_low,FD];
SE_low = [SE_low,SE];
LE_low = [LE_low,LE];
CD_low = [CD_low,CD];

[FD,SE,LE,CD] = nonlin_measures(RR7);
FD_low = [FD_low,FD];
SE_low = [SE_low,SE];
LE_low = [LE_low,LE];
CD_low = [CD_low,CD];

[FD,SE,LE,CD] = nonlin_measures(RR8);
FD_low = [FD_low,FD];
SE_low= [SE_low,SE];
LE_low = [LE_low,LE];
CD_low = [CD_low,CD];

[FD,SE,LE,CD] = nonlin_measures(RR9);
FD_low = [FD_low,FD];
SE_low = [SE_low,SE];
LE_low = [LE_low,LE];
CD_low = [CD_low,CD];

[FD,SE,LE,CD] = nonlin_measures(RR10);
FD_low = [FD_low,FD];
SE_low = [SE_low,SE];
LE_low = [LE_low,LE];
CD_low = [CD_low,CD];


