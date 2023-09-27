%% 1. EEG: Artifal removal

clear; clc
close all

% Required toolbox
addpath('/Users/julia/Documents/MATLAB/STBSP/mwf-artifact-removal-master/mwf')
addpath(genpath('/Users/julia/Documents/MATLAB/STBSP/eeglab2021.1'));


%% 1.2 Unsupervised artifact removal using CCA
% Load data
load('eegdata_artifacts.mat')

figure; eegplot_simple(eegdata, fs, '', 200);
title('EEG data', 'FontSize', 20);

[rows, columns] = size(eegdata);

t = linspace(0,columns/fs,columns);

d=1;
x = eegdata(1:rows,1:columns-d);
y = eegdata(1:rows,1+d:columns);

% --- canoncorr function ---
[Wx,Wy,r, U, V] = canoncorr(x',y');

S_est = Wx' * x;

figure; eegplot_simple(S_est, fs, '', 12);
title('Estimated sources', 'FontSize', 20);

figure; subplot(1,2,1); plot(r, 'o-')
title('Autocorrelation values', 'FontSize', 20)
axis tight

%% Muscle artifal removal

Z = zeros(size(S_est));
S_muscle = S_est;

for i =37:64    %Ch with muscle artifact
    S_muscle(i,(7000:14279)) = Z(i,(7000:14279));
end

figure, eegplot_simple(S_muscle, fs);
title('Estimated sources witouth muscle artifact', 'FontSize', 20);

X_clean = pinv(Wx)*S_muscle;

figure; eegplot_simple(eegdata, fs, '', 200);
title('EEG data', 'FontSize', 20);

figure; eegplot_simple(X_clean, fs, '', 200);
title('Clean EEG data', 'FontSize', 20);

figure;

plot(t, eegdata(40,:), 'g'); axis tight; ylim([-400 400])
hold on; plot(t(1:14279), X_clean(40,:), 'b');   axis tight; ylim([-400 400])
legend('EEG data', 'Clean EEG', 'FontSize', 15)
title('EEG data before and after muscle artifact removal', 'FontSize', 20)
subtitle('Channel 40', 'FontSize', 14)


%% 1.3 Supervised artifact removal using MWF --> EYE-BLINK

% EEG segmentation
mask = mwf_getmask(eegdata, fs);
clear ans

%% MWF artifact removal -- delay=0

[clean_EEG, artifact, W, SER, ARR]= mwf_process(eegdata, mask);

figure;
subplot(2,1,1);
plot(eegdata(33,(1:6694)), 'k'); hold on;
plot(artifact(33,(1:6694)), 'r');
legend('Original EEG signal', 'Eye-blink artifact', 'FontSize', 16)

subplot(2,1,2);     
plot(eegdata(33,(1:6694)), 'k'); hold on;
plot(clean_EEG(33,(1:6694)), 'b');
legend('Original EEG signal', 'Clean EEG signal', 'FontSize', 16)

txt1 = ['Ch33 without delay']; %Without SER and ARR
txt2 = ['Ch33 without delay - SER: ', num2str(SER), '. ARR: ', num2str(ARR)]; %With SER and ARR
sgtitle(txt2, 'FontSize', 20) %Choose txt
clear txt1; clear txt2

    
%% MWF artifact removal -- delay=3
[clean_EEG_3, artifact_3, W_3, SER_3, ARR_3]= mwf_process(eegdata, mask, 3);

figure;
subplot(2,1,1);     
plot(eegdata(33,(1:6694)), 'k'); hold on;        
plot(artifact_3(33,(1:6694)), 'r');
legend('Original EEG signal', 'Eye-blink artifact', 'FontSize', 16)

subplot(2,1,2);     
plot(eegdata(33,(1:6694)), 'k'); hold on;
plot(clean_EEG_3(33,(1:6694)), 'b');
legend('Original EEG signal', 'Clean EEG signal', 'FontSize', 16)

txt1 = ['Ch33 with delay=3']; %Without SER and ARR
txt2 = ['Ch33 with delay=3 - SER: ', num2str(SER_3), '. ARR: ', num2str(ARR_3)]; %With SER and ARR
sgtitle(txt2, 'FontSize', 20) %Choose txt
clear txt1; clear txt2

%% MWF artifact removal -- delay=0 -- different rank

p = mwf_params('rank', 'first', 'rankopt', 2);
[W_r, Lambda] = mwf_compute(eegdata, mask, p);
[clean_EEG_r, artifact_r] = mwf_apply(eegdata, W_r);

figure;
subplot(2,1,1);
plot(eegdata(33,(1:6694)), 'k'); hold on;
plot(artifact_r(33,(1:6694)), 'r');
legend('Original EEG signal', 'Eye-blink artifact', 'FontSize', 16)

subplot(2,1,2);     
plot(eegdata(33,(1:6694)), 'k'); hold on;
plot(clean_EEG_r(33,(1:6694)), 'b');
legend('Original EEG signal', 'Clean EEG signal', 'FontSize', 16)

txt1 = ['Ch33 without delay. Rank = 2']; %Without SER and ARR
txt2 = ['Ch33 without delay - SER: ', num2str(SER), '. ARR: ', num2str(ARR)]; %With SER and ARR
sgtitle(txt1, 'FontSize', 20) %Choose txt
clear txt1; clear txt2


%% 1.3 Supervised artifact removal using MWF --> MUSCLE

load('eegdata_artifacts.mat')
[channels, samples] = size(eegdata);

% EEG segmentation
mask_muscle = mwf_getmask(eegdata, fs);
clear ans
%% MWF artifact removal -- delay=3
[clean_EEG_3, artifact_3, W_3, SER_3, ARR_3]= mwf_process(eegdata, mask_muscle, 3);

figure;
subplot(2,1,1);     
plot(eegdata(33,(6694:14280)), 'k'); hold on;        
plot(artifact_3(33,(6694:14280)), 'r');
legend('Original EEG signal', 'Muscle artifact', 'FontSize', 16)

subplot(2,1,2);     
plot(eegdata(33,(6694:14280)), 'k'); hold on;
plot(clean_EEG_3(33,(6694:14280)), 'b');
legend('Original EEG signal', 'Clean EEG signal', 'FontSize', 16)

txt1 = ['Ch33 with delay=3']; %Without SER and ARR
txt2 = ['Ch33 with delay=3 - SER: ', num2str(SER_3), '. ARR: ', num2str(ARR_3)]; %With SER and ARR
sgtitle(txt1, 'FontSize', 20) %Choose txt
clear txt1; clear txt2
