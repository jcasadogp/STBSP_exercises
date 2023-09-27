%% 1. INITIALIZE
disp('1: INITIALIZE')

% IMPORTANT: DO NOT REMOVE OR MODIFY ANY OF THE GIVEN CODE, UNLESS
% EXPLICITLY SPECIFIED

clear; clc;close all;

% load the training data
load('trainingData.mat');

% data related parameters
fs = 25000; % sampling frequency
spike_window = 25; % samples spike window
nbChannels = size(trainingData,2);

figure; plotMultiChannel(trainingData((1:1000), :), 2.5, '');
title('Training Data (1-1000)', 'FontSize', 20);

%% 2. CALCULATE TEMPLATE
disp('2: CALCULATE TEMPLATE')

% calculate template
template = calcTemplate(trainingData, trainingLabels, spike_window);
% vectorize template for use in filter design
vecTemplate = mat2stacked(template');

% TODO: Visualize the calculated template using plotMultiChannel
% visualize template

figure; plotMultiChannel(template, 2.5, '');
title('Template', 'FontSize', 20);
%% 3. Template as a filter
disp('3: Template as a filter')

% TODO: complete applyMultiChannelFilter function
%  make sure to transform the filter coefficients back to a matrix
%  using stacked2mat.

templateFilter = template;
templateFilterOut = applyMultiChannelFilter(trainingData, templateFilter);

% calculate output power
templateFilterOut = templateFilterOut.^2;

figure; plot(templateFilterOut);
title('Template Filter Output', 'FontSize', 20);
%% 4. CHOOSE F1-OPTIMAL THRESHOLD FOR template filter
disp('4: CHOOSE F1-OPTIMAL THRESHOLD FOR template filter')

% shift labels to match output
outLabels = trainingLabels + floor(spike_window/2);

% try different thresholds on the filter output
outStd = std(templateFilterOut);
C = 10;
max_perf = 0;
SNR_senss = [];
SNR_precs = [];
SNR_thrs = [];
while 1
    % threshold segments
    thr = C*outStd;
    detections = templateFilterOut > thr;
    cuttedDetections = cutMask(detections, spike_window);
    
    % output training segments
    labels = zeros(size(detections));
    labels(outLabels) = 1;
    cuttedLabels = cutMask(labels, spike_window);
    
    if length(cuttedDetections) == 0
        break;
    end
    
    % validate the detections
    [sens, prec] = validateDetections(cuttedDetections, cuttedLabels);
    
    SNR_thrs = [SNR_thrs; thr];
    SNR_senss = [SNR_senss; sens];
    SNR_precs = [SNR_precs; prec];

    C = C + 3;
end
%% 5. PLOT P-R CURVE FOR template filter AND CHOOSE A THRESHOLD
disp('5: PLOT P-R CURVE FOR template filter AND CHOOSE A THRESHOLD')

% TODO: plot a P-R curve using SNR_senss and SNR_precs

figure; plot(SNR_senss, SNR_precs, '.-'); axis tight
xlabel('Recall', 'FontSize', 15)
ylabel('Precision', 'FontSize', 15)
title('P-R curve', 'FontSize', 20)

% TODO: based on the plotted P-R curve choose a threshold

f_score = zeros(size(SNR_precs));
for idx=1:length(SNR_precs)
    f_score(idx)=2*SNR_precs(idx)*SNR_senss(idx)/(SNR_precs(idx)+SNR_senss(idx));
    if isnan(f_score(idx))
        f_score(idx)=0;
    end
end

max_f = max(f_score);
idx = find(f_score == max_f);
max_threshold = SNR_thrs(idx)

x = 1:1:length(f_score);
figure; 
subplot(2,1,1); plot(x, f_score); xline(idx, 'r--', {'Maximum F-score'}); title('F-score', 'FontSize', 20)
subplot(2,1,2); plot(x, SNR_thrs); xline(idx, 'r--', {'Threshold for Maximum F-score'}); title('Threshold', 'FontSize', 20)
clear x
%% 6. VALIDATE template filter ON TESTING DATA
disp('6: VALIDATE template filter ON TESTING DATA')

% load the testing data
load('testingData.mat');

% calculate filter output
testingtemplateFilterOut = applyMultiChannelFilter(testingData, templateFilter);
testingtemplateFilterOut = testingtemplateFilterOut.^2;

figure; plot(testingtemplateFilterOut);
title('Testing - Template Filter Output', 'FontSize', 20);

% shift label to match output delay
testingOutLabels = testingLabels + floor(spike_window/2);

% threshold segments
detections = testingtemplateFilterOut > max_threshold;
cuttedDetections = cutMask(detections, spike_window);

% output testing segments
labels = zeros(size(detections));
labels(testingOutLabels) = 1;
cuttedLabels = cutMask(labels, spike_window);

% validate the detections
[sens_templateFilter, prec_templateFilter] = validateDetections(cuttedDetections, cuttedLabels);
fprintf('template-filter: for your threshold: recall: %.3f, precision: %.3f\n\n', sens_templateFilter, prec_templateFilter);

% visualize Matched filter output power
figure; hold on;
plot(testingtemplateFilterOut, 'DisplayName', 'Matched filter');
plot(testingOutLabels, testingtemplateFilterOut(testingOutLabels), 'g*', 'DisplayName', 'Testing labels');
title('Template filter output on testing data', 'FontSize', 20)
xlabel('Discrete time [samples]', 'FontSize', 15)
ylabel('Output power [arb. unit]', 'FontSize', 15)
legend('show', 'FontSize', 16)

%% 7. DESIGN Matched-filter
disp('7: DESIGN Matched-filter')

% find noise segments
noiseSegments = findNoiseSegments(trainingData, spike_window);

% calculate noise covariance matrix

% TODO: Write a function which computes the noise covariance matrix using only the
% data in noise segments.

% Let the estimated noise covariance be this variable.
tempNoiseCovariance =  calcCovMatrix2(trainingData, noiseSegments, spike_window);
figure; imagesc(tempNoiseCovariance); colorbar
title('Lagged Noise Covariance Matrix', 'FontSize', 20);

ch_3 = (51:1:75);
figure; imagesc(tempNoiseCovariance(ch_3,ch_3)); colorbar
title('Lagged Noise Covariance Matrix - Channel 3', 'FontSize', 20);
clear ch_3

% regularize the noise covariance matrix
noiseCov = regularizeCov(tempNoiseCovariance,1);

% TODO: calculate the Matched filter using noiseCov and the vectorized
% template, make sure to transform the filter coefficients back to a matrix
% Store the matched filter in maxSNR

maxSNR = pinv(noiseCov)*vecTemplate;
maxSNR = stacked2mat(maxSNR', nbChannels);

figure; plotMultiChannel(maxSNR, 2.5, '');
title('Matched Filter', 'FontSize', 20);

% TODO: complete applyMultiChannelFilter function
maxSNROut = applyMultiChannelFilter(trainingData, maxSNR);

% calculate output power
maxSNROut = maxSNROut.^2;

figure; plot(maxSNROut);
title('Matched Filter Output', 'FontSize', 20);
%% 8. CHOOSE F1-OPTIMAL THRESHOLD FOR Matched filter
disp('8: CHOOSE F1-OPTIMAL THRESHOLD FOR Matched filter')

% shift labels to match output
outLabels = trainingLabels + floor(spike_window/2);

% try different thresholds on the filter output
outStd = std(maxSNROut);
C = 10;
max_perf = 0;
SNR_senss = [];
SNR_precs = [];
SNR_thrs = [];
while 1
    % threshold segments
    thr = C*outStd;
    detections = maxSNROut > thr;
    cuttedDetections = cutMask(detections, spike_window);
    
    % output training segments
    labels = zeros(size(detections));
    labels(outLabels) = 1;
    cuttedLabels = cutMask(labels, spike_window);
    
    if length(cuttedDetections) == 0
        break;
    end
    
    % validate the detections
    [sens, prec] = validateDetections(cuttedDetections, cuttedLabels);
    
    SNR_thrs = [SNR_thrs; thr];
    SNR_senss = [SNR_senss; sens];
    SNR_precs = [SNR_precs; prec];

    C = C + 3;
end

%% 9. PLOT P-R CURVE FOR Matched filter AND CHOOSE A THRESHOLD
disp('9: PLOT P-R CURVE FOR Matched filter AND CHOOSE A THRESHOLD')

% TODO: plot a P-R curve using SNR_senss and SNR_precs
figure; plot(SNR_senss, SNR_precs, '.-'); axis tight
xlabel('Recall', 'FontSize', 15)
ylabel('Precision', 'FontSize', 15)
title('P-R curve', 'FontSize', 20)

% TODO: based on the plotted P-R curve choose a threshold

f_score = zeros(size(SNR_precs));
for idx=1:length(SNR_precs)
    f_score(idx)=2*SNR_precs(idx)*SNR_senss(idx)/(SNR_precs(idx)+SNR_senss(idx));
    if isnan(f_score(idx))
        f_score(idx)=0;
    end
end

max_f = max(f_score);
idx = find(f_score == max_f);
max_threshold = SNR_thrs(idx);

x = 1:1:length(f_score);
figure; 
subplot(2,1,1); plot(x, f_score); xline(idx, 'r--', {'Maximum F-score'}); title('F-score', 'FontSize', 20)
subplot(2,1,2); plot(x, SNR_thrs); xline(idx, 'r--', {'Threshold for Maximum F-score'}); title('Threshold', 'FontSize', 20)
clear x
%% 10. VALIDATE Matched filter ON TESTING DATA
disp('10: VALIDATE Matched filter ON TESTING DATA')

% load the testing data
load('testingData.mat');

% calculate filter output
testingMaxSNROut = applyMultiChannelFilter(testingData, maxSNR);
testingMaxSNROut = testingMaxSNROut.^2;

% shift label to match output delay
testingOutLabels = testingLabels + floor(spike_window/2);

% threshold segments
detections = testingMaxSNROut > max_threshold;
cuttedDetections = cutMask(detections, spike_window);

% output testing segments
labels = zeros(size(detections));
labels(testingOutLabels) = 1;
cuttedLabels = cutMask(labels, spike_window);

% validate the detections
[sens_SNR, prec_SNR] = validateDetections(cuttedDetections, cuttedLabels);
fprintf('matched-filter: for your threshold: recall: %.3f, precision: %.3f\n\n', sens_SNR, prec_SNR);

% visualize Matched filter output power
figure; hold on;
plot(testingMaxSNROut, 'DisplayName', 'Matched filter');
plot(testingOutLabels, testingMaxSNROut(testingOutLabels), 'g*', 'DisplayName', 'Testing labels');
title('Matched filter output on testing data', 'FontSize', 20)
xlabel('Discrete time [samples]', 'FontSize', 15)
ylabel('Output power [arb. unit]', 'FontSize', 15)
legend('show', 'FontSize', 16)

%% 11. DESIGN MAX-SPIR FILTER
disp('11: DESIGN MAX-SPIR FILTER')

% find interference segments
intSegments = findInterferenceSegments(trainingData, maxSNR, outLabels, spike_window);

% visualize interference segments
figure; hold on;
plot(maxSNROut, 'DisplayName', 'Matched filter');
plot(intSegments, maxSNROut(intSegments), 'r*', 'DisplayName', 'Detected interference segments');
plot(outLabels, maxSNROut(outLabels), 'g*', 'DisplayName', 'Training labels');
title('Matched filter output on training data', 'FontSize', 20)
xlabel('Discrete time [samples]', 'FontSize', 15)
ylabel('Output power [arb. unit]', 'FontSize', 15)
legend('show', 'FontSize', 16)

% calculate interference covariance matrix and store it in tempIntCov
% TODO: Re-use the function used to calculate noise covariance here.
tempIntCov = calcCovMatrix2(trainingData, intSegments, spike_window);

figure; imagesc(tempIntCov); colorbar
title('Lagged Interference Covariance Matrix', 'FontSize', 20);

% regularize the interference covariance matrix
intCov = regularizeCov(tempIntCov,0.01);

% TODO: calculate the max-SPIR filter using intCov and the vectorized
% template, make sure to transform the filter coefficients back to a matrix
% using stacked2mat.

maxSPIR = pinv(intCov)*vecTemplate;
maxSPIR = stacked2mat(maxSPIR', nbChannels);

figure; plotMultiChannel(maxSPIR, 2.5, '');
title('Max-SPIR Filter', 'FontSize', 20);

% calculate filter output
maxSPIROut = applyMultiChannelFilter(trainingData, maxSPIR);
maxSPIROut = maxSPIROut.^2;

figure; plotMultiChannel(maxSPIROut, 2.5, '');
title('Max-SPIR Filter Output', 'FontSize', 20);

%% 12. CHOOSE F1-OPTIMAL THRESHOLD FOR MAX-SPIR
disp('12: CHOOSE F1-OPTIMAL THRESHOLD FOR MAX-SPIR')

% try different thresholds on the filter output
outStd = std(maxSPIROut);
C = 10;
max_perf = 0;
SPIR_senss = [];
SPIR_precs = [];
SPIR_thrs = [];
while 1
    % threshold segments
    thr = C*outStd;
    detections = maxSPIROut > thr;
    cuttedDetections = cutMask(detections, spike_window);
    
    % output training segments
    labels = zeros(size(detections));
    labels(outLabels) = 1;
    cuttedLabels = cutMask(labels, spike_window);
    
    if length(cuttedDetections) == 0
        break;
    end
    
    % validate the detections
    [sens, prec] = validateDetections(cuttedDetections, cuttedLabels);
    
    SPIR_thrs = [SPIR_thrs; thr];
    SPIR_senss = [SPIR_senss; sens];
    SPIR_precs = [SPIR_precs; prec];

    C = C + 3;
end

%% 13. PLOT P-R CURVE FOR MAX-SPIR AND CHOOSE A THRESHOLD
disp('13: PLOT P-R CURVE FOR MAX-SPIR AND CHOOSE A THRESHOLD')

% TODO: plot a P-R curve using SPIR_senss and SPIR_precs

figure; plot(SPIR_senss, SPIR_precs, '.-'); axis tight
xlabel('Recall', 'FontSize', 15)
ylabel('Precision', 'FontSize', 15)
title('P-R curve', 'FontSize', 20)

% TODO: based on the plotted P-R curve choose a threshold

f_score = zeros(size(SPIR_precs));
for idx=1:length(SPIR_precs)
    f_score(idx)=2*SPIR_precs(idx)*SPIR_senss(idx)/(SPIR_precs(idx)+SPIR_senss(idx));
    if isnan(f_score(idx))
        f_score(idx)=0;
    end
end

max_f = max(f_score);
idx = find(f_score == max_f);
max_threshold_SPIR = SPIR_thrs(idx);

x = 1:1:length(f_score);
figure; 
subplot(2,1,1); plot(x, f_score); xline(idx, 'r--', {'Maximum F-score'}); title('F-score', 'FontSize', 20)
subplot(2,1,2); plot(x, SPIR_thrs); xline(idx, 'r--', {'Threshold for Maximum F-score'}); title('Threshold', 'FontSize', 20)
clear x

%% 14. VALIDATE MAX-SPIR FILTER ON TESTING DATA
disp('14: VALIDATE MAX-SPIR FILTER ON TESTING DATA')

testingMaxSPIROut = applyMultiChannelFilter(testingData, maxSPIR);
testingMaxSPIROut = testingMaxSPIROut.^2;

testingOutLabels = testingLabels + floor(spike_window/2);

detections = testingMaxSPIROut > max_threshold_SPIR;
cuttedDetections = cutMask(detections, spike_window);

labels = zeros(size(detections));
labels(testingOutLabels) = 1;
cuttedLabels = cutMask(labels, spike_window);

% validate the detections
[sens_SPIR, prec_SPIR] = validateDetections(cuttedDetections, cuttedLabels);
fprintf('max-SPIR: for the maximum F1-score threshold: recall: %.3f, precision: %.3f\n', sens_SPIR, prec_SPIR);

%% 15. COMPARE BOTH FILTER OUTPUTS FROM TESTING DATA
disp('15: COMPARE BOTH FILTER OUTPUTS FROM TESTING DATA')

% normalize filter outputs for plotting purposes
[normSNR, normSPIR] = normalizeOutputs(testingMaxSNROut, testingMaxSPIROut, testingOutLabels, spike_window);

% plot normalized filter outputs on testing data
figure; hold on;
plot(normSPIR, 'DisplayName', 'normalized max-SPIR');
plot(normSNR, 'DisplayName', 'normalized Matched filter');
plot(testingOutLabels, normSPIR(testingOutLabels), 'g*', 'DisplayName', 'Testing labels');
title('Filter outputs on testing data', 'FontSize', 20)
xlabel('Discrete time [samples]', 'FontSize', 15);
ylabel('Output power [arb. unit]', 'FontSize', 15);
legend('show', 'FontSize', 16)
