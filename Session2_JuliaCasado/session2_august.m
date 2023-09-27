clear; clc; close all

addpath('/Users/julia/Documents/MATLAB/tensorlab4.0beta');

%% 2.1 MLPCA
load('ex1data.mat');

T1 = tens2mat(Tn, 1); T2 = tens2mat(Tn, 2); T3 = tens2mat(Tn, 3);

[Ue,Se,sve] = mlsvd(Tn);

% Mode-1
figure(1)
plot3(T1(1,:),T1(2,:),T1(3,:),'x'), grid    % mode-1 vectors
title('tensor mode-1 vectors', 'FontSize', 15)
disp('tensor mode-1 singular values:'), disp(sve{1}(1:9)')
figure(4); plot3(T1(4,:),T1(5,:),T1(6,:),'x'), grid    % mode-1 vectors
title('tensor mode-1 vectors', 'FontSize', 15)
figure(5); plot3(T1(7,:),T1(8,:),T1(9,:),'x'), grid    % mode-1 vectors
title('tensor mode-1 vectors', 'FontSize', 15)

% Mode-2
figure(2)
plot3(T2(1,:),T2(2,:),T2(3,:),'x'), grid    % mode-2 vectors
title('tensor mode-2 vectors', 'FontSize', 15)
disp('tensor mode-2 singular values:'), disp(sve{2}(1:9)')

% Mode-3
figure(3)
plot3(T3(1,:),T3(2,:),T3(3,:),'x'), grid    % mode-3 vectors
title('tensor mode-3 vectors', 'FontSize', 15)
disp('tensor mode-3 singular values:'), disp(sve{3}(1:9)')
axis([-1 1 -0.6 1 -3 3])
% figure; plot3(T3(4,:),T3(5,:),T3(6,:),'x'), grid    % mode-3 vectors
% figure; plot3(T3(7,:),T3(8,:),T3(9,:),'x'), grid    % mode-3 vectors

% Truncation - low rank approximation
Ttrunc1 = Ue{1}(:,1:3) * Ue{1}(:,1:3)' * T1;
Ttrunc1_bis = Ue{1}(:,1:2) * Ue{1}(:,1:2)' * T1;
Ttrunc2 = Ue{2}(:,1:2) * Ue{2}(:,1:2)' * T2;
Ttrunc3 = Ue{3}(:,1:3) * Ue{3}(:,1:3)' * T3;

figure(1), hold, plot3(Ttrunc1(1,:),Ttrunc1(2,:),Ttrunc1(3,:),'rx')
figure(4); hold, plot3(Ttrunc1(4,:),Ttrunc1(5,:),Ttrunc1(6,:),'rx'), grid    % mode-1 vectors
figure(5); hold, plot3(Ttrunc1(7,:),Ttrunc1(8,:),Ttrunc1(9,:),'rx'), grid    % mode-1 vectors

figure(1); plot3(Ttrunc1_bis(1,:),Ttrunc1_bis(2,:),Ttrunc1_bis(3,:),'gx')
figure(4); plot3(Ttrunc1_bis(4,:),Ttrunc1_bis(5,:),Ttrunc1_bis(6,:),'gx'), grid    % mode-1 vectors
figure(5); plot3(Ttrunc1_bis(7,:),Ttrunc1_bis(8,:),Ttrunc1_bis(9,:),'gx'), grid    % mode-1 vectors

figure(2), hold, plot3(Ttrunc2(1,:),Ttrunc2(2,:),Ttrunc2(3,:),'rx')
figure(3), hold, plot3(Ttrunc3(1,:),Ttrunc3(2,:),Ttrunc3(3,:),'rx')

%% 2.2 Source separation with ICA and PCA

clear; clc; close all

load('ex2data.mat');
[rows, columns] = size(A);

N = 500;
iterations = 100;
SNR = 0:5:50;

qual_ICA = zeros(iterations, length(SNR));
qual_PCA = zeros(iterations, length(SNR));

for it = 1:iterations
    S = zeros(rows,N);
    
    %New sources
    for j = 1:rows
        S(j,:) = -0.25 + (0.25+0.25)*rand(1,N);
    end
    
    AS = A*S;
    
    for snr_index = 1:length(SNR)
        X = noisy(AS, SNR(snr_index));
        
        %Estimate the sources S from the observations X using ICA (aci.m)
        [A_est_ICA,delta]=aci(X);
        
        %Try to separate the sources by means of PCA
        A_est_PCA = pca(X');
        
        % SIR: [sir,P,D] = sir(A,Aest). sir = quality.
        
        [sir_ICA,P_ICA,D_ICA] = sir(A,A_est_ICA);
        [sir_PCA,P_PCA,D_PCA] = sir(A,A_est_PCA);
        
        qual_ICA(it, snr_index) = sir_ICA;
        qual_PCA(it, snr_index) = sir_PCA;
    end
end



mean_qual_ICA = mean(qual_ICA);
mean_qual_PCA = mean(qual_PCA);

figure; plot(SNR, qual_ICA, '--','Color','#85C1E9');
hold on; plot(SNR, mean_qual_ICA, 'r-', 'LineWidth', 4)
title('Separation quality with ICA', 'FontSize', 16);
xlabel('SNR', 'FontSize', 13); ylabel('SIR', 'FontSize', 13)

figure; plot(SNR, qual_PCA, '--','Color','#85C1E9');
hold on; plot(SNR, mean_qual_PCA, 'r-', 'LineWidth', 4)
title('Separation quality with PCA', 'FontSize', 16)
xlabel('SNR', 'FontSize', 13); ylabel('SIR', 'FontSize', 13)

%% 2.3 Synthetic CP

clear; clc; close all

load('ex3data.mat')

%% 2.3.1 PCA
% First part: Try to estimate the components from the noiseless slice 
% T(:,:,4) by means of PCA and explain.

T_4 = T(:,:,4);

% True values of T4
figure;
subplot(3,1,1); 
plot(1:80, A(:,1), 'LineWidth', 2); hold on; xlim([1,80]);
plot(1:80, A(:,2), 'LineWidth', 2); hold on; xlim([1,80]);
plot(1:80, A(:,3), 'LineWidth', 2); hold on; xlim([1,80]);
title('Components of A', 'FontSize', 16)
subplot(3,1,2); 
plot(1:90, B(:,1), 'LineWidth', 2); hold on; xlim([1,90]); 
plot(1:90, B(:,2), 'LineWidth', 2); hold on; xlim([1,90]); 
plot(1:90, B(:,3), 'LineWidth', 2); hold on; xlim([1,90]); 
title('Components of B', 'FontSize', 16)
subplot(3,1,3); 
plot(1:4, C(:,1), 'LineWidth', 2); hold on; xlim([1,4]); 
plot(1:4, C(:,2), 'LineWidth', 2); hold on; xlim([1,4]); 
plot(1:4, C(:,3), 'LineWidth', 2); hold on; xlim([1,4]); 
title('Components of C', 'FontSize', 16)


% Estimation of T4 factor matrices
[U,S,V] = svd(T_4); %3 sv --> 3 sources --> T_: STHx3 -->
figure; stem(diag(S));

coeff_pca = pca(T_4');
mixing_mat = coeff_pca(:,1:3);
sources = T_4'*mixing_mat;

Aest = mixing_mat(:,1)*sources(:,1)';
Best = mixing_mat(:,2)*sources(:,2)';
Cest = mixing_mat(:,3)*sources(:,3)';

figure;
plot(1:80, A(:,1)', 'r', 'LineWidth', 2); hold on;
plot(1:80, A(:,2)', 'b', 'LineWidth', 2); hold on;
plot(1:80, A(:,3)', 'g', 'LineWidth', 2); hold on;

plot(1:80, Aest(:,1)', 'r--', 'LineWidth', 2); hold on;
plot(1:80, Aest(:,2)', 'b--', 'LineWidth', 2); hold on;
plot(1:80, Aest(:,3)', 'g--', 'LineWidth', 2); hold on;

frontal_slide_est = mixing_mat*sources'; % check if correct
%% 2.3.2 CPD
% Second part: Compute the CPD of the noisy tensor for different numbers of 
% rank-1 terms --> different R.

T_noisy = noisy(T,15); %snr=15;
[Ue,Se,sve] = mlsvd(T_noisy); %This can also be seen in the gui_cpd

% Plot singular values
figure;
subplot(3,1,1)
plot(1:length(sve{1}), sve{1}, '.', 'MarkerSize', 15);
title('Mode-1 singular values', 'FontSize', 16)
subplot(3,1,2)
plot(1:length(sve{2}), sve{2}, '.', 'MarkerSize', 15);
title('Mode-2 singular values', 'FontSize', 16)
subplot(3,1,3)
plot(1:length(sve{3}), sve{3}, '.', 'MarkerSize', 15);
title('Mode-3 singular values', 'FontSize', 16)

% Basis of the algorithm
% ----------------------
% R = ; U = ; T = ;
% Uest = cpd(T,R,'Display',true)
% Test = cpdgen(Uest);
% [relerrU,P,D,Uest] = cpderr(U,Uest)
% P: permutation matrix
% D: scaling matrices
% Uest: permuted and scaled factor matrices: 

U_original = [{A},{B},{C}];

clear options
% compression default mlsvd_rsi OK
options.Initialization = @cpd_gevd; %If R≤2nd largest dim of T_noisy(80) OK

options.Algorithm = @cpd_nls; %default
options.AlgorithmOptions.Display = 1;
options.AlgorithmOptions.MaxIter = 100;      % Default 200
options.AlgorithmOptions.TolFun = eps^2;     % Tol for step size. Default 1e-12
options.AlgorithmOptions.TolX = eps;         % Tol for relative change in objective function. Default 1e-6

%R=1
R = 1;
disp('-----');disp('R = 1');disp('-----');
[Uest_1, output_nls_1] = cpd(T_noisy,R,'Display',true,options);
Test_1 = cpdgen(Uest_1);
[relerrU_1,P_1,D_1,U_per_sca_1] = cpderr(U_original,Uest_1);

%R=2
R = 2;
disp('-----');disp('R = 2');disp('-----');
[Uest_2, output_nls_2] = cpd(T_noisy,R,'Display',true,options);
Test_2 = cpdgen(Uest_2);
[relerrU_2,P_2,D_2,U_per_sca_2] = cpderr(U_original,Uest_2);

%R=3
R = 3;
disp('-----');disp('R = 3');disp('-----');
[Uest_3, output_nls_3] = cpd(T_noisy,R,'Display',true,options);
Test_3 = cpdgen(Uest_3);
[relerrU_3,P_3,D_3,U_per_sca_3] = cpderr(U_original,Uest_3);

%R=4
R = 4;
disp('-----');disp('R = 4');disp('-----');
[Uest_4, output_nls_4] = cpd(T_noisy,R,'Display',true,options);
Test_4 = cpdgen(Uest_4);
[relerrU_4,P_4,D_4,U_per_sca_4] = cpderr(U_original,Uest_4);

%R=5
R = 5;
disp('-----');disp('R = 5');disp('-----');
[Uest_5, output_nls_5] = cpd(T_noisy,R,'Display',true,options);
Test_5 = cpdgen(Uest_5);
[relerrU_5,P_5,D_5,U_per_sca_5] = cpderr(U_original,Uest_5);

iterations = [20 16 7 100 100];
rel_err = [0.522241 0.266521 0.173257 0.172209 0.171028];

figure;
subplot(3,1,1); plot(1:5, iterations, '.-', 'LineWidth', 2);
subplot(3,1,2); plot(1:5, rel_err, '.-', 'LineWidth', 2);
subplot(3,1,3); plot(1:5, rel_err.*iterations, '.-', 'LineWidth', 2);
xlabel('Rank', 'FontSize', 15);

figure;
semilogy(output_nls_1.Algorithm.fval); hold all;
semilogy(output_nls_2.Algorithm.fval); hold all;
semilogy(output_nls_3.Algorithm.fval); hold all;
semilogy(output_nls_4.Algorithm.fval); hold all;
semilogy(output_nls_5.Algorithm.fval); hold all;
ylabel('Objective function'); xlabel('Iteration');
title('Convergence plot');

Aest_cpd = Uest_3{1};
Best_cpd = Uest_3{2};
Cest_cpd = Uest_3{3};

figure;
subplot(3,1,1); 
plot(1:80, A(:,1), 'b', 'LineWidth', 1); hold on;
plot(1:80, A(:,2), 'r', 'LineWidth', 1); hold on;
plot(1:80, A(:,3), 'g', 'LineWidth', 1); hold on;

plot(1:80, Aest_cpd(:,1), 'b--', 'LineWidth', 1); hold on;
plot(1:80, Aest_cpd(:,2), 'r--', 'LineWidth', 1); hold on;
plot(1:80, Aest_cpd(:,3), 'g--', 'LineWidth', 1); hold on;

xlim([1,80]);

title('Components of A', 'FontSize', 16)
subplot(3,1,2); 
plot(1:90, B(:,1), 'b', 'LineWidth', 1); hold on;
plot(1:90, B(:,2), 'r', 'LineWidth', 1); hold on;
plot(1:90, B(:,3), 'g', 'LineWidth', 1); hold on;

plot(1:90, Best_cpd(:,1), 'b--', 'LineWidth', 1); hold on;
plot(1:90, Best_cpd(:,2), 'r--', 'LineWidth', 1); hold on;
plot(1:90, Best_cpd(:,3), 'g--', 'LineWidth', 1); hold on;

xlim([1,90]);

title('Components of B', 'FontSize', 16)
subplot(3,1,3);

plot(1:4, C(:,1), 'b', 'LineWidth', 1); hold on;
plot(1:4, C(:,2), 'r', 'LineWidth', 1); hold on;
plot(1:4, C(:,3), 'g', 'LineWidth', 1); hold on;

plot(1:4, Cest_cpd(:,1), 'b--', 'LineWidth', 1); hold on;
plot(1:4, Cest_cpd(:,2), 'r--', 'LineWidth', 1); hold on;
plot(1:4, Cest_cpd(:,3), 'g--', 'LineWidth', 1); hold on;

xlim([1,4]);

title('Components of C', 'FontSize', 16)
%% 2.4 Harmonic retrieval

clear; clc; close all;

load('ex4data.mat');
fs = 100; %Hz

H_mat = hankelize(x);

[U,S,V]=svd(H_mat);

figure; stem(diag(S)); xlim([0 20])
title('Sinular values of the Hankel matrix', 'FontSize', 20);

%ESPRIT --------------------------------------------------------------------
% ----------------------------------------------------------------------------------

    clearvars -except x H_mat H_ten poles_ESPRIT poles_LMLRA poles_CPD

[Ue,Se,Ve]=svd(H_mat);
Ueup = Ue(1:end-1,:);
Uedown = Ue(2:end,:);

% U_no_down = M * U_no_up --> Uup = M*Udown;
Me = pinv(Uedown)*Ueup;
[Fe,De] = eig(Me);
poles_ESPRIT = diag(De(1:5,1:5));

theta = 0:0.02:2*pi;
figure; plot(cos(theta), sin(theta), 'k')
for p = 1:length(poles_ESPRIT)
    hold on; plot(real(poles_ESPRIT(p)), imag(poles_ESPRIT(p)), 'bx', ...
        'LineWidth', 2)
end
axis('equal');
xlabel('Real part', 'FontSize', 15); ylabel('Imaginary part', 'FontSize', 15);
title('Poles', 'FontSize', 20);

% LMLRA-based harmonic retrieval ---------------------------------------------------
% ----------------------------------------------------------------------------------

    clearvars -except x H_mat H_ten poles_ESPRIT poles_LMLRA poles_CPD

H_ten = hankelize(x,'order',3);
[Ul,Sl,svl] = mlsvd(H_ten);

% Plot singular values
figure;
subplot(3,1,1)
plot(1:length(svl{1}), svl{1}, '.', 'MarkerSize', 15);
title('Mode-1 singular values', 'FontSize', 16)
subplot(3,1,2)
plot(1:length(svl{2}), svl{2}, '.', 'MarkerSize', 15);
title('Mode-2 singular values', 'FontSize', 16)
subplot(3,1,3)
plot(1:length(svl{3}), svl{3}, '.', 'MarkerSize', 15);
title('Mode-3 singular values', 'FontSize', 16)

R = [5 5 5];
Ulest = lmlra(H_ten,R);

Ul1 = Ulest{1};

Ul1up = Ul1(1:end-1, :);
Ul1down = Ul1(2:end, :);

% U_no_down = M * U_no_up --> Uup = M*Udown;
Ml = pinv(Ul1down)*Ul1up;
[Fl,Dl] = eig(Ml);
poles_LMLRA = diag(Dl);

theta = 0:0.02:2*pi;
figure; plot(cos(theta), sin(theta), 'k')
for p = 1:length(poles_LMLRA)
    hold on; plot(real(poles_LMLRA(p)), imag(poles_LMLRA(p)), 'rx', ...
        'LineWidth', 2)
end
axis('equal');
xlabel('Real part', 'FontSize', 15); ylabel('Imaginary part', 'FontSize', 15);
title('Poles', 'FontSize', 20);

%CPD-based harmonic retrieval ---------------------------------------------------
%----------------------------------------------------------------------------------

    clearvars -except x H_mat H_ten poles_ESPRIT poles_LMLRA poles_CPD

R = 5;
size0 = [size(H_ten,1) size(H_ten, 2) size(H_ten, 2)];
Uc0 = cpd_rnd(size0, R, 'imag', @rand);
Uc = cpd(H_ten,Uc0);
Uc1 = Uc{1};

Uc1up = Uc1(1:end-1, :);
Uc1down = Uc1(2:end, :);

% U_no_down = M * U_no_up --> Uup = M*Udown;
Mc = pinv(Uc1down)*Uc1up;
[Fc,Dc] = eig(Mc);
poles_CPD = diag(Dc);

theta = 0:0.02:2*pi;
figure; plot(cos(theta), sin(theta), 'k')
for p = 1:length(poles_CPD)
    hold on; plot(real(poles_CPD(p)), imag(poles_CPD(p)), 'gx', ...
        'LineWidth', 2)
end
axis('equal');
xlabel('Real part', 'FontSize', 15); ylabel('Imaginary part', 'FontSize', 15);
title('Poles', 'FontSize', 20);

% Plot all poles together
theta = 0:0.02:2*pi;
figure; plot(cos(theta), sin(theta), 'k'); hold on;
for p = 1:length(poles_CPD)
    hold on; pe = plot(real(poles_ESPRIT(p)), imag(poles_ESPRIT(p)), 'rx', ...
        'LineWidth', 2, 'MarkerSize', 10);
    hold on; pl = plot(real(poles_LMLRA(p)), imag(poles_LMLRA(p)), 'bx', ...
        'LineWidth', 2, 'MarkerSize', 10);
    hold on; pc = plot(real(poles_CPD(p)), imag(poles_CPD(p)), 'gx', ...
        'LineWidth', 2, 'MarkerSize', 10);
end
axis('equal');
legend([pe pl pc], [{'ESPRIT'},{'LMLRA'},{'CPD'}], 'FontSize', 15)
xlabel('Real part', 'FontSize', 15); ylabel('Imaginary part', 'FontSize', 15);
title('Poles', 'FontSize', 20);

%% 2.5 EEG

load('demosignal3_963.mat');
eegplot(demosignal3_963);

fs = 250;
nbChannels = 21;
t = 52;

tstart = t-10;
tend = t+10;

[data_3D,m,s] = normalise_3D(demosignal3_963,tstart,tend,5:70);



