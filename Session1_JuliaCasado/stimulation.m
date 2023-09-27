%% INITIALIZE (section1)
disp('Section 1: INITIALIZE')

close all;
clear; clc;

% load the training data
load('artifactData.mat');

lb = labels(1,1);
lb2 = labels(2,1);

maximo = max(max(data));
minimo = min(min(data));

v1 = [180 minimo; 180 maximo; 220 maximo; 220 minimo];
v2 = [379 minimo; 379 maximo; 419 maximo; 419 minimo];
faces = [1 2 3 4];

figure; plot(data(4,(lb-200:lb2+200)));
xline(200, 'r', {'Neural Spike'}, 'LineWidth', 2, 'FontSize', 20, 'LabelOrientation', 'horizontal');
xline(200+(lb2-lb), 'r', {'Neural Spike'}, 'LineWidth', 2, 'FontSize', 20, 'LabelOrientation', 'horizontal');
title('Channel 4 - Neural recording with stimulation artifacts', 'FontSize', 20);

figure; plot(data(4,(lb-200:lb2+200)));
patch('Faces',faces,'Vertices',v1,'FaceColor','red', 'FaceAlpha',.2)
patch('Faces',faces,'Vertices',v2,'FaceColor','red', 'FaceAlpha',.2)
axis tight
title('Channel 4 - Neural recording with stimulation artifacts', 'FontSize', 20);

%% BUILD ARTIFACT DATA MATRIX (section2) 

nbChannels = size(data,1);

artifact_times = [];
for idx=1:2:length(events)
    start_art_t = events(idx);
    stop_art_t = events(idx+1);
    times = [start_art_t:stop_art_t];
    artifact_times = [artifact_times, times];
end
clearvars start_art_t stop_art_t idx times;

%% A\b (section3) --> with real data

L = 3;
channelsVec = 1:nbChannels;
header = zeros(1, L-1); %For later

tm = size(data,2);

T=zeros(1,nbChannels);

cleanData=zeros(32, 300000);

for k = 1:nbChannels
    tic
    fprintf('<strong>k = %d</strong>\n',k);
    
    f = waitbar(0, 'Initiating...', 'Name', 'Calculating X-k matrix.');
    
    % The b vector on the LS: A\b --> data of channel k.
    b = data(k,:);
    
    % Find the no-adjacent channels for the k-th channel.
    no_adj = setdiff(channelsVec, adjacent_channels{k});

        fprintf('\n\tnbChannels = %d\n',nbChannels); 
        fprintf('\tlength(no_adj) = %d\n',length(no_adj));
        fprintf('\tL = %d\n',L);
        fprintf('\tnbChannels*L = %d\n\n',nbChannels*L);
    
    % The A matrix on the LS: A\b --> matrix X-k of the assignment (7).
    A = zeros(tm,L*length(no_adj));
    Z = zeros(1,L*length(no_adj));
    
    % xk = zeros(1, L*nbChannels); %Each row of A.
    for t=1:tm
        % t - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        resto = mod(t,100);
        if resto == 0
            waitbar(t/tm, f, sprintf('CHANNEL = %d', k));
        end
        
        if ismember(t, artifact_times)
            % xk = [];
            xk = zeros(1,L*length(no_adj));
            i = 1;
            for ch=1:nbChannels
                if ismember(ch, no_adj)
                    if t<L
                        amplified_data_ch = [header, data(ch, :)]; %Zeros+data
                        cell = amplified_data_ch(:,(t:t+L-1));
                    else
                        cell = data(ch,(t-L+1:t));
                    end
                    %xk = [xk, cell];
                    xk(L*i-L+1:L*i) = cell;
                    i=i+1;
                end
            end
            A(t,:) = xk;
        else
            A(t,:) = Z;
        end
        
        % t - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    end
    
    % LS problem for k-th channel.
    [X_row, X_col]=size(A);
    [data_row, data_col]=size(b');
    fprintf('X-k:\t\t%d x %d\nData :\t%d x %d\n',X_row, X_col, data_row, data_col);
    clear X_row; clear X_col; clear data_row; clear data_col;
    
    wk=A\b';

    [wk_row, wk_col]=size(wk);
    fprintf('Wk:<strong>\t%d x %d</strong>\n',wk_row, wk_col);
    fprintf('----------------------------\n'); fprintf('\n'); fprintf('\n'); fprintf('\n');fprintf('----------------------------\n');
    clear wk_row; clear wk_col;
    
    delete(f); clear f;
        
    T(k)=toc;
    
    art = A*wk;
    clean_data = b'-art;
    
    cleanData(k,:)=clean_data;
end

%% Artifact of the channel 4

load('ch4_art.mat')
figure; plot(art(lb-200:lb2+200)); axis tight;
patch('Faces',faces,'Vertices',v1,'FaceColor','red', 'FaceAlpha',.2)
patch('Faces',faces,'Vertices',v2,'FaceColor','red', 'FaceAlpha',.2)
title('Channel 4 - Artifact', 'FontSize', 20);

%% Plot clean sigal

lb = labels(1,1);
lb2 = labels(2,1);

maximo = max(max(cleanData));
minimo = min(min(cleanData));

v1 = [180 minimo; 180 maximo; 220 maximo; 220 minimo];
v2 = [379 minimo; 379 maximo; 419 maximo; 419 minimo];

figure; plot(cleanData(4,(lb-200:lb2+200))); axis tight;
patch('Faces',faces,'Vertices',v1,'FaceColor','red', 'FaceAlpha',.2)
patch('Faces',faces,'Vertices',v2,'FaceColor','red', 'FaceAlpha',.2)
title('Channel 4 - Neural recording without stimulation artifacts', 'FontSize', 20);

clear lb; clear lb2;
