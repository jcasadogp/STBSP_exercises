function R = calcCovMatrix2(data,segments, L)
%calcNoiseCovMatrix This function calculates the noise covariance matix
%   INPUTS:
%       data
%       noiseSegments
%   OUTPUT:
%       noiseCovMatrix: lagged noise covariance matrix

nbChannels=size(data, 2);

for i = segments
    data_segments = data(i,:);
end

matrix=[];
for k=1:nbChannels
    fprintf('<strong>k = %d</strong>\n',k);
    
    toe_col=data_segments(:,k);
    toe_row=zeros(1,L);
    toe_row(1,1)=toe_col(1,1);
    
    mat_k=toeplitz(toe_col,toe_row);
    matrix=[matrix,mat_k];
end

fprintf('<strong>Ya ha terminado todos los canales.</strong>\n');
R = cov(matrix); %Matrix is 475x1195145
end

