function [ output ] = applyMultiChannelFilter( data, filter )
%APPLYMULTICHANNELFILTER Calculate multi-channel FIR filter output on given data.
%
% inputs:
%   data : multi-channel data (numberSamples, numberChannels)
%   filter : multi-channel FIR filter (numberTaps, numberChannels)
%
% output:
%   output : filter output
%

% extract the temporal window length from the give matrix filter
% coefficients
L = size(filter, 1);
% initialize the filter output
output = zeros(size(data, 1), 1);

nbChannels = size(data, 2);

% calculate output
for idx=L:length(output)
    
    % TODO: Calculate the filter output(idx). At time idx consider only idx
    % and its past samples (i.e., we implement a causal filter).
    
    data_idx = data([idx-L+1:idx],:);
    
    for m=1:nbChannels
        data_m = data_idx(:,m);
        filter_m = filter(:,m);
        prod_m = dot(data_m, filter_m);
        output(idx) = output(idx)+prod_m;
    end
end

end

