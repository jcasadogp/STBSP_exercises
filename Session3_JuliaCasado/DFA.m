function [a1, a2, Fall] = DFA(RR)

K = length(RR);
y = [];
y_ = [];
Fall = [];

winlength = [4  : 2 : 70];    %matrix which contains the length of the different windows

for i = 1 : length(winlength())
    n = winlength(i);
    
    %Step 1 - calculation of the integrated time series
    m = mean(RR);
    for k = 1 : K
        sum_ = 0;
        for i = 1 : k
            sum_ = sum_ + RR(i);
        end
        y = [y, sum_- m];
    end
    
    %Step 2- division of the integrated signal into windows of length n
    newy = [];
    x = [1 : n];
    for i = 1 : (K/n)
        y_ = y(((i-1)*n+1) : (i* n)) ;
        %calculation of a least-squares line of the data representing the trend in our window
        coefs = polyfit(x, y_,1);
        z = coefs(1)*x+coefs(2);
        %plot(x,z);
        newy = [newy,z];
    end
    
    %Step 3 - calculation of the root-mean-square fluctuation
    summ = 0;
    for i = 1 : length(newy)
        summ = summ + (y(i)-newy(i))^2;
        
    end
    
    F1 = sqrt(1/K*summ);
    Fall =[Fall,F1];
end

%Step 4- Plot of the logF-logn Diagram
x = log(winlength);
y = log(Fall);
figure
plot(log10(winlength),log10(Fall));
title('logF - logn Diagram' )
xlabel('Length of window n ')
ylabel('Root-mean-square fluctuation F')

%Step 5 - Calculation of the short-term fractal scaling exponent α1 and
%the long-term scaling exponent α2
x1 = log10(4); x2 = log10(16); x3 = log10(64);
a1 = (y(7)-y(1)) / (x2-x1);
a2 = (y(31)-y(7)) / (x3-x2);
