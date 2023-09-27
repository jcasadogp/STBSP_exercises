function [sl] = slope(RR)
Fs = 0.02;
[Pxx,f] = pwelch(RR,500,300,500,Fs); % [Pxx,F] = pwelch(X,WINDOW,NOVERLAP,F,Fs)
plot(log10(f),log10(Pxx))
xlabel('Frequency ')
ylabel('PSD ')
title ('PSD in log - log scale')

ax = gca;
h = findobj(gca,'Type','line');
x = h.XData;
y = h.YData;
%I will keep only the range -4 : -2
x_ = x(3:length(x));
y_ = y(3:length(y));

hold on; plot(x_,y_,'r');
%slope of th regression line slope of the log(power) versus logf relation
%in the 10^4 to 10^−2 Hz frequency range
%sl = (y_(length(y_))- y_(1)) /( x_(length(x_))- x_(1));
coefs = polyfit(x_, y_,1);
sl = coefs(1);
hold on
z = coefs(1)*x_+coefs(2);
plot(x_,z);
