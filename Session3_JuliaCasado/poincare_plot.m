function [SD1,SD2,SDRR] = pointcare_plot(RR)
xp = RR;
xp(end) = [];
xm = RR;
xm(1) = [];
%SD1
SD1 = std(xp-xm)/sqrt(2);
%set(handles.sd1, 'string', num2str(SD1*1000));
%SD2
SD2 = std(xp+xm)/sqrt(2);
%SDRR
SDRR=sqrt(SD1^2+SD2^2)/sqrt(2);
%Poincare plot
plot(xm,xp,'.')
title('Poincaré plot')

