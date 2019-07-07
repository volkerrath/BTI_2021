clear all, close all;
A=[ 0.13   0.75   0.64  1.18   0.73    0.7];
B=[ 1073    705   807    474   1293    770];

T=[10:1:300];

figure;
ll=[];
for i=1:1:length(A)
[f,l]=lamt(T,10.,A(i),B(i));
LL=[LL l];
end

plot1(T,LL);hold on;
title('\lambda = f(T)');grid on;
xlabel('T (C)');ylabel('\lambda (W K^{-1} m^{-1})')
legend('limestone','metamorphics',...
       'acidic rocks','basic rocks','ultra-basic rocks')
