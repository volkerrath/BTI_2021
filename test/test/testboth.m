clear all, close all;
L=[ 2.3    2.5    3.0    2.6   4.4     3.0];
A=[ 0.13   0.75   0.64  1.18   0.73    0.7];
B=[ 1073    705   807    474   1293    770];

T=[1:1:160];

LL1=[];
for i=1:1:length(A)
[l]=lamtl(L(i),T,A(i),B(i));
LL1=[LL1 l'];
end
LL2=[];
for i=1:1:length(A)
[l]=lamth(L(i),T,A(i),B(i));
LL2=[LL2 l'];
end

opts = struct('Format','psc2','Color','rgb','Resolution',600);
figure;
plot(T,LL1);hold on;
title('\lambda = f(T) (Lehmann, 1998)');grid on;
xlabel('T (C)');ylabel('\lambda (W K^{-1} m^{-1})')
legend('limestone','metamorphics',...
       'acidic rocks','basic rocks','ultra-basic rocks',' mean crust')
exportfig(gcf,'lambda_lehmann.eps',opts)
figure;
plot(T,LL2);hold on;
title('\lambda = f(T) (Haenel, 1988)');grid on;
xlabel('T (C)');ylabel('\lambda (W K^{-1} m^{-1})')
legend('limestone','metamorphics',...
       'acidic rocks','basic rocks','ultra-basic rocks',' mean crust')
exportfig(gcf,'lambda_haenel.eps',opts)
