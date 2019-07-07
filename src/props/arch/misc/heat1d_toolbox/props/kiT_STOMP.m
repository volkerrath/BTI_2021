function ki=kiT(T)
% calculate ice thermal conductivity [W/(m*K)]
% if nargin <2, nonlinear='yes'; end;
[n1 n2]=size(T);if n1==1, T=T'; end 
T2=T.*T;
%ki=7.39519-2.86936e-2*T+3.54452e-5*T2;
ki=0.39519-2.86936e-2*T+3.54452e-5*T2;