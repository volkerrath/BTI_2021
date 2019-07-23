function [cpmT]=cpmT(cp,T,T0,T1)
% calculates heat capacity of rocks as function of temperature
% based on Vosteen and Schellschmidt 2003 - figure 4
% T           =  temperature in C  acording to flag unit
% vr april 25, 2004
if nargin < 3, T0=25.;T1=800; end
if nargin < 3, T0=25.;T1=800; end

[n1 n2]=size(T);if n1==1, T=T'; end 
[n1 n2]=size(cp);if n1==1, cp=cp'; end

 a = -0.0015;
 b = 1.0367;
 Ts = T-T0;         % T ist in situ, T0 Labor (25 Grad)
 
 cpmT = a*Ts.*Ts. + b*Ts+cp;
 
