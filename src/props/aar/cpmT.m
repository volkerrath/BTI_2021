function [cpmT]=cpmT(cpm,T,P,T0,T1)
% calculates heat capacity of rocks as function of temperature
% based on Vosteen and Schellschmidt 2003 - figure 4
% T           =  temperature in C  acording to flag unit
% vr april 25, 2004
if nargin < 4, T0=25.;T1=800; end

T=T(:);
cpmT=cpm(:);

%  a = -0.0015*ones(size(T));
%  b = 1.0367*ones(size(T));
%  Ts = T-T0;         % T ist in situ, T0 Labor (25 Grad)
%  
%  cpmT = a.*Ts.*Ts + b.*Ts + cpm;
 
