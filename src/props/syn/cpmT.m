function [cpmT]=cpmT(cpm,T,P,T0,T1)
% calculates heat capacity of rocks as function of temperature
% based on Vosteen and Schellschmidt 2003 - figure 4
% T           =  temperature in C  acording to flag unit
% vr april 25, 2004
if nargin < 4, T0=25.;T1=800; end

T=T(:);
cpm=cpm(:);

 cpmT = cpm;
 
