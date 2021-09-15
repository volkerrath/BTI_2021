function [cpmT]=cpmT(S,T)
% calculates heat capacity of rocks as function of temperature
% based on Mottaghy(2005,2007)
T=T(:);
cpmT=(1.e3./(S(1)*T+S(2)))';