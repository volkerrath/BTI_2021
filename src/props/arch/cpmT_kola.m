function [cpmT]=cpmT(cp,T)
% calculates heat capacity of rocks as function of temperature
% based on Mottaghy(2005,2007)
T=T(:);
cpmT=1.e3./(0.0044*T+2.134);


