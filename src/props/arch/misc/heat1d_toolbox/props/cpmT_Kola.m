function [cpmT]=cpmT(cp,T)
% calculates heat capacity of rocks as function of temperature
% based on Kola


[n1 n2]=size(T);if n1==1, T=T'; end 

cpmT=1.e3./(0.0044*T+2.134);


