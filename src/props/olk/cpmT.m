function [cpmT]=cpmT(cpm0,T)
% calculates heat capacity of rocks as function of temperature
% based on Mottaghy(2005,2007), fitted to OLK data
% vr sep 13, 2014
A =[-0.00253, 1.46];cpbas=712;

f=cpm0./cpbas;

cpm=(1.e3./(A(1)*T+A(2)))';
cpmT=f.*cpm(:);

