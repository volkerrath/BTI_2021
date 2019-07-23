function [rcmT]=rcmT(rcm,cpm,T)
% calculates heat capacity of rocks as function of temperature
% rm ist rho von z

global rm

T=T(:);rm=rm(:); 

rcmT=rm.*cpmT(cpm,T);




