function [rcmT]=rcmT(rcm,cpm,T)
% calculates heat capacity of rocks as function of temperature
% rm ist rho von z

global rm

T=T(:);rm=rm(:);  cpm=cpm(:); 

rcmT=rm.*cpmT(cpm,T);





