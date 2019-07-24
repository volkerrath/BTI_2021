function [rcmT]=rcmT(rm,cpm,T,P)
% calculates heat capacity of rocks as function of temperature
% rm ist rho von z

T=T(:);rm=rm(:);  cpm=cpm(:); 

rcmT=rm.*cpmT(cpm,T,P);





