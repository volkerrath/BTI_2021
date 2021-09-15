function [rcmT]=rcmT(rcm,cpm,T,P)
% calculates heat capacity of rocks as function of temperature

T=T(:);rcm=rcm(:); 

rcmT=rcm.*cpmT(cpm,T)./cpm;




