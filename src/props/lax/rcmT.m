function [rcmT]=rcmT(rcm,cpm,T)
% calculates heat capacity of rocks as function of temperature

T=T(:);rcm=rcm(:); cpm=cpm(:);

rcmT=rcm.*cpmT(cpm,T)./cpm;



