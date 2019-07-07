function [rcmT]=rcmT(rcm,cpm,T)
% calculates heat capacity of rocks as function of temperature
% based on Kola data

T=T(:);rcm=rcm(:); cpm=cpm(:);

rcm0=rcm./cpm;
rcmT=rcm0*1./(0.0044*T+2.134);




