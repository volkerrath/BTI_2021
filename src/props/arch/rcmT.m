function [rcmT]=rcmT(rcm,T)
% calculates heat capacity of rocks as function of temperature
% based on Kola data

T=T(:);rcm=rcm(:); 

rcmT=rcm.*1.e3./(0.0044*T+2.134);




