function [cpmT]=cpmT(cp,T,T0,T1,A,B,C)
% calculates heat capacity of rocks 
% as function of temperature (dummy)
% VR RWTH Aachen University,   April 25, 2004
[n1,n2]=size(P); if n1~=1, T=T'; end 
cpmT=cp;
