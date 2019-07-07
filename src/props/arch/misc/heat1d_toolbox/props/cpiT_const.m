function [cpi]=cpiT(T)
% calculate ice isobaric heat capacity [J/(kg*K]
% as function of temperature (dummy)
% VR RWTH Aachen University,   April 25, 2004
[n1,n2]=size(P); if n1~=1, T=T'; end 
cpi=2110*ones(size(T));
