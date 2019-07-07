function [rho]=rhofT(T,P)
% [rho]=rhofT(T,P) calculate the density i(in kg/m^3) of pure water, 
% as function of temperature (dummy)
% VR RWTH Aachen University,   April 25, 2004
[n1,n2]=size(P); if n1~=1, T=T'; end 
rho=999.84*ones(size(T));
