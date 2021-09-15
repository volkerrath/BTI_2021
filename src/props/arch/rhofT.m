function [rho]=rhofT(T,P)
% [rho]=rhofT(T,P) calculate the density i(in kg/m^3) of pure water, 
% as function of temperature (dummy)
% VR RWTH Aachen University,   April 25, 2004
T=T(:); 
if nargin<2, P=0.1*ones(size(T)); end
rho=999.84*ones(size(T));
