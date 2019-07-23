function [kmT]=kmT(km0,T)
% set thermal conductivity as function of temperature (dummy)
% v. r. oct. 23, 2001
if nargin < 2, T=20.;
[n1,n2]=size(P); if n1~=1, T=T'; end 
kmT=km0;
