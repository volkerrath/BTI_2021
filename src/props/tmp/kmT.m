function [kmT]=kmT(k,T,kA,kB)
% set thermal conductivity as function of temperature (dummy)
% v. r. oct. 23, 2001
global P

if nargin < 2, T=20.;end
[n1,n2]=size(T); if n2~=1, T=T'; end 
[n1,n2]=size(k); if n2~=1, k=k'; end
kmT=k.*ones(size(T));

