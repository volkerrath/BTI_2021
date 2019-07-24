function [kmT]=kmT(k,T,kA,kB,P)
% set thermal conductivity as function of temperature (dummy)
% v. r. oct. 23, 2001

if nargin < 2, T=20.;end
T=T(:); P=P(:);k=k(:); kA=kA(:);kB=kB(:);

kmT=k.*ones(size(T));

