function [kmT]=kmT(k,T,kA,kB)
% set thermal conductivity as function of temperature (dummy)
% v. r. oct. 23, 2001
T=T(:); k=k(:);
kmT=k.*ones(size(T));

