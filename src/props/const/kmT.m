function [kmT]=kmT(k,T,kA,kB,P)
% set thermal conductivity as function of temperature (dummy)
% last change: vr,  Sep 2019
if nargin < 4, P =  0.;end
if nargin < 2, T = 20.;end
[n1,n2]=size(T); if n2~=1, T=T'; end 
[n1,n2]=size(k); if n2~=1, k=k'; end
kmT=k.*ones(size(T));

