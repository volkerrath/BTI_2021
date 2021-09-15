function [C]= Ce(s,L,thresh)
% calculates markovian (exponential) covariance C
if nargin<3,thresh=0;end 
n=length(s);

s2=s(:).*s(:);

V=[0:n-1];
A=toeplitz(V);
C=exp(-abs(A)/L);

C=C*spdiags(s2,0,n,n);


if abs(thresh)>0, S=thresh*eps*ones(n,n);C=C+S; end
