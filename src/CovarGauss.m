
function [C]= Cg(s,L,thresh)
% calculates gaussian covariance C
if nargin<3,thresh=1e3*eps;end 
n=length(s);

s2=s(:).*s(:);L2=2*L^2;

V=[0:n-1];
A=toeplitz(V);
C=exp(-abs(A).^2./L2);
C=C*spdiags(s2,0,n,n);
[A,D] = eig(C); 
C = A * diag(max(diag(D),thresh)) / A;

