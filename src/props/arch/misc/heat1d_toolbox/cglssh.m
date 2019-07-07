function [x,nrm,active] = cglssh(A,b,tau,maxit,thresh,debug)
%CGLS Conjugate gradient algorithm applied implicitly to the normal equations.
%
% [X,rho,eta] = cgls(A,b,tau,maxit,thresh)
%
% Performs maxit steps of the conjugate gradient algorithm applied
% implicitly to the normal equations (A'*A+ lambda*I)*x = A'*b for nl
% values of lambda.
%
% The routine returns all nl solutions, stored as columns of
% the matrix X.  The corresponding solution and residual norms
% are returned in the vectors eta and rho, respectively.
%
%
% References: 
% J. van den Eshof & G L. G. Sleijpen, "Accuate Conjugate 
%       Gradient methods for shifted systems", Tech. Rep. 1265, 2003
% A. Bjorck, "Numerical Methods for Least Squares Problems",
%       SIAM, Philadelphia, 1996.
% A. Frommer and P. Maass: Fast CG-based Methods for Tikhonov-Phillips 
%       regularization, Bergische Universität GH Wuppertal,
%       Fachbereich Mathematik,Preprint BUGHW-SC 96-10  
%
% Volker Rath   RWTH   May 19, 2004 
 
% Initialization.

[m,n] = size(A);
nl = length(tau);

if (nargin < 4), maxit = 16; end; if (maxit < 1), error('Number of steps <maxit> must be positive'), end
if (nargin < 5), thresh=1.e-6; end; 
if (nargin < 6), debug=0; end; 


% Prepare for CG iteration.

z = b;
r =(z'*A)';    % A'*z;
x = zeros(n,nl);p=repmat(r,[1,nl]);
u=r; 

rho = norm(r);
gamma=ones(1,nl); 
t=tau;
tol=thresh*rho;rho=rho*rho;

k=1;
nrm(k,:)=sqrt(rho)./abs(gamma);
active=1:nl;

while k<=maxit,
 y=A*u; 
 sigma=y'*y;
 alpha=rho/sigma;
 z=z-alpha*y;
 r =(z'*A)';    % A'*z;
 sigma=rho;
 rho=r'*r;
 beta=rho/sigma;
 u=r+beta*u;
 
 for j=active
     l=1+alpha*t(j);
     t(j)=tau(j)+(beta/l)*t(j);
     gamma(j)=l*gamma(j);
     x(:,j)=x(:,j)+p(:,j)*alpha/gamma(j);
     p(:,j)=p(:,j)*beta/l + r;
 end
 
 normr=sqrt(rho)./abs(gamma(active));
 nrm(k+1,:)=nrm(k,:);
 
 k=k+1;
 
 nrm(k,active)=normr;
 
 active=find(rho>tol);
 if isempty(active), break; end
 
end

