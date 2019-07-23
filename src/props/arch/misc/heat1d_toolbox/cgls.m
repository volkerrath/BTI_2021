function [x,rho,eta,diagT,diag1T,xout] = cgls(A,b,tol,maxit,debug)
% CGLS applies conjugate gradient algorithm 
% implicitly to the normal equations
%
% [x,rho,eta,diagT,diag1T] = cgls(A,b,maxit,tol,debug)
% performs maxit steps of the conjugate gradient algorithm applied
% implicitly to the normal equations A'*A*x = A'*b.
% Iteration stops when either maxit steps have been made or the residual
% has been reduced by a factor of tol.
%
% The solution norm and residual norm are returned
% in eta and rho, respectively.
%
% If the argunents out are 5 then it calculates the triagonal matrix T such that
% T=R'AR, the eigenvalues of T are the Ritz values and can be used as 
% good approximations for the eigenvalues of A

% References: A. Bjorck, "Least Squares Methods", in P. G.
% Ciarlet & J. L Lions (Eds.), "Handbook of  Numerical Analysis,
% Vol. I", Elsevier, Amsterdam, 1990; p. 560.
% C. R. Vogel, "Solving ill-conditioned linear systems using the
% conjugate gradient method", Report, Dept. of Mathematical
% Sciences, Montana State University, 1987.

% based on cgls routines by Per Christian Hansen & Eldad Haber
% last modification  V. R.  October 3, 2006 
reorth=1;
[m,n] = size(A);
t0=cputime;

% Initialization.

if (nargin < 3), maxit = 16; end; if (maxit < 1), error('Number of steps <maxit> must be positive'), end
if (nargin < 4), tol=1.e-6; end; 
if (nargin < 5), debug=0; end; 
if (reorth==1), ATr = zeros(n,maxit); end
if (nargout > 1) eta = zeros(maxit,1); rho = eta;end
if (nargout > 3), Alfa=zeros(maxit,1);Beta=Alfa;diagT=Alfa;diag1T=zeros(maxit-1,1);end


if debug>0,
    disp([' CGLS iteration: ' ])
    disp(['   stop tol:   ', num2str(tol) ])
    disp(['   max iter:   ', num2str(maxit) ])
%    disp(['   reorth:     ', num2str(reorth) ])
end

% Prepare for CG iteration.
x = zeros(n,1);
d = (b'*A)';    % faster than A'*b;
r = b;
normr2 = d'*d;

if (reorth==1), ATr(:,1) = d/norm(d); end

% Iterate.
for j=1:maxit
  % Update x and r vectors.
  Ad = A*d; 
  alpha = normr2/(Ad'*Ad);
  x  = x + alpha*d;
  r  = r - alpha*Ad;
  s  =  (r'*A)';   % faster than A'*r; 
  
  % Reorthogonalize s to previous s-vectors, if required.
  if (reorth==1)
    for i=1:j-1, s = s - (ATr(:,i)'*s)*ATr(:,i); end
    ATr(:,j) = s/norm(s);
  end

  % Update d vector. 
  normr2_new = s'*s;
  beta = normr2_new/normr2;
  normr2 = normr2_new;
  d = s + beta*d;
  
  if (nargout>1), rho(j) = norm(r); end
  if (nargout>2), eta(j) = norm(x); end
  if (nargout>3), Alfa(j)=alpha;Beta(j)=beta;end;
  if (nargout==6),xout(:,j) = x;end;
  
  if (debug>1),
   if (j==1),
     normres=norm(r);normres0=normres;
     disp([ 'CGLS residual norm for iteration ',...
              num2str(j),' = ',num2str(normres)])
    else
      normres=norm(r);
      test=normres/normres0;
      disp([ 'CGLS residual norm for iteration ',...
      num2str(j),' = ',num2str(normres),'  rfact = ', num2str(test)])
      if (test <= tol),break;end
    end;
  end 
end

if(nargout>3),
   for i=2:maxit,   diagT(i)=1/Alfa(i)+Beta(i)/Alfa(i-1);end; 
   diagT(1)=1/Alfa(1);
   for i=1:maxit-1, diag1T(i)=-Beta(i+1)^0.5/Alfa(i); end;
end;

if (debug>0),
  disp(['   ',num2str(j) ' CGLS iteration used ' num2str(cputime-t0) 's' ])  
end


