function [x,rho,eta,diagT,diag1T,xall] = ...
                         cgls(A,b,maxit,tol,stoprule,debug)
% CGLS applies conjugate gradient algorithm 
% implicitly to the normal equations
%
% [x,rho,eta,diagT,diag1T] = cgls(A,b,maxit,tol,stoprule,debug)
% performs maxit steps of the conjugate gradient algorithm applied
% implicitly to the normal equations A'*A*x = A'*b.
% Iteration stops when either maxit steps have been made (stoprule=0),  the residual
% has been reduced by a factor of tol (stoprule=1) or ACB limit is reached 
% (stoprule=2).
%
% The solution norm and residual norm are returned
% in eta and rho, respectively.
%
% If the argunents out are 5 then it calculates the triagonal matrix T such that
% T=R'AR, the eigenvalues of T are the Ritz values and can be used as 
% good approximations for the eigenvalues of A

% References: 
% A. Bjorck, "Least Squares Methods", in P. G.
%     Ciarlet & J. L Lions (Eds.), "Handbook of  Numerical Analysis,
%     Vol. I", Elsevier, Amsterdam, 1990; p. 560.
% C. R. Vogel, "Solving ill-conditioned linear systems using the
%     conjugate gradient method", Report, Dept. of Mathematical
%     Sciences, Montana State University, 1987.
% A. C. Berglund,"Nonlinear Regularization -wirh applicatios in geophysics"
%     PhD thesis DTU, IMM LYNGBY,2002
%
% based on cgls routines by Per Christian Hansen & Eldad Haber
% last modification  V. R.  October 3, 2006 

reorth=1;
[m,n] = size(A);
t0=cputime;

% Initialization.
if (nargin<6), debug=0; end; 
if (nargin<5), stoprule=0; tol=0.0001; end; if (stoprule==2), tol=1.; end
if (nargin<4), tol=1.e-6; end
if (nargin<3), maxit=16; end
if (maxit < 1), error('Number of steps maxit must be positive'), end
if (reorth==1), ATr = zeros(n,maxit); end

if (nargout>1), eta = zeros(maxit,1); rho = eta;end
if (nargout>3), Alfa=zeros(maxit,1);Beta=Alfa;diagT=Alfa;diag1T=zeros(maxit-1,1);end

if debug>0,
    disp([' CGLS iteration: ' ])
    disp(['   stop tol:   ', num2str(tol) ])
    disp(['   max iter:   ', num2str(maxit) ])
    disp(['   stop rule:  ', num2str(stoprule) ])
%    disp(['   reorth:     ', num2str(reorth) ])
end

aprs=zeros(maxit,1);fourc=zeros(maxit,1);               % ACB

% Prepare for CG iteration.
x = zeros(n,1);
d = (b'*A)';    % faster than A'*b;
r = b;
normr2 = d'*d;

if (reorth==1), ATr(:,1) = d/norm(d); end

% begin of iteration

for j=1:maxit
    
    % Update x and r vectors.
    Ad = A*d;                       
    alpha = normr2/(Ad'*Ad);  
    aprs(j) = alpha;      
    x  = x + alpha*d;
    rold = r; 
    r  = r - alpha*Ad;
    fourc(j) = abs(norm(r - rold)/sqrt(alpha));
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
    
    
    
    xall(:,j) = x;
    
    normres=norm(r);
    if (j==1),
        normres0=normres;test=0;
        %disp([ '    starting CGLS iterations '])
    else 
        test=normres/normres0;
        if(stoprule==1),
            % norm/norm0
            
            if (debug>1),
                disp([ '      CGLS residual norm reduction for iteration ',...
                        num2str(j),' = ',num2str(normres),'  norm/normold = ', num2str(test)])
            end
            if (test <= tol),break;end
            
        else
            % ACB
            if(j>2) 
	           testACB= (fourc(j)*sqrt(aprs(j)))/(fourc(j-1)*sqrt(aprs(j-1)));
               if (debug>1),
                  disp([ '      CGLS ACB-value for iteration ', num2str(j),' = ', num2str(testACB)...
                        '   norm(r) = ',num2str(normres) '   norm(r)/norm = ',num2str(test)])
               end
               if (testACB >= tol),break;end
            end
        end 
    end
    
end;     % end of iteration

if (debug>0),
  disp(['   ',num2str(j) ' CGLS iteration used ' num2str(cputime-t0) 's' ])  
end


if(nargout>3),
    for i=2:maxit,   diagT(i)=1/Alfa(i)+Beta(i)/Alfa(i-1);end; 
    diagT(1)=1/Alfa(1);
    for i=1:maxit-1, diag1T(i)=-Beta(i+1)^0.5/Alfa(i); end;
end;

