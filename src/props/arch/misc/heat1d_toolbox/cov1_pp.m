function [Cpp,sqCpp]=cov1_pp(p,par,typ)
% COV1_PP calculates a gaussian covariance a-priori matrix
% 
% [Cpp,sqCpp]=cov1_pp(p,par,typ) calculates
% covariance matrix (and optionally its square root)
% of type type, where typ for now can only take the value
% 'gauss'. in this case:
%        p   is a vector of np 1D coordinates
%        par consists of sigma sig and correlation length L
% V. R., Oct. 15, 2001 

if nargin < 3, typ='gauss'; end

np=length(p);

switch typ
  case('gauss')
     sig=par(1); L=par(2);
      for i=1:1:np
         for j=1:1:np
           Cpp(j,i)=sig^2*exp(-0.5*((p(i)-p(j))/(L))^2);
         end 
      end
   otherwise
      error(['typ ',typ,' not defined']);
      break;
   end
if nargout> 1, sqCpp=sqrtm(Cpp); end
