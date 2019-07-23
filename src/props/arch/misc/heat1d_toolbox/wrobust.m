function [weights]=wrobust(res,method,param)
% WROBUST calculates IRLS weights for robust regression 
% 
% [weights]=wrobust(res,method,param) calculates weights 
% based on the formulas of Huber, Hampel, Tukey, or  for l1 estimation.
%
% Arguments:
%               res:    vector of residuals
%               method: type of funtion ('lsq','l1',huber','hampel','tukey')
%               param:  parameter of the chosen function;
%                       default value depends on method.
% Scale is calculated as  robust estimate MAD
%
% V. R. 8/01                 last change: Aug.18,2001 
% (using code from P.B. Stark  stark@stat.berkeley.edu)

if nargin==3
  if (param >= 0)
      error('parameter must be positive');
  end
end

mad=median(res);  
scale=1.483*mad;


switch lower(method)
    
  case 'lsq'
    % dummy least squares weights      
    weights=ones(size(res));
    
    
  case 'l1'
    % computes the weight function for minimum l_1 regression
    % for iteratively reweighted least squares: scale/|res|         
    % no param. 
    thresh = eps^(1/3); 
    unit = ones(size(res))*thresh;
    res(abs(res) <= thresh ) = unit(abs(res) <= thresh);
    weights = (abs(res)).^(-1);     
    

  case 'huber'
     % Huber's weight function for robust regression:
     % min( 1, param/(|res|/scale) )    
     p = 2.5;               
     if (nargin == 3), 
        p = param;
     end
     unit = ones(size(res));
     res(res == 0 ) = p*scale*unit(res==0); 
     weights = min(unit, scale*p*(abs(res)).^(-1));
     
     
  case 'hampel'
     % campel's weight function for robust regression:
     % hampel(x) =   { 1,                     |x| < param(1),
     %               { a/|x|,                 a <= |x| < b,
     %               { a/|x|*(c-|x|)/(c-b),   b <= |x| <  c,
     %               { 0,                     |x| >= c.
     % must have 0 <= a <= b <= c. 
     % param is a vector [a, b, c]. A standard choice is [2 4 8], which is the
     % default. The resulting influence function is nearly identical to the 
     % tukey (biweight) function with  parameter 8.
     if ((param(1) < 0) | (param(2) < param(1)) | (param(3) < param(2))),
        error([' illegal choice of parameters in Hampel: ' ...
                num2str(param) ]')
     end

     if (nargin == 3),
        a = param(1);b = param(2);c = param(3);
     else
        a = 2; b = 4;c = 8;
     end
     weights = ones(size(res));
     weights(abs(res) >= a) = a*(abs(res(abs(res) >= a))).^(-1);
     weights(abs(res) >= b) = weights(abs(res) >= b).*(-abs(res(abs(res) >= b)) +c) ...
                 /(c-b);
    z = zeros(size(res));
    weights(abs(res) > c) = z(abs(res) > c);
case 'tukey'
    % Tukey's biweightseight function for robust regression
    % biweightsght(x) =   (1-x^2/param^2)^2 , |x| < param; 0, |x| >= param.
    p = 8;
    if (nargin == 3),
        p = param;
    end
    res = res/scale;
    weights = (ones(size(res)) - res.^2/p^2).^2 .* (abs(res) < p);
otherwise
   % dummy least squares weights      
   weights=ones(size(res));

end
