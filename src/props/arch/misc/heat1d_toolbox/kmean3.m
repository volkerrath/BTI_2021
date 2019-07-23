 function [keff]=keff3(por1,k1,por2,k2,por3,k3,method)
% calculates the mean vector  k of a three-phase medium by various methods
% methods:
%  'arithmetic','ari','a'      arithmetic mean 
%  'geometric','geo','g'       geometric mean (default)
%  'harmonic','har','h'        harmonic mean
%  'sqrmean','sqr','s'         square root mean
%
%  vr,     April 19, 2004

 if nargin < 7, method='g', end;
 
 switch lower(method)
        case {'arithmetic','ari','a'}
            keff  = por1.*k1+  por2.*k2 + por3.*k3;
        case {'geometric','geo','g'}
            keff= exp(log(k1).*por1+log(k2).*por2+log(k3).*por3);
        case {'harmonic','har','h'}
            keff= 1./ (por1./kf +por2./ki + por3./k3);
        case {'sqrmean','sqr','s'}
            keff=(por1*.sqrt(k1)+por2*.sqrt(k2)+por3.*sqrt(k3)).^2;
        otherwise
            disp(['KMEAN3: mode set to arithmetic, >', mean, '<  not defined'])
            keff  = por1.*k1 +  por2.*ki +  por3.*k3;
    end
    
