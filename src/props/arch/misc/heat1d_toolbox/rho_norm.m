function [T]=rho_norm(v,par,flag)
% [T]=T_norm(v,par,flag) calculates rho for
% ITLS solution of norm_p minimisation
% 
% References:
%
% Borsic, A.; A. Adler, B. M. G. & Lionheart, W. T. B. 
%   In--vivo impedance imaging with total variation regularization 
%   IEEE Trans. Med. Imaging, 2010, 29(1), 44-54
% Farquharson, C. G. 
%   Constructing piecewise-constant models in
%   multidimensional minimum-structure inversions 
%   Geophysics, 2008, 73(1), K1-K9
% Farquharson, C. G. & Oldenburg, D. W. 
%   Non-linear inversion using general 
%   measures of data misfit and model structure 
%   Geoph. J. Int., 1998, 134, 213-227
% Gersztenkorn, A.; Bednar, J. B. & Lines, L. R. 
%   Robust iterative inversion for the one-dimensional 
%   acoustic wave equation 
%   Geophysics, 1986, 51, 357-368
% 
% vr Feb 13, 2010
% 
if nargin<3, flag='x'; end
if nargin<2, par=[]; end

n=size(v);

switch lower(flag)
    case {'1' 'p' 'l_p'}
        % Norm p, see Gersztenkorn 1986
        p=par(1);
        Ti=abs(v).^p;
        T=sum(Ti);
    case {'2' 'e' 'ekblom'}
        % Norm p approximated, see Farquharson 1998,2008
        p=par(1);ep2=par(2)^2;
        Ti=(v.^2+ep2).^(p/2);
        T=sum(Ti);
    case {'3' 's'  'minsup'}
        % Minimum support, see Farquharson 1998,2008
        c=par(1);cc=c^2;
        vv=v.^2;
        T=sum(vv./(vv+cc));
    case {'4' 'm' 'r' 'robust'}
        % Huber Norm, see Farquharson 1998,2008
        c=par(1);cc=c^2;
        Ti(1:n)=0;
        ind=abs(v) <= c;Ti(ind)=v(ind).*v(ind);
        ind=abs(v) >  c;Ti(ind)=2*c*abs(v(ind))-cc;
        T=sum(Ti);
    otherwise
        T=sum(v.*v);
end
