function [R]=drho_norm(v,par,flag)
% [R]=R_d(Obs,funflag) calculates R Matrix for
% IRLS solution of norm_p minimisation
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

n=length(v);


switch lower(flag)
    case {'1' 'p' 'l_p'}
        % Norm p, see Gersztenkorn 1986
        p=par(1);gam=par(2);
        Rii(1:n)=p*gam^(p-2);
        Rii(abs(v) > gam)=p*v(abs(v) > gam).^(p-2);
        R=spdiags(Rii,0,n,n);
    case {'2' 'e' 'ekblom'}
        % Norm p, see Farquharson 1998,2008
        p=par(1);ep2=par(2)^2;
        Rii=p*(v.^2+ep2).^(p/2-2);
        R=spdiags(Rii,0,n,n);
    case {'3' 's'  'minsup'}
        % Minimum support, see Farquharson 1998,2008
        pp=par(1)^2;
        vv=v.^2;
        Rii=2*pp/(vv+pp).^2;
        R=spdiags(Rii,0,n,n);
    case {'4' 'm' 'r' 'robust'}
        % Norm p, see Farquharson 1998,2008
        c=par(1);
        Rii(1:n)=2;
        Rii(abs(v) > c)=2*c/abs(v);
        R=spdiags(Rii,0,n,n);
    otherwise
        R=2*speye(n);
end
