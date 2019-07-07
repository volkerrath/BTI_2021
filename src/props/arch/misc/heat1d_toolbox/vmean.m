function [mean]=vmean(v,w,flag)
% [MEAN]=HMEAN(V,W) calculates harmonic mean of
% vector V with weights W.
% vr Feb 12, 2010
if nargin < 2, W=ones(size(v)); end
if nargin < 3, flag='a'; end
d=sum(w);
switch lower(flag)
    case {'a' 'arithmetic'}
        mean=sum(v.*(w/d));
    case {'g' 'geometric'}
        w=w/d;
        mean=exp(sum(log(v).*w));
    case {'h' 'harmonic'}
        mean=d/sum(w./v);
    otherwise
        disp(['VMEAN: unknown average type (default=a)'])
        mean=sum(v.*(w/d));
end

