function [avg]=vavg(v,w,flag)
% [avg]=vavg(V,W) calculates harmonic avg of
% vector V with weights W.
% vr Feb 12, 2010

if nargin < 2, w=ones(size(v)); end
if nargin < 3, flag='a'; end

v=v(:);w=w(:);

d=sum(w);w=w./d;

switch lower(flag)
    case {'a' 'arithmetic'}
        avg=sum(v.*w);
    case {'g' 'geometric'}
        avg=exp(sum(log(v).*w));
    case {'h' 'harmonic'}
        avg=1./sum(w./v);
    case {'s' 'sqrt'}
        avg=sum(w.*sqrt(v)).^2;
    otherwise
        disp(['VAVG: unknown average type (default=a)'])
        avg=sum(v.*w);
end

