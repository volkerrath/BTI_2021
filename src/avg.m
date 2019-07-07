function [avg]=avg(w,v,mode)
% calculates the wwigthed arithmetic, geometric,harmonic, or 
% sqare-root mean from matrices w and v, acording to parameter mode 
% vr Aug 4, 2012
if nargin<3 | isempty(mode), mode='a'; end
sw= sum(w,2);wnorm=sw;
if sum(w,2)~=1, warning('avg: sum of volumes not 1!'); end

switch lower(mode)
case {'a' 'ari' 'arithmetic'}
   avg= sum(w.*v,2);avg=avg./wnorm;
case {'g' 'geo' 'geometric'}
   avg= exp(sum(log(v).*w,2)./wnorm);
case {'h' 'har' 'harmonic'}
   avg= wnorm/ sum(w./v,2);
case {'s' 'sqr' 'square'}
   avg = (sum(w.*sqrt(v),2)./wnorm).^2;
otherwise
   warning(['avg: mode <',mode,' not implemented, arithmetic assumed!'])
   avg= sum(w.*v,2);avg=avg./wnorm;
end
