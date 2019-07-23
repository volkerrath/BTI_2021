function [km]=ke2km(ke,kf,phi,mode)
% matrix contistivity using an  arithmetic, geometric,harmonic, or 
% sqare-root mixing law 
% vr Apr 3, 2014
if nargin<4 | isempty(mode), mode='a'; end

switch lower(mode)
case {'a' 'ari' 'arithmetic'}
   km=(ke(:)-phi(:).*kf(:))./(1-phi(:));             %(arithmetic mean)
case {'g' 'geo' 'geometric'}
   km=exp((log(ke(:))-phi(:).*log(kf(:)))./(1-phi(:))); %(geometric mean)
case {'h' 'har' 'harmonic'}
   km=(1-phi(:))./(1./ke(:)-phi(:)./kf(:));
case {'s' 'sqr' 'square'}
   km=((sqrt(ke(:))-phi(:).*sqrt(kf(:)))./(1-phi(:))).^2;
otherwise
   warning(['avg: mode <',mode,' not implemented, arithmetic assumed!'])
   km=(ke(:)-phi(:).*kf(:))./(1-phi(:));             %(arithmetic mean)
end


