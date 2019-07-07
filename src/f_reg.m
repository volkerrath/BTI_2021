function [J,Q,S]=f_reg(Dm)
% Calculates the values of the weights matrices used in the
% edge preserving inversion 
%
% Youzwishen & Sacchi
% Edge preserving Imaging
% Jour. Seismic Exploration
% 15(4), 45-56, 2006
%

% vr August 2012
Dm=Dm(:);n=length(Dm);
Dm2=abs(Dm).^2;
J=sum(Dm2./(1+Dm2));
if nargout > 1,
    d=1./(1+Dm2).^2;
    Q=spdiags([d(:)],[0],n,n);
end
if nargout > 2,
    S=spdiags([sqrt(d(:))],[0],n,n);
end

end
