function [S]=f_scal(J,mode)
% calulates scale for Jacobian
% see:
% Mehanee & Zhdanov
% Two-dimensional magnetotelluric inversion of
% blocky geoelectric structures
% J. Geophys. Res., 2002, 107
% doi:10.1029/2001JB000191

% vr NOV 2014
if nargin < 2, mode='col'; end

[n,m]=size(J);
switch lower(mode)
    case {'col'}
        V=sum((J).^2,1);
    case {'row'}
        V=sum((J).^2,2);
end
V=sqrt(V);
S=spdiags(V(:),[0],m,m)
end
