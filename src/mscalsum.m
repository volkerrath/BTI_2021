function [S,JS]=scaljac(J,mode)
% calulates scale for Jacobian and optionally scales it
% vr July 2014
if nargin < 2, mode='col'; end

[n,m]=size(J);
switch lower(mode)
    case {'col'}
        V=1./sum(J,2); S=spdiags(V(:),[0],n,n);
        whos
        if nargout == 2,
            JS=S*J;
        end
    case {'row'}
        V=1./sum(J,1); S=spdiags(V(:),[0],m,m)
        if nargout == 2,
            JS=J*S;
        end
      
end ;

end