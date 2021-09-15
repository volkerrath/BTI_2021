function [S,JS]=scaljac(J,mode)
% calulates scale for Jacobian and optionally scales it
% vr July 2014
if vargin < 2, mode='col'; end

[n,m]=size(J);
switch lower(mode)
    case {'col'}
        V=sum(J,2); S=spdiags(V(:),[0],n,n)
        if vargout = 2,
            JS=S*J;
        end
    case {'row'}
        V=sum(J,1); S=spdiags(V(:),[0],m,m)
        if vargout = 2,
            JS=J*S;
        end
      
end ;

end