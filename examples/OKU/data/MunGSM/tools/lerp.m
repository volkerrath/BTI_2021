function yi = lerp(x,y,xi,extrapspan)
%LERP Quick 1-D linear interpolation with bounded extrapolation region.
%   F=LERP(X,Y,XI) returns the value of the 1-D function Y at the points
%   of column vector XI using linear interpolation. Length(F)=length(XI).
%   The column vector X specifies the coordinates of the underlying
%   interval.
% 
%   Vectorization allows LERP to be much faster that built-in MATLAB
%   interpolation functions for some applications that would otherwise
%   require many looped function calls.
%   
%   If Y is a matrix and XI is a column vector, then the interpolation is
%   performed for each column of Y in which case F is
%   length(XI)-by-size(Y,2).
% 
%   If Y is a column vector and XI is an ND array, F is size(XI).
% 
%   If Y is a matrix and XI is an array with N dimensions, F will have N+1
%   dimensions. The first N dimensions correspond to size(XI), and the N+1
%   dimension is size(Y,2). If XI is a row vector it is treated as a matrix
%   with 2 dimensions.
% 
%   F=LERP(X,Y,XI,EXTRAP) will expand the function domain (defined in
%   coordinates X) with linear extrapolation by a distance EXTRAP in either
%   direction. EXTRAP may be a 2-element vector with the first element
%   defining the extrapolation distance below the beginning of X and the
%   second element defining the extrapolation distance beyond the final
%   value of X. Elements of EXTRAP must be non-negative. 
% 
%   NaN's are returned for values of XI outside the interval defined by X
%   and EXTRAP (if provided).
% 
%   LERP does no input checking. For LERP to work properly:
%       X must be a monotonically increasing column vector.
%       Y must be a column vector or matrix with length(X) rows.
%
%   Example: interpolation for log and log10 functions with plot comparing
%   extrapolation to reality
%       X  = (linspace(.01, 2, 7)').^2;
%       Y  = [log(X) log10(X)];
%       XI = magic(7)*.2 - .5;
%       F  = lerp(X, Y, XI, [0 4])
%       x = linspace(4,8);
%       plot(X,Y,'-k',x,log(x),'r:',x,log10(x),'r:',...
%            XI,F(:,:,1),'bo',XI,F(:,:,2),'ko')   
% 
%   See also INTERP1Q, INTERP1.
% 
%   F=LERP(X,Y,XI,EXTRAP)

%   Author: Sky Sartorius
%       http://www.mathworks.com/matlabcentral/fileexchange/authors/101715

%% Set up extrapolation regions
if nargin == 4 %if extrapolating, alter inputs
    if any(extrapspan<0)
        error('Extrapolation region(s) must be non-negative')
    end
    if numel(extrapspan)>2
        error('Extrapolation definitions has too many elements')
    end
    if extrapspan(1)
        yl = y(1,:)-(y(2,:)-y(1,:))/(x(2)-x(1))*extrapspan(1);
        xl = x(1)-extrapspan(1);
    else yl = []; xl = []; 
    end
    if extrapspan(end)
        yu = y(end,:)+(y(end,:)-y(end-1,:))/...
            (x(end)-x(end-1))*extrapspan(end);
        xu = x(end)+extrapspan(end);
    else yu = []; xu = [];
    end
    y =[yl; y; yu];
    x =[xl; x; xu];
end

siz = size(xi);
if length(xi)~=1
   xi = xi(:);
   [xxi,k] = sort(xi);
   [~,j]=sort([x;xxi]);
   n = length(x);
   ni = length(xxi);
   r(j)=1:(n+ni);
   r=r(n+1:end)-(1:ni);
   r(k)=r;
   r(xi==x(end))=n-1;
   ind=find((r>0) & (r<n));
   ind = ind(:);
   yi = NaN(ni,size(y,2), superiorfloat(x,y,xi));
   rind = r(ind);
   u = (xi(ind)-x(rind))./(x(rind+1)-x(rind));
   yi(ind,:)=y(rind,:)+(y(rind+1,:)-y(rind,:)).*u(:,ones(1,size(y,2)));
else % Special scalar xi case exactly from interp1q
   r = find(x <= xi, 1, 'last' );
   r(xi==x(end)) = length(x)-1;
   if isempty(r) || (r<=0) || (r>=length(x))
      yi = NaN(1,size(y,2),superiorfloat(x,y,xi));
   else
      u = (xi-x(r))./(x(r+1)-x(r));
      yi=y(r,:)+(y(r+1,:)-y(r,:)).*u(:,ones(1,size(y,2)));
   end
end

if numel(siz) == 2
    if siz(2) == 1
        siz = siz(1);
    end
end
yi = reshape(yi,[siz size(y,2)]);
end

% Revision history:
% V1.0 9 August 2010
% V1.1 couple of logic and coding fixes; some edits to help and example
