function [xi,yi,zi] = clinspace(x,y,z,n)

%CLINSPACE Curvilinearly spaced points.
%   [XI,YI] = CLINSPACE(X,Y,N) generates N equally curvilinearly spaced
%   points interpolating the curve represented by the series (X,Y).
%
%   [XI,YI,ZI] = CLINSPACE(X,Y,Z,N) generates N equally curvilinearly
%   spaced points interpolating the curve represented by the series
%   (X,Y,Z).
%
%   If N is omitted, 100 points are generated.
%
%   Examples:
%   --------
%   x = 0:0.01:10;
%   y = 5*sin(x.^3/100).^2;
%   [xi,yi] = clinspace(x,y,10);
%   plot(x,y,xi,yi,'o'), axis equal
%
%   x = -2*pi:0.1:2*pi;
%   y = 10*sin(x);
%   z = linspace(0,10,length(x)).^1.5;
%   [xi,yi,zi] = clinspace(x,y,z,10);
%   plot3(x,y,z,xi,yi,zi,'o')
%   box on, axis equal
%
%   See also LINSPACE.
%
%   -- Damien Garcia -- 2008/01, revised 2009/02


%% Check input arguments
error(nargchk(2,4,nargin))

if nargin==2 % clinspace(x,y)
    z = zeros(size(x));
    n = 100;
elseif nargin==3
    if isscalar(z) % clinspace(x,y,n)
        n = z;
        z = zeros(size(x));
    else % clinspace(x,y,z)
        n = 100;
    end
end
        
if ~isscalar(n)
    error('MATLAB:clinspace:BadNArgument','N must be a scalar')
end

x = x(:); y = y(:); z = z(:);
N = length(x);

if length(y)~=N || length(z)~=N
    error('MATLAB:clinspace:NumelMismatch',...
        'X, Y (and Z) must have same number of elements')
end

n = floor(n);
if n<1
    error('MATLAB:clinspace:BadNArgument','N must be >=1')
end

%% Particular cases: N=1 or N=2
if N==1
    error('MATLAB:clinspace:OnlyOnePoint',...
        'At least two points are required')
elseif N==2
    xi = linspace(x(1),x(2),n);
    yi = linspace(y(1),y(2),n);
    zi = linspace(z(1),z(2),n);
    return
end

%% Test if CUMSIMPS exists
if ~exist('cumsimps','file')
    error('MATLAB:clinspace:MissingFunction',...
        ['CUMSIMPS is required. Download <a href="matlab:web(''',...
        'http://www.biomecardio.com/matlab/cumsimps.html'')',...
        '">CUMSIMPS</a>.'])
end

%% Curvilinearly spaced points

[dx,dy,dz] = deal(zeros(size(x)));

% [x y z]'(t1)
dx(1) = sum(x(1:3).*[-3 4 -1]')/2;
dy(1) = sum(y(1:3).*[-3 4 -1]')/2;
dz(1) = sum(z(1:3).*[-3 4 -1]')/2;

% [x y z]'(tN)
dx(end) = sum(x(N-2:N).*[1 -4 3]')/2;
dy(end) = sum(y(N-2:N).*[1 -4 3]')/2;
dz(end) = sum(z(N-2:N).*[1 -4 3]')/2;

% [x y z]'(ti) (i>2 & i<(N-1))
Dc = [1 0 -1]';
tmp = conv(Dc,x)/2; dx(2:N-1) = tmp(3:N,:);
tmp = conv(Dc,y)/2; dy(2:N-1) = tmp(3:N,:);
tmp = conv(Dc,z)/2; dz(2:N-1) = tmp(3:N,:);

L = cumsimps(hypot(hypot(dx,dy),dz));
L = L + linspace(0,L(end),N)'*eps; % to avoid equal X values in INTERP1
Li = linspace(0,L(end),n);

xi = interp1(L,x,Li,'cubic');
yi = interp1(L,y,Li,'cubic');
zi = interp1(L,z,Li,'cubic');

