function [x,dx]=set_mesh(x0,x1,nn,type,direction,debug)
% [X,DX]=T_MESH(X0,END,NC,TYPE,DIRECTION) generates a lin/log mesh
%      x0,end           x0 and end value of mesh (in seconds!)
%      nn               number of nodes, i.e. nn-1 cells
%      type	        meshtype ('lin', 'log'), default='linear'
%      direction        1 if from ts=>te, -1,if from -ts=>te  
% V. R.,Oct 13, 2002 


if nargin < 6, debug=0;end
if nargin < 5, direction=1;end
if nargin < 4, type='lin';end
if nargin < 3, nn=128;end

switch lower(type)
 case {'lin' 'linear'}
     x=direction*linspace(x0,x1,nn);
 case {'log' 'logarithmic'}
     nn=nn-1;
     x=direction*logspace(log10(x0),log10(x1),nn);
     if direction < 0, x=[x 0]; else  x=[0 x]; end
 otherwise
     x=direction*linspace(x0,x1,nn);
end

nx=length(x);dx=diff(x);

if debug==1,
 disp(['    mesh: ', num2str(nx),   ' nodes,  x:  ',...
                     num2str(x(1)),' - ', ...
                     num2str(x(nx)), ' units ' ] );
 disp(['          ', num2str(nx-1), ' cells, dx: ',...
                     num2str(dx(1)),' - ', ...
                     num2str(dx(nx-1)), ' units ' ] );
end	       
