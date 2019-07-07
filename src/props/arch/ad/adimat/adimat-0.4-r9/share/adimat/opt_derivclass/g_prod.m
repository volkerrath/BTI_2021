function g_y= g_prod(g_x, x, cr)
% G_PROD -- Derivative computation for prod.
% 
% Computes the derivative of prod observing the dimension along which 
% prod was used.
%
% Support routine of ADiMat. Used by differentiated code.
%
% Copyright 2005 Andre Vehreschild, Institute for Scientific Computing
%           RWTH Aachen University.
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!

if nargin<3
  cr=1;
end

num= size(x, cr);

if (cr== 1)
  g_y= g_x(1,:).* prod(x(2: num, :), 1);  
  for i= 2: (num- 1)
    g_y=g_y+ prod(x(1: (i-1), :), 1).* g_x(i, :).* ...
	  prod(x((i+1): num, :), 1);
  end
  g_y= g_y+ prod(x(1: (num-1), :), 1).* g_x(num, :);
elseif (cr== 2)
  g_y= g_x(:,1).* prod(x(:, 2: num), 2);  
  for i= 2: (num- 1)
    g_y=g_y+ prod(x(:, 1: (i-1)), 2).* g_x(: ,i).* ...
	  prod(x(:, (i+1): num), 2);
  end
  g_y= g_y+ prod(x(:, 1: (num-1)), 2).* g_x(:, num);

else
  error 'Value of cr illegal. Has to be 1 or 2!';
end

