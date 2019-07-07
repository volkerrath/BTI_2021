function g_r= g_abs(g_p, p)
% G_ABS -- Compute the derivative g_p according to the value of p
% g_r is the negative of g_p if p is less than zero,
% g_r is g_p if p is greater than zero, and
% if p is zero then an error is issued, because abs is not differentiable
% at zero.
%
% Copyright 2001-2004 Andre Vehreschild, Institute for Scientific Computing   
%                     RWTH Aachen University
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!

sig_p= sign(p);
if (any(find(sig_p==0.0)))
   error('g_abs(g_p, p) not defined for p==0.0');
end

g_r= sig_p.* g_p;

% vim:sts=3:

