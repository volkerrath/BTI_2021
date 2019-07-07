function g_res= ls_mexp(g_s1, s1, g_s2, s2, res)
%ADDERIV/LS_MEXP Execute the exponentiation rule for two operands.
%
% Copyright 2001-2004 Andre Vehreschild, Institute for Scientific Computing   
%                     RWTH Aachen University
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!

g_res= adderivsp(g_s1);

tmp= log(s1);
for i= 1: g_res.ndd
  g_res.deriv{i}= cond_sparse((g_s2.deriv{i}* tmp+ s2* (g_s1.deriv{i}./ s1))* res);
end

