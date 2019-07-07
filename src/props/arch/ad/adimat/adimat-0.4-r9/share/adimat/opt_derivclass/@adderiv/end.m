function res=end(g, k, n)
% ADDERIV/END(G, K, N) hands the index-expression down to the objects stored
% in the derivative.
%
% Copyright 2001-2004 Andre Vehreschild, Institute for Scientific Computing   
%                     RWTH Aachen University
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!

res= size(g.deriv{1}, k);

