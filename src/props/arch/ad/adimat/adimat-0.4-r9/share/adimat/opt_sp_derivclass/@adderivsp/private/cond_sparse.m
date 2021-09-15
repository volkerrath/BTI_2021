function res= cond_sparse(arg)
% ADDERIV/PRIVATE/COND_SPARSE -- Convert the matrix to a sparse one, if
%   the number of nonzeros in the matrix is less then 1/3 of the matrix
%   size.
%
% Copyright 2003, 2004 Andre Vehreschild, Inst. f. Scientific Computing
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!

if (issparse(arg))
  res= arg;
else
  cap= prod(size(arg));
  tmp= sparse(arg);
  if (nnz(tmp)< cap/3)
    res= tmp;
  else
    res= arg;
  end
end

