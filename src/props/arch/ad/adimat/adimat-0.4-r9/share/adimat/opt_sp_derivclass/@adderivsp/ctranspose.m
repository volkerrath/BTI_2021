function res=ctranspose(g)
%ADDERIV/CTRANSPOSE Complex conjugate transpose the derivative.
%
% Copyright 2001-2004 Andre Vehreschild, Institute for Scientific Computing   
%                     RWTH Aachen University
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!


res=adderivsp(g);

if g.dims==1
  for i= 1: g.ndd(1)
    res.deriv{i}= g.deriv{i}';
  end
else
  for i= 1: g.ndd(1)
    for j= 1: g.ndd(2)
      res.deriv{i,j}= g.deriv{i,j}';
    end
  end
end

