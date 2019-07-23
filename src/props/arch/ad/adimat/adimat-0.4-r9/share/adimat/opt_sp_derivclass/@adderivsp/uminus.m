function res=uminus(g)
%ADDERIV/UMINUS Negate the derivative
%
% Copyright 2001-2004 Andre Vehreschild, Institute for Scientific Computing   
%                     RWTH Aachen University
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!

res= adderivsp(g);

if res.dims==1
   for i= 1: res.ndd
     res.deriv{i}= -g.deriv{i};
   end
else
   for i= 1: res.ndd(1)
     for j= 1: res.ndd(2)
        res.deriv{i,j}= -g.deriv{i,j};
     end
   end
end

