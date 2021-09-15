function g_r= g_norm(r, g_x, x, p)
% G_NORM -- Compute the derivative of the norm specified by p of x.
%
% Copyright 2001-2005 Andre Vehreschild, Institute for Scientific Computing   
%                     RWTH Aachen University
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!


if nargin<4
   p= 2; % Ensure default case
end

if ischar(p) 
   if strcmp(lower(p), 'fro')
      if all(size(x)>1)
         g_r= call(@sum, call(@diag, g_x'*x+ x'*g_x))./(2*r);
         warning('Derivative of Frobenius-norm not tested yet.');
      else
         error('Frobenius-norm specified for matrix only.');
      end
   else
      error('Only "fro" is a valid string for p-norm computation currently.');
   end
else
   switch p
      case inf
         if all(size(x)>1)
            error('Derivative of matrix-Inf-norm not implemented.');
         else
            [val, ind]= max(abs(x));
            if (val== r)
               g_r= g_x(ind);
            else
               error('The value of norm function and the value on the index position do not correspond.');
            end
         end
      case -inf
         if all(size(x)>1)
            error('Derivative of matrix-minimum-norm not implemented.');
         else
            [val, ind]= min(abs(x));
            if (val== r)
               g_r= g_x(ind);
            else
               error('The value of norm function and the value on the index position do not correspond.');
            end
         end
      case 2
         if all(size(x)>1)
            error('Derivative of matrix-2-norm not implemented.');
         else
            g_r=  call(@sum, abs(x).* g_abs(g_x, x))./ r; 
         end
      otherwise
         if all(size(x)>1)
            error('Derivatives of matrix-p-norm not implemented yet.');
         else
            g_r= r.^(1-p).* call(@sum, abs(x).^(p-1).* g_abs(g_x, x));
         end
   end
end

% vim:sts=3:sw=3:ts=3:

