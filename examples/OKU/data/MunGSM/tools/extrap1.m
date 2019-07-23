function y_out = extrap1(x_in,y_in,x_out,method)
%function y_out = extrap1(x_in,y_in,x_out,method)
%
%This function is analogous to interp1.m but tests for NaN's
%associated with interp1.m
%Dan Feldman (C) 2004
%California Institute of Technology
%
%%%%%%%%%%%%%%%%%%%%%
%%Input variable(s)%%
%%%%%%%%%%%%%%%%%%%%%
%x_in   = independent variable for input vector
%y_in   = dependent variable for input vector
%x_out  = independent variable for output vector
%method = interp1 methos
%        'nearest'  - nearest neighbor interpolation
%        'linear'   - linear interpolation
%        'spline'   - piecewise cubic spline interpolation (SPLINE)
%        'pchip'    - piecewise cubic Hermite interpolation (PCHIP)
%        'cubic'    - same as 'pchip'
%        'v5cubic'  - the cubic interpolation from MATLAB 5, which does not
%                   extrapolate and uses 'spline' if X is not equally spaced.
%
%%%%%%%%%%%%%%%%%%%%%%
%%Output variable(s)%%
%%%%%%%%%%%%%%%%%%%%%%
%y_out  = dependent variable for output vector

s_identify = 'extrap1.m';

y_out = interp1(x_in,y_in,x_out,method);

for i=1:length(y_out)
  if isnan(y_out(i)) | isinf(y_out(i))
    if strcmp(method,'nearest') | strcmp(method,'linear')
      %get nearest y_in(j)
      a = [1:length(y_in)];
      a_prime = zeros(size(a));
      index = 1;
      for j=1:length(a)
	if ~isnan(y_in(j))
	  a_prime(index) = a(j);
	  index = index + 1;
	end
      end
      a_prime = a_prime(1:index-1);
      %y_out(i) = y_out(a_prime(find(min(abs(x_out(a_prime)-x_out(i))))));
      diff = 1e9;
      for j=1:length(a_prime)
	diff_comp = abs(x_in(a_prime(j))-x_out(i));
	if diff_comp<diff
	  diff = diff_comp;
	  index = a_prime(j);
	end
      end
      y_out(i) = y_in(index);
    elseif strcmp(method,'pchip')
      y_prime = interp1(x_in,y_in,x_out,'cubic');
      y_out(i) = y_prime(i);
    end
  end
end