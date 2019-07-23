function res= vertcat(varargin)
%ADDERIV/VERTCAT Concatenate derivatives vertically
%
% Copyright 2001-2004 Andre Vehreschild, Institute for Scientific Computing
%                     RWTH Aachen University
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!

if nargin==1
   res= varargin{1};
elseif nargin==2 % Faster, because no reallocation occurs.
   if isempty(varargin{1})
      res= varargin{2};
      return;
   end
   if isempty(varargin{2})
      res= varargin{1};
      return;
   end
   if ~(isa(varargin{1}, 'adderiv') && isa(varargin{2}, 'adderiv'))
      error('Both arguments of a vertical matrix concatenation have to be adderivs.');
   end;

   res= adderiv(varargin{1});

   if res.dims==1
      for i= 1: res.ndd
         res.deriv{i}= [varargin{1}.deriv{i}; varargin{2}.deriv{i}];
      end
   else
      for i= 1: res.ndd(1)
         for j= 1: res.ndd(2)
            res.deriv{i,j}= [varargin{1}.deriv{i,j}; varargin{2}.deriv{i,j}];
         end
      end
   end
else
   c=1;
   while (c<=nargin) && isempty(varargin{c})
      c= c+1;
   end

   if c==nargin
      if isempty(varargin{c})
         res= g_dummy;
      else
         res= varargin{c};
      end
      return;
   end
   
   if ~isa(varargin{c}, 'adderiv')
      error(sprintf('All arguments of a vertical matrix concatenation have to be adderivs. %d is a %s', c, isa(varargin{c})));
   end;

   res= varargin{c}; % Intended deep-copy !!!
   c= c+ 1;
   while c<= nargin
      if isempty(varargin{c})
         c= c+ 1;
         continue;
      end
      if ~isa(varargin{c}, 'adderiv')
         error(sprintf('All arguments of a vertical matrix concatenation have to be adderivs. %d is a %s', c, isa(varargin{c})));
      end;

      if res.dims==1
         for i= 1: res.ndd
            res.deriv{i}= [res.deriv{i}; varargin{c}.deriv{i}];
         end;
      else
         for i= 1: res.ndd(1)
            for j= 1: res.ndd(2)
               res.deriv{i,j}= [res.deriv{i,j}; varargin{c}.deriv{i,j}];
            end
         end
      end
      c= c+1;
   end
end

% vim:sts=3:
