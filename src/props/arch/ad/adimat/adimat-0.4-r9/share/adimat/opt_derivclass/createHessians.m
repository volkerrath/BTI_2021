function [varargout]= createHessians(g_dirs, varargin)
%CREATEHESSIANS Create Hessians.
%   [hessians]= createHessians(g_dirs, objects OR gradients)
%   Create a Hessian-matrix for every object OR gradient supplying
%   enough slots for g_dirs number of directional derivatives.
%   Specify either the "normal" object or the gradients, but
%   NEVER both.
%   The hessians will be sparse matrices using the MATLAB 
%   sparse datatype, if applicable.
%   The number of input objects ("normal" objects or gradients)
%   has to be equal to the number of desired Hessians.
%
%   g_dirs is scalar. The Hessian will g_dirs x g_dirs. If g_dirs
%   is a vector than g_dirs(1)==g_dirs(2) has to be true!
%
%   Example call:
%   [h_A, h_b, h_c]= createHessians(7, A, g_b, c);
%
% Copyright 2001-2005 Andre Vehreschild, Institute for Scientific Computing   
%                     RWTH Aachen University
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!


if nargin-1~= nargout
   error('The number of input objects and the number of outputs has to be equal.');
end

% Ensure that g_dirs specifies a matrix.
if isempty(g_dirs)
   g_dirs= get(adderiv([], [], 'empty'), 'NumberOfDirectionalDerivatives');
   if g_dirs(1)==1
      g_dirs= [g_dirs(2) g_dirs(2)];
   end
elseif length(g_dirs)==1
   g_dirs=[g_dirs g_dirs];
elseif ~ (g_dirs(1)==g_dirs(2))
   error('The number of directional derivatives of a Hessian has to be n x n.');
end

for i=1: nargout
   if isa(varargin{i}, 'adderiv')
      adobj= varargin{i};
      if getDims(adobj)==1
         sz= size(adobj{1});
      else
         sz= size(adobj{1,1});
      end
   else
      sz= size(varargin{i});
   end
   
   if isstruct(varargin{i})
      res= deepcopy(varargin{i});
      varargout{i}= adderiv(g_dirs, res, 'object');
   else
      if length(sz)>2
         varargout{i}= adderiv(g_dirs, sz, 'zeros');
      else
         varargout{i}= adderiv(g_dirs, sz, 'sparse');
      end
   end
end

set(varargout{1}, 'NumberOfDirectionalDerivatives', g_dirs);

% vim:sts=3:

