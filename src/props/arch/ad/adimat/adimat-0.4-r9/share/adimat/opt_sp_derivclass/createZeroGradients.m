function [varargout]= createZeroGradients(g_dirs, varargin)
%CREATEZEROGRADIENTS Create the gradient shells.
%   [gradients]= createEmptyGradients(g_dirs, realobjects)
%   Create for every realobject a gradient object supplying
%   enough slots for g_dirs number of directional derivatives.
%   All gradient objects are zero matrices. 
%   The number of the realobjects has to be equal to the
%   number of gradients.
%
%   Example:
%   [g_A, g_b, g_c]= createZeroGradients(7, A, b, c);
%
% Copyright 2001-2004 Andre Vehreschild, Institute for Scientific Computing   
%                     RWTH Aachen University
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!


if nargin-1~= nargout
   error('The number of realobjects and gradients have to be equal.');
end

for i=1: nargout
   if isstruct(varargin{i})
      res= deepcopy(varargin{i});
      varargout{i}= adderivsp(g_dirs, res, 'object');
   else
      varargout{i}= adderivsp(g_dirs, size(varargin{i}), 'zeros');
   end
end

% Set the global option.
set(varargout{1}, 'NumberOfDirectionalDerivatives', g_dirs);

% vim:sts=3:

