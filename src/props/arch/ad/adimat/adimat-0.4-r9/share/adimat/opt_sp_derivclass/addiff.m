function [success, msg]=addiff(fname, indeps, deps, flags)
% ADDIFF Differentiate fname using ADiMat
%
% Use automatic differentiation to get the derivatives of fname with respect to
% deps and indeps. Both indeps and deps enumerate variable names of the 
% parameter and result list of fname, respectively. If indeps and/or deps 
% is empty, all variables of the associated list are selected.
%
% Returns (logical) true on success, false otherwise.
%
% flags gives the command-line flags as described in the html-manual of ADiMat
% provided with this installation.
%
% The function fname and all user-defined functions it depends on either have to
% be in the current directory for ADiMat being able to find them or a searchpath
% (flag -I <PATH>) has to be specified.
%
% Examples: Suppose f is function of x and c and returns y
%
%     function y= f(x,c)
%
% 1. To get the derivative of f with respect to both inputs issue:
%
%	     addiff(@f);
%
%    If no errors occured the function g_f is generated in the current
%    directory. To evaluate it do something like:
%
%        [g_x, g_c]= createFullGradients(x,c);
%        [g_r, r]= g_f(g_x, x, g_c, c);
%
% 2. To compute the derivative of f with respect to x only, do:
%
%       addiff(@f, 'x');
%
% 3. Supplying a flag to ADiMat is possible, too:
%
%       addiff(@f, 'x', [], '-v20 -f');
%
%    This statement issues a lot of information (-v20) and applies the
%    forward mode of AD (-f). Which is the only one available currently.
%
% 4. Return true or false, whether ADiMat succeeded or failed.
%
%       r=addiff(...anything...)
% 
%    r contains a logical value now, indicating true, if ADiMat succeeded.
%    Version information, error and warning messages are printed nevertheless.
%
% 5. Like 4., but put messages into the string msg.
%
%       [r, msg]= addiff(...anything...)
%
%    No output it made. The messages are stored in msg and r contains the 
%    status of the command as in 4.
%
% Copyright 2003- 2005 Andre Vehreschild, Institute for Scientific Computing
%           RWTH Aachen University.
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!
%

% Ensure, that all inputs are in a valid format.
if (nargin<4) , flags= ''; end
if (nargin<3 | isempty(deps)) , deps=''; else deps= ['-d', deps]; end
if (nargin<2 | isempty(indeps)) , indeps=''; else indeps= ['-i', indeps]; end
if (nargin<1) ; help addiff; return; end
if (isa(fname, 'function_handle')) , fname= func2str(fname); end

adimat_home= get(g_dummy, 'ADiMatHome');
[stat, message]= system(['"', adimat_home, '/bin/adimat" ', flags, ...
                         ' ', deps, ' ', indeps, ' ', fname]);
if nargout<2
   if (stat)
      disp(message);
   else
      disp(message);
   end
end

if nargout>0
	success= ~stat;
   if nargout>1
      msg=message;
   end
end

