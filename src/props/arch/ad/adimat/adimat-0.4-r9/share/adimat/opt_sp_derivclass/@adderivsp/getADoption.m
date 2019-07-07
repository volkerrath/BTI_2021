function res= getADoption(adobj, cmd)
%ADDERIV/getADOption -- Gets the value of an ADiMat-option.
%   The adobj is just a dummy needed for correct class dispatching.
%   This method is not accessable, if no adderiv-object is supplied.
%   The option set is nevertheless global to the current MATLAB session.
%   Internal use only!
%
% Currently no valid options are available.
%
% Copyright 2001-2004 Andre Vehreschild, Institute for Scientific Computing   
%                     RWTH Aachen University
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!

warning('ADiMat:deprecation', 'Deprecation warning: getADoption() will be deprecated in future versions of ADiMat. Use get() instead.');
res= option(cmd);

