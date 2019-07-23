function setADoption(adobj, cmd, val)
%ADDERIV/setADOption -- Sets an ADiMat-option to the desired value.
%   The adobj is just a dummy needed for correct class dispatching.
%   This method is not accessable, if no adderivsp-object is supplied.
%   The option set is nevertheless global to the current MATLAB session.
%   Internal use only!
%
% Copyright 2001-2004 Andre Vehreschild, Institute for Scientific Computing   
%                     RWTH Aachen University
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!

warning('ADiMat:deprecation', 'Deprecation warning: setADoption() will be deprecated in future versions of ADiMat. Use set() instead.');
option(cmd, val);

