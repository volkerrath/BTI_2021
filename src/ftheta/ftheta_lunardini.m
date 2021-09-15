function [Theta,dTheta]=ftheta(T,Tf,w) 
% [theta,dtheta]=ftheta(T,Tf,w) calculates the fluid/ice  
% partition function used for the apparent heat capacity approach
% to phase change.
%
% On Input: 
% T            =     temperature
% Tf           =     phase boundary temperature
%                   (DEFAULT: 0)
% w            =     scaling factor for "mushy" region 
%                   (DEFAULT: 1)
%
% On Ouput: 
% Theta            =    value of partition function
%                       (0<Theta<1)
% dTheta            =   derivative of partition function
%                       with respect to temperature
%
if nargin < 3, w=.5;            end
if nargin < 2, Tf=0;            end


Theta=exp(-((T-Tf)./w).^2);
b=find(T>Tf);Theta(b)=1.;

if nargout>1,
    a=-2/(w*w);
    dTheta=a*(T-Tf).*Theta;
    dTheta(b)=0.;
end