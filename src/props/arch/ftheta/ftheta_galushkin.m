function [Theta,dTheta]=ftheta_galushkin(T,Tf) 
% [theta,dtheta]=ftheta(T,Tf,) calculates the fluid/ice  
% partition function used for the apparent heat capacity approach
% to phase change.
%
% On Input: 
% T            =     temperature
% Tf           =     phase boundary temperature
%                   (DEFAULT: 0)
% % On Ouput: 
% Theta            =    value of partition function
%                       (0<Theta<1)
% dTheta            =   derivative of partition function
%                       with respect to temperature
%
if nargin < 2, Tf=0.;            end


a0=   1.0;
a1=   0.60152823763179;
a2=   0.23218232347212;
a3=   0.04669792788297;
a4=   0.00535597924776;
a5=   0.00036415588418;
a6=   0.00001450956751;
a7=   0.00000031279149;
a8=   0.00000000281502;

Theta=a0+a1*T+a2*T.^2+a3*T.^3+a4*T.^4+a5*T.^5+a6*T.^6 +  a7*T.^7 + a8*T.^8;
b=find(T>Tf);Theta(b)=1.;

if nargout>1,
    dTheta=a1+2*a2*T+3*a3*T.^2+4*a4*T.^3+5*a5*T.^4+6*a6*T.^5 +  7*a7*T.^6 + 8*a8*T.^7;
    dTheta(b)=0.;
end
