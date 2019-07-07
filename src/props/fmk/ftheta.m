function [Theta,dTheta]=ftheta(T,Tf,w,rfl) 
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
if nargin < 4, rfl=.025;            end
if nargin < 3, w=1.;            end
if nargin < 2, Tf=0;            end


Theta=exp(-((T-Tf)./w).^2);
b=find(T>Tf);Theta(b)=1.;
a=find(Theta<rfl);Theta(a)=rfl;
if nargout>1,
    c=-2/(w*w);
    dTheta=c*(T-Tf).*Theta;
    dTheta(b)=0.;
    dTheta(a)=0.;
end
