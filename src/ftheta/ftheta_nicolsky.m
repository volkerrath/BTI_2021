function [Theta,dTheta]=ftheta_nicolsky(T,Ts,p) 
% [theta,dtheta]=ftheta(T,Ts,m) calculates the fluid/ice  
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
if nargin < 2, Ts=1.;            end
if nargin < 3,  w=1.;            end

b=find(T<Ts);Theta=ones(size(T));
Theta(b)=abs(Ts)^p * abs(T(b)).^(-p) ;

if nargout>1,
    dTheta=zeros(size(T));
    dTheta(b)=abs(Ts)^p*(-p)*abs(T(b)).^(-p-1).*sign(T(b));
end
