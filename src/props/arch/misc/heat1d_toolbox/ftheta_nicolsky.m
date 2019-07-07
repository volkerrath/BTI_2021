function [Theta,dTheta]=ftheta_nicolsky(T,Ts,p) 
% [theta,dtheta]=ftheta(T,Ts,m) calculates the fluid/ice  
% partition function used for the apparent heat capacity approach
% to phase change.
%
% On Input: 
% T            =      temperature
% Ts          =       phase boundary temperature
%                    (DEFAULT: 0)  
%
% On Ouput: 
% Theta            =    value of partition function
%                       (0<Theta<1)
% dTheta            =   derivative of partition function
%                       with respect to temperature
%
if nargin < 3,   p= 0.5; end
if nargin < 2,  Ts=-0.01;end

if p >= 0.5 & p <= 0.8 & Ts < 0,
        b=find(T<Ts);Theta=ones(size(T));
        Theta(b)=abs(Ts)^p * abs(T(b)).^(-p) ;
        
    if nargout>1,
        dTheta=zeros(size(T));
        dTheta(b)=abs(Ts)^p*(-p)*abs(T(b)).^(-p-1).*sign(T(b));  
    end
else
    
    error ('ftheta_nicolsky: incorrect values for parameters p/Ts ')
    
end
