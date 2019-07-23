function [Theta,dTheta]=ftheta_galushkin(T,Tb,p) 
% [theta,dtheta]=ftheta(T,Tf,w) calculates the fluid/ice  
% partition function used for the apparent heat capacity approach
% to phase change.
%
% On Input: 
% T            =     temperature
% Tb           =    1 solidus temperature
%                     (DEFAULT: -2)
%                   2 liquidus temperature
%                     (DEFAULT: 0)
% p            =     scaling factor for "mushy" region 
%                   (DEFAULT: 10)
%
% On Ouput: 
% Theta            =    value of partition function
%                       (0<Theta<1)
% dTheta            =   derivative of partition function
%                       with respect to temperature
%
if nargin < 3, p=10;            end
if nargin < 2, Tb=[-2 -0.001];            end

b=find(T>=Tb(2))

f=p/(Tb(2)-Tb(1));
Theta=(exp(-f*T)-exp(-f*Tb(2)))./(exp(-f*Tb(1))-exp(-f*Tb(2)));
Theta(b)=1.;

if nargout>1,
    fp=f*p;
    dTheta = (fp*exp(T*fp))/(exp(Tb(1)*fp) - exp(Tb(2)*fp));
    dTheta(b)=0.;
end
