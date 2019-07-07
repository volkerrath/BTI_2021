function [kmT]=kmT(k,T,Tref)
% calcuates thermal conductivity as function of temperature
% based on the formula given by Sass et al.(1992, JGR 97, p5017ff)
% (see also Clauser, AGU Bookshelf, 1995).
% k           =  thermal conductivity at Tref
% T           =  temperature in C
% Tref        =  temperature for k measurement (default=25C)
% v. r. Nov 3. , 2001
if nargin < 3, Tref=25.; end
T=T(:);k=k(:);

k0=k.*(1.007+Tref.*(0.0037-(0.0074/k)));
kd=    1.007+   T.*(0.0036-(0.0072/k0));

kmT=k0./kd;

