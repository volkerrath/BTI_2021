function [kmT,dmT]=kmT(k,T,A,B)
% calcuates thermal conductivity as function of temperature
% based on the formula given by 
% Vosteen et al.,Phys. Chem. Earth, 28 (2003), 499-509
%
%  k        =  thermal conductivity at 25C 
%  T        =  temperature in C 
%  A, B     =  coefficients: 0.0030 0.0042 (crystalline)
%                            0.0034 0.0039 (sedimentary)
%              default is sediment
% vr April 3, 2005
if nargin < 3, A=0.0034; B=0.0039; end
[n1,n2]=size(k); if n1~=1, k=k'; end 
[n1,n2]=size(T); if n1~=1, T=T'; end 
[n1,n2]=size(A); if n1~=1, A=A'; end 
[n1,n2]=size(B); if n1~=1, B=B'; end 

km0=0.53*k + 0.5*sqrt(1.13*k.^2-0.42*k);
kmT= km0./(0.99+T*(A-B/km0));
if nargout>1, dmT=0.45*kmT;end
