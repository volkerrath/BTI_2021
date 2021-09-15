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
T=T(:);k=k(:);A=A(:);B=B(:);


km0=0.53*k + 0.5*sqrt(1.13*k.^2-0.42*k);
kmT= km0./(0.99+T.*(A-B./km0));

if nargout>1, dmT=0.45*kmT;end
