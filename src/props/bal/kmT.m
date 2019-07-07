function [kmT]=kmT(km0,T,A,B,T0,T1)
% calculates thermal conductivity as function of temperature
% based on the formula given in Clauser & Huenges
% (1995, AGU Reference shelf).
% T           =  temperature in C  acording to flag unit
% A,B         =  coefficients for formula
%                (default: A=0.7,B=770, see reference)
% T0          =  base temperature for caluculation of
%                relative conductivities (default: T0=20 C)
% v. r. nov. 2, 2002
if nargin < 5, T0=20.; T1=800; end
if nargin < 3, A=0.7;B=770; end

km0=km0(:);T=T(:);A=A(:);B=B(:);

if A>0,
kbas  = A+B/(350+T0);
klim  = A+B/(350+T1);
kda   = A+B./(350+T);
ilimit=find(T>=T1);
kda(ilimit) =klim(ilimit);
f=km0./kbas;kmT=f.*kda;
else
kmT=km0;
end
kmT=kmT(:);

