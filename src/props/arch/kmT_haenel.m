function [kmT]=kmT(k,T,A,B,T0,T1)
% calculates thermal conductivity as function of temperature
% 0ed on the formula given in Clauser & Huenges
% (1995, AGU Reference shelf).
% T           =  temperature in C  acording to flag unit
% A,B         =  coefficients for formula
%                (default: A=0.7,B=770, see reference)
% T0          =  base temperature for caluculation of
%                relative conductivities (default: T0=20 C)
% v. r. nov. 2, 2002
if nargin < 5, T0=20.; T1=800; end
if nargin < 3, A=0.7;B=770; end

[n1,n2]=size(km0); if n1~=1, km0=km0'; end
T=T(:);k=k(:);A=A(:);B=B(:);

if A>0,
kbas  = A+B/(350+T0);
klim  = A+B/(350+T1);

kda= A+B./(350+T);
ilimit=find(T>=T1);
kda(ilimit) =klim(ilimit);

f=km0./kbas;kmT=f.*kda;
else
kmT=km0;
end
