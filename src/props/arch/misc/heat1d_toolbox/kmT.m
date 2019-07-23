function [kmT]=kmT(lamb0,T,A,B,T0,T1)
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

[n1,n2]=size(T); if n1==1, T=T'; end 

if A>0,
labas  = A+B/(350+T0); 
lalim  = A+B/(350+T1);

lambda                  = A+B./(350+T);
ilimit=find(T>=T1);
lambda(ilimit) =lalim(ilimit);

f=lamb0./labas;kmT=f.*lambda;
else
kmT=lamb0;
end
