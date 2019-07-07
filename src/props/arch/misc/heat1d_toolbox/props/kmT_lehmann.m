function [kmT]=kmT_lehmann(km0,T,A,B,T0,T1)
% calculates thermal conductivity as function of temperature
% based on the formula given in Clauser & Huenges (1995) and 
% the modification by Lehmann (1998) for high temperatures.
%
% km0       =  thermal conductivity at base Temperature T0
% T           =  temperature (C) 
% A,B         =  coefficients for formula A+B/(350+T)
%                (default: A=0.7,B=770, see reference)
% T0          =  base temperature for calculation of 
%                temperature-dependent conductivities (default: T0=20 C)
% Tlimit          =  base temperature for calculation of 
%                temperature-dependent conductivities (default: T0=20 C)
% v. r. july 13, 2004

if nargin < 3, A=0.7;B=770; end
if nargin < 5, T0=0.;T1=400; end

[n1,n2]=size(km0); if n1~=1, k0=k0'; end 
[n1,n2]=size(T); if n1~=1, T=T'; end 
[n1,n2]=size(A); if n1~=1, A=A'; end 
[n1,n2]=size(B); if n1~=1, B=B'; end 

if A>0,
 kmdf = A+B./(350+T);kmdf0= A+B/(350+T0);kmdf1= A+B/(350+T1);
 c0=km0./kmdf0;fac=c0-(c0-1).*(T-T0)./(T1-T0);
 kmT=kmdf.*fac;
% ilimit0=find(T<T0);
% if ~isempty(ilimit0)
%     kmT(ilimit0)=kmdf0(ilimit0);
% end
% ilimit1=find(T>T1); 
% if ~isempty(ilimit1)
%     kmT(ilimit1)=kmdf1(ilimit1);
% end

else
 kmT=km0;
end
