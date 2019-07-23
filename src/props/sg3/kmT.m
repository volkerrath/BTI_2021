function [kmT,dmT]=kmT(k,T,A,B)
% calcuates thermal conductivity as function of temperature
% based on the formula developed for the Kola SG3 area 
% given by 
% Vosteen et al.,Phys. Chem. Earth, 28 (2003), 499-509  and \n
% Mottaghy et al,Int J Earth Sci (Geol Rundsch) (2008) 97,435–442\n
%  k        =  thermal conductivity at 25C 
%  T        =  temperature in C 
%  A, B     =  coefficients: 	0.0030 0.0042 (crystalline)
%                            	0.0034 0.0039 (sedimentary)
%				0.0013 0.0029 (Kola)	
%              default is sediment
% vr March 2013
Tlim=800;
%if nargin < 3, A=0.0034; B=0.0039; end
if nargin < 3, A=0.0013; B=0.0029; end
k=k(:);T=T(:);A=A(:);B=B(:);

km0=0.53*k + 0.5*sqrt(1.13*k.^2-0.42*k);

i1=T<=Tlim;i2=T>Tlim;
kmT(i1)= km0(i1)./(0.99+T(i1).*(A(i1)-B(i1)./km0(i1)));
kmT(i2) = km0(i2)./(0.99+Tlim.*(A(i2)-B(i2)./km0(i2)));
kmT=kmT(:);
if nargout>1, dmT=0.45*kmT;end
