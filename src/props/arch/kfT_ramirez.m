function [kf,kfw]=kfT(T,S)
% [kf,kfw]=kfT(T,S,method) calculate the thermal conductivity kf in W/(m*K) of
% formation water, given temperature in C, and salinity in mass fraction
% (g/g)of NaCl.
%
% Range of validity:  unknown
%
% Code  verification:
%		INPUT:	TEMP = 298.15K	P =0.1013 MPa  S = 0.25 g/g	 OUTPUT: TC_BR2 = 0.587 W/(m*oC)
%		INPUT:	TEMP = 393.15K	P =   30 MPa   S = 0.10 g/g	 OUTPUT: TC_BR2 = 0.688 W/(m*oC)
%
%
% References: Ramirez, MLV and CAN deCastro 1994. Thermal conductivity of aqueous
%     sodium chloride solutions. Journal of Chemistry and Engineering Data, 39, p. 186-190.
%
% Fortran source written written by JJAadams for Alberta Geological Survey
% May. 2001; Matlab code V. Rath, RWTH Aachen University, April 2004
T=T(:);
if nargin < 2, S=0; end
if length(S)==1;S=S*ones(size(T));

T2=T.*T;
f0=0.5621 + 0.00199*T -8.6e-6*T2;
kfw=f0;kf=kfw;


if any(S),
	 C = S/58.443*1.d3; C2=C.*C;% Salinity mol/kg
         f1 = -0.01394 + 0.000294*T - 2.3e-6*T2;
  	 f2 =  0.00177 - 6.3e-5*T     + 4.5e-7*T2;
	 kf = (f0 + f1*C + f2*C2);
end

