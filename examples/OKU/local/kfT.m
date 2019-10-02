function [kf,kfw]=kfT(T,S)
% [kf,kfw]=kfT(T,S,method) calculate the thermal conductivity kf in W/(m*K) of
% formation water, given temperature in C, and salinity in mass fraction
% (g/g)of NaCl. Thermal conductivity of freshwater, kfw is calculated using
% the Phillips (1981) formulation. 
% Range of validity:  20 to 330°C and up to 4 molal NaCl
%
% Code  verification:
%		INPUT:	TEMP = 298.15K	P =0.1013 MPa  S = 0.25 g/g	 
%     OUTPUT: TC_BR2 = 0.587 W/(m*K)
%		INPUT:	TEMP = 393.15K	P =   30 MPa   S = 0.10 g/g	
%     OUTPUT: TC_BR2 = 0.688 W/(m*K)
%
%  References: Phillips SL, Igbene A, Fair JA, Ozbek H, Tavana M (1981) A
%         technical databook for geothermal energy utilization. Lawrence Berkeley
%         Laboratory Report 12810, 46 p.
%
% Fortran source written written by JJAadams for Alberta Geological Survey
% May. 2001; Matlab code V. Rath, RWTH Aachen University, April 2004
%
if nargin < 2, S=0. ; end
% normalized temperature Tr 
Tr=(T+273.15)/273.15; 
Tr2=Tr.*Tr; 
Tr3=Tr2.*Tr;
Tr4=Tr3.*Tr; 
kfw = (-0.92247d0 + 2.8395d0*Tr - 1.8007d0*Tr2 + 0.52577d0*Tr3 ...
                                                - 0.07344d0*Tr4);
kf=kfw;

if S > 0.0000001,
    T2=T.*T;
    C = S./(1 + S)*1.d2;C2=C.*C;    % C=salinity in mol/kg
    kf = kfw.*(1.d0 - (2.3434d-3 - 7.924d-6*T + 3.924d-8*T2).*C ...
                         + (1.06d-5 - 2.d-8*T - 1.2d-10*T2).*C2)
end