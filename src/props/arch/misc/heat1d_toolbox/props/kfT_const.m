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
[n1,n2]=size(P); if n1~=1, T=T'; end 
kf=0.561*ones(size(T));
