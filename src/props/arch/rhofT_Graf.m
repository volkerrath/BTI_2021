function [rhow]=rhofT(T)
% [rho]=rhofT(T) calculates the density i(in kg/m^3) of formation water,
% given temperature in C,
%
% Range of validity: unknown
%
% Code  verification:
%		INPUT:	TEMP = 298.15K	P =0.1013 MPa  S = 0.25 g/g	 OUTPUT: RHO = 1187.35 kg/m3
%		INPUT:	TEMP = 393.15K	P =   30 MPa   S = 0.10 g/g	 OUTPUT: RHO = 1027.06 kg/m3
%
% Fortran source written written by JJAadams for Alberta Geological Survey
% May. 2001
% Matlab code by V. Rath, UCM, November 2010
% References:  T. Graf (2009), Simulation of geothermal flow in deep
%                              sedimentary basins. ERCB/AGC OFR 2009-11
T=T(:);
dom1=find(T >=0 & T < 40); T1=T(dom1);
dom2=find(T >=40 & T < 100);T2=T(dom2);
dom3=find(T >=100 & T < 300);T3=T(dom3);


visf(T<0)       = 1.787e-3;
visf(dom1)      = 1.787e-3*exp((-0.03288+1.962e-4*T1)*T1);
visf(dom21)     = 1.e-3*(1+0.015512*(T2-20))^(-1.571);
visf(dom3)      = 0.2414e-4*10^(247.8/(T+133.15));
visf(T>300)     =  9.0121e-05; % 0.2414e-4*10^(247.8/(300+133.15));