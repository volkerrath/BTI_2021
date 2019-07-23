function [rho,rhow]=rhofT(T,P,S)
% [rho]=rhofT(T,P,S) calculate the density i(in kg/m^3) of formation water, 
% given temperature in C, pressure in MPa, and salinity in mass fraction
% (g/g).
%
% Range of validity: Pressures 5-100 MPa, Temperature 20-350°C, Salinity <=320000 mg/L
%
% Code  verification:
%		INPUT:	TEMP = 298.15K	P =0.1013 MPa  S = 0.25 g/g	 OUTPUT: RHO = 1187.35 kg/m3	 
%		INPUT:	TEMP = 393.15K	P =   30 MPa   S = 0.10 g/g	 OUTPUT: RHO = 1027.06 kg/m3	 
% 
% Fortran source written written by JJAadams for Alberta Geological Survey
% May. 2001
% Matlab code by A. Schmidt and V. Rath, RWTH Aachen University, April 2004
% References: Batzle & Wang, Seismic properties of pore fluids. 
%                    Geophysics, Nov. 1992, v.57, no.11, pp. 1396-1408, 1992.
%             Adams & Bachu, EOS for basin geofluids:algorithm review and
%                   intercomaperison for brines, Geofluids,2,257-271,2002.

if nargin < 3, S=0. ; end
if nargin < 2, P=0.1; end
[n1,n2]=size(T); if n1~=1, T=T'; end 
[n1,n2]=size(P); if n1~=1, P=P'; end 
[n1,n2]=size(S); if n1~=1, S=S'; end 

% freshwater 
P1=P/1.e6;P2=P1.*P1;T2=T.*T;T3=T2.*T;
rhow = 1+1.e-6* ...
       (-80*T-3.3*T2+0.00175*T3+489*P1-2*T.*P1+0.016*T2.*P1-1.3e-5*T3.*P1-(0.333-0.002*T).*P2);
rhow=rhow*1000.;rho = rhow;

% brine 
if S > 0.0000001,
   rho = (rhow/1.e3 + S.*(0.668 + 0.44*S + 1.0e-6 .* (3.e2*P1 - 2.4e3 ...
        .* P1 .* S + T.* (80 + 3*T - 3.3e3*S -13*P1 + 47*P1.*S))))*1.e3 ;   
end
    
