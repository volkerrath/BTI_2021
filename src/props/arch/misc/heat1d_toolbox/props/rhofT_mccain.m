function [rho,rhow]=rhofT(T,P,S)
% [rho]=rhofT(T,P,S) calculate the density i(in kg/m^3) of formation water, 
% given temperature in C, pressure in Pa, and salinity in mass fraction
% (g/g).
%
% Written by JJA for Alberta Geological Survey May. 2000 for NT
%
% References: McCain: Reservoir-Fluid Property Correlations - State of the Art
%                   SPE Reservoir Engineering, May 1991, v.6, no.2, pp. 266-272, 1991
%             Adams & Bachu: EOS for basin geofluids:algorithm review and
%                   intercomaperison for brines, Geofluids,2,257-271,2002.
%
%  Range of Validity: Pressure          0.69-69 MPa, 
%                     Temperature of       <=127 °C, 
%                     Salinity          up to <=450000 
%
%  Code verification:
%     INPUT:  TEMP = 298.15K  P =0.1013 MPa  S = 0.25 g/g      OUTPUT: RHO = 1186.52 kg/m3     
%     INPUT:  TEMP = 393.15K  P =   30 MPa   S = 0.10 g/g      OUTPUT: RHO = 1023.06 kg/m3     
%
% Fortran source written written by JJAadams for Alberta Geological Survey
% May. 2001
% Matlab code by A. Schmidt and V. Rath, RWTH Aachen University, April 2004
%  
if nargin < 3 S=0.; end
if nargin < 2 P=0.1; end
[n1,n2]=size(T); if n1~=1, T=T'; end 
[n1,n2]=size(P); if n1~=1, P=P'; end 
[n1,n2]=size(S); if n1~=1, S=S'; end 

Sm = S*100;
Tm = 1.8*T+32;
Pm = P*1.45038e-4;

rhos = 62.368 + 0.438603*Sm + 1.60074e-3*Sm.*Sm;
dvt = -1.0001e-02 + 1.33391e-04*Tm+5.50654e-07*Tm.*Tm;
dvp = Pm.*((-1.95301e-09 - 1.72834e-13*Pm).*Tm - 3.58922e-07 - 2.25341e-10*Pm);

bw = (1 + dvp).*(1 + dvt);
rho = 16.01846*rhos./bw;
