function [visw]=visfT(T)
% [visf]=visfT(T) calculates the density i(in kg/m^3) of formation water,
% given temperature in C, pressure in MPa, and salinity in mass fraction
% (g/g).
%%
% Fortran source written written by JJAadams for Alberta Geological Survey
% May. 2001
% Matlab code by V. Rath, UCM, November 2010
% References:  T. Graf (2009), Simulation of geothermal flow in deep
%                              sedimentary basins. ERCB/AGC OFR 2009-11
T=T(:);
% freshwater
rhow = 1.e3*(1-((T-3.9863)^2/508929.2)*(T+288.9414)/(T+68.12963));
