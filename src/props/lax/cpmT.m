function [cpmT]=cpmT(cpm0,T)
% calculates heat capacity of rocks as function of T  
% after Sundberg 2009 , Rath et al 2018
% vr oct 2018
cpm0=cpm0(:);T=T(:);
f=.25e-2; 
cpmT=cpm0  +cpm0*f.*(T-25.);

