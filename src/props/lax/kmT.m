function [kmT]=kmT(km0,T,A,B,T0,T1)
% calculates thermal conductivity as function of temperature
% after Sundberg 2009 , Rath et al 2018
% vr oct 2018
km0=km0(:);T=T(:);
f=0.03e-2*km0; 
kmT=km0-f.*(T-25.);

kmT=kmT(:);

