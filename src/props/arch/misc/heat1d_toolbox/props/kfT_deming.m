function [kf]=kfT(T)
% [kf,kfw]=kfT(T,S,method) calculate the thermal conductivity kf in W/(m*K) of
% formation water, given temperature in C.
% Range of validity:  0 to 300°C 
%
%  References: Deming and Chapmann (1988) Heat flow from BHT data,
%         JGR93(11),13657ff
%
%  V. Rath, RWTH Aachen University, April 2004
%
[n1 n2]=size(T);if n1==1, T=T'; end 
a=find(T<137):b=find(T>137)
kf(a)=0.5648+1.878e-3*T(a)-7.231e-6*T(a).*T(a);
kf(b)=0.6020+1.309e-3*T(a)-5140e-6*T(a).*T(a);
