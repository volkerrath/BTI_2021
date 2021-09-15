function [cpi]=cpiT(T)
% calculate ice isobaric heat capacity [J/(kg*K]
% regression on published data
T=T(:);
 cpi=2110+7.7*T;
