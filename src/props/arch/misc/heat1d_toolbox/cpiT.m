function [cpi]=cpiT(T)
% calculate ice isobaric heat capacity [J/(kg*K]
% regression on published data
 cpi=2110+7.7*T;
% cpi=2110