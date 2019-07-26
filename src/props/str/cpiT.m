function [cpi]=cpiT(T)
% calculate ice isobaric heat capacity [J/(kg*K]
% regression on published data
[n1 n2]=size(T);if n1==1, T=T'; end 
cpi=0.00001*ones(size(T));
cpi(T<0)=2110+7.7*T(T<0);
