function [cpi]=cpiT(T)
% calculate ice isobaric heat capacity [J/(kg*K]
% regression on published data
% formulation according to Fukusako (1990) 

[n1 n2]=size(T);if n1==1, T=T'; end 

cpi=185.+6.89*(273.15+T)
