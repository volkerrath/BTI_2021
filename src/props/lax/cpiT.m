function [cpi]=cpiT(T)
% calculate ice isobaric heat capacity [J/(kg*K]
% After:
% Fukusako, S.:
% Thermophysical Properties of Ice, Snow,and Sea Ice
% International Journal of Thermophysics, 1990, 11, 353-372,
% p. 356;
T=T(:);
Tk=T+273.15;
cpi=185+6.89*Tk;
