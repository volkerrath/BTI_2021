function rhoi=rhoiT(T)
% calculate ice density in [kg/m**3]
% After:
% Fukusako, S.:
% Thermophysical Properties of Ice, Snow,and Sea Ice
% International Journal of Thermophysics, 1990, 11, 353-372,
% p. 357;
T=T(:);
rhoi=917. - 0.1073*T;

