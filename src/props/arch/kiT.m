function ki=kiT(T)
% calculate ice thermal conductivity [mW/(m*K)]
% After:
% Fukusako, S.:
% Thermophysical Properties of Ice, Snow,and Sea Ice
% International Journal of Thermophysics, 1990, 11, 353-372,
% p. 355;
T=T(:);
ki=2.24.*ones(size(T));

