function rhoi=rhoiT(T)
% calculate ice density in [kg/m**3]
T=T(:);
rhoi=917. - 0.151*T;
