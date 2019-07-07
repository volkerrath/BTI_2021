function rhoi=rhoiT(T)
% calculate ice density in [kg/m**3]
% as function of temperature (dummy)
% VR RWTH Aachen University,   April 25, 2004
T=T(:);
% if nargin<2, P=0.1*ones(size(T));
rhoi=917.*ones(size(T));
