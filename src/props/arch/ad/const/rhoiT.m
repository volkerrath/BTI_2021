function rhoi=rhoiT(T)
% calculate ice density in [kg/m**3] 
% as function of temperature (dummy)
% VR RWTH Aachen University,   April 25, 2004
[n1,n2]=size(T); if n2~=1, T=T'; end 
rhoi=917.*ones(size(T));
