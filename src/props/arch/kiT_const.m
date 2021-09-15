function ki=kiT(T)
% calculate ice thermal conductivity [mW/(m*K)]
% as function of temperature (dummy)
% VR RWTH Aachen University,   April 25, 2004
% regression on published data
T=T(:);
ki=2.164*ones(size(T));
