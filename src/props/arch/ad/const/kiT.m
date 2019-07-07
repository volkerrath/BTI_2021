function ki=kiT(T)
% calculate ice thermal conductivity [mW/(m*K)]
% as function of temperature (dummy)
% VR RWTH Aachen University,   April 25, 2004
% regression on published data
[n1,n2]=size(T); if n2~=1, T=T'; end 
 ki=2.164*ones(size(T));
