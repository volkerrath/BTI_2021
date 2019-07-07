function Lh=LathT(T)
% calculate ice thermal conductivity [mW/(m*K)]
% After:
% Ling, F. & Zhang (2004):
% 	A Numerical Modelfor surface energy balance and the thermal 
% 	regime of the active layer and permafrost containing 
% 	unfrozen water or brine
% 	Cold Regions Science & Technology, 38, 1-15,
% Osterkamp (1987):
%	 Freezing and thawing of soils  and permafrost containing 
% 	unfrozen water or brine, WRR, 23(12),2279ff. 
% vr 6/2011
T=T(:);
Lh=rhofT(T)*(333.2+4.995*T+0.2987*T.^2);
