function [kf,kfw]=kfT(T,S)
% [kf,kfw]=kfT(T,S,method) calculate the thermal conductivity kf in W/(m*K) of
% formation water, given temperature in C, and salinity in mass fraction
% (g/g)of NaCl. Thermal conductivity of freshwater, kfw is calculated using
% the Phillips (1981) formulation.
% Ling, F. & Zhang (2004):
% 	A Numerical Modelfor surface energy balance and the thermal 
% 	regime of the active layer and permafrost containing 
% 	unfrozen water or brine
% 	Cold Regions Science & Technology, 38, 1-15,
% Osterkamp (1987):
%	 Freezing and thawing of soils  and permafrost containing 
% 	unfrozen water or brine, WRR, 23(12),2279ff. 

T=T(:);
if nargin < 2, S=0; end
if length(S)==1;S=S*ones(size(T));

% normalized temperature Tr
Tr=(T+273.15)/273.15; Tr2=Tr.*Tr;Tr3=Tr2.*Tr;Tr4=Tr3.*Tr;
kfw = (-0.92247 + 2.8395*Tr - 1.8007*Tr2 + 0.52577*Tr3 ...
                                                - 0.07344*Tr4);
kf=kfw;

if T<-0.001
kfw = 0.22455+1.6318e-3.*(T+273.15d0)
end 

if any(S),
    T2=T.*T;
    C = S./(1 + S)*1.d2;C2=C.*C;    % C=salinity in mol/kg
    kf = kfw.*(1.d0 - (2.3434d-3 - 7.924d-6*T + 3.924d-8*T2).*C ...
                         + (1.06d-5 - 2.d-8*T - 1.2d-10*T2).*C2)
end
