function l = kfT_kaye(T)
% kwT_kaye Thermal conductivity (lambda) of water
%
% k = kfT_kaye(T) computes the thermal conductivity of water at temperature
% T for fresh water.
%
% Units: 
% T: [°C]
% K: [W (m K)^-1]
%
% After Griffiths (1992) and Kaye & Laby (1968).

% aha, 31-dec-2003

l = 0.56 + 0.002*T - 1.01e-5*T.^2 + 6.71e-9*T.^3;
