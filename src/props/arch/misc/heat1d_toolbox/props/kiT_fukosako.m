function ki=kiT(T)
% calculate ice thermal conductivity [mW/(m*K)]
% formulation according to Fukusako (1990) 

if nargin <2, nonlinear='yes'; end;

 ki = 0.4685 +488.19e0./(T+273.15);
% ki = 2.2156 -0.0100456*T +3.4452e-5*T.*T;
