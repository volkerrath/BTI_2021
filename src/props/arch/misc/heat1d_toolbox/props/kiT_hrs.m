function ki=kiT(T)
% calculate ice thermal conductivity [mW/(m*K)]
% formulation according to HRS Manual 

if nargin <2, nonlinear='yes'; end;
ki = 2.2156 -0.0100456*T +3.4452e-5*T.*T;
