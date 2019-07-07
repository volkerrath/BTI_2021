function ki=kiT(T)
% calculate ice thermal conductivity [mW/(m*K)]
if nargin <2, nonlinear='yes'; end;

% regression on published data
        ki=2.164- 0.0112*T;
%        ki=2.164
