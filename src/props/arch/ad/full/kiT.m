function ki=kiT(T)
% calculate ice thermal conductivity [mW/(m*K)]
% regression on published dataif nargin <2, nonlinear='yes'; end;
ki=0.00001*ones(size(T));
ki(T<0)=2.164- 0.0112*T(T<0);
