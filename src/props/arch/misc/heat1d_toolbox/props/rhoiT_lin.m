function rhoi=rhoiT(T)
% calculate ice density in [kg/m**3]
[n1 n2]=size(T);if n1==1, T=T'; end 
       rhoi=917. - 0.151*T;
%        rhoi=917.;
