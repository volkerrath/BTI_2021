function rhoi=rhoiT(T)
% calculate ice density in [kg/m**3]
[n1 n2]=size(T);if n1==1, T=T'; end 
rhoi=0.00001*ones(size(T));
rhoi(T<0)=917. - 0.151*T(T<0);
