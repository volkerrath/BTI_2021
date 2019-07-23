
function [kmT]=kmT(A,T)
% calculates thermal conductivity as function of temperature
% based on the formula given in Clauser & Huenges
T0=25;
km0=2.91;T=T(:);
kbas  = A(1)+A(2)/(350+T0);
kda   = A(1)+A(2)./(350+T);
f=km0./kbas;kmT=f.*kda;
kmT=kmT';
%whos
end
