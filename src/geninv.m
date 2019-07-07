function [G] = GenInv(Jw,reg)
% calculates Generalized Inverse G
[nd,np]=size(Jw);
L0 = reg1d(np,'l0');L1 = reg1d(np,'l1');L2 = reg1d(np,'l2');
A =[ Jw;sqrt(reg(3))*L2 ;sqrt(reg(2))*L1 ;sqrt(reg(1))*L0 ];
G=inv(A'*A)*Jw';
end