function [rho]=rhofT(T,P)
% [rho]=rhofT(T,P) calculate the density i(in kg/m^3) of pure water, 
% given temperature (T, in C), and optionally pressure (P,in Pa).
%
% Derived from the Formulation given in:
% Zylkovskij et al: Models and Methods Summary for the FEHMN Application,
%       ECD 22, LA-UR-94-3787, Los Alamos NL, 1994.
% 
% Range of validity: Pressures   0.01 - 110 MPa, 
%                    Temperature   15 - 350 °C
%
% VR RWTH Aachen University,   April 25, 2004

if nargin < 2, P=0.1; end
[n1 n2]=size(T);if n1==1, T=T'; end 
[n1 n2]=size(P);if n1==1, P=P'; end 
% Pressure in MPa 
P=P/1.e6;P2=P.*P;P3=P2.*P;P4=P3.*P;
         T2=T.*T;
         T3=T2.*T;TP=P.*T;T2P=T2.*P;TP2=T.*P2;

c=[0.10000000e+01  0.17472599e-01 -0.20443098e-04 -0.17442012e-06 ...
   0.49564109e-02 -0.40757664e-04  0.50676664e-07  0.50330978e-04 ...
   0.33914814e-06 -0.18383009e-06 ...
   0.10009476e-02  0.16812589e-04 -0.24582622e-07 -0.17014984e-09 ...
   0.48841156e-05 -0.32967985e-07  0.28619380e-10  0.53249055e-07 ...
   0.30456698e-09 -0.12221899e-09];
%    liquid density
%     
a=c(1)+c(2)*P+c(3)*P2+c(4)*P3+c(5)*T+c(6)*T2+c(7)*T3+c(8)*TP+c(10)*T2P+c(9)*TP2;
b=c(11)+c(12)*P+c(13)*P2+c(14)*P3+c(15)*T+c(16)*T2+c(17)*T3+c(18)*TP+c(20)*T2P+c(19)*TP2;
rho=a./b;
