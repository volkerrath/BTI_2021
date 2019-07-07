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
T=T(:); P=P(:);


ac=[0.10000000e+01  0.17472599e-01 -0.20443098e-04 -0.17442012e-06 ...
   0.49564109e-02 -0.40757664e-04  0.50676664e-07  0.50330978e-04 ...
   0.33914814e-06 -0.18383009e-06 ...
   0.10009476e-02  0.16812589e-04 -0.24582622e-07 -0.17014984e-09 ...
   0.48841156e-05 -0.32967985e-07  0.28619380e-10  0.53249055e-07 ...
   0.30456698e-09 -0.12221899e-09];

% Pressure in MPa
P=P/1.e6;P2=P.*P;P3=P2.*P;P4=P3.*P;
         T2=T.*T;
         T3=T2.*T;TP=P.*T;T2P=T2.*P;TP2=T.*P2;
%    liquid density
%
a=ac(1)+ac(2)*P+ac(3)*P2+ac(4)*P3+ac(5)*T+ac(6)*T2+ac(7)*T3+ac(8)*TP+ac(10)*T2P+ac(9)*TP2;
b=ac(11)+ac(12)*P+ac(13)*P2+ac(14)*P3+ac(15)*T+ac(16)*T2+ac(17)*T3+ac(18)*TP+ac(20)*T2P+ac(19)*TP2;
rho=a./b;

% after Speedy (1987) for T < 0 to -46 C
bc=[ 0.999195706402050*901.5328593 -0.0011761652 ...
    0.0038442382 -0.0157270761  0.0744064614 -0.1406432653];
ic=T<0.;
T(T<-45.)=-45;
r=(T(ic)+273.15-227.15)/227.15;
r2=r.*r;r3=r2.*r;r4=r3.*r;
rho(ic)=bc(1)*exp(-227.15*( bc(3)*r + ...
                             0.5*bc(4)*r2 + ...
                           0.333333*bc(5)*r3 + ...
                            0.25*bc(6)*r4 + ...
			       2*bc(2)*sqrt(r)));


