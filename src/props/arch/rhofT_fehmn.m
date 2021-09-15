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
%     new: after Speedy (1987) for T < 0 to -46 C
% VR RWTH Aachen University,   April 25, 2004
T=T(:);
if nargin < 2, P=1.e5*ones(size(T)); end
rho=zeros(size(T);

% Pressure in MPa
P=P/1.e6;P2=P.*P;P3=P2.*P;P4=P3.*P;
         T2=T.*T;T3=T2.*T;TP=P.*T;T2P=T2.*P;TP2=T.*P2;

a=[0.10000000e+01  0.17472599e-01 -0.20443098e-04 -0.17442012e-06 ...
   0.49564109e-02 -0.40757664e-04  0.50676664e-07  0.50330978e-04 ...
   0.33914814e-06 -0.18383009e-06 ...
   0.10009476e-02  0.16812589e-04 -0.24582622e-07 -0.17014984e-09 ...
   0.48841156e-05 -0.32967985e-07  0.28619380e-10  0.53249055e-07 ...
   0.30456698e-09 -0.12221899e-09];
b=[901.5328593 -0.0011761652 0.0038442382 -0.0157270761,...
   0.0744064614, -0.1406432653];


%    liquid density

p1=a(1)+a(2)*P+a(3)*P2+a(4)*P3+a(5)*T+a(6)*T2+a(7)*T3+a(8)*TP+a(10)*T2P+a(9)*TP2;
p2=a(11)+a(12)*P+a(13)*P2+a(14)*P3+a(15)*T+a(16)*T2+a(17)*T3+a(18)*TP+a(20)*T2P+a(19)*TP2;
rho=p1./p2;
rho(T>=0.001)=rho;

fac = 0.999195706402050
Tr=(T(T<0.001)+273.15-227.15)./227.15;
rho(T<0.001) = fac*b(1)*exp(-227.15*(b(3)*Tred+0.5D0*b(4)* ...
            Tr.*Tr+0.333*b(5)*Tr.*Tr.*Tr+0.25*b(6)* ...ç
            Tr.*Tr.*Tr.*+2.*b(2)*sqrt(Tr)));
