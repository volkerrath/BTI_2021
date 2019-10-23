function [cpf]=cpfT(T,P)
% calculates pure water heat capacity depending on
% pressure (P,in Pa), and temperature (T, in C). 
% 
% derived from the Formulation given in:
% Zylkovskij et al: Models and Methods Summary for the FEHMN Application,
%       ECD 22, LA-UR-94-3787, Los Alamos NL, 1994.
% 
% Range of validity: Pressures   0.01-110 MPa, 
%                    Temperature 15-350�C
%
% VR RWTH Aachen Univerity,   April 25, 2004

if nargin < 2, P=0.1*ones(size(T)); end
T=T(:); P=P(:);

% Pressure in MPa 
P=P/1.e6;P2=P.*P;P3=P2.*P;P4=P3.*P;
         T2=T.*T;T3=T2.*T;TP=P.*T;T2P=T2.*P;TP2=T.*P2;

pmin=0.001;pmax=110.;tmin=15.;tmax =360.;
% enthalpy of liquid ****
% numerator coefficients 1:10; denominator coefficients 11:20
ac =[ 0.25623465e-3  0.10184405e-2  0.22554970e-4  0.34836663e-7 ...
      0.41769866e-2 -0.21244879e-4  0.25493516e-7  0.89557885e-4 ...
      0.10855046e-6 -0.21720560e-6  ...
      0.10000000e+1  0.23513278e-1  0.48716386e-4 -0.19935046e-8 ...
     -0.50770309e-2  0.57780287e-5  0.90972916e-9 -0.58981537e-4 ...
     -0.12990752e-7  0.45872518e-8];
          
a=ac(1)+ac(2)*P+ac(3)*P2+ac(4)*P3+ac(5)*T+ac(6)*T2+ac(7)*T3+ac(8)*TP+ac(10)*T2P+ac(9)*TP2;
b=ac(11)+ac(12)*P+ac(13)*P2+ac(14)*P3+ac(15)*T+ac(16)*T2+ac(17)*T3+ac(18)*TP+ac(20)*T2P+ac(19)*TP2;
da=ac(5) +2*ac(6) *T+3*ac(7) *T2+ac(8) *P+2*ac(10)*TP+ac(9)* P2;
db=ac(15)+2*ac(16)*T+3*ac(17)*T2+ac(18)*P+2*ac(20)*TP+ac(19)*P2;
b2=b.*b;cpf=da./b-a.*db./b2;
cpf=cpf*1.e6;

% if nargout > 1; enthalp=a./b; end
 
% after Speedy (1987) for T < 0 to -46 C
%      Cx    B0     B1      B2      B3        B4
bc = [14.2 25.952 128.281 -221.405 196.894 -64.812];
ic=T<0.;
T(T<-45.)=-45;
r=(T(ic)+273.15-227.15)/227.15;
r2=r.*r;r3=r2.*r;r4=r3.*r;
cpf(ic)=bc(1)./sqrt(r)+ bc(2) + bc(3)*r + bc(4)*r2 + bc(5)*r3 + bc(6)*r4;
cpf(ic)=cpf(ic)*0.99048992406520*1000/18.;
        
        
