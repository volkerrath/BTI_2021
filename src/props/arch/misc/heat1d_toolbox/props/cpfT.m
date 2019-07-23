function [cpf,enthalp]=cpfT(T,P)
% calculates pure water heat capacity depending on
% pressure (P,in Pa), and temperature (T, in C). 
% 
% derived from the Formulation given in:
% Zylkovskij et al: Models and Methods Summary for the FEHMN Application,
%       ECD 22, LA-UR-94-3787, Los Alamos NL, 1994.
% 
% Range of validity: Pressures   0.01-110 MPa, 
%                    Temperature 15-350?C
%
% VR RWTH Aachen Univerity,   April 25, 2004

if nargin < 2, P=0.1*ones(size(T)); end
[n1 n2]=size(T);if n1==1, T=T'; end 
[n1 n2]=size(P);if n1==1, P=P'; end 

% Pressure in MPa 
P=P/1.e6;P2=P.*P;P3=P2.*P;P4=P3.*P;
T2=T.*T;T3=T2.*T;TP=P.*T;T2P=T2.*P;TP2=T.*P2;

pmin=0.001;pmax=110.;tmin=15.;tmax =360.;

% enthalpy of liquid ****
% numerator coefficients 1:10; denominator coefficients 11:20
c =[ 0.25623465e-3  0.10184405e-2  0.22554970e-4  0.34836663e-7 ...
      0.41769866e-2 -0.21244879e-4  0.25493516e-7  0.89557885e-4 ...
      0.10855046e-6 -0.21720560e-6  ...
      0.10000000e+1  0.23513278e-1  0.48716386e-4 -0.19935046e-8 ...
     -0.50770309e-2  0.57780287e-5  0.90972916e-9 -0.58981537e-4 ...
     -0.12990752e-7  0.45872518e-8];
          
a=c(1)+c(2)*P+c(3)*P2+c(4)*P3+c(5)*T+c(6)*T2+c(7)*T3+c(8)*TP+c(10)*T2P+c(9)*TP2;
b=c(11)+c(12)*P+c(13)*P2+c(14)*P3+c(15)*T+c(16)*T2+c(17)*T3+c(18)*TP+c(20)*T2P+c(19)*TP2;
da=c(5) +2*c(6) *T+3*c(7) *T2+c(8) *P+2*c(10)*TP+c(9)* P2;
db=c(15)+2*c(16)*T+3*c(17)*T2+c(18)*P+2*c(20)*TP+c(19)*P2;
b2=b.*b;cpf=da./b-a.*db./b2;
cpf=cpf*1.e6;
if nargout > 1; enthalp=a./b; end
    
        
        
