clc, clear, close all

% Bedingungen synthetisches Profil: GSTH nach Balling 1981, 
% gradT = 30
% surfaceT = 6
% kappa       = 0.9301*10^(-6);           % mean(TD)-Templin  % thermal diffusivity in m^2/s
% rho         = 10^6;
% c           = 2.5;	                    % rho*c = 2.5 10^6 nach Balling 81: "normally values of rho*c range between 2 and 3 *10^6 J m-3 K-1, and we apply the constant value of 2.5"
% lambda      = kappa*rho*c;  
dT_TA=importdata('DeltaTBallingA.csv'); % Tiefe,  disturbance, T-profile, T-profile minus disturbance
dT_TB=importdata('DeltaTBallingB.csv'); % Tiefe,  disturbance, T-profile, T-profile minus disturbance
GSTH_A=importdata('GSTHBallingA.csv');  % rel, abs
GSTH_B=importdata('GSTHBallingB.csv');  % rel, abs