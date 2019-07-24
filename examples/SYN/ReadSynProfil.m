clc, clear, close all

% Bedingungen synthetisches Profil: GSTH nach Balling 1981, 
gradT = 30
surfaceT = 6
kappa       = 0.9301*10^(-6);           % mean(TD)-Templin  % thermal diffusivity in m^2/s
rho         = 10^6;
c           = 2.5;	                    % rho*c = 2.5 10^6 nach Balling 81: "normally values of rho*c range between 2 and 3 *10^6 J m-3 K-1, and we apply the constant value of 2.5"
lambda      = kappa*rho*c;  
dT_TA=importdata('DeltaTBallingA.csv'); % Tiefe,  disturbance, T-profile, T-profile minus disturbance
dT_TB=importdata('DeltaTBallingB.csv'); % Tiefe,  disturbance, T-profile, T-profile minus disturbance
GSTH_A=importdata('GSTHBallingA.csv');  % rel, abs
GSTH_B=importdata('GSTHBallingB.csv');  % rel, abs


N = length(dT_TA(:,1));

L=3;
Err = 1.;
Cov=CovarGauss(Err*ones(N,1),L);

C    = chol(Cov);


err_nor =  0.3*randn(N,1);
err_cor =  err_nor'*C;

figure 
plot(1:N,[err_nor(:)], '-b'); hold on
plot(1:N,[err_cor(:)], '-r','LineWidth',2); hold on

Tprofile = dT_TA(:,4)+err_cor;

err_nor =  0.3*randn(N,1);
err_cor =  err_nor'*C;
Tprofile = dT_TB(:,4)+err_cor;
