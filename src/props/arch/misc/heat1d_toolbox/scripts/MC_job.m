% PARAMETERS FOR INVERSION
clear
close all
clc
%borehole={'3200','3209','3356', '3359'};
% borehole={'SG3_97','2400','2915','3200','3209','3356','3359'};



% PARAMETERS FOR INVERSION


tolrms= 0.1;
tolmc= 1;% tolerance for rms in Monte Carlo
maxiter=10000;             % maximal number of iterations

nsteps=24; % steps for inversion grid
delta_m=8; % max. Intervall for model generation
maxdelta=10; % absolute value of maximum temperature interval
middle=-5; % Mittelwerte des Temperaturintervalls
q_range=0.03;
qmin=0.02;
gt_range=2;
gtmin=5;

freeze='no';



%load Kola_apriori

%m_apr=mod;
borehole={'Elk1'};
number=borehole{1};
NAME=strcat(['MC-' number '-' num2str(maxiter) '-'  num2str(nsteps) '-' freeze]);

MC_Poland;

freeze='yes';
borehole={'Elk1'};
number=borehole{1};
NAME=strcat(['MC-' number '-' num2str(maxiter) '-'  num2str(nsteps) '-' freeze]);

MC_Poland;