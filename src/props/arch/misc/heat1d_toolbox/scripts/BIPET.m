% BAYESVINV_PETRO performs petrophysical inversion in parameter space
% V. R., March 2, 2002 
close all;clear all
disp(['Nonlinear Bayesian Inversion for thermal parameters'] )
disp(['V. R., ',date] )
ALG='BIPET'

% GENERATE WEAKLY GRADED SPATIAL MESH 
dz=logspace(1,1.76338,185);z=[0,cumsum(dz)];nz=length(z);
disp(['  ']);disp([' ...set up spatial grid   ']);
disp(['spatial mesh : ', num2str(nz),   ' nodes,  z:  ',...
                         num2str(z(1)),' - ', ...
                         num2str(z(nz)), ' m ' ] )
disp(['               ', num2str(nz-1), ' cells, dz: ',...
                         num2str(dz(1)),' - ', ...
                         num2str(dz(nz-1)), ' m ' ] )

% GENERATE TEMPORAL MESH
year=31557600;
tstart=150000*year;tend=10*year;
t=-logspace(log10(tstart),log10(tend),150);
nt=length(t);dt=diff(t);
disp([' ']);disp([' ...set up temporal grid   ']);
disp(['temporal mesh: ', num2str(nz),   ' nodes,  t:  ',...
               num2str(t(1)/year),' - ', ...
               num2str(t(nt)/year), ' a ' ] )
disp(['               ', num2str(nz-1), ' cells, dt: ',...
               num2str(dt(1)/year),' - ', ...
               num2str(dt(nt-1)/year), ' a ' ] )


% GENERATE MODEL PARAMETER 
disp(['   ']);disp([' ...set up model   ']);
[ic]   =setmodel('testmodel_synthetic',z);  
set_properties_molasse;units=length(pc);
%plot_true_model 

% SET DISTURBED_LOG
disp(['   ']);disp([' ...set up disturbed log   ']);
set_disturbed_log_no_paleo


%plot_true_climate
%plot_disturbed_log


% INVERSION SETUP 
disp(['   ']);disp([' ...set up inversion parameters  ']);
disp([' ']);
% define parameters for inversion 
tolrms= 0.5;            % tolerance for rms
maxiter=5 ;           % maximal number of iterations
%beta = 0.001;		% zeroth regularization (marquardt beta)
tau  = .0;		% weighting for paramer covariance
k    = 6;              % number of CGLS iterations

% define parameters for sensitivity calculations
dp=0.01;                % fixed perturbation for sensitivity calculation
maxitdp=5;                 % maximal number of iterations
toldp=1e-7;                % maximal number of iterations

% define parameters for forward modeling
tolnl=.001;                % fixed perturbation for sensitivity calculation
maxitnl=2;               % maximal number of iterations

% SET APRIORI MODEL & ERRORS

mpc_apr=[ log(3*ones(1,units))];	epc_apr=[ 0.1*ones(1,units)];
%mpc_apr=[ log(pc)];	                epc_apr=[ 0.001*ones(1,units)];

mph_apr=[ log(ph)];	                eph_apr=[ 0.001*ones(1,units)];
%mph_apr=[ log(1.e-6*ones(1,units))];	eph_apr=[ 0.5*ones(1,units)];

mpr_apr=[ log(pr)];	                epr_apr=[ 0.001*ones(1,units)]; 
%mpr_apr=[ log(2.e6*ones(1,units))];	epr_apr=[ 0.6*ones(1,units)]; 

mqb_apr=[ log(85e-3)];			eqb_apr=[ 0.0001];

m_apr =[mpc_apr mph_apr mpr_apr mqb_apr pcA pcB];
e_apr =[epc_apr eph_apr epr_apr eqb_apr];
w     =[1*ones(1,3*units+1)];

% plot_apriori_model

% SETUP DATA COVARIANCE A PRIORI ANS WEIGHTS
nd=length(Tobs(id));
Vd=ErrT^2*ones(1,nd);
Wd =diag(1./sqrt(Vd),0);
Cdi=diag(1./Vd,0);

% SETUP PARAMETER COVARIANCE MATRIX A PRIORI ANS WEIGHTS
np=3*units+1;
Ep=e_apr';
Vp=Ep.^2;
Wp =diag(1./sqrt(Vp),0);
Cpi=diag(1./Vp,0);



vcp=[];vhp=[];vrp=[];vgt=[];     % used for plotting

% important index arrays for indirekt addresses:
% id        points to nodes defining data
% ic        associates parameter values of diffusivity and production to calls
% ip        points to the parameters associated with the present inversion 

% ITERATION
m  = m_apr;
for iter=1:maxiter 
%  calculate residual
   p=m2p(m,units); 
   Tcalc=heat1dnt_nl(p,units,ic,dz,dt,Ti,Ts,theta,maxitnl,tolnl);
   resid=Tobs(id)-Tcalc(id,nt);r=resid/ErrT;
   rms=norm(r)/sqrt(length(r));
   disp([ 'rms for iteration ',num2str(iter-1),' = ',num2str(rms)])
   if rms < tolrms, break; end
%     [J,ip]=sensfdt_pet(m,units,ic,dz,dt,T0,Ts,...
%                         theta,maxitnl,tolnl,dp);
     [J,ip]=sensfdt_pet0(m,units,ic,dz,dt,T0,Ts,...
                         theta,maxitnl,tolnl,dp,maxitdp,toldp);
     S=J(id,ip)*Wd(ip,ip);
     S_Aug=[     S;...
                 sqrt(tau) *Wp(ip,ip)];
     resid_Aug=[ r;...
                -sqrt(tau) *Wp(ip,ip)*(m(ip)'-m_apr(ip)')];
     [delta_m] = cgls(S_Aug,resid_Aug,k,0.001,1);        
     m(ip) = m(ip) + delta_m' ;
     m_iter(iter,:)=m(1,:);
     rms_iter(iter)=rms;
     d_iter(iter,:)=Tcalc(:,nt)'; 
% used for plotting
vcp=[exp(m(        1:  units))' vcp];     % used for plotting
vhp=[exp(m(1*units+1:2*units))' vhp];     % used for plotting
vrp=[exp(m(2*units+1:3*units))' vrp];     % used for plotting
vgt=[exp(m(3*units+1)) vgt];

plot_iterate

end

disp(['   ']);disp([' ...calculate aposteriori quantities ']);
% CALCULATE COVARIANCES A POSTERIORI
S_Aug =[S ;sqrt(tau)*Wp(ip) ];
S_Augt=[S' sqrt(tau)*Wp(ip)'];
% GENERALIZED INVERSE
GT_Aug=inv(S_Augt*S_Aug)*S';
% PARAMETER & DATA COVARIANCE MATRIX A POSTERIORI
% (FROM GENERALIZED INVERSE NOLET(99))
Cmm=GT_Aug *GT_Aug'; 
Cdd=GT_Aug'*GT_Aug;


% PLOTS 
disp(['   ']);disp([' ...plot results for ',ALG]);
plot_aposteriori_covar
plot_aposteriori_model
plot_model_response

% SAVE RESULTS
disp(['   ']);disp([' ...save results  ']);
filename=strcat('BIPET_Results.mat');
save(filename,'m_iter','d_iter','rms_iter' ) 
epc=exp(m(        1:  units));delpc=pc-epc;dpc=norm(delpc);
eph=exp(m(  units+1:2*units));delph=ph-eph;dph=norm(delph);
epr=exp(m(2*units+1:3*units));delpr=pr-epr;dpr=norm(delpr);
eqb=exp(m(3*units+1));delqb=qb-eqb;dqb=norm(delqb);
disp([ '|delta m|   :  ',num2str(dpc),'  ', ...
                         num2str(dph),'  ', ...
			 num2str(dpr),'  ', ...
			 num2str(dqb)])
