% v.r. 8/01                    last change: Aug.18,2001 
close all;
% clear all
disp(['Test for loglenght: ', num2str(logend),' - ', ...
      num2str(logstart), ' m ' ] )

% generate weakly graded mesh 
dz=logspace(1,1.5,161);z=[0,cumsum(dz)];nz=length(z);
% generate constant mesh 
% z=[0:10:3000];dz=diff(z);nz=length(z);
% generate log  
[ip]   =setmodel('./testmodel_synthetic',z);  
ic=ip;
% define  beginninmg and length of log
%logstart=30;logend=500;
id=z > logstart & z < logend;

molasse_prop

np=length(pc);

% calculate inital model  and get data

par=[pc ph qb gt];

%[Tobs]=    setdata('./logs/Moenchsrot20',z);   
%[ip  ]   =setmodel('./logs/Moenchsrot20',z);  
%ic=ip;id=~isnan(Tobs);

T0=heat1dns(par,ic,dz); 


% set  parameter for time-dependent calculation
clear par;
par=[pc ph pr qb];


% set paloeclimate
year=31557600;
tstart=30000;tend=1;
t=-logspace(log10(tstart),log10(tend),101);
   % t=[-20000:100:0];
nt=length(t);t=t*year;dt=diff(t);
[pt,it]=paleo_haenel(t,gt);

figure;
stairs(t/year,pt(it))
set(gca,'xscal','log')
xlabel('time (a)');ylabel('\Delta T');
title(' Paleoclimate from Haenel (1988)')
grid on;
print( '-depsc','-r600',strcat('paleoclimate_haenel88.eps'))
close;

Ts=pt(it);
delT=Ts(1); % -T0(1);
Ti=T0;



theta=0.5*ones(1,nt);
theta(1:10)=1.;
Tobstrue=heat1dnt(par,ic,dz,dt,Ti,Ts,theta);

ErrT=0.25;
% randn('state',sum(100*clock));
Tobs=Tobstrue(:,nt)+ErrT*randn(nz,1);


figure;
plot(Tobs(id),-z(id),'-g');hold on;
xlim([0 120]);ylim([-2500 0]);xlabel('T (^\circ C)');
legend('Observed ',3)
ylabel('z (m)');grid on
close;

% define parameters for inversion loop

dp=0.01;
tolrms= 0.5;           % tolerance for rms
maxiter=16;            % maximal number of iterations
tau0 = 0.1;		% intial/center regularization parameters
%tau1 =1.0;		% Gradient
tau2 =1.0;	        % Lapalcian
k = 12;                 % number of CGLS iterations

% define regularization matrix 
np= length(pt);
L0 = reg1d(np,'l0'); 
L1 = reg1d(np,'l1');
%L2 = reg1d(np,'l2');

% ;
nd=length(Tobs);
ACov=ErrT^2*ones(1,nd);
Wdfi=diag(1./sqrt(ACov),0);
Cdfi=diag(1./ACov,0);

% a-priori model

m_apr=gt*ones(1,np);
m_0  = m_apr;
m    = m_0;

Ts=m(it);
Tcalc=heat1dnt(par,ic,dz,dt,Ti,Ts,theta);
resid=Tobs(id)-Tcalc(id,nt);r=resid/ErrT;

for iter=1:maxiter 
 
%  calculate residual  
    
  Ts=m(it);
  
  Tcalc=heat1dnt(par,ic,dz,dt,Ti,Ts,theta);

  resid=Tobs(id)-Tcalc(id,nt);r=resid/ErrT;
  rms=norm(r)/sqrt(length(r));
  
%  plot(r,-z(id),':r');hold on;
  disp([ 'rms for iteration ',num2str(iter), ...
       ' = ',num2str(rms)])
       
  if rms < tolrms, break; end


  [J]=sensfdt(par,ic,dz,dt,T0,m,it,dp,theta);
  
  S=Wdfi(id,id)*J(id,:);
  S_Aug=[              S;...
           sqrt(tau1)*L1;
           sqrt(tau0)*L0 ];
  resid_Aug=[                       r;...
             -sqrt(tau1)*L1*(m'-m_apr');...
             -sqrt(tau0)*L0*(m'-m_apr')];
             [delta_m] = cgls(S_Aug,resid_Aug,k,0.001,1);        
  m = m + delta_m' ;
  mall(iter,:)=m(1,:);
  rmsall(iter)=rms;
  tall(iter,:)=Tcalc(:,nt)'; 
 end

filename=strcat('test_',num2str(logstart),'-',num2str(logend),'_',num2str(tau1),'.mat');
save(filename,'mall','tall','rmsall' ) 

 
figure;
stairs(t/year,m(it))
set(gca,'xscal','log')
xlabel('time (a)');ylabel('\Delta T');ylim([-6 6]);
title(' Regularized Inversion : final model')
text(-1000,-0.5,['temperature log: ',num2str(logstart),' - ',num2str(logend),' m'])
text(-1000,-1.5,['\tau_{l1} = ',num2str(tau1)])
text(-1000,-2.5,['Error_T  = ',num2str(ErrT),' K '])
text(-1000,-3.5,['rms  = ',num2str(rms),' K '])
grid on;
print( '-depsc','-r600',strcat('test_',num2str(logstart),'-',num2str(logend),'_',num2str(tau1),'_model.eps'))
close
%convergence 
figure;
subplot(1,2,1)
plot(rmsall); grid on;
xlabel('iteration');ylabel('rms');
subplot(1,2,2) 
semilogy(rmsall); grid on;
text(10,5.,['rms  = ',num2str(rms)])
xlabel('iteration');ylabel('rms');
suptitle(' Regularized Inversion : final model')
print( '-depsc','-r600',strcat('test_',num2str(logstart),'-',num2str(logend),'_',num2str(tau1),'_iter.eps'))
close;

% residuals
figure;
plot(r,-z(id),'-r');
hold on;
xlim([-5 5]);ylim([-2500 0]);
xlabel('\Delta T (K)');
ylabel('z (m)');grid on
title(' Regularized Inversion : final model')
print( '-depsc','-r600',strcat('test_',num2str(logstart),'-',num2str(logend),'_',num2str(tau1),'_resid.eps'))
close;

% data covariance apriori     
%Cd = inv(Wd'*Wd');
Cd = inv(Cdfi(id,id));
% [Cmm_post,R_mm,R_dd]=function resol(S,Cmm_prior,Cdd_prior)
S_Aug =[S ;sqrt(tau1)*L1 ;sqrt(tau0)*L0];
S_Augt=[S',sqrt(tau1)*L1',sqrt(tau0)*L0'];
GT_Aug=inv(S_Augt*S_Aug)*S';
% resolution matrix covariance matrix I 
%(from generalized inverse Nolet(99))
%Cmm=GT_Aug*Cd*GT_Aug';
Cmm=GT_Aug*GT_Aug';
figure;imagesc(Cmm);colorbar;title('C_{pp}^{aposteriori} based on Generalized Inverse');  
print( '-depsc','-r600',strcat('test_',num2str(logstart),'-',num2str(logend),'_',num2str(tau1),'_cpp.eps'))



Errap=sqrt(diag(Cmm));
ty=t/year;
mod=m(it);
err=Errap(it);
filename=strcat('Mod_',num2str(logstart),'-',num2str(logend),'_',num2str(tau1),'.mat');
save(filename,'ty','mod','err') 

figure;
stairs(t/year,m(it)); hold on;
stairs(t/year,m(it)+2*Errap(it)',':g')
stairs(t/year,m(it)-2*Errap(it)',':g')
set(gca,'xscal','log')
xlabel('time (a)');ylabel('\Delta T');ylim([-6 6]);
title(' Regularized Inversion : final model')
text(-1000,-0.5,['temperature log: ',num2str(logstart),' - ',num2str(logend),' m'])
text(-1000,-1.5,['\tau_{l1} = ',num2str(tau1)])
text(-1000,-2.5,['Error_T  = ',num2str(ErrT),' K '])
text(-1000,-3.5,['rms  = ',num2str(rms),' K '])
grid on;
print( '-depsc','-r600',strcat('test_',num2str(logstart),'-',num2str(logend),'_',num2str(tau1),'_model.eps'))
    
figure; stairs(t/year,Errap(it),'-g');
set(gca,'xscal','log');
xlabel('time (a)');ylabel('\Delta T');grid on;
title(' Regularized Inversion : Error based on C_{ii}^{aposteriori}')
print( '-depsc','-r600',strcat('test_',num2str(logstart),'-',num2str(logend),'_',num2str(tau1),'_error.eps'))
