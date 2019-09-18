% v.r. 8/01                    last change: Aug.18,2001 
clear all;close all;

disp('Test transient heat conduction')
clear all;close all;
% generate weakly graded mesh 
dz=logspace(1,1.5,161);z=[0,cumsum(dz)];nz=length(z);
% generate constant mesh 
% z=[0:10:3000];dz=diff(z);nz=length(z);
% generate log  
[ip]   =setmodel('./testmodel_synthetic',z);  
ic=ip;
% define  beginninmg and length of log
logstart=30;
logend=2000;
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

Ts=pt(it);
delT=Ts(1); % -T0(1);
Ti=T0;



theta=0.5*ones(1,nt);
theta(1:10)=1.;
Tobstrue=heat1dnt(par,ic,dz,dt,Ti,Ts,theta);

ErrT=0.1;
% randn('state',sum(100*clock));
Tobs=Tobstrue(:,nt)+ErrT*randn(nz,1);


figure;
plot(Tobs(id),-z(id),'-g');hold on;
xlim([0 120]);ylim([-2500 0]);xlabel('T (^\circ C)');
legend('Observed ',3)
ylabel('z (m)');grid on


% define parameters for inversion loop

dp=0.01;
tolrms= 0.5;           % tolerance for rms
maxiter=20;            % maximal number of iterations
tau0 = 0.1;		% intial/center regularization parameters
tau1 =2.0;		% Gradient
tau2 =1.0;	        % Lapalcian
k = 12;                 % number of CGLS iterations

fct=[0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000];
num=(length(fct)-1)/2+1;

% define regularization matrix 
np= length(pt);
L0 = reg1d(np,'l0'); 
L1 = reg1d(np,'l1');
L2 = reg1d(np,'l2');

% ;
Wd=diag(1./(ErrT));

% a-priori model

m_apr=zeros(1,np);
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
   S=Wd*J(id,:);

          if length(fct)>0,    
           for itau =1:length(fct)   
             S_Aug=[S;...
                   sqrt(tau1*fct(itau))*L1;
                   sqrt(tau0)          *L0];
             resid_Aug=[r;...
                   -sqrt(tau1*fct(itau))*L1*(m'-m_apr');...
                   -sqrt(tau0          )*L0*(m'-m_apr')];
             [delta_m] = cgls(S_Aug,resid_Aug,k,0.001,1);        
             
             tau_lc(itau,iter)= tau1*fct(itau);
	     Ts=m(it)+delta_m(it)';
             T_tau=heat1dnt(par,ic,dz,dt,Ti,Ts,theta);
             rho_lc(itau,iter)= norm((Tobs(id)-T_tau(id,nt))/ErrT);
%             eta_lc(itau,iter) = norm( L1'*((m'+delta_m)-m_apr'));
             eta_lc(itau,iter) = norm( L1'*((m'+delta_m)));
          end;
          
          %plot L-curve
          
          
          % now choose tau/k?
           
         end; 

 
  
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
   
 end
 
figure;
stairs(t/year,m(it))
set(gca,'xscal','log')
xlabel('time (a)');ylabel('\Delta T');
title(' Regularized Inversion : final model')
text(-1000,-2.0,['temperature log: ',num2str(logstart),' - ',num2str(logend),' m'])
text(-1000,-3.0,['\tau_{l1} = ',num2str(tau1)])
text(-1000,-4.0,['Error_T  = ',num2str(ErrT),' K '])
text(-1000,-5.0,['rms  = ',num2str(rms),' K '])
grid on;
print( '-depsc','-r600',strcat('test',num2str(logend),'_model.eps'))
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

print( '-depsc','-r600',strcat('test',num2str(logend),'_iter.eps'))
% residuals
figure;
plot(r,-z(id),'-r');
hold on;
xlim([-5 5]);ylim([-2500 0]);
xlabel('\Delta T (K)');
ylabel('z (m)');grid on
title(' Regularized Inversion : final model')
print( '-depsc','-r600',strcat('test',num2str(logend),'_resid.eps'))

figure;
loglog(rho_lc,eta_lc,'-g');hold on;
loglog(rho_lc(num,:),eta_lc(num,:),'ro')
loglog(rho_lc(1,:),eta_lc(1,:),'bv')
loglog(rho_lc(length(fct),:),eta_lc(length(fct),:),'mv')
loglog(rho_lc(:,maxiter),eta_lc(:,maxiter),'-r')
%xlim([0.0001 1000.]);ylim[0.0001,10000.]);     
ylabel('| L*(m-m_{apr} |^2 (parameter seminorm)');
xlabel('| r |^2 (residual norm)' );
title('L-curve'); % legend('1',4);
grid on;
print( '-depsc','-r600',strcat('test',num2str(logend),'_lcurve.eps'))
         
