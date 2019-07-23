% TESTHEAT1D
%
% Script um HEAT1DN zu testen
% v.r. 8/01                    last change: Aug.16,2001 
clear all;close all;

disp('Test transient heat conduction')

% physical parameters
z=[0:10:5000];       			
dz=diff(z);
nz=length(z);
gt= 10.  ; 
qb= 0e-3;               				
pc =  [2];  
ph =  0. ;
pr =  2.e6;
np=length(pc);
ic=ones(1,length(dz)); 
par=[pc ph pr qb];




% initial values
T0= zeros(nz,1);T0(1,1)=gt;

% time steps &time dependent surface temperatures
year=31557600;
t=[0:100:30000];nt=length(t);
t=t*year;dt=diff(t);
[Ts]=gt*ones(1,nt);


beta=.667;

%
tic
T=heat1dnt(par,ic,dz,dt,T0,Ts,beta);
toc
% solution bear discontinuity 
figure
plot(T(2:10,1:50)');
title(' Temperature node 2-10');grid on;
ylabel('Temperature [K]');xlabel('time step');
legend(' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9','10',4)  ;
text(30,1, strcat('Theta = ',num2str(beta))) ;  

% analytical solution 
kappa=pc/pr;
scal(1)=1;
scal(2:nt)=0.5./sqrt(kappa*t(2:nt));
Tc(:,1)=T0;Tc(1,1)=gt;    
for i=2:nt
 Tc(:,i)=gt*erfc(z'*scal(i));
end

% plot results and analytical solution 
figure;
plot(T(:,[1:50:300]),-z);hold on;
title(' Temperature (FD)');grid on;
xlabel('Temperature [K]');ylabel('Depth [m]');
plot(Tc(1:10:nz,[1:50:300]),-z(1:10:nz),'o')
legend('t =  0 ka','  =  5 ka','  =  10 ka','  =  15 ka','  =  20 ka', ...
       '  =  25 ka',4)  

text(3,-2500, ...
 'Analytical Solution: Turcotte & Schubert (1979)')
text(3,-3000, strcat('Theta = ',num2str(beta))) ;  

% plot residuals  
Tres= T(:,[1:50:300])-Tc(:,[1:50:300]);figure;plot(Tres,-z); grid on;    
legend('t =  0 ka','  =  5 ka','  =  10 ka','  =  15 ka','  =  20 ka', ...
       '  =  25 ka',3)  
text(-0.0095,-3200, strcat('Theta = ',num2str(beta))) ; 
 
figure;
subplot(1,2,1);
TT=T-Tc;imagesc(TT);colorbar   
title(' T-T_a');xlabel('timestep');ylabel('Depth [nodes]');
subplot(1,2,2);
TT=log10(abs(T-Tc)+eps);imagesc(TT);colorbar   
title(' log (|T-T_a|)');xlabel('timestep');
suptitle(strcat('Theta = ',num2str(beta)));

