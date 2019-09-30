%% Skript to compare Ballings results with Westaways results and the numerical scheeme
clc, clear, close all
%% Model auswählen
Model = 'B'

% Rock properties
kappa       = 0.9301*10^(-6);           % mean(TD)-Templin  % thermal diffusivity in m^2/s
rho         = 10^6;
c           = 2.5;	                    % rho*c = 2.5 10^6 nach Balling 81: "normally values of rho*c range between 2 and 3 *10^6 J m-3 K-1, and we apply the constant value of 2.5"
lambda      = kappa*rho*c;              % kappa*rho*c       % thermal conductivity in W/(m*K)
% depth 
z_min       = 0;                        % min depth in m
z_max       = 7000;                     % max depth in m
z_interv    = 50;
z           = (z_min:z_interv:z_max)';

if Model == 'A'
    %Temperature variation set up - Model A
    t_obs       = 0;                    % time of observation
    t_change_y    = [70000 10000 0];    % time when change in T happens
    T_change    = [-10 10 0];           % change in T at time t_change
elseif Model == 'A1'
    %Temperature variation set up - Model A
    t_obs       = 0;                    % time of observation
    t_change_y    = [10000 0];          % time when change in T happens
    T_change    = [10 0];               % change in T at time t_change
elseif Model == 'B'
    %Temperature variation set up - Model B
    t_obs       = 0;                    % time of observation
    t_change_y    = [190000 130000 70000 10000 8000 5000 1600 1300 1000 525 350 175 0];    % time when change in T happens
    T_change    = [-10 10 -10 10 2 -1 -1 0.5 1 -1 -1 0.5 0.5 ];% change in T at time t_change
elseif Model == 'B3'
    %Temperature variation set up - Model B3
    t_obs       = 0;                     % time of observation
    t_change_y    = [1000 525 350 175 0];% time when change in T happens
    T_change    = [1 -1 -1 0.5 0.5 0];   % change in T at time t_change
elseif Model == 'C'
    %Temperature variation set up - Model C
    t_obs       = 0;                     % time of observation
    t_change_y    = [65000 35000 10000 7000 2000 1000 600 500 125 100 0];    % time when change in T happens
    T_change    = [9 -3 12 1.5 -1.5 0.75 -1.5 -1.5 1.25 1.25 0]; % change in T at time t_change
end
%% Balling - without linear part - linear parts being steptified in half
t_change = t_change_y*365.25*24*60*60; % time in seconds

DeltaT_B      = zeros(length(z),length(t_change)); % matrix for Temperature disturbance over depth and time
DeltaG_B      = DeltaT_B;                          % matrix for disturbance on temperature gradient
HF_B          = DeltaT_B;                          % matrix for heat flow
DeltaT_cum_B  = DeltaT_B;                          % cumulative DeltaT
DeltaG_cum_B  = DeltaT_B;                          % cumulative DeltaG

for i = 1:length(t_change)-1
    for n = 1:length(z)
        DeltaT_B(n,i) = T_change(i)*erfc(z(n)/...  % step equation
            (4*kappa*(t_change(i)-t_obs))^(1/2));
        
        DeltaG_B(n,i) = (-T_change(i)/(pi*kappa*...% step equation
            (t_change(i)-t_obs))^(1/2)*exp(-z(n)^2/...
            (4*kappa*(t_change(i)-t_obs))));
        
        HF_B(n,i)     = -lambda*DeltaG_B(n,i);
        
        DeltaT_cum_B(n,i+1) = sum(DeltaT_B(n,1:i));
        DeltaG_cum_B(n,i+1) = sum(DeltaG_B(n,1:i));
    end
    
end


DeltaT_cum_B(:,1) = z;        % Tiefen in erste Spalte schreiben, zum besseren Lesen in der Matrix
DeltaG_cum_B(:,1) = z;        % Tiefen in erste Spalte schreiben, zum besseren Lesen in der Matrix

% close(figure(2))
% figure(2)
% subplot(1,2,1)
% titlesss = strcat('Model ',Model,' - DeltaTcum - Balling');
% title(titlesss);
% hold on
% for i = 3:length(DeltaT_cum_B(1,:)-1)
%     h2 = plot(DeltaT_cum_B(:,i),z,'--');
%     h2.Color(4)=0.5;
% end
% plot(DeltaT_cum_B(:,end),z,'k');
% set(gca,'YDir','reverse');
% grid on, box on 
% hold off
% 
% subplot(1,2,2)
% titlesss = strcat('Model ',Model,' - DeltaGcum - Balling');
% title(titlesss);
% hold on
% for i = 3:length(DeltaG_cum_B(1,:)-1)
%     h2 = plot(DeltaG_cum_B(:,i),z,'--');
%     h2.Color(4)=0.5;
% end
% plot(DeltaG_cum_B(:,end),z,'k');
% set(gca,'YDir','reverse');
% grid on, box on 
% % legend('location','southeast')
% hold off


%% Westaway

t_change_y  = fliplr(t_change_y);
T_change    = fliplr(T_change);
t_change    = t_change_y*365.25*24*60*60; % time in seconds
for i = 1:length(T_change)
    T_change_sum(i) = sum(T_change(i:end));
end


% output matrices
DeltaT_W      = zeros(length(z),length(t_change_y));  % matrix for Temperature disturbance over depth and time
DeltaG_W      = DeltaT_W;                             % matrix for disturbance on temperature gradient
HF_W          = DeltaT_W;                             % matrix for heat flow
DeltaT_cum_W  = DeltaT_W;                             % cumulative DeltaT
DeltaG_cum_W  = DeltaT_W;                             % cumulative DeltaG

for i = 2:length(t_change_y)
    for n = 1:length(z)
        DeltaT_W(n,i) = T_change_sum(i)*(erf(z(n)/(4*kappa*t_change(i-1))^(1/2))...
            -erf(z(n)/(4*kappa*t_change(i))^(1/2)));
        
        if i == 2
            DeltaG_W(n,i) = -1/(pi^(1/2))*T_change_sum(i)*(1/(kappa*t_change(i))^(1/2)...
                *exp(-z(n)^2/(4*kappa*t_change(i))));
        else
            DeltaG_W(n,i) = 1/(pi^(1/2))*T_change_sum(i)*(1/(kappa*t_change(i-1))^(1/2)...
                *exp(-z(n)^2/(4*kappa*t_change(i-1)))-1/(kappa*t_change(i))^(1/2)...
                *exp(-z(n)^2/(4*kappa*t_change(i))));
        end
        HF_W(n,i)     = -lambda*DeltaG_W(n,i);
        
        DeltaT_cum_W(n,i+1) = sum(DeltaT_W(n,1:i));
        DeltaG_cum_W(n,i+1) = sum(DeltaG_W(n,1:i));
    end
    
end


DeltaT_cum_W(:,1) = z;        % Tiefen in erste Spalte schreiben, zum besseren Lesen in der Matrix
DeltaG_cum_W(:,1) = z;        % Tiefen in erste Spalte schreiben, zum besseren Lesen in der Matrix

% close(figure(1))
% figure(1)
% subplot(1,2,1)
% titlesss = strcat('Model ',Model,' - DeltaTcum - Westaway');
% title(titlesss);
% hold on
% for i = 3:length(DeltaT_cum_W(1,:)-1)
%     h2 = plot(DeltaT_cum_W(:,i),z,'--');
%     h2.Color(4)=0.5;
% end
% plot(DeltaT_cum_W(:,end),z,'k');
%  set(gca,'YDir','reverse');
%  grid on, box on 
%  hold off
% 
% subplot(1,2,2)
% titlesss = strcat('Model ',Model,' - DeltaGcum - Westaway');
% title(titlesss);
% hold on
% for i = 3:length(DeltaG_cum_W(1,:)-1)
%     h2 = plot(DeltaG_cum_W(:,i),z,'--');
%     h2.Color(4)=0.5;
% end
% plot(DeltaG_cum_W(:,end),z,'k');
% set(gca,'YDir','reverse');
% grid on, box on 
% % legend('location','southeast')
% hold off


%% Numerical part

t_change_y_num  = max(t_change_y)-fliplr(t_change_y);
T_change_num    = fliplr(T_change_sum);
t_change_num    = t_change_y_num*365.25*24*60*60; % time in seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numerical parameters
nz      =   z_max/z_interv+1;        %   number of gridpoints in z-direction
endtime =   max(t_change_y_num);     %   total modeled time in a
nt      =   25000;                   %   number of timesteps to compute
day     =   3600*24;	             %   # seconds per day
year    =   day*365;	             %   # days per year
dt      =   endtime/(nt-1)*year;     %   time step  [s]
time    =   0:dt/year:endtime;

%initial temperature parameters
L       =   z_max;                   %   thickness of modeled domain  [m]
Tsurf   =   0;                       %   surface temperature          [K]
gradT   =   25;                      %   "normal" temperature gradient in K/km
Tbot    =   Tsurf+L*gradT/1000;      %   initial bottom temperature    [K] - normal thermal gradient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Spatial and temporal grids %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup grid
dz      =   L/(nz-1);	             %   grid spacing 
z       =   0:dz:L;                  %   grid

% Setup initial temperature profile and surface temperature history
DeltaT_cum_num       =   zeros(nz,nt);
for i = 2:length(t_change_y_num)
    [~,t_start] = min(abs(time-t_change_y_num(i-1)));
    [~,t_stop] = min(abs(time-t_change_y_num(i)));
    DeltaT_cum_num(1,t_start+1:t_stop) = T_change_num(i-1);
end

% Setup layer model with physical parameters
layers   = [z_min z_max];            %   depths of the bottom surface of the different layers
kappa_ini= kappa;                    %   thermal diffusivity of rocks [m2/s] according to layers
cp_ini   = c;                        %   heat capacity               [J/kg/K]
H_ini    = 0;%[0.5]*10^(-9);         %   crustal radioactive heat production [W/kg]

% setup arrays with the layer parameters, that are the same length as the
% depth array
kappa    = zeros(size(z));
cp       = kappa;
H        = kappa;
kappa(1) = kappa_ini(1);
cp(1)    = cp_ini(1);
H(1)     = H_ini(1);
for i = 2:length(layers)
    [~,l_top] = min(abs(z-layers(i-1)));
    [~,l_bottom] = min(abs(z-layers(i)));
    kappa(l_top+1:l_bottom) = kappa_ini(i-1);
    cp(l_top+1:l_bottom) = cp_ini(i-1);
    H(l_top+1:l_bottom) = H_ini(i-1);
end
% Setup time array to plot the thermal gradient
DeltaG_cum_num   =   zeros(nz,nt);	% numerical solution

% Stability
disp(['Stability: ',num2str(kappa_ini*dt/dz^2),' should be < 0.5.'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(figure(3))
tic
for n=2:nt	% Timestep loop

    % Compute new temperature 
    for i=2:nz-1
        DeltaT_cum_num(i,n) =  DeltaT_cum_num(i,n-1) + kappa(i)*dt*(DeltaT_cum_num(i+1,n-1)-2*DeltaT_cum_num(i,n-1)+DeltaT_cum_num(i-1,n-1))/(dz^2) + dt*H(i)/cp(i);
        DeltaG_cum_num(i,n)    =   (DeltaT_cum_num(i,n)-DeltaT_cum_num(i-1,n))/dz; % compute thermal gradient in K/km
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if mod(n,1000)==0
        % Plot solution
%         figure(3)%, clf
%         %set(subplot(2,2,[3 4]),'position',[0.05 0.05 0.8 0.2]);
% 
%         subplot(2,2,1)
%         hold on
%         if n == nt
%             plot(DeltaT_cum_num(:,n),z)
%         else
%             plot(DeltaT_cum_num(:,n),z,':')
%         end
%         axis ij
%         grid on, box on 
%         xlabel('T [^oC]')
%         ylabel('z [m]')
%         title('temperature T(z)')
%         for f = 2:length(layers)
%             plot([-max(DeltaT_cum_num(:,n))*1.5 max(DeltaT_cum_num(:,n))*1.5],[layers(f) layers(f)],'--k')
%         end
%         %xlim([-10 10])
%         
%         subplot(2,2,2)
%         hold on
%         grid on, box on  
%         if n == nt
%             plot(DeltaG_cum_num(:,n),z)
%         else
%             plot(DeltaG_cum_num(:,n),z,':')
%         end
%         axis ij
%         ylabel('z [m]')
%         xlabel('\DeltaT/\Deltaz [K/km]')
%         title('geothermal gradient')
%         for f = 2:length(layers)
%             plot([-max(DeltaG_cum_num(:,n))*1.5 max(DeltaG_cum_num(:,n))*1.5],[layers(f) layers(f)],'--k')
%         end
%         
%         
%         subplot(2,2,[3 4]);
%         hold on
%         plot(max(time)/1000-time/1000,DeltaT_cum_num(1,:))
%         plot(max(time)/1000-time(n)/1000,DeltaT_cum_num(1,n),'*')
%         ylim([min(DeltaT_cum_num(1,:))-5 max(DeltaT_cum_num(1,:))+5])
%         title('Surface temperature history')
%         xlabel('time [ka]');
%         ylabel('T_{surface} [^oC]')
%         set(gca,'XDir','reverse');
%         grid on, box on 
    end
    
%     drawnow
end
toc

%% comparison of results: Balling-Westaway and Balling-numerical
for i = 1:length(DeltaT_cum_B(:,end))
    Diff_DT_BW(i,1) = DeltaT_cum_B(i,end)-DeltaT_cum_W(i,end);
    Diff_DT_Bnum(i,1) = DeltaT_cum_B(i,end)-DeltaT_cum_num(i,end);
    Diff_DT_Wnum(i,1) = DeltaT_cum_W(i,end)-DeltaT_cum_num(i,end);
    
    Diff_DG_BW(i,1) = DeltaG_cum_B(i,end)-DeltaG_cum_W(i,end);
    Diff_DG_Bnum(i,1) = DeltaG_cum_B(i,end)-DeltaG_cum_num(i,end);
    Diff_DG_Wnum(i,1) = DeltaG_cum_W(i,end)-DeltaG_cum_num(i,end);
end

close(figure(5))
fig5 = figure(5);
set(fig5,'Position',[10 10 800 1100])
subplot(3,3,1,'position',[0.05 0.7 0.27 0.26])
title({'Difference in \DeltaT'})
hold on
zeroline = zeros(length(Diff_DT_BW),1);
fill([Diff_DT_BW';zeroline'],[z ;z],[0.5 0.5 0.5],'handlevisibility','off')
plot(DeltaT_cum_B(:,end),z,'r','DisplayName','Balling (1981)')
plot(DeltaT_cum_W(:,end),z,'b--','DisplayName','Westaway (2013)')
set(gca,'YDir','reverse');
ylim([0 6000])
xlim([-7 3])              % [-7 3] fuer A und B, [-20 20] fuer C
xlabel('[K]')
grid on, box on, legend('location','southwest')

subplot(3,3,4,'position',[0.05 0.35 0.27 0.26])
title({'Difference in \DeltaG'})
hold on
fill([Diff_DG_BW';zeroline'],[z ;z],[0.5 0.5 0.5])
plot(DeltaG_cum_B(:,end),z,'r')
plot(DeltaG_cum_W(:,end),z,'b--')
set(gca,'YDir','reverse');
ylim([0 6000])
xlim([-7 3]*10^-3)        % [-7 3] fuer A und B, [-20 20] fuer C
xlabel('[K/m]')
grid on, box on

subplot(3,3,2,'position',[0.37 0.7 0.27 0.26])
title({'Difference in \DeltaT'})
hold on
fill([Diff_DT_Bnum';zeroline'],[z ;z],[0.5 0.5 0.5],'handlevisibility','off')
plot(DeltaT_cum_B(:,end),z,'r','DisplayName','Balling (1981)')
plot(DeltaT_cum_num(:,end),z,'g--','DisplayName','numerical FTCS')
set(gca,'YDir','reverse');
ylim([0 6000])
xlim([-7 3])              % [-7 3] fuer A und B, [-20 20] fuer C
grid on, box on, legend('location','southwest')
xlabel('[K]')
set(gca,'yticklabels','')

subplot(3,3,5,'position',[0.37 0.35 0.27 0.26])
title({'Difference in \DeltaG'})
hold on
fill([Diff_DG_Bnum(2:end)';zeroline(2:end)'],[z(2:end) ;z(2:end)],[0.5 0.5 0.5])
plot(DeltaG_cum_B(:,end),z,'r')
plot(DeltaG_cum_num(2:end,end),z(2:end),'g--')
set(gca,'YDir','reverse');
ylim([0 6000])
xlim([-7 3]*10^-3)        % [-7 3] fuer A und B, [-20 20] fuer C
grid on, box on
xlabel('[K/m]')
set(gca,'yticklabels','')

subplot(3,3,3,'position',[0.69 0.7 0.27 0.26])
title({'Difference in \DeltaT'})
hold on
fill([Diff_DT_Wnum';zeroline'],[z ;z],[0.5 0.5 0.5],'handlevisibility','off')
plot(DeltaT_cum_W(:,end),z,'b','DisplayName','Westaway (2013)')
plot(DeltaT_cum_num(:,end),z,'g--','DisplayName','numerical FTCS')
set(gca,'YDir','reverse');
ylim([0 6000])
xlim([-7 3])              % [-7 3] fuer A und B, [-20 20] fuer C
grid on, box on, legend('location','southwest')
xlabel('[K]')
set(gca,'yticklabels','')

subplot(3,3,6,'position',[0.69 0.35 0.27 0.26])
title({'Difference in \DeltaG'})
hold on
fill([Diff_DG_Wnum(2:end)';zeroline(2:end)'],[z(2:end) ;z(2:end)],[0.5 0.5 0.5])
plot(DeltaG_cum_W(:,end),z,'b')
plot(DeltaG_cum_num(2:end,end),z(2:end),'g--')
set(gca,'YDir','reverse');
ylim([0 6000])
xlim([-7 3]*10^-3)        % [-7 3] fuer A und B, [-20 20] fuer C
grid on, box on
xlabel('[K/m]')
set(gca,'yticklabels','')

subplot(3,3,[7:9],'position',[0.1 0.05 0.8 0.2])
hold on
plot(max(time)/1000-time/1000,DeltaT_cum_num(1,:),'k')
% plot(max(time)/1000-time(n)/1000,DeltaT_cum_num(1,n),'*')
ylim([min(DeltaT_cum_num(1,:))-5 max(DeltaT_cum_num(1,:))+5])
Titeltext = strcat('Surface temperature model: ',Model,' from Balling (1981)');
title(Titeltext)
xlabel('time [ka]');
ylabel('T_{surface} [^oC]')
set(gca,'XDir','reverse');
grid on, box on

savePic = strcat('TdisturbanceComparison_Model',Model,'_new');
print(savePic,'-dpdf')
