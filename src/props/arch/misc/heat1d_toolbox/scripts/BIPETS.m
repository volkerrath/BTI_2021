% BAYESVINV_PETRO performs petrophysical inversion in parameter space


% important index arrays for indirekt addresses:
% id        points to nodes defining data
% ip        associates parameter values of petrophysical to cells
% 
% structures used:
% unit      information associated with geological units (layers)
% model     information defining model for given borehole
% obs       information cncerning the observations in borehole 
% 
%
% V. R. May 6, 2003

close all;clear all
disp(['Nonlinear Bayesian Inversion for Thermal Parameters - Steady-State'] )
disp(['V. R., ',date] ); 
NAME=strcat(['test']);

disp([' ']);


% READ MODEL PARAMETER & DATA 
borehole={'test_1', 'test_2' 'test_3'};
% borehole={'synthetic_1'};
nb=length(borehole);

% INVERSION SETUP 
disp(['   ']);disp([' ...set up inversion parameters  ']);
disp([' ']);

% define parameters for inversion 
tolrms  =   0.1;                  % tolerance for rms
maxiter =   300 ;                 % maximal number of iterations

% define parameters for CGLS 
kc      =   30;                   % number of CGLS iterations
stoprule=    2;                   % 0=none 1=norm/normold 2=acb 
reorth  =    1;
tol     =    1;
debug   =    1;
% define parameters for sensitivity calculations
dp=0.01;                    % fixed perturbation for sensitivity calculation
maxitdp=1;                   % maximal number of iterations
toldp=1e-4;                  % maximal number of iterations

% define parameters for nonlinear forward modeling
tolnl=.001;                % fixed perturbation for nonlinear iterationb
maxitnl=5;                 % maximal number of iterations

% SET APPROXIMATE FLUID/ICE PROPERTIES
rhof=1000;cpf=4217.;kf=0.561;
rhoi=917 ;cpi=2112.;ki=2.14;
freeze='no';
Err=0.2;

disp([' ']);
% GENERATE WEAKLY GRADED SPATIAL MESH 
zstart=10;zend=5000;nt=251;type='logarithmic';dir=1;
disp([' ']);disp([' ...set up ' type ' spatial mesh ']); 
[z,dz]=set_mesh(zstart,zend,nt,type,dir,1);


% READAPRIORI PROPERTIES FOR BOREHOLE
filename=strcat([NAME '_apriori.prp']);
disp([' ']);disp([' ...read geological units from ' filename]); 
[unit] = get_props(filename);units=length(unit);
unit_a=unit;



ncount=0;
for ib=1:nb
    bore=borehole{ib};
    disp([' ...reading data and model for borehole ' bore]); 
    
    % READ MODEL GEOMETRY FOR BOREHOLE
    filename=strcat(bore,'.mod');
    [model]=get_model(filename,z);
    well(ib).model=model;
    
    % READ DATA
    filename=strcat(bore,'.dat');
    [data]=get_data(filename,z,1);
    well(ib).data=data;id=data.id;ncount=ncount+length(id);
    
    
    well(ib).unit=unit;
end




% DIMENSIONS
% SET APRIORI COVARIANCE AND MODEL VECTOR

mpara=2*units+2*nb;ndata=ncount;
kl=[unit.k];hl=[unit.h];kAl=[unit.kA];kBl=[unit.kB];porl=[unit.p];
dkl=[unit.dk];dhl=[unit.dh];

p=[log(kl) log(hl)];perr=[log(1+dkl./kl) log(1+dhl./hl)];

for ib=1:nb
    model=well(ib).model;dQb=model.dqb;dTs=model.dgt;Qb=model.qb;Ts=model.gt;ip=model.ip;
    perr=[perr log(1+dQb./Qb) dTs];
    p=[p log(Qb) Ts];
    ipb(1,ib)=length(ip);
    ibp(2:length(ip)+1,ib)=ip';
end
pvar=perr.*perr;Cp=diag(pvar);CpI=diag(1./pvar);Wp=diag(1./perr);
p_apr=p;p_0=p;


% ITERATION
disp(['   ']);disp([' ...start iteration   ']);

p=p_0;
for iter=1:maxiter 
    
    r=[];
    for ib=1:nb
        
        % FORWARD MODEL
        
        kl=exp(p(1:units));hl=exp(p(units+1:2*units));qb=exp(p(2*units+(ib-1)*2+1));gt=p(2*units+(ib-1)*2+2);
        model=well(ib).model;ip=model.ip;
        Tcalc=heat1dns(kl,kAl,kBl,hl,porl,kf,ki,qb,gt,ip,dz,maxitnl,tolnl,freeze);
        
        % CALCULATE RESIDUAL   
        data=well(ib).data;
        Tobs=data.T;id=data.id;derr=data.Err;dcov=data.Cov;       
        W= diag(1./sqrt(dcov),0);
        resid=Tobs(id)'-Tcalc(id);rb=W*resid;
        rms=norm(rb)/sqrt(length(rb));r=[r; rb];
        disp([ 'rms for iteration ',num2str(iter), ...
                ' at borehole ', borehole{ib}, ' = ',num2str(rms)])
        
        result.Tobs=Tobs(id);result.Tcalc=Tcalc(id);
        result.z = z(id); result.r= rb;result.rms=rms;
        
        well(ib).result=result;
        
        
    end
    
    rms=norm(r)/sqrt(length(r));
    disp([ ' *** total rms for iteration ',num2str(iter), ...
            ' = ',num2str(rms)])
    if rms < tolrms, break; end
    
    disp([' ... calculate sensitivities   ']);
   
    S=zeros(ndata,mpara);    %l=zeros(ndata,2*units);Sb=zeros(ndata,2*nb);
    nc=0;mc=2*units;
    
    
    for ib=1:nb
        
        kl=exp(p(1:units));hl=exp(p(units+1:2*units));qb=exp(p(2*units+(ib-1)*2+1));gt=p(2*units+(ib-1)*2+2);
        model=well(ib).model;ip=model.ip;
        [Jk,Jh,JQb,JTs]=sensfds_pet1(kl,kAl,kBl,hl,porl,kf,ki,qb,gt,ip,dz, maxitnl,tolnl,dp,maxitdp,toldp,freeze,'no');
        
        
        data=well(ib).data;id=[data.id];nd=length(id);
        dcov=[data.Cov];Wd= diag(1./sqrt(dcov),0);
        % Sl0=Wd*[Jk(id,:) Jh(id,:) ]; 
        S(nc+1:nc+nd,1:2*units)          =   Wd*[Jk(id,:) Jh(id,:) ]; 
        % Sb0=Wd*[JQb(id)   JTs(id) ]; 
        S(nc+1:nc+nd,mc+1:mc+2)          =   Wd*[JQb(id)  JTs(id) ]; 
        
%         opts = struct('Format','psc2','Color','rgb','Resolution',600);
%         filename=strcat(borehole{ib},'-J-',num2str(j),'.eps');
%         %       Jplot=[S(nc+1:nc+nd,1:2*units) S(nc+1:nc+nd,mc+1:mc+2)];
%         figure;
%         subplot(1,4,1);imagesc(Wd*Jk(id,:));colorbar;title(strcat('\lambda/',borehole{ib},'/',num2str(j)));  
%         subplot(1,4,2);imagesc(Wd*Jh(id,:));colorbar;title(strcat('h/',borehole{ib},'/',num2str(j)));  
%         subplot(1,4,3);imagesc(Wd*JQb(id,:));colorbar;title(strcat('Q_b/',borehole{ib},'/',num2str(j)));  
%         subplot(1,4,4);imagesc(Wd*JTs(id,:));colorbar;title(strcat('T_s/',borehole{ib},'/',num2str(j)));  
%         exportfig(gcf,filename,opts);close(gcf) 
%         %       filename=strcat(borehole{ib},'-J-',num2str(iter));save filename iter, Jplot
       
        
        
        nc=nc+nd;mc=mc+2;
        
    end
    
    
    disp([' ... calculate parameter increment ']);
    
    S_Aug=[ S;Wp];resid_Aug=[ r;-Wp*(p'-p_apr')];
    
    [delta] = cglsACB(S_Aug,resid_Aug,kc,reorth,stoprule,tol,debug);        
    p = p + delta' ;
    p_iter(:,iter)=p';rms_iter(iter)=rms;
    % RESTORE PARAMETERS TO STRUCTS
    for ib=1:nb
        kl=exp(p(1:units));hl=exp(p(units+1:2*units));qb=exp(p(2*units+(ib-1)*2+1));gt=p(2*units+(ib-1)*2+2);
        unit=well(ib).unit; 
        units=length(unit);
        for i=1:units
            unit(i).k=kl(i); unit(i).h=hl(i);
        end
        well(ib).unit=unit;
        model=well(ib).model; 
        model.qb=exp(p(2*units+(ib-1)*2+1));
        model.gt=p(2*units+(ib-1)*2+2); well(ib).model=model;
    end
    
end

disp(['   ']);disp([' ...calculate aposteriori quantities ']);
% CALCULATE COVARIANCES A POSTERIORI BY GENERALIZED INVERSE
GT_Aug=inv(S_Aug'*S_Aug)*S';
% PARAMETER & DATA COVARIANCE MATRIX A POSTERIORI
% (FROM GENERALIZED INVERSE NOLET(99))
Cmm=GT_Aug *GT_Aug'; 
Cdd=GT_Aug'*GT_Aug;


% PLOTS 
disp(['   ']);disp([' ...plot results for ',NAME]);
opts = struct('Format','psc2','Color','rgb','Resolution',600);

% PLOT APOSTERIORI COVARIANCE MATRICES
filename=strcat('BIPETS_cpp_aposteriori.eps');
figure;imagesc(Cmm);colorbar;title('C_{pp}^{aposteriori} based on Generalized Inverse');  
exportfig(gcf,filename,opts);%close(gcf)
filename=strcat('BIPETS__cdd_aposteriori.eps');
figure;imagesc(Cdd);colorbar;title('C_{dd}^{aposteriori} based on Generalized Inverse');  
exportfig(gcf,filename,opts);%close(gcf) 

% PLOT APOSTERIORI MODELS
for ib=1:nb
    data=well(ib).data;
    Tobs=data.T;id=data.id;derr=data.Err;dcov=data.Cov;   
    model=well(ib).model;ip=model.ip;Qb=model.qb;Ts=model.gt;
    unit=well(ib).unit; 
    kl=[unit.k];hl=[unit.h];rl=[unit.r];cl=[unit.c];kAl=[unit.kA];kBl=[unit.kB]; porl=[unit.p];rcl=rl.*cl; 
    dkl=[unit.dk];dhl=[unit.dh];
    Tc=heat1dns(kl,kAl,kBl,hl,porl,kf,ki,qb,gt,ip,dz,maxitnl,tolnl,freeze);
    filename=strcat(borehole{ib},'.prp');
    [unit] = get_props(filename);units=length(unit);
    klt=[unit.k];hlt=[unit.h];rlt=[unit.r];clt=[unit.c];kAlt=[unit.kA];kBlt=[unit.kB]; porlt=[unit.p];rclt=rl.*cl; 
    dklt=[unit.dk];dhlt=[unit.dh];
    
    
    filename=strcat(borehole{ib},'_aposteriori.eps');

    figure
    subplot(1,4,1);
    kp=[kl(1) kl(ip)];dkp=[dkl(1) dkl(ip)];
    kpt=[klt(1) klt(ip)];dkpt=[dklt(1) dklt(ip)];
    stairs(kp',z', '-r');hold on;grid on;
    stairs(kpt',z', '-b');
    stairs(kpt-dkpt,z, ':g');hold on
    stairs(kpt+dkpt,z, ':g');hold on
    xlim([1  5]);ylim([0 3000]);
    xlabel('\lambda (W m^{-1} K^{-1})');ylabel('z (m)');
    set(gca,'ydir','reverse')
    title(['T_s = ' num2str(Ts) ' mW/m^2'])
    
    subplot(1,4,2);
    hp=[hl(1) hl(ip)]*1.e6;dhp=[dhl(1) dhl(ip)]*1.e6;
    hpt=[hlt(1) hlt(ip)]*1.e6;dhpt=[dhlt(1) dhlt(ip)]*1.e6;
    stairs(hp',z', '-r');hold on;grid on;    
    stairs(hpt',z', '-b');
    stairs(hpt-dhpt,z, ':g');hold on
    stairs(hpt+dhpt,z, ':g');hold on
    xlim([0  2]);ylim([0 3000]);
    xlabel('h (\mu W m^{-3})');ylabel('z (m)');
    set(gca,'ydir','reverse','yticklabel',{})
    
    
    subplot(1,4,3);
    plot([Tc(id) Tobs(id)'],z(id));hold on;grid on;
    xlim([0 120]);ylim([0 3000]);
    xlabel('T (C)');ylabel('z (m)');
    set(gca,'ydir','reverse','yticklabel',{})
    % legend('calculated',' with random error')
    title(['Q_b = ' num2str(Qb*1000) ' mW/m^2'])
    
    subplot(1,4,4);
    plot(Tobs(id)'-Tc(id),z(id), '-r')
    xlim([-2 2]);ylim([0 3000]);;hold on;grid on;
    xlabel('\Delta T (C)');ylabel('z (m)');
    set(gca,'ydir','reverse','yticklabel',{})
    
       
    exportfig(gcf,filename,opts);%close(gcf)
end


figure;
filename=strcat('Covergence_rms.eps');
Conv=[rms_iter];grid on;
plot(Conv);xlabel('Iteration');ylabel('ln(total rms)');
exportfig(gcf,filename,opts);%close(gcf)

figure;
filename=strcat('Covergence_k.eps');
Conv=[exp(p_iter(1:units,:))'];grid on;
plot(Conv);xlabel('Iteration');ylabel('\lambda (W m^{-1} K^{-1})');
legend('1','2','3','4');
exportfig(gcf,filename,opts);%close(gcf)

figure;
filename=strcat('Covergence_h.eps');
Conv=[exp(p_iter(units+1:2*units,:))'];grid on;
plot(Conv);xlabel('Iteration');ylabel('h (\mu W m^{-3})');
legend('1','2','3','4');
exportfig(gcf,filename,opts);%close(gcf)

figure;
filename=strcat('Covergence_Qb.eps');
Conv=[exp(p_iter(2*units+1:2:2*units+2*nb,:))'*1000];grid on;
plot(Conv);xlabel('Iteration');ylabel('Q_b (m W m^{-2})')
legend('test_1','test_2','test_3');
exportfig(gcf,filename,opts);%close(gcf)

figure;
filename=strcat('Covergence_Ts.eps');
Conv=[p_iter(2*units+2:2:2*units+2*nb,:)'];grid on;
plot(Conv);xlabel('Iteration');ylabel('T_s (C)')
legend('test_1','test_2','test_3');
exportfig(gcf,filename,opts);%close(gcf)


% SAVE RESULTS
disp(['   ']);disp([' ...save results  ']);
filename=strcat('BIPET_Results.mat');
save(filename,'p_iter','rms_iter' ) 