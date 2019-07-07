% TIPAL_MULTI performs paleoclimatic inversion in parameter space for
% multiple boreholes

% important index arrays for indirekt addresses:
% id        points to nodes defining data
% ip        associates parameter values of diffusivity and production to cells
% it        accociates paleotemperatures to temporal grid cells
%
% structures used:
% unit      information associated with geological units (layers)
% model     information defining model for given borehole
% obs       information cncerning the observations in borehole
%
%
% V. R. April 2, 2003
%close all;clear all
close all;
keep kc stoprule tolcgls reorth tolrms maxiter tau0 tau1 tau2 borehole freeze NAME m_apr POR debug;
% keep m_apr;?

disp(['Nonlinear Tikhonov Inversion for Paleoclimate from Multiple Boreholes'])
disp(['- Permafrost -'])

disp(['V. R., ',date] )

disp([' ']);
% READ MODEL PARAMETER & DATA
nb=length(borehole);

% GENERATE SPATIAL MESH
zstart=10;zend=5000;nz=251;type='logarithmic';dir=1;
disp([' ']);disp([' ...set up ' type ' spatial mesh ']);
[z,dz]=set_mesh(zstart,zend,nz,type,dir,1);

% GENERATE TEMPORAL MESH
year2sec=31557600;
tstart=100000*year2sec;tend=10*year2sec;nt=257;type='logarithmic';dir=-1;
disp([' ']);disp([' ...set up ' type ' temporal mesh ']);
[t,dt]=set_mesh(tstart,tend,nt,type,dir,1);

% DEFINE LOGARITHMIC INVERSION GRID
disp(['   ']);disp([' ...set up parametrization for paleoclimate inversion  ']);
nsteps=24; base=0.;
[pt,it]=set_paleo_grid(t,base,tstart,tend,nsteps);
np= length(pt);

% INVERSION SETUP
disp(['   ']);disp([' ...set up inversion control parameters  ']);
% PARAMETER FOR FWD CALCULATIONS
theta=0.5*ones(1,nt);
theta(1:10)=1.;
% NONLINEAR ITERATION PARAMETERS
maxitnl=2;
tolnl=0.0001;
% SENSITVITY ITERATION PARAMETERS
dp=0.1;
maxitdp=1;
toldp=0.01;
% PARAMETERS FOR INVERSION
% tolrms= 0.5;           % tolerance for rms
% maxiter=16;           % maximal number of iterations
% tau0 = 0.001;		% intial/center regularization parameters
% tau1 = 5;		 % Gradient
% tau2 = 1.0;	        % Laplacian
% kc = 12;                 % number of CGLS iterations



% SETUP REGULARIZATION MATRIX

L0 = reg1d(np,'l0');
L1 = reg1d(np,'l1');
L2 = reg1d(np,'l2');

disp([' ']);
for ib=1:nb
    bore=borehole{ib};
    disp([' ...reading data and model for borehole ' bore]);

%     % READ MEASUREMENTS AND DATA FOR BOREHOLE
%     file=strcat(['TS_' bore '.mat']);
%     load(file);clear ip
%     % store into structure model
%     nz=length(z);
%     ip(1:nz-1)=[1:nz-1];nip=length(ip);nones=ones(nip,1)';
%     zd=Z;Pd=P;Ld=L;
%     Pm=sum(Pd)/length(Pd);
%     Lm=sum(Ld)/length(Ld);
% 
%     p=interp1(zd,Pd,z,'linear');
%     k=interp1(zd,Ld,z,'linear');
%     indexnan=find(isnan(p));
%     p(indexnan)=Pm;
%     indexnan=find(isnan(k));
%     k(indexnan)=Lm;
% 
%     model.ip=ip;model.z=z;model.gt=gts;model.dgt=0.;model.qb=qbs;model.dqb=0.;
%     model.k=k;
%     model.por=p;
%     model.rho=2650*nones; model.cp =850*nones;
%     model.h=1.e-6.*nones;model.kA=0.7*nones;model.kB=770*nones;
%     % store into structure data
%     index=find(zd>=top & zd<=bot);
%     zd=zd(index);Td=T(index);
%     Tm=interp1(zd,Td,z,'linear');
%     id=find(~ isnan(Tm));nd=length(id);
%     data.T=Tm;data.id=id;data.nd=nd;data.z=zd;
%     err=0.3;data.Err=err;data.Cov=err^2*ones(1,nd);
% 
%     well(ib).data=data;
%     well(ib).model=model;
% 
%     file=strcat(['TS_' bore '_out.mat']);
%     save(file,'data','model');
    
    
%     ylimits=[0 1500];
%     tlimits=[0 60];
%     climits=[0 6];
%     hlimits=[10 70];
%     figure;
%     subplot(1,2,1);
%     plot(Tm(id),z(id),'-g','LineWidth',2);grid on;hold on;
%     ylim(ylimits);xlim(tlimits);
%     ylabel('z (m)');xlabel('T (?C)');
%     set(gca,'YDir','reverse')
%     subplot(1,2,2);
%     plot(k,z,'ob');hold on; grid on;
%     ylim(ylimits);xlim(climits);
%     ylabel('');xlabel('\lambda (J m^{-1}K^{-1})');
%     set(gca,'YDir','reverse')
%     suptitle(['Borehole:' bore])
%     file=strcat(['TS-' bore '.ps']);
%     saveas(gcf,file,'epsc2');
%     close(gcf)
    file=strcat(['TS_' bore '_out.mat']);
    load(file)
    well(ib).data=data;
    well(ib).model=model;


end
disp([' ']);

% SETUP A-PRIORI MODEL
disp([' ']); disp([' ...setup apriori model ' ]);
% [m_apr]=set_gst_prior('constant',struct('n',np,'gt', 10) )
%gtm=mean(gtb);qbm=mean(qbb);

m_apr=zeros(1,np);

m_0 = m_apr; m = m_0;

% ITERATION
disp(['   ']);disp([' ...start iteration   ']);
for iter=1:maxiter

    Ts=m(it);
    r_tot=[];
    for ib=1:nb

        % FORWARD MODEL
        model=well(ib).model;
        kl=[model.k];hl=[model.h];kAl=[model.kA];kBl=[model.kB];porl=[model.por];
        cpml=[model.cp];rhoml=[model.rho];qb=model.qb;gt=model.gt;ip=model.ip;
        POM=m(1);
        Ti(:,ib)=heat1dns(kl, kAl, kBl,hl,porl,qb,POM+gt,ip,dz,maxitnl,tolnl,freeze);
        Tini=Ti(:,ib);         Ts=m(it)+gt;
        [Tcalc,zout,tout,N,k_eff,rc_eff,ipor,lheat,rci]=heat1dnt(kl,kAl,kBl,hl,rhoml,cpml,porl,qb,...
            ip,dz,dt,Tini,Ts,theta,maxitnl,tolnl,freeze);
        % CALCULATE RESIDUAL
        data=well(ib).data;
        Tobs=data.T;id=data.id;err=data.Err;
        resid=Tobs(id)'-Tcalc(id,nt);res=resid/err;
        rms=norm(res)/length(res);
        r_tot=[r_tot; res];
        disp([ 'rms for iteration ',num2str(iter), ...
            ' at borehole ', borehole{ib}, ' = ',num2str(rms)])

        result.Tobs=Tobs(id);result.Tcalc=Tcalc(id,nt);
        result.z = z(id); result.resid= resid;result.errT=err;
        result.rms=rms;result.m=m;
        well(ib).result=result;
        file=strcat(['TS_' bore '_out.mat']);
        save(file,'data','model','result');


    end

    rms_tot=norm(r_tot)/length(r_tot);
    disp([ ' *** total rms for iteration ',num2str(iter), ...
        ' = ',num2str(rms_tot)])
    if rms_tot < tolrms, break; end

    disp([' ... calculate sensitivities   ']);
    S=[];
    for ib=1:nb
        model=well(ib).model;
        kl=[model.k];hl=[model.h];kAl=[model.kA];kBl=[model.kB];porl=[model.por];cpml=[model.cp];rhoml=[model.rho];
        qb=model.qb;gt=model.gt;ip=model.ip;Tini=Ti(:,ib);mib=m+gt;


        data=well(ib).data;
        Tobs=data.T;id=data.id;err=data.Err;cov=data.Cov;
        W=diag(1./sqrt(cov),0);



        [J]=sensfdt_pal(kl,kAl,kBl,hl,rhoml,cpml,porl,qb,...
            ip,dz,dt,Tini,mib,it,theta,...
            maxitnl,tolnl,dp,maxitdp,toldp,freeze);

        S=[S;W*J(id,:)];

    end

    disp([' ... calculate parameter increment ']);
    S_Aug=[              S;...
        sqrt(tau2)*L2; ...
        sqrt(tau1)*L1; ...
        sqrt(tau0)*L0 ];
    resid_Aug=[                     r_tot;...
        -sqrt(tau2)*L2*(m'-m_apr');...
        -sqrt(tau1)*L1*(m'-m_apr');...
        -sqrt(tau0)*L0*(m'-m_apr')];
    %     [delta_m] = cglsACB(S_Aug,resid_Aug,reorth,kc,tolcgls,stoprule,1);
    %     m = m + delta_m' ;
    [delta_m] = cgls(S_Aug,resid_Aug,kc,reorth,tolcgls,debug);
    m = m + delta_m' ;

    m_iter(iter,:)=m(1,:);rms_iter(iter)=rms;t_iter(iter,:)=Tcalc(:,nt)';
end

disp(['   ']);disp([' ...calculate aposteriori quantities ']);
% CALCULATE COVARIANCES A POSTERIORI
S_Aug =[ S;sqrt(tau2)*L2;sqrt(tau1)*L1 ;sqrt(tau0)*L0];
S_Augt=[S',sqrt(tau2)*L2',sqrt(tau1)*L1',sqrt(tau0)*L0'];
% GENERALIZED INVERSE
GT_Aug=inv(S_Augt*S_Aug)*S';
% PARAMETER & DATA COVARIANCE MATRIX A POSTERIORI
% (FROM GENERALIZED INVERSE NOLET(99))
Cmm=GT_Aug*GT_Aug'; Cdd=GT_Aug'*GT_Aug;




% SAVE DATA
disp(['   ']);disp([' ...save results  ']);
ty=t/year2sec;errall=sqrt(diag(Cmm));err=errall(it);mod=m;
filename=strcat(NAME,'_results.mat');
save(filename,'ty','mod','it','err','tau0','tau1','tau2','rms','well','rci','J','S')

% PLOTS
disp(['   ']);disp([' ...plot results for ',NAME]);

opts = struct('Format','psc2','Color','rgb','Resolution',600);
lc={'b','r','g','c','y','m','k','b','r','g','c','y','m','k'};
ls={'-','-','-','-','-','-','-' '--','--','--','--','--','--','--'};
bc=borehole;

figure;
filename=strcat(NAME,'_climate_posteriori.ps');

[X,Y]=stairs(-ty,mod(it));
plot(X,Y,'LineWidth',2,'Color','b','LineStyle','-');hold on;
[X,Y]=stairs(-ty,mod(it)+2*err');
plot(X,Y,'LineWidth',1,'Color','r','LineStyle','--');hold on;
[X,Y]=stairs(-ty,mod(it)-2*err');
plot(X,Y,'LineWidth',1,'Color','r','LineStyle','--');hold on;

set(gca,'XScale','log','XDir','reverse')
xlabel('time b. p.(a)','FontSize',14);ylabel('\Delta T','FontSize',14);
ylim([-10 10]);
%text(1000,-3.00,['L_{1}  = ',num2str(tau1)],'FontSize',14)
%text(1000,-4.00,['rms  = ',num2str(rms),' K '],'FontSize',14)
title([NAME ': final model'],'FontSize',14)
grid on;
exportfig(gcf,filename,opts)
%close(gcf)


% plot Covariances
filename=strcat(NAME,'_cpp_aposteriori.eps');
figure;imagesc(Cmm);colorbar;title('C_{pp}^{aposteriori} based on Generalized Inverse');
exportfig(gcf,filename,opts);close(gcf)
filename=strcat(NAME,'_cdd_aposteriori.eps');
figure;imagesc(Cdd);colorbar;title('C_{dd}^{aposteriori} based on Generalized Inverse');
exportfig(gcf,filename,opts);close(gcf)
%close(gcf)


% plot single fits
for ib=1:length(well);
    bore=borehole{ib};

    result=well(ib).result;
    Tobs=[result.Tobs]';
    Tcalc=[result.Tcalc];
    z=[result.z];resid=[result.resid];
    col=lc{ib}; namb=cellstr(bc{ib});
    figure;
    subplot(1,2,1)
    plot(resid,z, 'LineWidth',2,'Color',col); hold on;
    set(gca,'YDir','reverse');grid on;
    xlabel(' Residual \Delta T (^\circ C)'):ylabel ('z (m)');
    grid on;legend(namb,'Location','northeast');

    subplot(1,2,2)
    plot(Tobs,z, 'LineWidth',2,'Color',col,'LineStyle','-'); hold on;
    plot(Tcalc,z, 'LineWidth',2,'Color',col,'LineStyle','--');
    set(gca,'YDir','reverse');grid on;
    xlabel(' Residual \Delta T (^\circ C)'):ylabel ('z (m)');
    suptitle([NAME ': Residuals'])
    grid on;legend('Observed','Calculated','Location','northeast');
    filename=strcat([NAME '-' bore '_aposteriori_log.ps']);
    exportfig(gcf,filename,opts)
    %close(gcf)

end

% plot total fit
figure;
for ib=1:length(well);
    result=well(ib).result;
    Tobs=[result.Tobs]';
    Tcalc=[result.Tcalc];
    z=[result.z];resid=[result.resid];
    col=lc{ib};
    plot(resid,z, 'LineWidth',2,'Color',col); hold on;

end
set(gca,'YDir','reverse');grid on;
xlabel(' Residual \Delta T (^\circ C)'):ylabel ('z (m)');
%text(-0.4,2500,['L_{2}  = ',num2str(tau1)],'FontSize',14)
%text(-0.4,3000,['rms  = ',num2str(rms),' K '],'FontSize',14)
title([NAME ': Residuals'],'FontSize',14)
grid on;legend(cellstr(bc),'Location','northeast');
filename=strcat(NAME,'_aposteriori_log.ps');
exportfig(gcf,filename,opts)
%close(gcf)
