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
% V. R. may 16, 2004

close all;clear all



disp(['Nonlinear Tikhonov Inversion for Paleoclimate'])
disp(['Minimum (gradient) support'])
disp(['V. R., ',date] )
NAME='FOPAL-MGS';

disp([' ']);
borehole={'borehole0'}
nb=length(borehole);


% GENERATE SPATIAL MESH 
zstart=0;zend=3000;nz=301;methz='linear';dirz=1;
disp([' ']);disp([' ...set up ' methz ' spatial mesh ']); 
[z,dz]=set_mesh(zstart,zend,nz,methz,dirz,0);

% GENERATE TEMPORAL MESH
year2sec=31557600;
tstart=100000*year2sec;tend=100*year2sec;nt=301;metht='logarithmic';dirt=-1;
disp([' ']);disp([' ...set up ' metht ' temporal mesh ']); 
[t,dt]=set_mesh(tstart,tend,nt,metht,dirt,0);

% DEFINE LOGARITHMIC INVERSION GRID
disp(['   ']);disp([' ...set up parametrization for paleoclimate inversion  ']);
nsteps=36; base=0.;
[pt,it]=set_paleo_grid(t,base,tstart,tend,nsteps);


% PARAMETER FOR FWD CALCULATIONS

freeze='yes';           % include permafrost
theta=1.*ones(1,nt);
maxitnl=2;
tolnl=0.001;



% INVERSION SETUP 
disp([' ']); disp([' ...setup inversion parameters ' ]); 

% PARAMETERS FOR INVERSION 
tolrms= 0.1;            % tolerance for rms
kc = 32;                % number of CGLS iterations
stoprule=2;             % use 1/ norm reduction 2/ACB-criterium
tolcgls=0.001;          % 1 for  ACB-criterium
reorth=1;               % reorthogonalization
tau0 = 0.1;		% weight for identity matrix  
tau1 = 3;               % weight for Gradient
tauf = 0.1;             % weight for focusing algorithm
epsl  = 1.e-3;          % eps for focusing algorithm

% SENSITVITY CALCULATION 
dp=0.1;
maxitdp=1;
toldp=0.001;

debug=1;

regmeth=[1 1 1 3 1 3 1];

maxiterG   = length(find(regmeth==1));                % number of CGLS iterations
maxiterMS  = length(find(regmeth==2));
maxiterMGS = length(find(regmeth==3));

maxiter   = length(regmeth);
disp([' ']); 
disp(['    number of gradient smoothing iterations       : ' num2str(maxiterG) ]); 
disp(['    number of minimum support iterations          : ' num2str(maxiterMS) ]); 
disp(['    number of minimum gradient support iterations : ' num2str(maxiterMGS) ]); 

disp(['    tolat number of  iterations   : ' num2str(maxiter) ]); 

% SETUP REGULARIZATION MATRICES 
np= length(pt);
[L0 W0] = l_1d(np,'l0'); 
[L1 W1] = l_1d(np,'l1');
% [L2 W2] = l_1d(np,'l2');

tau_ratio=sqrt(tau0/tau1);
L1(np,np)=1;LS=L1'+tau_ratio*L0;


% SETUP A-PRIORI MODEL
disp([' ']); disp([' ...setup apriori model ' ]); 
% load SYNTH_Apriori_Zero.mat; m_apr =mod;
% [m_apr]=set_gst_prior('constant',struct('n',np,'gt', 10) )
% gtm=mean(gtb);qbm=mean(qbb);
m_apr=zeros(1,np); 
m_0 = m_apr; m = m_0;



% READ AND ORGANIZE DATA

disp([' ']);
for ib=1:nb
    bore=borehole{ib};
    disp([' ...reading data and model for borehole ' bore]); 
    file=strcat([ bore '.mat']);load(file);
    % store into structure data 
    Tm=interp1(zp,T,z,'linear');id=find(z > min(zp) & z < max(zp));
    nd=length(id); data.T=Tm(id);data.id=id;data.nd=nd;data.z=z(id);
    err=0.5; data.Err=err;data.Cov=err^2*ones(1,nd);
    well(ib).data=data;
    
    % store into structure model 
    lamc=set_cell(zp,abs(lam),z,'linear'); 
    porc=set_cell(zp,abs(por),z,'linear');
    ip(1:nz-1)=[1:nz-1];repv=ones(1,nz-1); gts=Tm(1);
    model.ip=ip;model.z=z;model.gt=gts;model.dgt=0.;model.qb=qbs;model.dqb=0.;
    model.k=lamc;model.por=porc; model.rho=2321*repv; model.cp =840*repv;
    model.h=0.*repv;model.kA=0.7*repv;model.kB=770*repv;
    well(ib).model=model;
    
    % save data and model
    filename=strcat([bore  '_IN.mat' ]);
    disp(['    borehole <' bore '> read and saved  to: '  filename]); 
    save(filename,'model','data') 
    
end
disp([' ']); 


% SETUP INITIAL VALUES FOR TEMPERATURE AT TSTART
disp([' ']); disp([' ...setup initial values ' ]); 

for ib=1:nb
    
    kl=[model.k];hl=[model.h];kAl=[model.kA];kBl=[model.kB];porl=[model.por];
    qb=model.qb;gt=model.gt;
    
    Ti(:,ib)=heat1dns(kl, kAl, kBl,hl,porl,qb,gt,ip,dz,maxitnl,tolnl,freeze); 
end
disp([' ']);


% ITERATION
disp(['   ']);disp([' ...start iteration   ']);

for iter=1:maxiter 
    disp([' ']);
    disp([' ']);
    disp([ '    iteration ',num2str(iter-1)])

    time0=cputime;
    
    Ts=m(it);
    r_tot=[];
    for ib=1:nb
        
        % FORWARD MODEL
        unit=well(ib).model;
        kl=[model.k];hl=[model.h];kAl=[model.kA];kBl=[model.kB];porl=[model.por];cpml=[model.cp];rhoml=[model.rho];
        qb=model.qb;gt=model.gt;ip=model.ip;
        Tini=Ti(:,ib);Ts=m(it)+gt;
        
        [Tcalc,Td,N]=heat1dnt(kl,kAl,kBl,hl,rhoml,cpml,porl,qb,...
            ip,dz,dt,Tini,Ts,theta,maxitnl,tolnl,freeze);
        
%         if debug==1,
%             warn=find(N>0.5);numw=length(warn);
%             disp(['    neumann criterium failed at ' num2str(numw) ' time steps !']); 
%         end
        
        % CALCULATE RESIDUAL   
        data=well(ib).data;
        Tobs=data.T;id=data.id;err=data.Err;cov=data.Cov;Wd=diag(1./sqrt(cov),0);        
        
        resid=Tobs'-Tcalc(id,nt);
        rms=norm(Wd*resid)/sqrt(length(resid));
        
        disp([ '       rms for iteration ',num2str(iter-1), ...
                  ' at borehole <', borehole{ib}, '> = ',num2str(rms)])
        r_tot=[r_tot; Wd*resid];
        
        result.Tobs=Tobs;result.Tcalc=Tcalc(id,nt);
        result.z = z(id); result.resid= resid;result.errT=err;
        result.rms=rms;
        well(ib).result=result;
                
        
    end
    % 
  
    n_tot=length(r_tot);
    norm_data=norm(r_tot);
    rms_tot=norm_data/sqrt(n_tot);
    disp([ '    total rms   = ',num2str(rms_tot)])
    phi_tot=0;norm_para=0;
    if iter > 1,
        norm_para=norm(Wp*(m-m_apr+gt)');
        phi_tot=norm_data+norm_para;
        disp([ '    norm objf   = ',num2str(phi_tot)])
        disp([' ']); 
    end
    mod_iter(iter,:)=m(1,:);   
    rms_iter(iter)=rms_tot;   
    phi_iter(iter)=phi_tot;   
    nrd_iter(iter)=norm_data;
    nrp_iter(iter)=norm_para;

    
    
    if rms_tot < tolrms, break; end
    
    %    if iter > 1,
    %       if rms_tot/rms_tot_old > 1 break; end
    %    end
    %    rms_tot_old=rms_tot;
    
    
    disp([' ...calculate sensitivities   ']);
    S=[];
    for ib=1:nb
        model=well(ib).model;
        kl=[model.k];hl=[model.h];kAl=[model.kA];kBl=[model.kB];porl=[model.por];cpml=[model.cp];rhoml=[model.rho];
        qb=model.qb;gt=model.gt;ip=model.ip;Tini=Ti(:,ib);Ts=m(it)+gt;
        
        
        data=well(ib).data;
        Tobs=data.T;id=data.id;err=data.Err;cov=data.Cov; 
        Wd=diag(1./sqrt(cov),0);        
        
        
        
        [J]=sensfdt_pal(kl,kAl,kBl,hl,rhoml,cpml,porl,qb,...
            ip,dz,dt,Tini,m,it,theta,...
            maxitnl,tolnl,dp,maxitdp,toldp,freeze);
        
        S=[S;Wd*J(id,:)]; 
        
  
    switch regmeth(iter) 
        case{1}
	    disp([ '    gradient smoothing regularization'])
            Wp= sqrt(tau1)*LS; 
        case{2}
	    disp([ '    minimum support regularization'])
            Wp=tauf*wp_ms(m,m_apr+gt,epsl,1);
        case{3}
	    disp([ '    minimum gradient support regularization'])
            Wp=tauf*wp_ms(m,m_apr+gt,epsl,2);
    otherwise
            Wp= sqrt(tau1)*LS; 
    end
    
 
   
    
    disp([' ...calculate parameter increment ']);
    S_Aug=    [              S;...
                           Wp];
    r_Aug=[              r_tot;...
            -Wp*(m'-m_apr'+gt)];
    
    
    [delta_m] = cglsACB(S_Aug,r_Aug,reorth,kc,tolcgls,stoprule,0);        
    m = m + delta_m' ;
    
    disp([ '    cpu time for iteration ',num2str(cputime-time0),' s '])  
    
end

disp(['   ']);disp([' ...calculate aposteriori quantities ']);
% CALCULATE COVARIANCES A POSTERIORI
S_Aug =[ S;Wp];S_Augt=[S',Wp'];
% GENERALIZED INVERSE
GT_Aug=inv(S_Augt*S_Aug)*S';
% PARAMETER & DATA COVARIANCE MATRIX A POSTERIORI
% (FROM GENERALIZED INVERSE NOLET(99))
Cmm=GT_Aug*GT_Aug'; Cdd=GT_Aug'*GT_Aug;




% SAVE DATA
disp(['   ']);disp([' ...save results  ']);
ty=t/year2sec;errall=sqrt(diag(Cmm));err=errall(it);mod=m;rms=rms_tot;
normd=norm_data;normp=norm_para;
filename=strcat(NAME,'_results.mat');
save(filename,'ty','mod','it','err','normd','normp','tau0','tau1','epsl','rms','well','J','S','Wp') 




% PLOTS
disp(['   ']);disp([' ...plot results for ',NAME]);

opts = struct('Format','psc2','Color','rgb','Resolution',600);
lc={'b','r','g','c','y','m','k'};
ls={'-','-','-','-','-','-','-'};
% bc={'2271' ,'2731', '2908' ,'3200' ,'3209' ,'3356' ,'3359'};
bc={'2271' ,'2731', '2908' ,'3200' ,'3356' ,'3359'};

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
ylim([-12 2]);xlim([5 2e5]);
% text(1000,-3.00,['L_{1}  = ',num2str(tau1)],'FontSize',14)
text(1000,-3.50,['rms  = ',num2str(rms),' K '],'FontSize',14)
title([NAME ': final model'],'FontSize',14)
grid on;
exportfig(gcf,filename,opts)
%close(gcf)    

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

% text(-0.4,1700,[NAME],'FontSize',10)
% text(-0.4,1650,['L_{0}  = ',num2str(tau0)],'FontSize',10)
% text(-0.4,1625,['L_{1}  = ',num2str(tau1)],'FontSize',10)
% text(-0.4,1600,['L_{2}  = ',num2str(tau2)],'FontSize',10)
text(-0.4,1500,['rms  = ',num2str(rms),' K '],'FontSize',14)
title([NAME ': Residuals'],'FontSize',14)
grid on;
%legend(cellstr(bc),4);
filename=strcat(NAME,'_aposteriori_log.ps');
exportfig(gcf,filename,opts)
%close(gcf)
filename=strcat(NAME,'_cpp_aposteriori.ps');
figure;imagesc(Cmm);colorbar;title('C_{pp}^{aposteriori} based on Generalized Inverse');  
exportfig(gcf,filename,opts);close(gcf)
filename=strcat(NAME,'_cdd_aposteriori.ps');
figure;imagesc(Cdd);colorbar;title('C_{dd}^{aposteriori} based on Generalized Inverse');  
exportfig(gcf,filename,opts);close(gcf)
%close(gcf)
