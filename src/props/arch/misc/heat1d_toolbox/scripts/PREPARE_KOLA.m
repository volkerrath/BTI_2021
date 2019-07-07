clear all; close all;
%Number ={'2271','2731','2908','3200','3209','3356','3359'};
Number ={'2400','2915','3200','3209','3356','3359','3396','2908','2271','2731'};
Long   ={0,0,0,0,0,0,0,0,0,0};
Lati   ={0,0,0,0,0,0,0,0,0,0};
Elev   ={     9999,     9999,     9999,    9999,    9999,   9999,  9999,    9999,   9999,  9999};
Top    ={       50,       50,       50,      50,      50,    50,  50,      50,      50,     50};
Bottom ={     9999,     9999,     9999,    9999,    9999,   9999,  9999,    9999,   9999,  9999};

%SurfT  ={      2.5,      2.5,      2.5,     2.5,     2.5,    2.5 };
SurfT  ={      1.5,       1.0,      1.5,       1.5,        1.0,       2.5,       2.5,      2.5,        2.5,       2.5,  };
HFD    ={   40.e-3,   41.e-3,   37.e-3,  38.e-3,  38.e-3,   37.e-3, 32.e-3, 38.e-3,   41.e-3, 40.e-3, };
Lam    ={      3.41,      3.44,      3.08,    3.14,     3.18,    3.2,    2.82,     3.19,    3.16,    3.34  };
Poro    ={     .01,      .01,      .01,     .01,     .01,    .01,    0.1,     .01,    .01,    0.1 };
nSite=10;


site=2;from=site;to=site;

%borehole={'2908','3200','3209','3356'};

opts = struct('Format','psc2','Color','rgb','Resolution',600);
ylimits=[0 1500];
tlimits=[0 60];
climits=[0 3];
hlimits=[10 70];
plimits=[0 1];
rlimits=[-1 1];

for iSite=from:to


    num=Number{iSite};
    top=Top{iSite};
    bot=Bottom{iSite};
    gts=SurfT{iSite};
    qbs=HFD{iSite};
    lat=Lati{iSite};
    lon=Long{iSite};
    elv=Elev{iSite};
    lam=Lam{iSite};
    lamb=Lam{iSite};
    por=Poro{iSite};
    hprod=1e-6;
    %   gr=Grad{iSite};
    %     disp([' ...reading data and model for borehole ' num]);
    disp([' ...Borehole ' num ' :  LAM=' num2str(lam) '  ' num2str(lamb) ...
        '   HDF=' num2str(qbs) '   POR=' num2str(por)]);

    
    file=strcat(['T' num '.xls']);    Tin=xlsread(file);
    file=strcat(['L' num '.xls']);    Lin=xlsread(file);

    ZT=-Tin(:,1);
    ZL=-Lin(:,1);
    
    T=Tin(:,2);
    L=Lin(:,2);
    

    index=find(ZT>=top);ZT=ZT(index,:);T=T(index,:);
    index=find(ZT<=bot);ZT=ZT(index,:);T=T(index,:);
% ----------
    find(diff(ZL)<=0);
    Lmov=movav(ZL,L,5,50); 
    L=Lmov(:,2);
    ZL=-Lmov(:,1);
    indexnan=find(isnan(L));
    L(indexnan)=lam;
    index=find(ZT<max(ZL)&ZT>min(ZL)); % index leer?
    Ln(index)=interp1(ZL,L,ZT(index));
    
    Lm=sum(Ln(index))/length(index);
    Ln(index)=interp1(ZL,L,ZT(index),'linear',Lm);
% ----------
    L=ones(size(T))*lam;
    %L(index)=Ln(index);
    Z=ZT;
    

    P=ones(size(L))*por;
%     figure;
%     subplot(1,3,1);
%     plot(T,Z,':b','LineWidth',1);hold on;
%     plot(T(ok),Z(ok),'-b','LineWidth',2);grid on;hold on;
%     ylim(ylimits);xlim(tlimits);
%     ylabel('z (m)');xlabel('T (C)');
%     set(gca,'YDir','reverse')
%     subplot(1,3,2);
%     plot(P,Z,'-b','LineWidth',2);hold on;grid on
%     ylim(ylimits);xlim(plimits);
%     ylabel('z (m)');xlabel('\phi (-)');
%     set(gca,'YDir','reverse')
%     subplot(1,3,3);
%     plot(Leff,Z,'-b','LineWidth',2);hold on;grid on
%     plot(L,Z,'-g','LineWidth',2);hold on;grid on
%     plot(kfT(T),Z,'-m','LineWidth',2);hold on;grid on
%     ylim(ylimits);xlim(climits);
%     ylabel('z (m)');xlabel('\lambda (J m^{-1}K^{-1})');
%     set(gca,'YDir','reverse');legend('e','m','f','Location','southeast');
%     file=strcat(['TS_' num '_in.ps']);
%     exportfig(gcf,file,opts);close(gcf)

    file=strcat(['Kola-',num '.mat']);
    %    save(file,'T','L','S','gts','qbs','lat','lon','elv','lam','por','lamb','gr')
    save(file,'Z','T','L','P','gts','qbs','lat','lon','elv','lam',...
        'por','lamb','hprod','top','bot');
end
disp([' ']);

% GENERATE SPATIAL MESH
zstart=10;zend=5000;nz=251;type='logarithmic';dir=1;
%disp([' ']);disp([' ...set up ' type ' spatial mesh ']);
[z,dz]=set_mesh(zstart,zend,nz,type,dir,0);
dz=dz';
% GENERATE TEMPORAL MESH
year2sec=31557600;
tstart=100000*year2sec;tend=10*year2sec;nt=257;type='logarithmic';dir=-1;
%disp([' ']);disp([' ...set up ' type ' temporal mesh ']);
[t,dt]=set_mesh(tstart,tend,nt,type,dir,0);
% DEFINE LOGARITHMIC INVERSION GRID
%disp(['   ']);disp([' ...set up parametrization for paleoclimate inversion  ']);
%nsteps=24; base=0.;
%[pt,it]=set_paleo_grid(t,base,tstart,tend,nsteps);

% PARAMETER FOR FWD CALCULATIONS
theta=1.*ones(1,nt);
theta(1:10)=1.;
maxitnl=2;
tolnl=0.0001;
freeze='no';

for iSite=from:to

%    disp([' ... modelling Borehole ' num ]);
    %        % READ MEASUREMENTS AND DATA FOR BOREHOLE
    file=strcat(['Kola-' num '.mat']);
    load(file);clear ip
    
    % store into structure model
    nz=length(z);ip(1:nz-1)=[1:nz-1];nip=length(ip);nones=ones(nip,1)';
    zd=Z;Pd=P;Ld=L;Pm=sum(Pd)/length(Pd);Lm=sum(Ld)/length(Ld);
    zm=0.5*(z(1:nz-1)+z(2:nz));
    
    p=interp1(zd,Pd,zm,'linear');
    k=interp1(zd,Ld,zm,'linear');
    indexnan=find(isnan(p));p(indexnan)=Pm;
    indexnan=find(isnan(k));k(indexnan)=Lm;
    model.ip=ip';model.z=z;model.gt=gts;model.dgt=0.;model.qb=qbs;model.dqb=0.;
    model.k=k';
    model.por=p';
    model.rho=2300*nones';
    model.cp =1000*nones';
    model.rc=model.rho.*model.cp;
    model.h=1.e-6.*nones';model.kA=0.*nones';model.kB=0*nones';

    % store into structure data
    index=find(zd>=top & zd<=bot);
    zd=zd(index);Td=T(index);
    Tm=interp1(zd,Td,z,'linear');
    id=find(~ isnan(Tm));nd=length(id);
    data.T=Tm;data.id=id;data.nd=nd;data.z=zd;
    err=0.3;data.Err=err;data.Cov=err^2*ones(1,nd);
    % save
    file=strcat(['Kola-' num '_out.mat']);
    save(file,'data','model');

    % modeling

    kl=[model.k];hl=[model.h];kAl=[model.kA];kBl=[model.kB];porl=[model.por];
    cpml=[model.cp];rhoml=[model.rho];qb=model.qb;gt=model.gt;ip=model.ip;

      
    GTemp=[-5 -10 0];GTime=[-100000 -10000]*year2sec;
    [Ts]=paleo_boxcar_smooth(t,GTemp,GTime,8);np= length(Ts);
    
  
    
    
    
%whos kl kAl kBl hl rhoml cpml porl qb ip dz dt Tini Ts
    POM=Ts(np);
    % inital values
    Ti(:,iSite)=heat1dns(kl, kAl, kBl,hl,porl,qb,POM+gt,ip,dz,maxitnl,tolnl,freeze);
    % transients
    Tini=Ti(:,iSite);Ts=Ts+gt;
 %   whos kl kAl kBl hl rhoml cpml porl qb ip dz dt Tini Ts
        
    [Tcalc,zout,tout,N,k_eff,rc_eff,ipor,lheat,rci]=heat1dnt(kl,kAl,kBl,hl,rhoml,cpml,porl,qb,...
        ip,dz,dt,Tini,Ts,theta,maxitnl,tolnl,freeze);
    % plot
    Tobs=data.T;id=data.id;err=data.Err;T1=Tobs(id);T2=Tcalc(id,nt);RES=(T1'-T2);
    figure;
    
    
    plot(T1,z(id), 'LineWidth',2,'Color','r','LineStyle','-');
    hold on
    plot(T2,z(id), 'LineWidth',2,'Color','b','LineStyle','-');
    set(gca,'YDir','reverse');grid on;ylim(ylimits);xlim(tlimits);
    xlabel('T (°C)'):ylabel ('z (m)');
    title(['Site: ' num])
    grid on;
    %legend('observed','calculated','Location','northeast');
    filename=strcat(['Kola-' num '_residual.ps']);
    exportfig(gcf,filename,opts)
    %close(gcf)
    
%     [p,S]=polyfit(z(id),T1,1);T3=polyval(p,z(id));lamnew1=qbs/p(1);
%     disp(['Coeffizients are: ',num2str(p(1)) '  ' num2str(p(2)) ])  
%     disp(['lambda (NEW) = ',num2str(lamnew1) ]) 
    
    
%     [p,S]=polyfit(T2',T1,1); T3=polyval(p,T2');
%     qbsnew=p(1)*qbs;
%     gtsnew=p(1)*gts+p(2);
%     hprodnew=hprod*p(1);
%     
%     gts=gtsnew;qbs=qbsnew;hprod=hprodnew;
%     file=strcat(['Kola-',num '_new.mat']);
%     %    save(file,'T','L','S','gts','qbs','lat','lon','elv','lam','por','lamb','gr')
%     save(file,'Z','T','L','P','gts','qbs','lat','lon','elv','lam',...
%         'por','lamb','hprod','top','bot','ok');
%     model.gts=gts;model.qbs=qbs;model.h=hprod;
%         % save
%     file=strcat(['Kola-' num '_newout.mat']);
%     save(file,'data','model');
% 
%     
% %     disp(['Coeffizients are: ',num2str(p(1)) '  ' num2str(p(2)) ])  
% %     disp(['Borehole ' num '     lambda (NEW) = ',num2str(lamnew) '   T0 (NEW)= ',num2str(gtsnew) ]) 
%       disp(['Borehole ' num ]) 
% 
% 
%     figure;
%     plot(T1 ,z(id), 'LineWidth',2,'Color','b','LineStyle','--'); hold on;
%     plot(T2 ,z(id), 'LineWidth',2,'Color','r','LineStyle','-.');
%     plot(T3 ,z(id), 'LineWidth',2,'Color','g','LineStyle','-');
%     set(gca,'YDir','reverse');grid on;ylim(ylimits);xlim(tlimits);
%     xlabel('T (^\circ C)'):ylabel ('z (m)');
%     title(['Site: ' num]);
%     text(5,1200,strcat(['H (NEW)= ',num2str(hprod)]))
%     text(5,1300,strcat(['Q_{b} (NEW)= ',num2str(qbs)]))
%     text(5,1400,strcat(['T_{0} (NEW)= ',num2str(gts)]))
% 
% 
%     grid on;legend('observed','calculated','polyfit','Location','northeast');
%     filename=strcat(['Kola-' num '_logs.ps']);
%     exportfig(gcf,filename,opts)
%     close(gcf)
%     
%     figure;
%     %plot(T1'-T2,z(id), 'LineWidth',2,'Color','r','LineStyle','-');
%     plot(T1'-T3',z(id), 'LineWidth',2,'Color','r','LineStyle','--');
%     set(gca,'YDir','reverse');grid on;ylim(ylimits);xlim(rlimits);
%     xlabel('\Delta T = T_{obs} -T_{calc} (K)'):ylabel ('z (m)');
%     title(['Site: ' num]);
%     text(-.9,1200,strcat(['H (NEW)= ',num2str(hprod)]))
%     text(-.9,1300,strcat(['Q_{b} (NEW)= ',num2str(qbs)]))
%     text(-.9,1400,strcat(['T_{0} (NEW)= ',num2str(gts)]))
% 
%     grid on;%legend('orig','polyfit','Location','northeast');
%     filename=strcat(['Kola-' num '_polyfit.ps']);
%     exportfig(gcf,filename,opts)
%     close(gcf)
% 
%     % inital values
%     Ti(:,iSite)=heat1dns(kl, kAl, kBl,hl,porl,qb,POM+gts,ip,dz,maxitnl,tolnl,freeze);
%     % transients
%     Tini=Ti(:,iSite);Ts=pt(it)+gts;
%     [Tcalc,zout,tout,N,k_eff,rc_eff,ipor,lheat,rci]=heat1dnt(kl,kAl,kBl,hl,rhoml,cpml,porl,qb,...
%         ip,dz,dt,Tini,Ts,theta,maxitnl,tolnl,freeze);
%     % plot
%     Tobs=data.T;id=data.id;err=data.Err;T1=Tobs(id);T2=Tcalc(id,nt);RES=(T1'-T2);
%        
%     
%     figure;
%     plot(T1 ,z(id), 'LineWidth',2,'Color','b','LineStyle','--'); hold on;
%     plot(T2 ,z(id), 'LineWidth',2,'Color','r','LineStyle','-.');
%     set(gca,'YDir','reverse');grid on;ylim(ylimits);xlim(tlimits);
%     xlabel('T (^\circ C)'):ylabel ('z (m)');
%     title(['Site: ' num '(Control Run)']);
% %     text(5,1200,strcat(['\lambda (NEW)= ',num2str(lamnew)]))
% %     text(5,1300,strcat(['T_0 (NEW)= ',num2str(gtsnew)]))
% 
%     grid on;legend('observed','calculated','Location','northeast');
%     filename=strcat(['Kola-' num '_logs.ps']);
%     exportfig(gcf,filename,opts)
%     close(gcf)

end

