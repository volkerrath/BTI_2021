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
% V. R. April 17, 2005
%close all;clear all
close all;
keep borehole nonlinear timestep freeze inverse jacpal jacpet debug NAME


disp(['Nonlinear Tikhonov Inversion for Paleoclimate from Multiple Boreholes'])
disp(['- Permafrost -'])
disp(['V. R., ',date] )

disp([' ']);
% READ MODEL PARAMETER & DATA
nb=length(borehole);
POM=m(1);

% SETUP REGULARIZATION MATRIX
np= length(pt);
L0 = reg1d(np,'l0');
L1 = reg1d(np,'l1');
L2 = reg1d(np,'l2');

disp([' ']);
for ib=1:nb
    bore=borehole{ib};
    disp([' ...reading data and model for borehole ' bore]);

    % READ MEASUREMENTS AND DATA FOR BOREHOLE
    file=strcat([ bore '.mat']);
    load(file);clear ip
    % store into structure model
    nz=length(z);
    ip(1:nz-1)=[1:nz-1];nip=length(ip);nones=ones(nip,1)';
    zd=Z;Pd=P;Ld=L;
    Pm=sum(Pd)/length(Pd);
    Lm=sum(Ld)/length(Ld);
    p=interp1(zd,Pd,z,'linear');
    k=interp1(zd,Ld,z,'linear');
    indexnan=find(isnan(p));
    p(indexnan)=Pm;
    indexnan=find(isnan(k));
    k(indexnan)=Lm;

    model.ip=ip;
    model.z=z;
    model.gt=gts;
    model.dgt=0.;
    model.qb=qbs;
    model.dqb=0.;
    model.k=k;
    model.kA=0.7*nones;
    model.kB=770*nones;
    model.por=p;
    model.rho=2650*nones;
    model.cp =850*nones;
    model.h=0.*nones;
    % store into structure data
    index=find(zd>=top & zd<=bot);
    zd=zd(index);Td=T(index);
    Tm=interp1(zd,Td,z,'linear');
    id=find(~ isnan(Tm));
    nd=length(id);
    data.T=Tm;data.id=id;data.nd=nd;data.z=zd;
    % under construction ....
    err=0.5;
    data.Err=err;
    data.Cov=err^2*ones(1,nd);
    well(ib).data=data;
    well(ib).model=model;

    file=strcat([ bore '-out.mat']);
    save(file,'data','model');
    %    ylimits=[0 1800];
    %    tlimits=[0 24];
    %    climits=[0 6];
    %    hlimits=[10 70];
    %    figure;
    %    subplot(1,2,1);
    %    plot(Tm(id),z(id),':g','LineWidth',2);grid on;hold on;
    %    ylim(ylimits);xlim(tlimits);
    %    ylabel('z (m)');xlabel('T (?C)');
    %    set(gca,'YDir','reverse')
    %    subplot(1,2,2);
    %    plot(k,z,'ob');hold on; grid on;
    %    ylim(ylimits);xlim(climits);
    %    ylabel('');xlabel('\lambda (J m^{-1}K^{-1})');
    %    set(gca,'YDir','reverse')
    %    suptitle(['Borehole:' bore])
    %    file=strcat(['Kola_' bore '.ps']);
    %    saveas(gcf,file,'epsc2');


end
disp([' ']);


m_0 = m_apr; m = m_0;inverse.m=m;
% ITERATION
disp(['   ']);disp([' ...start iteration   ']);
for iter=1:maxiter

    Ts=m(it);
    r_tot=[];
    for ib=1:nb

        % FORWARD MODEL
        model=well(ib).model;
        gt=model.gt;
        model.gt=model.gt+m(1);
        Tini=heat1ds(model,nonlinear,freeze,out);
        model.gt=gt;model.Ts=m+gt;
        Tcalc=(model,timestep,nonlinear,freeze,out)

        % CALCULATE SCALED RESIDUAL
        data=well(ib).data;Tobs=data.T;
        id=data.id;err=data.Err;
        resid=Tobs(id)'-Tcalc(id,nt);res=resid/err;
        rms=norm(res)/sqrt(length(res));r_tot=[r_tot; res];
        disp([ 'rms for iteration ',num2str(iter), ...
            ' at borehole ', borehole{ib}, ' = ',num2str(rms)])

        % STORE RESULTS FOR THIS SITE AND THIS ITERATE
        result.Tobs=Tobs(id);result.Tcalc=Tcalc(id,nt);
        result.z = z(id); result.resid= resid;result.errT=err;
        result.rms=rms;
        well(ib).result=result;
    end
    %SUM UP AND COMPARE TO TOLERANCWE
    rms_tot=norm(r_tot)/sqrt(length(r_tot));
    disp([ ' *** total rms for iteration ',num2str(iter), ...
        ' = ',num2str(rms_tot)])
    if rms_tot < inverse.tolrms, break; end

    % SETUP JACOBIANS AND BUKDUP WEIGHTED MATRIX
    disp([' ... calculate sensitivities   ']);
    S=[];
    for ib=1:nb
        data=well(ib).data;
        Tobs=data.T;id=data.id;err=data.Err;cov=data.Cov;
        W=diag(1./sqrt(cov),0);
        model=well(ib).model;
        model.Ts=m+model.gt;
        [J]=jacfdt_pal(jacpal,model,timestep,nonlinear,freeze,out)
        S=[S;W*J(id,:)];
    end
    % CALCULATE MODEL UPDATE
    disp([' ... calculate parameter increment ']);
    S_Aug=[              S;...
        sqrt(tau2)*L2; ...
        sqrt(tau1)*L1; ...
        sqrt(tau0)*L0 ];
    r_Aug=[                     r_tot;...
        -sqrt(tau2)*L2*(m'-m_apr');...
        -sqrt(tau1)*L1*(m'-m_apr');...
        -sqrt(tau0)*L0*(m'-m_apr')];
    [delta_m] = cgls(S_Aug,r_Aug,kc,reorth,tolcgls);
    % UPDATE MODEL
    m = m + delta_m' ;inverse.m=m;
    
%     % STORE FOR PLOTTING
%     m_iter(iter,:)=m(1,:);
%     rms_iter(iter)=rms;
%     t_iter(iter,:)=
%     Tcalc(:,nt)';
end
filename=strcat(NAME,'_results.mat');
save;

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
filename=strcat(NAME,'_results.mat');
save(filname);


ty=t/year2sec;errall=sqrt(diag(Cmm));err=errall(it);mod=m;

% PLOTS
disp(['   ']);disp([' ...plot results for ',NAME]);

opts = struct('Format','psc2','Color','rgb','Resolution',600);
lc={'b','r','g','c','y','m','k'};
ls={'-','-','-','-','-','-','-'};
bc=borehole;
% bc={'2271' ,'2731', '2908' ,'3200' ,'3356' ,'3359'};

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
ylim([-6 2]);
%text(1000,-3.00,['L_{1}  = ',num2str(tau1)],'FontSize',14)
%text(1000,-4.00,['rms  = ',num2str(rms),' K '],'FontSize',14)
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
%text(-0.4,2500,['L_{2}  = ',num2str(tau1)],'FontSize',14)
%text(-0.4,3000,['rms  = ',num2str(rms),' K '],'FontSize',14)
title([NAME ': Residuals'],'FontSize',14)
grid on;legend(cellstr(bc),4);
filename=strcat(NAME,'_aposteriori_log.ps');
exportfig(gcf,filename,opts)
%close(gcf)
filename=strcat(NAME,'_cpp_aposteriori.eps');
figure;imagesc(Cmm);colorbar;title('C_{pp}^{aposteriori} based on Generalized Inverse');
exportfig(gcf,filename,opts);close(gcf)
filename=strcat(NAME,'_cdd_aposteriori.eps');
figure;imagesc(Cdd);colorbar;title('C_{dd}^{aposteriori} based on Generalized Inverse');
exportfig(gcf,filename,opts);close(gcf)
%close(gcf)
