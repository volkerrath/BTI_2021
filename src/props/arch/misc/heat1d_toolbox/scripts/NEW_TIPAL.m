% TIPA performs paleoclimatic inversion in parameter space for
% multiple boreholes
%
% VR July 25, 2005
% close all;%
clc
keep  NAME debug iterout borehole model linsolve timestep nonlinear freeze inverse jacpal

disp([' ']);
disp(['====================================== ']);
disp(['Nonlinear Tikhonov Inversion for GSTH'])
disp(['- New - VR, ',date])
disp(['====================================== ']);
disp([' ']);



%==========================================================================
%  Prepocessing
%==========================================================================


% set constants
nb=length(borehole);
np=inverse.np;
tolrms =inverse.tolrms;
% setup regularization matrices

L0 = reg1d(np,'l0');
L1 = reg1d(np,'l1');
L2 = reg1d(np,'l2');


disp([' ']);
% load specific model parameter and data
% preprocess both and store into structure 'well'
% for each borehole
for ib=1:nb
    bore=borehole{ib};
    disp([' ...reading data and model for borehole ' bore]);
    file=strcat([bore '.mat']); load(file);
    disp([' ... modelling Borehole ' bore ]);

    % interpolate parameter to modelling mesh cells and
    % store into structure 'model'
    z=model.z;nz=length(z);nc=nz-1;zm=0.5*(z(1:nz-1)+z(2:nz));
    ip=model.ip;ni=length(ip);nones=ones(ni,1)';
    zd=Z;Pd=P;Ld=L;
    p=interp1(zd,Pd,zm,'linear');
    k=interp1(zd,Ld,zm,'linear');
    index=find(~isnan(p));topindex=min(index);botindex=max(index);
    p(1:topindex-1)=p(topindex);p(botindex+1:length(zm))=p(botindex);
    index=find(~isnan(k));topindex=min(index);botindex=max(index);
    k(1:topindex-1)=k(topindex);k(botindex+1:length(zm))=k(botindex);
    model.gt=gts;model.dgt=0.1;
    model.qb=qbs;model.dqb=0.005;
    model.k=k';
    model.p=p';
    model.r=2300*nones';
    model.c =840*nones';
    model.rc=model.r.*model.c;
    model.h=hprod*nones';
    model.kA=0.7*nones';
    model.kB=770*nones';

    % interpolate data to modelling mesh nodes
    %  and store into structure data
    index=find(zd>=top & zd<=bot);zd=zd(index);Td=T(index);
    Tm=interp1(zd,Td,z,'linear');
    id=find(~ isnan(Tm));nd=length(id);

    data.T=Tm;data.id=id;data.nd=nd;data.z=zd;
    err=0.25;data.Err=err;data.Cov=err^2*ones(1,nd);
    well(ib).name=bore;
    well(ib).data=data;
    well(ib).model=model;


    % save
    file=strcat([ bore '_out.mat']);
    save(file,'data','model');
end
disp([' ']);

file=strcat([ NAME '_all-boreholes.mat']);
save(file,'well');


%==========================================================================
%  Inversion
%==========================================================================


% setup prior
disp([' ']); disp([' ...setup apriori model ' ]);
m_apr=inverse.prior; m = inverse.mod;

% start inverse iteration
disp(['   ']);disp([' ...start iteration   ']);
maxiter=inverse.maxiter;
for iter=1:maxiter
    m
    r_tot=[];
    for ib=1:nb
        % forward modeling
        model=well(ib).model;
        data=well(ib).data;


        % initial values
        POM=m(1);gt=model.gt;
        model.Ts=POM+gt;
        Ti(:,ib)=heat1ds(model,linsolve,nonlinear,freeze,debug);

        % forward model
        timestep.T0 =Ti(:,ib);
        model.Ts    =m(it)+gt;
        Tcalc=heat1dt(model,timestep,linsolve,nonlinear,freeze,debug);

        % calculate residual
        Tobs=data.T;id=data.id;err=data.Err;
        resid=Tobs(id)'-Tcalc(id);res=resid/err;
        rms=sqrt(sum(res.*res)/length(res));
        r_tot=[r_tot; res];
        disp([ 'rms for iteration ',num2str(iter), ...
            ' at borehole ', borehole{ib}, ' = ',num2str(rms)])

        %         data.Tcalc(iter)(:)=Tcalc(id);
        %         data.resid(iter)=resid;
        data.rms(iter)=rms;
        well(ib).data=data;
    end
    rms_tot=sqrt(sum(r_tot.*r_tot)/length(r_tot));
    disp([ ' *** total rms for iteration ',num2str(iter), ...
        ' = ',num2str(rms_tot)])

    rnorm(3)= tau2*norm(L2*(m'-m_apr'));
    rnorm(2)= tau1*norm(L1*(m'-m_apr'));
    rnorm(1)= tau0*norm(L0*(m'-m_apr'));

    if iterout ~= 0,
        % store model
        file=strcat([ NAME '_iter-' num2str(iter) '.mat']);
        save(file,'inverse','well','timestep','m','rms_tot','rnorm1','rnorm2','rnorm3');
    end

    % convergence achieved?

    if rms_tot < tolrms, break; end


    disp([' ... calculate sensitivities   ']);
    S=[];
    for ib=1:nb
        model=well(ib).model;
        timestep.T0 =Ti(:,ib);
        m=inverse.mod;
        model.Ts=m(it)+model.gt;

        data=well(ib).data;cov=data.Cov;
        W=diag(1./sqrt(cov),0);

        inverse.mod=m+gt;
        [J]=jacfdt_pal(inverse,model,timestep,linsolve,nonlinear,freeze,debug);
        inverse.mod=m;

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
    
    % solve equations by CGLS
    stoprule   = inverse.stoprule;
    maxcgls    = inverse.maxcgls;
    tolcgls    = inverse.tolcgls;
    reorth     = inverse.reorth;        
    %[delta_m] = cglsACB(S_Aug,resid_Aug,reorth,maxcgls,tolcgls,stoprule,debug);
    [delta_m] = cgls(S_Aug,resid_Aug,maxcgls,reorth,tolcgls,debug);

    % update model 
    m = m + delta_m' ;
    inverse.mod = m;

end


%==========================================================================
%  Postpocessing
%==========================================================================
% parameter and data covariance matrices a posteriori
% (from generalized inverse, see NOLET 1999).
disp(['   ']);disp([' ...calculate aposteriori quantities ']);
S_Aug =[ S;sqrt(tau2)*L2;sqrt(tau1)*L1 ;sqrt(tau0)*L0];
S_Augt=[S',sqrt(tau2)*L2',sqrt(tau1)*L1',sqrt(tau0)*L0'];
GT_Aug=inv(S_Augt*S_Aug)*S';
Cmm=GT_Aug*GT_Aug'; Cdd=GT_Aug'*GT_Aug;



%==========================================================================
% save all
%==========================================================================

disp(['   ']);disp([' ...save results  ']);
filename=strcat(NAME,'_results.mat');
save(filename,'linsolve','timestep','nonlinear','freeze','inverse','jacpal','S')


%==========================================================================
% Plots
%==========================================================================
disp(['   ']);disp([' ...plot results for ',NAME]);

opts = struct('Format','psc2','Color','rgb','Resolution',600);
lc={'b','r','g','c','y','m','k', 'b', 'r', 'g', 'c', 'y','m','k'};
ls={'-','-','-','-','-','-','-','--','--','--','--','--','--','--'};
bc=borehole;

ty=t/year2sec;errall=sqrt(diag(Cmm));err=errall(it);mod=m;

%--------------------------------------------------------------------------
figure;
filename=strcat(NAME,'_GSTH_aposteriori.ps');

[X,Y]=stairs(-ty,mod(it));
plot(X,Y,'LineWidth',2,'Color','b','LineStyle','-');hold on;
[X,Y]=stairs(-ty,mod(it)+2*err');
plot(X,Y,'LineWidth',1,'Color','r','LineStyle','--');hold on;
[X,Y]=stairs(-ty,mod(it)-2*err');
plot(X,Y,'LineWidth',1,'Color','r','LineStyle','--');hold on;

set(gca,'XScale','log','XDir','reverse')
xlabel('time b. p.(a)','FontSize',14);ylabel('\Delta T','FontSize',14);
ylim([-15 15]);
title([NAME ': final model'],'FontSize',14)
grid on;
exportfig(gcf,filename,opts)
%close(gcf)

%--------------------------------------------------------------------------
% plot covariances
filename=strcat(NAME,'_cpp_aposteriori.eps');
figure;imagesc(Cmm);colorbar;title('C_{pp}^{aposteriori} based on Generalized Inverse');
exportfig(gcf,filename,opts);close(gcf)

filename=strcat(NAME,'_cdd_aposteriori.eps');
figure;imagesc(Cdd);colorbar;title('C_{dd}^{aposteriori} based on Generalized Inverse');
exportfig(gcf,filename,opts);close(gcf)
%close(gcf)

%--------------------------------------------------------------------------
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
    filename=strcat([NAME '-' bore '-residuals_aposteriori.ps']);
    exportfig(gcf,filename,opts)
    %close(gcf)

end


%--------------------------------------------------------------------------
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
filename=strcat(NAME,'_logs_aposteriori.ps');
exportfig(gcf,filename,opts)
%close(gcf)
