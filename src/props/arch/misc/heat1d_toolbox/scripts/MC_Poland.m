% MC performs paleoclimatic monte carlo inversion in parameter space for
% multiple boreholes

% important index arrays for indirect addresses:
% id        points to nodes defining data
% ip        associates parameter values of diffusivity and production to cells
% it        accociates paleotemperatures to temporal grid cells
%
% structures used:
% unit      information associated with geological units (layers)
% model     information defining model for given borehole
% obs       information concerning the observations in borehole
%
%
% V. R. April 2, 2003
%close all;clear all
close all;
keep kc stoprule tolcgls reorth tolrms maxiter tau0 tau1 tau2...
    borehole freeze NAME m_apr POR debug nsteps delta_m maxdelta q_range qmin middle...
    gt_range gtmin

disp(['Nonlinear Monte Carlo Inversion for Paleoclimate from Multiple Boreholes'])
disp(['- Permafrost -'])
disp(['V. R. and D. M., ',date] )

disp([' ']);
% READ MODEL PARAMETER & DATA
nb=length(borehole);

% GENERATE SPATIAL MESH
zstart=10;zend=4000;nz=251;type='logarithmic';dir=1;
disp([' ']);disp([' ...set up ' type ' spatial mesh ']);
[z,dz]=set_mesh(zstart,zend,nz,type,dir,1);

dz=dz';
% GENERATE TEMPORAL MESH
year2sec=31557600;
tstart=150000*year2sec;tend=10*year2sec;nt=257;type='logarithmic';dir=-1;
disp([' ']);disp([' ...set up ' type ' temporal mesh ']);
[t,dt]=set_mesh(tstart,tend,nt,type,dir,1);

% DEFINE LOGARITHMIC INVERSION GRID
disp(['   ']);disp([' ...set up parametrization for paleoclimate inversion  ']);
%nsteps=10;
ia=0; % Variable for counting accepted runs

base=0.;
[pt,it]=set_paleo_grid(t,base,tstart,tend,nsteps);
np= length(pt);

delta_m=4; % max. Intervall f?r Modellgenerierung

% INVERSION SETUP
disp(['   ']);disp([' ...set up inversion control parameters  ']);

% PARAMETER FOR FWD CALCULATIONS
theta=0.5*ones(1,nt);
theta(1:10)=1.;

% NONLINEAR ITERATION PARAMETERS
maxitnl=1;
tolnl=0.00021;

% SENSITVITY ITERATION PARAMETERS
dp=0.1;
maxitdp=1;
toldp=0.1;

% SETUP MC RESULT MATRIX
M=NaN(maxiter,nsteps+3); % number of time steps plus q, rms, S


disp([' ']);
for ib=1:nb
    bore=borehole{ib};
    disp([' ...reading data and model for borehole ' bore]);

    % READ MEASUREMENTS AND DATA FOR BOREHOLE
    file=strcat([bore '_prep.mat']);
    load(file);clear ip
    % store into structure model
    nz=length(z);
    ip(1:nz-1)=[1:nz-1];nip=length(ip);nones=ones(nip,1)';
    zd=Z;Pd=P;Ld=L;
    Pm=sum(Pd)/length(Pd);
    Lm=sum(Ld)/length(Ld);
    zm=0.5*(z(1:nz-1)+z(2:nz));

    p=interp1(zd,Pd,zm,'linear');
    k=interp1(zd,Ld,zm,'linear');
    indexnan=find(isnan(p));
    p(indexnan)=Pm;
    indexnan=find(isnan(k));
    k(indexnan)=Lm;

    model.ip=ip';model.z=z;model.gt=gts;model.dgt=0.;model.qb=qbs;model.dqb=0.;
    model.k=k';
    model.por=p';
    model.rho=1000*nones';
    model.cp =850*nones';
    model.h=0.*nones';model.kA=0.0013*nones';model.kB=0.0029*nones';  %hier Koeffizienten!!!
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
    data.Cov=err^2*ones(1,nd); % s^2 in Metropolis rule?
    well(ib).data=data;
    well(ib).model=model;

    file=strcat(['MC-' bore '-out.mat']);
    save(file,'data','model');

end
disp([' ']);

% end variable input
% SETUP INITIAL VALUES FOR TEMPERATURE AT TSTART
disp([' ']); disp([' ...setup initial values ' ]);

for ib=1:nb
    model=well(ib).model;
    kl=[model.k];hl=[model.h];kAl=[model.kA];kBl=[model.kB];porl=[model.por];
    qb=model.qb;gt=model.gt;

    Tini=heat1dns(kl, kAl, kBl,hl,porl,qb,gt,ip,dz,maxitnl,tolnl,freeze);
    Ti(:,ib)=Tini;
end
disp([' ']);


% SETU A-PRIORI MODEL
disp([' ']); disp([' ...setup apriori model ' ]);
% [m_apr]=set_gst_prior('constant',struct('n',np,'gt', 10) )
%gtm=mean(gtb);qbm=mean(qbb);

m_apr=zeros(1,np);
m_0 = m_apr;

m = m_0;
m_new=m_0;

% ITERATION
disp(['   ']);disp([' ...start iteration   ']);
for iter=1:maxiter

    % generate new model

    loc=randperm(nsteps+2);  % choose location of change or q or gt

    delta=delta_m*rand(1,1) - delta_m/2;% choose temperature shift.

    q_delta=q_range*rand(1,1) + qmin; % choose q-shift.

    gt_delta=gt_range*rand(1,1) + gtmin; % choose gt-shift.


    m_old=m;
    qb_old=qb;
    qmc=qb;
    gt_old=gt;
    gtmc=gt;

    if loc(1)<=nsteps
        
        m_new(loc(1))=delta; % new value at random entry
        
        if m_new(loc(1))<=-maxdelta+middle
            m_new(loc(1))=-maxdelta+middle;
        end

        if m_new(loc(1))>=maxdelta+middle
            m_new(loc(1))=maxdelta+middle;
        end

    end

    if loc(1)==nsteps+1
        qmc=q_delta;

    end
    if loc(1)==nsteps+2
        gtmc=gt_delta;

    end



    gt=gtmc;
    qb=qmc;
    m=m_new;

    Ts=m(it);
    r_tot=[];
    for ib=1:nb

        % FORWARD MODEL
        model=well(ib).model;
        kl=[model.k];
        hl=[model.h];
        kAl=[model.kA];kBl=[model.kB];
        porl=[model.por];
        cpml=[model.cp];
        rhoml=[model.rho];
        % qb=model.qb; in MC!
        % gt=model.gt; % auch gt in MC!
        ip=model.ip;


        Tini=Ti(:,ib);Ts=m(it)+gt; % hier neue GSTH anwenden...

        [Tcalc,zout,tout,N,k_eff,rc_eff,ipor,lheat,rci]=heat1dnt(kl,kAl,kBl,hl,rhoml,cpml,porl,qb,...
            ip,dz,dt,Tini,Ts,theta,maxitnl,tolnl,freeze);

        % CALCULATE RESIDUAL
        data=well(ib).data;
        Tobs=data.T;
        id=data.id;err=data.Err;
        resid=Tobs(id)'-Tcalc(id,nt);

        if iter==1
            S=1;
        end

        Sold=S;
        zgrad=z(id); tgrad_obs=Tobs(id); tgrad_calc=Tcalc(id,nt)';
        Aobs=vertcat(zgrad,tgrad_obs)';
        Acalc=vertcat(zgrad,tgrad_calc)';

        [zgradobs,gradobs]=tgrad(Aobs);
        [zgradcalc,gradcalc]=tgrad(Acalc);

        gradresid=gradobs/1000-gradcalc/1000;

        %


        S=0.5*(norm(gradresid))^2; % Delta S in Metropolis rule!

        %S=0.5*(norm(resid))^2; % Delta S in Metropolis rule!


        res=resid/err;
        rms=norm(res)/sqrt(length(res));r_tot=[r_tot; res];
        disp([ 'rms for iteration ',num2str(iter), ...
            ' at borehole ', borehole{ib}, ' = ',num2str(rms)])

        result.Tobs=Tobs(id);result.Tcalc=Tcalc(id,nt);
        result.z = z(id); result.resid= resid;result.errT=err;
        result.rms=rms;
        well(ib).result=result;


    end



    %     if iter > 1, rms_tot_old=rms_tot; end

    rms_tot=norm(r_tot)/sqrt(length(r_tot));


    disp([ ' *** total S and Sold for iteration ',num2str(iter), ...
        ' = ',num2str(S) '  ' num2str(Sold)])
    %     if rms_tot < tolrms, break; end

    % testing results

    %P=exp(-S/err^2); % Probability accepting model

    P=exp(-S/2.75e-3); % Probability accepting model (gradient)

    iP=round(1/P);

    if  iP<1000
        nP=randperm(iP);    % Probability of first entry equals P
    end

    af=0;  % Acception flag

    if  Sold >S
        if loc(1)<=nsteps
            m=m_new;af=1;
        end
        if loc(1)==nsteps+1
            qb=qmc;af=1;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          