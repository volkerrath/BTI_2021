function []=Tikh_gsth(name)
% GSTH_Tikh performs paleoclimatic inversion in parameter
% space for a single SITE

% important index arrays for indirekt addresses:
% id        points to nodes defining obs
% ip        associates parameter values of diffusivity and production to cells
% it        accociates paleotemperatures to temporal grid cells
%
% structures used:
% sitepar     information defining parameter and observations forgiven site
% fwdpar      numerical control for forward modeling
% invpar      numerical control for inversion
%
% V. R. March 2019

inv_debug=1;
year2sec=31557600;

load('common')

% PATHS
addpath([srcpath,'/src']);
addpath([srcpath,'/tools']);

dfmt=1;
%ffmt='.zip';
%archive(mfilename,strcat([mfilename '_' datestr(now,dfmt)]),ffmt);


disp([' ']);
disp(strcat([ '============================================']));
disp(['Nonlinear Tikhonov GSTH Inversion (Single Site)'])
disp([name,' --- ',datestr(now,dfmt)])
disp(strcat([ '============================================']));
if run_parallel
    disp(['Number of processes: ',num2str(parcors)])
    disp(strcat([ '============================================']));
end
disp([' ']);


if run_parallel
   if isempty(gcp('nocreate')), parpool(parcors); end
end 

% if exist('matlabpool')==0
% if isempty(gcp('nocreate')), parpool(parcors); end
% else
%      if matlabpool('size') == 0
%         matlabpool (parcors)
%      end
% end

% READ PARAMETERS AND DATA FOR SITE
disp([' ']); disp([' ...read & organize obs ' ]);
nobs=0;

F=strcat([name '_SitePar.mat']);
if exist(F,'file')
    load(F);
    disp([' >>>>> site obs read from: ' F]);
    disp([' ']);
    mstruct(sitepar);
else
    error([F ' does not exist in path!']);
end
nd=length(id);nobs=nobs+nd;

% SITE SPECIFIC PATHS
addpath([srcpath,'/src/props/',props]);
disp([' ']); disp([' ...local property model: ' props]);

% READ LOCAL MODELLING PARAMETERS
F=strcat([name,'_FwdPar.mat']);
if exist(F,'file')
    disp([' ']);disp([mfilename,': local inversion pars loaded from file ',F]);
    load(F);mstruct(fwdpar);
else
    error([F ' does not exist in path! STOP.']);
end
% READ LOCAL INVERSION PARAMETERS
F=strcat([name,'_InvPar.mat']);
if exist(F,'file')
    disp([' ']);disp([mfilename,': local inversion pars loaded from file ',F]);
    load(F);mstruct(invpar); it = pt;
else
    error([F ' does not exist in path! STOP.']);
end

% SETUP A-PRIORI AND INITIAL MODELS
disp([' ']); disp([' ...setup apriori par ' ]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iscell(m_apr_set)
    mm=m_apr_set{1};ms=m_apr_set{2};
    m_apr=mm+ms*randn(1,nsteps);
else
    m_apr=ones(1,nsteps);
    m_apr=m_apr*m_apr_set;
    disp([' ']);disp([' constant prior ',num2str(m_apr_set),' assumed. ']);
end

if iscell(m_ini_set)
    mm=m_ini_set{1};ms=m_ini_set{2};
    m_ini=mm+ms*randn(1,nsteps);
else
    m_ini=ones(1,nsteps);
    m_ini=m_ini*m_ini_set;
    disp([' ']);disp([' constant initial ',num2str(m_ini_set),' assumed. ']);
end


F=[name '_Initial.mat'];
if exist(F,'file')
    load(F)  ;
    disp([' ']);disp(['Tikh_gsth: initial conditions loaded from file ',F]);
    T0=Tinit;
else
    disp('Tikh_gsth:: no T initial conditions loaded. equilibrium assumend')
    T0=[];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = m_ini(:);m_apr=m_apr(:);npar=length(m);

disp([' ']);
disp([' ...... no of parameters = ',num2str(npar) ]);
disp([' ...... no of observations = ',num2str(nobs) ]);
disp([' ']);

% SETUP REGULARIZATION MATRIX
L0 = reg1d(npar,'l0');L1 = reg1d(npar,'l1');L2 = reg1d(npar,'l2');

% ITERATION
disp(['   ']);disp([' ... start iteration   ']);
theold = 1.e6;
for iter=1:maxiter_inv
    % SOLVE FORWARD PROBLEM
    
    dt=diff(t);nt=length(t);
    dz=diff(z);nz=length(z);

    [Tcalc,dTcalc,Qcalc,Tcon]=forward_solve(m,...
        k,kA,kB,h,p,rc,ip,dz,qb,gts,...
        it,dt,T0,theta,maxitnl,tolnl,freeze,0);
    % CALCULATE RESIDUAL
    r=Tobs(:)-Tcalc(id,nt);
    resid=r./Terr(:);
    rms=sqrt(sum(resid.*resid)/length(resid));
    mae=sum(abs(resid))/length(resid);
    disp([ ' ...... iteration ',num2str(iter-1), ...
        '  rms= ',num2str(rms),...
        '  mae= ',num2str(mae),...
        '   (',name,')' ]);
    
    % CALCULATE OBJECTIVE FUNCTION THETA
    m_iter(iter,:)=m;r_iter(iter,:)=r(:)';
    if iter > 1
        Wm=sqrt(regpar_iter(iter-1,1))*L0+ ...
            sqrt(regpar_iter(iter-1,2))*L1+ ...
            sqrt(regpar_iter(iter-1,3))*L2;
    else
        Wm=sqrt(regpar0(1))*L0+ ...
            sqrt(regpar0(2))*L1+ ...    % SOLVE FORWARD PROBLEM
            sqrt(regpar0(3))*L2;
    end
    
    theta_m_iter(iter) = norm(Wm*(m(:)-m_apr(:)))^2;
    Wd = spdiags(1./Terr,0,nobs,nobs);
    theta_d_iter(iter) = norm(Wd*r)^2;
    rms_iter(iter) = sqrt(norm(Wd*r)^2/length(r));
    theta_iter(iter)=(theta_d_iter(iter)+theta_m_iter(iter));
    
    if iter > 1
        disp([ ' ... theta = ',num2str(theta_iter(iter)), ...
            '  theta_d = ',num2str(theta_d_iter(iter)), ...
            '  theta_m = ',num2str(theta_m_iter(iter)), ...
            '  tau = (',num2str(regpar_iter(iter-1,1),'%5.2g'),...
            ' ',num2str(regpar_iter(iter-1,2),'%5.2g'), ...
            ' ',num2str(regpar_iter(iter-1,3),'%5.2g'),')']);
    end
    
    thenew = rms_iter(iter);
    if outsteps,
        if thenew < theold,
            filename=strcat([name,...
                '_iter',num2str(iter),...
                '_rms',num2str(rms_iter(iter)),'.mat']);
            save(filename,'m','theta_m_iter','theta_d_iter','rms_iter');
            
        else
            disp([ ' ... rms = ',num2str(rms_iter(iter)),...
                '  increased at iteration ',num2str(iter),' !']);
        end
    end
    
    
    if theold-thenew <= tol_inv(2), break; end
    if thenew < tol_inv(1), break; end
    
    
    theold=thenew;
    disp([' ... calculate sensitivities   ']);
    
    J=[];
    GST=m+gts;
    
    if isempty(T0)
        POM=GST(1);
        T0=heat1dns(k, kA, kB,h,p,qb,POM,dz,ip,maxitnl,tolnl,freeze,0);
    end
    
    if run_parallel
        [Jb]=sensfdt_palp(k,kA,kB,h,rc,p,qb,...
            dz,ip,dt,it,GST,T0,theta,maxitnl,tolnl,dp,freeze,0);
    else
        [Jb]=sensfdt_pals(k,kA,kB,h,rc,p,qb,...
            dz,ip,dt,it,GST,T0,theta,maxitnl,tolnl,dp,freeze,0);
    end
    J=[J;Jb(id,:)];[~,ms]=size(J);
    Jw = Wd*J(:,1:ms-1);
    %

    if iter <= start_regpar
        regpar=setregpar(regpar0,reg0par,reg1par,reg2par);
        regpar_loc = regpar;
        regpar_iter(iter,:) =  [0. 0. 0.];
        disp([' ...solving for new par']);
        A=[Jw;...
            sqrt(regpar_loc(3))*L2; ...
            sqrt(regpar_loc(2))*L1; ...
            sqrt(regpar_loc(1))*L0 ];
        rhs=[r_iter(iter,:)';...
            -sqrt(regpar_loc(3))*L2*(m(:)-m_apr(:));...
            -sqrt(regpar_loc(2))*L1*(m(:)-m_apr(:));...
            -sqrt(regpar_loc(1))*L0*(m(:)-m_apr(:))];
        %[delta_m] = cglsACB(A,rhs,maxiter_solve,tol_solve,stoprule,0);
        [delta_m] = lsqr(A,rhs,tol_solve,maxiter_solve);
        m = m + delta_m ;
        
    else if iter >= start_regpar && mod(iter,modul_regpar) == 0
            disp(['   ']);
            disp([' ...test for optimal regpar using ' reg_opt]);
            % SETUP ARRAY OF REGPAR VALUES
            regpar=setregpar(regpar0,reg0par,reg1par,reg2par);
            
                   
            [lregpar,~]=size(regpar); tol_loc=tol_solve;
            maxiter_loc=maxiter_solve;

            r_itr=r_iter(iter,:);reg_loc=reg_opt;
                        
            locpar=mstruct(k,kA,kB,h,p,c,r,rc,ip,dz,qb,gts,...
                it,dt,T0,theta,maxitnl,tolnl,freeze,Tobs,id,r_itr,Terr);
            
            parfor L=1:lregpar

                regpar_loc=regpar(L,:)
                
                [m_loc,T_loc,r_loc,theta_d,theta_m,r_norm,m_norm,gcv,upr]=...
                    solvereg(Jw,L0,L1,L2,m,m_apr,...
                    regpar_loc,tol_loc,maxiter_loc,locpar)
                
                mL(L,:) = m_loc ;
                theta_mL(L) = theta_m;
                theta_dL(L) = theta_d;
                rms_L(L)    = sqrt(r_norm^2/length(r_norm));
                theta_L(L)  = (theta_dL(L)+theta_mL(L));
                mL(L,:) = m_loc ;
                norm_r(L) = r_norm;
                norm_m(L) = m_norm;
                m_norm2(L) = norm(Wm*(m_loc-m_apr));
                
                GCV(L)=gcv;UPR(L)=upr;
                
                if inv_debug>1 disp([ ' >>>#',num2str(L),' ',reg_loc,...
                        '  regpar = ' num2str(regpar(L,:),'%5.3g %5.3g %5.3g'), ...
                        ' norm(res) = ',num2str(norm_r(L),'%5.3g'),...
                        ' norm(mod) = ',num2str(norm_m(L),'%5.3g'),...
                        ' GCV = ',num2str(GCV(L),'%5.3g'),...
                        ' UPR = ',num2str(UPR(L),'%5.3g'),...
                        ]); end
                
            end
            
            switch lower(reg_opt)
                case {'g' 'gcv'}
                    val=abs(GCV-min(GCV));
                case {'u','upr'}
                    val=abs(UPR-min(UPR));
                    
                otherwise
                    error(['GSTH_Tikh: method ' reg_opt ' not implemented!'])
            end
            index=find(val==min(val));
            if reg_shift~=0,
                index=index+reg_shift;
                index=min([index,length(val)]); index=max([index,1]);
            end
            m=mL(index,:);m=m(:);
            regpar_iter(iter,:)=regpar(index,:,:);
            disp([' ']);
            disp([' min(',upper(reg_opt),') = ',num2str(val(index)), ' at index: ',num2str(index)]);
            disp([' ...RegPar = ',num2str(regpar(index,:))]);
            
        end
        
        
    end
    %     if inv_debug
    %         filename=strcat([name '_Iteration' num2str(iter)]);
    %         save(filename);
    %     end
end


%==========================================================================
%  Postpocessing
%==========================================================================
% parameter and obs covariance matrices a posteriori
% (from generalized inverse, see NOLET 1999).
disp(['   ']);disp([' ...calculate aposteriori quantities ']);
regs=regpar_iter(iter-1,:);

A =[ Jw;sqrt(regs(3))*L2 ;sqrt(regs(2))*L1 ;sqrt(regs(1))*L0];At=A';
%At=[Jw',sqrt(regs(3))*L2',sqrt(regs(2))*L1',sqrt(regs(1))*L0'];
Jdag=inv(At*A)*Jw';
Rmm=Jdag*Jdag';Rdd=Jdag'*Jdag;Cmm=inv(At*A);

% if exist('matlabpool')==0
if ~isempty(gcp('nocreate')),delete(gcp); end
% else
%      if matlabpool('size') ~= 0
%         matlabpool ('close')
%      end
% end

% SAVE DATA
filename=strcat([name,'_results.mat']);
disp(['   ']);disp([' ...save results to: ',filename]);save(filename)
% SAVE MODEL for RESTART
[nit,nst]=size(m_iter);
m_ini=m_iter(nit,:);m_apr=m_ini;
filename=strcat([name,'_initial_out.mat']);
save(filename,'m_ini','it','dt','t')
filename=strcat([name,'_prior_out.mat']);
save(filename,'m_apr','it','dt','t')


% whos regpar_iter rms_iter theta_m_iter theta_d_iter theta_iter
% iter
S=fopen('INFO.dat','a+');
sline=strcat([...
    ' ',num2str(rms_iter(iter),'% 10.6g'),...
    ' ',num2str(theta_m_iter(iter),'%10.6g'),...
    ' ',num2str(theta_d_iter(iter),'%10.6g'),...
    ' ',num2str(theta_iter(iter),'%10.6g'),...
    ' ',num2str(iter,'%6i'),...
    ' ',num2str(regpar_iter(iter-1,:),'%10.6g %10.6g %10.6g'),...
    '  Tikh',reg_opt,' ' name]);

fprintf(S,'%s \n',sline);
fclose(S);


clear


disp(['   '])

