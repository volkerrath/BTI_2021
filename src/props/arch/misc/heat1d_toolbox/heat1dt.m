function [T,zo,to,N,k_eff,rc_eff,ipor,lheat,rci]=...
    heat1dt(model,timestep,linsolve,nonlinear,freeze,out)
% HEAT1DT  solves the nonlinear time-dependent heat equation
% for general 1-D grid with permafrost.
%       T = heat1dt(model,timestep,nonlinear,freeze,out)
% calculates  numerically (FD) stationary temperatures, given a model for
% thermalconductivity, heat production, basal heat flow, and surface
% temperature. % Thermal properties are assumed as nonlinear functions
% of temperature (and pressure). If nu is the number of units, and nc
% the number of cells, the input is controlled by four structures:
%
%  structure model:
%       model.k (1:nu)          = k
%       model.kA, kB(1:nu)      = thermal conductivity coefficient A and B
%       model.h(1:nu)		    = volumetric heat production
%       model.r(1:nu)		    = density of rock
%       model.r(1:nu)		    = heat capacity of rock
%       model.dz(1:nc)      	= cell size (m)
%
%       model.qb			    = basal heat flow
%       model.Ts                = time-dependent boundary temperatures
%                               at the top (e.g., paleoclimate)
%       model.it                = ip vector mapping model.Ts to
%                               actual time steps defined in tinestep.
%
%  structure timestep:
%       timestep.dt (1:nt)      = time step (s)
%       timestep.theta(1:nt-1) 	= time stepping control parameter,
%                      		      (.5 = Crank-Nicholson, 1=Backward Euler)
%       timestep.T0(1:nc+1)    	= initial temperatures
%       timestep.out         	= true if output for all time steps required
%
%  structure linsolve:
%       linsolve.solver         = 'direct','bicgstab','gmres'  
%       if not 'direct':
%           linsolv.liniter         = number of iterations if not 'direct'
%           linsolv.lintol          = stopping tollerance  if not 'direct'
%           linsolv.precon          = preconditioner (ILU only)
%           linsolv.pretol          = tollerance for ILU
%           linsolv.prefreq         = frequency of preconditioner updates
%           linsolv.restart         = restart (for 'gmres')
%
%  structure nonlinear:
%       nonlinear.mean      = n/a/g/s calculation of effective properties
%                             (fix, arithmetic, geometric, square-root)
%       nonlinear.maxiter   = maximal number of nonlinear Picard iterations
%       nonlinear.tol       = stopping tolerance for Picard iterations
%                             (based on max deviation from last iteration)
%       nonlinear.relax     = relaxation parameter for Picard iterations
%                             (usually 0 < relax <1)
%
%  structure freeze:
%       freeze.flag         = yes/no
%       freeze.Tf           = freezing temperature (default 0C)
%       freeze.w            = width of freezing interval (in K, ususlly 1)
%
%  out                      = debug output
%
% On output:
% T(1:nc+1)          		= temperatures at nodes zout
% z(1:nc+1)                 = depth of nodes
% t(1:nt)                   = times
% N(1:nt)                   = Neuman numbers
% k_eff(1:nc)               = effective thermal conductivity
% rc_eff(1:nc)              = effective volumetric heat capacity
% lheat(1:nc)               = latent heat
% ipor(1:nc)                = ice fraction in porosity
% rci(1;nc)                 = volumetric heat capacity od ice fraction
%
% vr, May 21, 2005

if nargin<6, out= 'no';end

if nargin<5,
    flg='n';Tf=0;w=1;Lh=333600;
else
    flg   = freeze.flag;
    Tf    = freeze.Tf;
    w     = freeze.w;
    Lh    = freeze.Lh;
end

if nargin<4,
    mean='g'; maxiter=10;tol=1.e-5;relax=0.5;
else
    mean    = nonlinear.mean;
    maxiter = nonlinear.maxiter;
    tol     = nonlinear.tol;
    relax   = nonlinear.relax;
end

if nargin<3,
    solver='d';
else
    solver  = linsolve.solver;
    switch lower(solver)
        case {'b', 'bicgstb'}
            liniter = linsolve.liniter;
            lintol  = linsolve.lintol;
            precon  = linsolve.precon;
            pretol  = linsolve.pretol;
            prefreq = linsolve.prefreq;
        case {'g', 'gmres'}
            liniter = linsolve.liniter;
            lintol  = linsolve.lintol;
            precon  = linsolve.precon;
            pretol  = linsolve.pretol;
            prefreq = linsolve.prefreq;
            restart = linsolve.restart
    end
end

if nargin<2,
    error([' STOP. No timestep given. ']); return;
else
    dt      = timestep.dt;
    theta   = timestep.theta;
    T0      = timestep.T0;
    stepsout= timestep.out;
    [n1 n2]=size(dt);    if n1==1, dt=dt';end
    nt=length(dt)+1;
end

% setup model
ip = model.ip;[n1 n2]=size(ip); if n1==1, ip=ip';end
k= model.k(ip);  		% thermal conductivity
kA =    model.kA(ip);         % thermal conductivity coefficient A
kB =    model.kB(ip);         % thermal conductivity coefficient B
h = model.h(ip);          % heat production
por =   model.p(ip);        % Porosity
rhos  = model.r(ip);        % density
cps   = model.c(ip);        % heat capacity
dz =    model.dz;
[n1 n2]=size(k);    if n1==1, k=k';end
[n1 n2]=size(kA);   if n1==1, kA=kA';end
[n1 n2]=size(kB);   if n1==1, kB=kB';end
[n1 n2]=size(h);    if n1==1, h=h';end
[n1 n2]=size(rhos); if n1==1, rhom=rhom';end
[n1 n2]=size(cps);  if n1==1, cp=cp';end
[n1 n2]=size(por);  if n1==1, por=por';end
[n1 n2]=size(dz);    if n1==1, dz=dz';end
nc=length(ip);nz=nc+1;


qb =    model.qb;
Ts =    model.Ts;


%
one=ones(size(ip));
zero=zeros(size(ip));

%INITIAL VALUES

% initialize time, depth and pressure
to=[0;cumsum(dt)];
zo=[0;cumsum(dz)];
zc=0.5*(zo(1:nz-1)+zo(2:nz));Pc=998.*9.81*zc;


Tlast=T0;
Titer=T0;
Tlast(1)=Ts(1);
Titer(1)=Ts(1);

dc= 0.5*(dz(2:nc)+dz(1:nc-1));

I= speye(nz,nz);

% START TIME STEPPING
for i = 1:nt-1

    %   start NONLINEAR ITERATION

    for  iter =1:maxiter
        %        build SYSTEM MATRIX
        %        initialize diagonals
        a(1:nz)=0.;b(1:nz)=0.;c(1:nz)=0.;
        %        define  coefficients for interior points

        if iter==1
            Tc=n2c(Tlast,dz);
        else
            Tlast=relax*Titer+(1-relax)*Tlast;
            Tc=n2c(Tlast,dz);
        end
        Pc=9.81*rhofT(Tc,Pc).*zc;
        
        rhoi=rhoiT(Tc);
        ki=kiT(Tc);  
        cpi=cpiT(Tc);
        
        rhof=rhofT(Tc,Pc); 
        kf=kfT(Tc);  
        cpf=cpfT(Tc,Pc);
        
        rhom=rhos;  
        km=kmT(k,Tc,kA,kB); 
        cpm=cpmT(cps,Tc); 

        %  permafrost

        switch lower(flg)
            case{'y' 'yes'},
                [gf dgf]=ftheta(Tc,Tf,w);
            otherwise
                gf=one;dgf=zero;
        end

        porm=one-por;
        porf= por.*gf;
        pori=por-porf;
        switch lower(mean)
            case {'a' ,'ari', 'arithmetic'}
                keff  = porf.*kf +  pori.*ki + porm.*km;
            case {'g','geo','geometric'}
                keff= exp(log(kf).*porf+log(ki).*pori+log(km).*porm);
            case {'h','har','harmonic'}
                keff= 1./ (porf./kf +pori./ki + porm./km);
            case {'s','sqr','sqrmean'}
                keff=(porf.*sqrt(kf)+pori.*sqrt(ki)+porm.*sqrt(km)).^2;
            otherwise
                disp(['WMEAN: mode set to arithmetic, >', mean, '<  not defined'])
                keff  = porf.*kf +  pori.*ki +  porm.*km;
        end
        rceff = porm.*rhom.*cpm + ...
            pori.*rhoi.*cpi +  ...
            porf.*rhof.*cpf + por.*rhof.*Lh.*dgf;

        if nargout > 4 && intermediate==1, k_eff(:,i)=keff'; end
        if nargout > 5 && intermediate==1, rc_eff(:,i)=rceff'; end
        if lower(flg)=='y'||lower(flg)=='yes',
            if nargout > 6 && intermediate==1, ipor(:,i)=(pori./por)'; end
            if nargout > 7 && intermediate==1, lheat(:,i)=(por.*rhof.*Lh.*dgf)'; end
            if nargout > 8 && intermediate==1, rci(:,i)=(rhoi.*cpi)'; end
        else
            disp(strcat(['Warning: No freezing parameters available!' ]));
        end

        dl = keff./dz;
        N(i)=max((keff*dt(i))./(rceff.*dz.^2));

        %  define matrix coefficients for interior points
        c(3:nz ) = dl(2:nc)   ./ dc;   % upper
        a(1:nz-2)= dl(1:nc-1) ./ dc ;  % lower
        b(2:nz-1)= -(a(1:nz-2)+c(3:nz));   % center
        % modify for dirichlet bc at top node
        b(1)=1.;
        % modify for neumann at bottom node
        acn=keff(nc)/(dz(nc)*dz(nc));
        a(nz-1)  = 2.*acn;
        b(nz)    = -a(nz-1);
        %       generate sparse system matrix from diagonals a,b,and c
        A = spdiags([a' b' c'], -1:1, nz, nz);
        %       build RIGHT HAND SIDE
        %       nodal heat sources
        qd=c2n(h,dz);
        rc=c2n(rceff,dz);
        rhs=qd';
        %       modify for neumann at bottom node
        rhs(nz)= rhs(nz) + 2*acn*dz(nc)*qb/keff(nc);
        % setup matrix L
        F= spdiags( 1./rc', 0, nz, nz);
        A=F*A; rhs=rhs./rc';
        L = I-dt(i)*theta(i)*A;
        %       right hand side
        r =(I+dt(i)*(1-theta(i))*A)*Tlast+dt(i)*rhs;
        %       modify for top (dirichlet) bc:
        r(1) = Ts(i+1); L(1,1)=1.;


        %       solve by LU or iterativ solver
        switch lower(solver)
            case {'d','direct'}
                Tnew = (L\r);
            case {'b','bicgstab'}
                if  (mod(iter,prefreq)==1), [LA UA]=luinc(A,pretol);end
                [Tn,flag,relres,niter] = bicgstab(A,rhs,lintol,liniter,LA,UA,Titer);
                Tnew = Tn;
            case {'g','gmres'}
                if  (mod(iter,prefreq)==1), [LA UA]=luinc(A,pretol);end
                [Tn,flag,relres,niter] = gmres(A,rhs,lintol,liniter,LA,UA,Titer);
                Tnew = Tn;
            otherwise
                Tnew = (L\r);
        end

        if iter>1;
            checktol=norm(abs(Tnew-Titer),inf);
            if strcmpi(out,'yes')==1,
                disp([' maximal deviation of fix point iteration at timestep ',...
                    num2str(i),' is ', num2str(checktol), ' K  (',...
                    num2str(iter),')']);end
            if checktol <= tol, break; end
            if iter >= maxiter, break; end
        end
        Titer=Tnew;
        %    end NONLINEAR ITERATION
    end
    Tlast=Tnew;
    if stepsout, T(:,i+1)=Tnew; end
    % end TIME STEPPING
end

if ~stepsout,
    T=Tnew;
    if nargout > 4,k_eff=keff';end
    if nargout > 5,rc_eff=rceff',end
    if lower(flg)=='y'||lower(flg)=='yes',
        if nargout > 6,  ipor=(pori./por)'; end
        if nargout > 7, lheat=(por.*rhof.*Lh.*dgf)'; end
        if nargout > 8, rci=(rhoi.*cpi)'; end
    end

else
    T(:,1)=T0;
end







