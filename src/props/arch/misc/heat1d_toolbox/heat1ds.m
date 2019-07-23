function [T,zout,k_eff,ipor] ...
    = heat1ds(model,linsolve,nonlinear,freeze,out)

% HEAT1DS solves nonlinear steady-state heat equation
% for general 1D grids (permafrost included).
%       T = heat1ds(model,nonlinear,freeze,out)
% calculates  numerically (FD) stationary temperatures,given a model for thermal
% conductivity, heat production, basal heat flow, and surface temperature.
% Thermal conductivity is assumed as nonlinear functions of temperature.
% If nu is the number of units, and nc the number of cells, the
% input is by three structures:
%
%  structure model:
%       model.k (1:nu)          = k
%       model.kA, kB(1:nu)      = thermal conductivity coefficient A and B
%       model.h(1:nu)		    = volumetric heat production
%       model.qb			    = basal heat flow
%       model.gt			    = boundary temperatures at the top
%       model.ip (1:nc)         = pointer to assign parameters to cells
%                   	   	  (nc gridsize in cells)
%       model.dz(1:nc)      	= cell size (m)
%                             (based on max deviation from last iteration)
%%  structure linsolve:
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
% z_out(1:nc+1)              = depth of nodes
% k_eff(1:nc)              = effective thermal conductivity
% ipor(1:nc)                = ice fraction in porosity
%
% vr, May 21, 2005

if nargin<5, out= 'no';end

if nargin<4,
    flg='y';Tf=0;w=1;
else
    flg   = freeze.flag;
    Tf    = freeze.Tf;
    w     = freeze.w;
end
if nargin<3,
    mean='g'; maxiter=10;tol=1.e-5;
else
    mean    = nonlinear.mean;
    maxiter = nonlinear.maxiter;
    tol     = nonlinear.tol;
    relax   = nonlinear.relax;
end
if nargin<2,
    solver='d';
else
    solver  = linsolve.solver;
    switch lower(solver)
        case {'b', 'bicgstb'}
            liniter = linsolve.liniter;
            lintol  = linsolve.lintol;
            precon  = linsolve.precon;
            pretol  = linsolve.pretol;
        case {'g', 'gmres'}
            liniter = linsolve.liniter;
            lintol  = linsolve.lintol;
            precon  = linsolve.precon;
            pretol  = linsolve.pretol;
            restart = linsolve.restart
    end
end


ip = model.ip;[n1 n2]=size(ip); if n1==1, ip=ip';end
k= model.k(ip);  		% thermal conductivity
kA =    model.kA(ip);           % thermal conductivity coefficient A
kB =    model.kB(ip);           % thermal conductivity coefficient B
h = model.h(ip);                % heat production
por =   model.p(ip);            % Porosity
dz =    model.dz;
[n1 n2]=size(k);    if n1==1, k=k';end
[n1 n2]=size(kA);   if n1==1, kA=kA';end
[n1 n2]=size(kB);   if n1==1, kB=kB';end
[n1 n2]=size(h);    if n1==1, h=h';end
[n1 n2]=size(por);  if n1==1, por=por';end
[n1 n2]=size(dz);   if n1==1, dz=dz';end
nc=length(ip);nz=nc+1;

qb=     model.qb;
Ts=     model.Ts;

z_out=[0;cumsum(dz)];% zc=0.5*(z(1:nz-1)+z(2:nz));Pc=998.*9.81*zc;
%
one=ones(size(ip));
zero=zeros(size(ip));


dc= 0.5 * (dz(2:nc)+dz(1:nc-1));


% NONLINEAR ITERATION LOOP
for iter=1:maxiter
    %      MATRIX A
    %      define  coefficients for interior points

    if iter==1
        km = k;ki=kiT(zeros(size(km)));
        kf=kfT(20*ones(size(km)));
        gf=one;
    else
        Tc=n2c(Titer,dz);
        km=kmT(k,Tc,kA,kB);
        ki=kiT(Tc);kf=kfT(Tc);
        % permafrost
        switch lower(flg)
            case{'y' 'yes'},
                [gf]=ftheta(Tc,Tf,w);
            otherwise
                gf=one;
        end
    end
    porm=(one-por);
    porf= por.*gf;
    pori= por-porf;
    switch lower(mean)
        case {'a' ,'ari', 'arithmetic'}
            keff  = porf.*kf +  pori.*ki + porm.*km;
        case {'g','geo','geometric'}
            keff= exp(log(kf).*porf+log(ki).*pori+log(km).*porm);
        case {'h','har','harmonic'}
            keff= 1./ (porf./kf +pori./ki + porm./km);
        case {'s','sqr','sqrmean'}
            keff=(porf.*sqrt(kf)+pori.*sqrt(ki)+porm.*sqrt(km)).^2;
        case {'n','none'}
            keff = km;
        otherwise
            disp(['mode set to arithmetic, >', mean, '<  not defined'])
            keff  = porf.*kf +  pori.*ki +  porm.*km;
    end


    if nargout > 2, k_eff=keff; end
    if lower(flg)=='y'||lower(flg)=='yes',
        if nargout > 3,ipor=pori;end
    else
        disp(strcat(['warning: no freezing parameters available!']));
    end

    dl = keff./dz;
    %      define matrix coefficients for interior points
    a(1:nz)=0.;b(1:nz)=0.;c(1:nz)=0.;q(1:nz)=0.;
    c(3:nz ) = dl(2:nc)   ./ dc;   % upper
    a(1:nz-2)= dl(1:nc-1) ./ dc;  % lower
    b(2:nz-1)= -(a(1:nz-2)+c(3:nz));   % center
    %      modify for dirichlet bc at top node
    b(1)=1.;
    %      modify for neumann at bottom node
    acn=keff(nc)/(dz(nc)*dz(nc));
    a(nz-1)  = 2.*acn;
    b(nz)    = -a(nz-1);
    %      generate sparse system matrix from diagonals a,b,and c
    A = spdiags([a' b' c'], -1:1, nz, nz);
    %      RIGHT HAND SIDE
    % nodal heat sources
    q=c2n(h,dz);
    rhs= -q';
    rhs(nz)= rhs(nz) - 2*acn*dz(nc)*qb/keff(nc);
    %  modify for dirichlet bc at top node
    rhs(1)=Ts;
    %      SOLVE by  LU od iterative methods
    switch lower(solver)
        case {'d','direct'}
            T = (A\rhs);
        case {'b','bicgstab'}
            [LA UA]=luinc(A,pretol);
            [Tn,flag,relres,itern] = bicgstab(A,rhs,lintol,liniter,LA,UA,Titer);
            T = Tn;
        case {'g','gmres'}
            [LA UA]=luinc(A,pretol);
            [Tn,flag,relres,itern] = gmres(A,rhs,lintol,liniter,LA,UA,Titer);
            T = Tn;
        otherwise
            T = (A\rhs);
    end
    if iter>1;

        checktol=max(abs(T-Titer));
        if strcmpi(lower(out),'yes') ~=0, disp([' maximal deviation from last iteration is ',...
                num2str(checktol), ' K']);end
        if checktol <= tol, break; end
        Titer=relax*T+(1-relax)*Titer;
    else
        Titer=T;
    end
end

