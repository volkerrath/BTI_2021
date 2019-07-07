function [T,zout,k_eff,ipor] ...
    = heat1ds(k,ka,kb,h,p,dz,ip,qb,Ts,nonlinear,freeze,out)

% HEAT1DS solves nonlinear steady-state heat equation
% for general 1D grids (permafrost included).
%       T = heat1ds(model,nonlinear,freeze,out)
% calculates  numerically (FD) stationary temperatures,given a model for thermal
% conductivity, heat production, basal heat flow, and surface temperature.
% Thermal conductivity is assumed as nonlinear functions of temperature.
% If nu is the number of units, and nc the number of cells, the
% input is by three structures:
%
%  parameters:
%  k (1:nu)         = lambda
%  kA, kB(1:nu)     = thermal conductivity coefficient A and B
%  h(1:nu)		    = volumetric heat production
%  p(1:nu)		    = porosity
%  dz(1:nc)      	= cell size (m)
%
%  qb			    = basal heat flow
%  Ts               = boundary temperature at the top 
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
% vr, April 3, 2005

if nargin<4, out= 'no';end

if nargin<3,
    flg='y';Tf=0;w=1;
else
    flg   = freeze.flag;
    Tf    = freeze.Tf;
    w     = freeze.w;
end
if nargin<2,
    mean='g'; maxiter=10;tol=1.e-5;
else
    mean    = nonlinear.mean;
    maxiter = nonlinear.maxiter;
    tol     = nonlinear.tol;
    relax   = nonlinear.relax;
end


lambda= k(ip);  			% thermal conductivity
kA =    ka(ip);           % thermal conductivity coefficient A
kB =    kb(ip);           % thermal conductivity coefficient B
hprod = h(ip);            % heat production
por =   p(ip);            % Porosity


z_out=[0 cumsum(dz)];% zc=0.5*(z(1:nz-1)+z(2:nz));Pc=998.*9.81*zc;
nc=length(ip); nz=nc+1;
one=ones(1,nc);
dc= 0.5 * (dz(2:nc)+dz(1:nc-1));


% NONLINEAR ITERATION LOOP
for iter=1:maxiter
    %      MATRIX A
    %      define  coefficients for interior points

    if iter==1
        km = lambda;ki=kiT(zeros(size(km)));
        kf=kfT(20*ones(size(km)));
        gf=one;
    else
        Tc=n2c(Titer,dz);
        km=kmT(lambda,Tc,kA,kB);
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
            disp(['WMEAN: mode set to arithmetic, >', mean, '<  not defined'])
            keff  = porf.*kf +  pori.*ki +  porm.*km;
    end
    ipor=pori';k_eff=keff';



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
    q=c2n(hprod,dz);
    rhs= -q';
    rhs(nz)= rhs(nz) - 2*acn*dz(nc)*qb/keff(nc);
    %  modify for dirichlet bc at top node
    rhs(1)=Ts;
    %      SOLVE by  LU
    T = (A\rhs);
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

