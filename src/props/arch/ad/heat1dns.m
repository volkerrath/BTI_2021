function [T,z,k_eff,ipor] ...
    = heat1dns(kl,kAl,kBl,hl,porl,qb,Ts,dz,ip,maxiter,tol,freeze,out)
% HEAT1DNS solves nonlinear stationary heat equation
% for general 1D grids (permafrost included).
%
% T = heat1dns(k,kA,kB,h,por,Qb,Ts,ip,dz,maxiter,tol,out) calculates
% numerically (FD) stationary temperatures, given a model for thermal
% conductivity, heat production, basal heat flow, and surface temperature.
% Thermal conductivity is assumed as nonlinear functions of temperature;
% Input :
% k (1:nu)              = lambda
% kA, kB(1:nu)	 	= thermal conductivity coefficient A and B
% h(1:nu)		= volumetric heat production
% qb			= basal heat flow
% gt			= boundary temperatures at the top
% ip (1:nc)     	= pointer to assign parameters to cells
%                 	   	(nc gridsize in cells)
% dz(1:nc)      	= cell size (m)
% dt(1:nt-1)    	= time step (s)
% maxiter,tol         	= maximal number of nonlinear Picard iterations and
%                    		   tolerance (max deviation from last iteration)
%
% Output:
% T(1:nc+1)  		= temperatures at given time steps dt
%
% V. R., May 18, 2005
solver='direct';
mean='g';Tf=0;w=1;


if nargin<11, maxiter=10;tol=1.e-5; end
if nargin<12, freeze='no';  end
if nargin<13, out= 'no';end


[n1 n2]=size(ip); if n1==1, ip=ip';end
[n1 n2]=size(dz);    if n1==1, dz=dz';end
z=[0 ;cumsum(dz)];
% zc=0.5*(z(1:nz-1)+z(2:nz));Pc=998.*9.81*zc;


nc=length(ip); nz=nc+1;

k=kl(ip(1:  nc));          % thermal conductivity
[n1 n2]=size(k);    if n1==1, k=k';end
kA =    kAl(ip(1:  nc));   % thermal conductivity coefficient A
[n1 n2]=size(kA);  if n1==1, kA=kA';end
kB =    kBl(ip(1:  nc));   % thermal conductivity coefficient B
[n1 n2]=size(kB);  if n1==1, kB=kB';end
h  =     hl(ip(1:  nc));   % heat production
[n1 n2]=size(h);   if n1==1, h=h';end
por   =  porl(ip(1:  nc));   % pcorosity
[n1 n2]=size(por); if n1==1, por=por';end
one=ones(size(ip));
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
        Tc=n2c(T,dz);

        km=kmT(k,Tc,kA,kB);
        ki=kiT(Tc);kf=kfT(Tc);
        % permafrost
        if strcmpi(lower(freeze),'yes')==1,
            [gf]=ftheta(Tc,Tf,w);
        else
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

    if nargout > 2, k_eff=keff; end
    if strcmpi(freeze,'yes')==1,
        if nargout > 3,ipor=pori;end
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
    T = (A\rhs);
    if iter>1;
        checktol=max(abs(T-Tlast));
        if strcmpi(lower(out),'yes') ~=0, disp([' maximal deviation from last iteration is ',...
                num2str(checktol), ' K']);end
        if checktol <= tol, break; end
    end
    Tlast=T;
end
