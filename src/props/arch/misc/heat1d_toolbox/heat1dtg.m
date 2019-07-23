function [Tx,T,zo,to,N,k_eff,rc_eff,ipor,lheat,rci]=...
         heat1dt(k,ka,kb,h,r,c,p,dz,ip,qb,Ts,it,dt,theta,T0,nonlinear,freeze,out)

% HEAT1DT  solves nonlinear time-dependent heat equation
% for general 1-D grid with permafrost.
%       T = heat1dt(model,timestep,nonlinear,freeze,out)
% calculates  numerically (FD) stationary temperatures, given a model for 
% thermalconductivity, heat production, basal heat flow, and surface 
% temperature. % Thermal properties are assumed as nonlinear functions 
% of temperature (and pressure). If nu is the number of units, and nc 
% the number of cells, the input is controlled by some input parameters
% including sime structures:
%
%  parameters:
%  k (1:nu)         = lambda
%  kA, kB(1:nu)     = thermal conductivity coefficient A and B
%  h(1:nu)		    = volumetric heat production
%  r(1:nu)		    = density of rock
%  c(1:nu)		    = heat capacity of rock
%  r(1:nu)		    = porosity of rock
%  dz(1:nc)      	= cell size (m)
%
%  qb			    = basal heat flow
%  Ts               = time-dependent boundary temperatures
%                      at the top (e.g., paleoclimate)
%  it               = index vector mapping model.Ts to
%                      actual time steps defined in tinestep.
%
%  dt (1:nt)        = time step (s)
%  theta(1:nt-1) 	= time stepping control parameter,
%                      		      (.5 = Crank-Nicholson, 1=Backward Euler)
%  T0(1:nc+1)    	= initial temperatures
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
% vr, April 3, 2005
if nargin<4, out= 'no';end

if nargin<4,
    flg='n';Tf=0;w=1;Lh=333600;
else
    flg   = freeze.flag;
    Tf    = freeze.Tf;
    w     = freeze.w;
    Lh    = freeze.Lh;
end
if nargin<3,
    mean='g'; maxiter=10;tol=1.e-5;relax=0.5;
else
    mean    = nonlinear.mean;
    maxiter = nonlinear.maxiter;
    tol     = nonlinear.tol;
    relax   = nonlinear.relax;
end

if nargin<2, error([' STOP. No timestep given. ']); return; end


lambda= k(ip);  		% thermal conductivity
kA =    ka(ip);         % thermal conductivity coefficient A
kB =    kb(ip);         % thermal conductivity coefficient B
hprod = h(ip);          % heat production
por =   p(ip);        % Porosity
rhom  = r(ip);        % density
cpm   = c(ip);        % heat capacity
Tsi =    Ts(it);
dz=dz;

nc=length(ip);nz=nc+1;
one=ones(nc,1);zero=zeros(nc,1);

% dt      = timestep.dt;
% theta   = timestep.theta;
% T0      = timestep.T0;
nt=length(dt)+1;

%INITIAL VALUES

% initialize time, depth and pressure
to=[0 cumsum(dt)]';
zo=[0 cumsum(dz')]';
zc=0.5*(zo(1:nz-1)+zo(2:nz));Pc=998.*9.81*zc;

T(:,1)=T0;T(1,1)=Tsi(1);
Tlast=T0;
dc= 0.5 * (dz(2:nc)+dz(1:nc-1));

% START TIME STEPPING
for i = 1:nt-1

    %   start NONLINEAR ITERATION

    for  iter =1:maxiter
        %        build SYSTEM MATRIX
        %        initialize diagonals
        a(1:nz)=0.;b(1:nz)=0.;c(1:nz)=0.;
        %        define  coefficients for interior points

        if iter==1
            Tc=n2c(Tlast,dz)';
        else
            Tlast=relax*Titer+(1-relax)*Tlast;
            Tc=n2c(Tlast,dz)';
        end
        Pc=9.81*rhofT(Tc,Pc).*zc;
        rhoi=rhoiT(Tc);    ki=kiT(Tc);  cpi=cpiT(Tc);
        rhof=rhofT(Tc,Pc); kf=kfT(Tc);  cpf=cpfT(Tc,Pc);
        cpm=cpmT(cpm,Tc);  km=kmT(lambda,Tc,kA,kB);

        %  permafrost

        switch lower(flg)
            case{'y' 'yes'},
                [gf dgf]=ftheta(Tc,Tf,w);
%                 gf=gf';dgf=dgf';
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

        ipor(:,i)=(pori./por)';k_eff(:,i)=keff';rc_eff(:,i)=rceff';
        lheat(:,i)=(por.*rhof.*Lh.*dgf)';rci(:,i)=(rhoi.*cpi)';

        dl = keff./dz;
        N(i)=max((keff*dt(i))./(rceff.*dz.^2));

        %        define matrix coefficients for interior points
        c(3:nz ) = dl(2:nc)   ./ dc;   % upper
        a(1:nz-2)= dl(1:nc-1) ./ dc ;  % lower
        b(2:nz-1)= -(a(1:nz-2)+c(3:nz));   % center
        %       modify for dirichlet bc at top node
        b(1)=1.;
        %       modify for neumann at bottom node
        acn=keff(nc)/(dz(nc)*dz(nc));
        a(nz-1)  = 2.*acn;
        b(nz)    = -a(nz-1);
        %       generate sparse system matrix from diagonals a,b,and c
        A = spdiags([a' b' c'], -1:1, nz, nz);
        %       build RIGHT HAND SIDE
        %       nodal heat sources
        qd=c2n(hprod,dz);
        rc=c2n(rceff,dz);
        rhs=qd';
        %       modify for neumann at bottom node
        rhs(nz)= rhs(nz) + 2*acn*dz(nc)*qb/keff(nc);

        I= spdiags( ones(nz,1), 0, nz, nz);
        F= spdiags( 1./rc', 0, nz, nz);
        A=F*A; rhs=rhs./rc';
        L = I-dt(i)*theta(i)*A;
        %       right hand side
        r =(I+dt(i)*(1-theta(i))*A)*T(:,i)+dt(i)*rhs;
        %       modify for top (dirichlet) bc:
        r(1) = Tsi(i+1); L(1,1)=1.;
        %       solve by LU
        T(:,i+1) = (L\r);

        if iter>1;
            checktol=norm(abs(T(:,i+1)-Titer),inf);
            if strcmpi(out,'yes')==1,
                    disp([' maximal deviation at timestep ',...
                    num2str(i),' is ', num2str(checktol), ' K  (',...
                    num2str(iter),')']);end
            if checktol <= tol, break; end
            if iter >= maxiter, break; end
        end
        Titer=T(:,i+1);
        %    end NONLINEAR ITERATION
    end
    Tlast=T(:,i+1);
 

    % end TIME STEPPING
end
Tx=T(:,nt);







