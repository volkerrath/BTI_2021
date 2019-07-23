function [T,z,t,N,k_eff,rc_eff,ipor,lheat,rci]=...
    heat1dnt(kl,kAl,kBl,hl,rhoml,cpml,porl,qb,...
    index,dz,dt,T0,Ts,theta,maxiter,tol,freeze,out)

% HEAT1DNTP  solves nonlinear time-dependent heat equation
% for general 1-D grid with permafrost.
%
% T = heat1dnt(par,nu,index,dz,dt,T0,Ts,methods) calculates
% numerically (FD) temperatures for a given set of timesteps,
% given a model for thermal conductivity, heat production and rho*c.
% transient heat conduction with time-dependent source term is assumed.
% thermal conductivity and heat capacity are assumed as nonlinear
% functions of temperature;
% Input :
% k (1:nu)              = lambda matrix
% kA, kB(1:nu)	 	= thermal conductivity coefficient A and B
% hl(1:nu)		= volumetric heat production matrix
% cml (1:nu)		= heat capacity c_p matrix
% rhoml 		= density c_p matrix
% porl			= porosity
% qb			= basal heat flow
% index (1:nc)     	= pointer to assign parameters to cells
%                 	   	(nc gridsize in cells)
% dz(1:nc)      	= cell size (m)
% dt(1:nt-1)    	= time step (s)
% T0(1:nc+1)        	= initial temperatures
% Ts(1:nt+1)        	= time-dependent boundary temperatures
%                     		at the top (e.g., paleoclimate)
% theta(1:nt-1)      	= methods(1:nt) is time stepping control parameter,
%                      		.5 =Crank-Nicholson, 1=Bacward Euler
% maxiter,tol         	= maximal number of nonlinear Picard iterations and
%                       	tolerance (max deviation from last iteration)
%
% Output:
% T(1:nc+1,1:nt )  = temperatures at given time steps
%
%
% V. R., July 20, 2005

solver='direct';
% solver='bicgstab',tolit=1.e-12;lutol=1.e-6;maxit=20;lufreq=5;restart=10;
% solver='gmres',tolit=1.e-12;lutol=1.e-6;maxit=20;lufreq=5;restart=10;
mean='g';relax=1.;
Lh=333600;Tf=0;w=1.;


if nargin<17, freeze='yes';  end
if nargin<18, out= 'no';end


[n1 n2]=size(index); if n1==1, index=index';end
[n1 n2]=size(dz);    if n1==1, dz=dz';end
[n1 n2]=size(dt);    if n1==1, dt=dt';end

nt=length(dt)+1;
nc=length(index);nz=nc+1;

k=kl(index(1:  nc));          % thermal conductivity
[n1 n2]=size(k);    if n1==1, k=k';end
kA =    kAl(index(1:  nc));   % thermal conductivity coefficient A
[n1 n2]=size(kA);   if n1==1, kA=kA';end
kB =    kBl(index(1:  nc));   % thermal conductivity coefficient B
[n1 n2]=size(kB);   if n1==1, kB=kB';end
h  =     hl(index(1:  nc));   % heat production
[n1 n2]=size(h);    if n1==1, h=h';end
rhom  = rhoml(index(1:  nc));   % density
[n1 n2]=size(rhom); if n1==1, rhom=rhom';end
cpm0   =  cpml(index(1:  nc));   % heat capacity
[n1 n2]=size(cpm0);  if n1==1, cpm0=cpm0';end
por   =  porl(index(1:  nc));   % pcorosity
[n1 n2]=size(por);  if n1==1, por=por';end
%
one=ones(size(index));
zero=zeros(size(index));
%INITIAL VALUES

% initialize time, depth and pressure
t=[0; cumsum(dt)];
z=[0 ;cumsum(dz)];zc=0.5*(z(1:nz-1)+z(2:nz));Pc=998.*9.81*zc;

Tlast=T0;Titer=T0;Tlast(1)=Ts(1);Titer(1)=Ts(1);

dc = 0.5 * (dz(2:nc)+dz(1:nc-1));

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
            Tc=(n2c(Tlast,dz));
        else
            Tc=(n2c(relax*Titer+(1-relax)*Tlast,dz));
        end

        Pc=9.81*rhofT(Tc,Pc).*zc;
        rhoi=rhoiT(Tc);    ki=kiT(Tc);  cpi=cpiT(Tc);
        rhof=rhofT(Tc,Pc); kf=kfT(Tc);  cpf=cpfT(Tc,Pc);
        if any(cpm0<0)
            cpm = -cpm0;
        else
            cpm=cpmT(cpm0,Tc);
        end
        km=kmT(k,Tc,kA,kB);

        %  permafrost
        if strcmpi(freeze,'yes')==1,
            [gf dgf]=ftheta(Tc,Tf,w);
            %            gf=gf';dgf=dgf';
        else
            gf=one;dgf=zero;
        end

        porm=one-por;
        porf=por.*gf;
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

        if nargout > 4 , k_eff(:,i)=keff'; end
        if nargout > 5 , rc_eff(:,i)=rceff'; end
        if nargout > 6 , ipor(:,i)=(pori./por)'; end
        if nargout > 7 , lheat(:,i)=(por.*rhof.*Lh.*dgf)'; end
        if nargout > 8 , rci(:,i)=(rhoi.*cpi)'; end


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
        qd=c2n(h,dz);
        rc=c2n(rceff,dz);
        rhs=qd';
        %       modify for neumann at bottom node
        rhs(nz)= rhs(nz) + 2*acn*dz(nc)*qb/keff(nc);

        F= spdiags( 1./rc', 0, nz, nz);
        A=F*A; rhs=rhs./rc';
        L = I-dt(i)*theta(i)*A;
        %       right hand side
        r =(I+dt(i)*(1-theta(i))*A)*Tlast+dt(i)*rhs;
        %       modify for top (dirichlet) bc:
        r(1) = Ts(i+1); L(1,1)=1.;
        %       solve by LU or iterative methods
        switch lower(solver)
            case {'d','direct'}
                Tnew = (L\r);
            case {'b','bicgstab'}
                if  (mod(iter,lufreq)==1), [LA UA]=luinc(A,lutol);end
                [Tn,flag,relres,iter] = bicgstab(A,rhs,tol,maxit,LA,UA,Titer);
                Tnew = Tn;
            case {'g','gmres'}
                if  (mod(iter,lufreq)==1), [LA UA]=luinc(A,lutol);end
                [Tn,flag,relres,iter] = gmres(A,rhs,tol,maxit,LA,UA,Titer);
                Tnew = Tn;
            otherwise
                Tnew = (L\r);
        end

        if iter>1;
            checktol=norm(abs(Tnew-Titer),inf);
            if strcmpi(lower(out),'yes') ~=0, disp([' maximal deviation of FPI at timestep ',...
                    num2str(i),' is ', num2str(checktol), ' K  (',...
                    num2str(iter),')']);end
            if checktol <= tol, break; end
            if iter >= maxiter, break; end
        end
        Titer=Tnew;
        %    end NONLINEAR ITERATION
    end
    Tlast=Tnew;
    % end TIME STEPPING
    if out ~=0, T(:,i+1)=Tnew; end
    % end TIME STEPPING
end

if out==0;
    T=Tnew;
else
    T(:,1)=T0;
end

if nargout > 4 , k_eff=keff'; end
if nargout > 5 , rc_eff=rceff'; end
if nargout > 6 , ipor=(pori./por)'; end
if nargout > 7 , lheat=(por.*rhof.*Lh.*dgf)'; end
if nargout > 8 , rci=(rhoi.*cpi)'; end








