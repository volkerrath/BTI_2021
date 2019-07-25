function [T,dT,Q,kbulk,ipor]=...
    heat1dnt(kl,kAl,kBl,hl,rl,cpl,rcl,porl,qb,...
        dz,ip,dt,it,GST,T0,theta,maxiter,tol,freeze,out)
%ADiMat BMFUNC $$=spdiags($1, $#) DIFFTO call(@spdiags, $@1, $#)
%
% HEAT1DNTP  solves nonlinear time-dependent heat equation
% for general 1-D grid with permafrost.
%
% T = heat1dnt(par,nu,ip,dz,dt,T0,GST,methods) calculates
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
% ip (1:nc)     	= pointer to assign parameters to cells
%                 	   	(nc gridsize in cells)
% dz(1:nc)      	= cell size (m)
% dt(1:nt-1)    	= time step (s)
% T0(1:nc+1)        	= initial temperatures
% GST(1:nt+1)        	= time-dependent boundary temperatures
%                     		at the top (e.g., paleoclimate)
% theta    	= methods(1:nt) is time stepping control parameter,
%                      		.5 =Crank-Nicholson, 1=Bacward Euler
% maxiter,tol         	= maximal number of nonlinear Picard iterations and
%                       	tolerance (max deviation from last iteration)
%
% Output:
% T(1:nc+1,1:nt )  = temperatures at given time steps
%
%
% V. R., July 20, 2005

rref=2650;
relax=1.;
Lh=333600;Tf=0;w=1.;
debug=1;

if nargin<17, freeze=1;  end
if nargin<18, out= 'no';end
if nargout>3,ipor=zeros(size(dz));end

ip=ip(:); dz=dz(:);z=[0 ;cumsum(dz)];nc=length(ip);nz=nc+1;
it=it(:); dt=dt(:);nt=length(dt)+1;

thetstep=theta*ones(1,nt-1);


k       =   kl(ip(1:  nc));     k=k(:);                      % thermal conductivity
kA      =   kAl(ip(1:  nc));   kA=kA(:);        % thermal conductivity coefficient A
kB      =   kBl(ip(1:  nc));   kB=kB(:);        % thermal conductivity coefficient B
rm       =   rl(ip(1:  nc));     rm=rm(:);         % rock density
cpm      =   cpl(ip(1:  nc));    cpm=cpm(:);       % rock cp
h       =   hl(ip(1:  nc));     h=h(:);         % heat production
por     =   porl(ip(1:  nc));   por=por(:);     % porosity


one=ones(size(ip));zero=zeros(size(ip));


%INITIAL VALUES

% initialize time, depth and pressure
t=[0; cumsum(dt)];
z=[0 ;cumsum(dz)];

%disp([mfilename '    spatial mesh: ',num2str(nz),' temporal mesh:',num2str(nz)]);

zc=0.5*(z(1:nz-1)+z(2:nz));
Pcl=[101325; 9.81*cumsum(dz.*rm)];Pcl=n2c(Pcl,dz);
Pch=[101325; 9.81*cumsum(dz.*998.)];Pch=n2c(Pch,dz);

GST=GST(:);T0=T0(:);
Tlast=T0(:);Tlast(1)=GST(1);
Titer=T0(:);Titer(1)=GST(1);

dc = 0.5 * (dz(2:nc)+dz(1:nc-1));

I = speye(nz,nz);

T=NaN*ones(nz,nt);
% nz
% nt
% whos T0

T(:,1)=T0(:);



% START TIME STEPPING
for itime = 1:nt-1

    %   start NONLINEAR ITERATION

    for  iter =1:maxiter
        %        build SYSTEM MATRIX
        %        initialize diagonals
        a(1:nz)=0.;b(1:nz)=0.;c(1:nz)=0.;
        %        define  coefficienGST for interior poinGST

      
            Titer=relax*Titer+(1-relax)*Tlast;
            Tc=n2c(Titer,dz);

% FLUID PROPS WITH HYDROSTATIC PRESSURE 
        rhof=rhofT(Tc,Pch);
        rcf=rhof.*cpfT(Tc,Pch);
        kf=kfT(Tc,Pch);
% ICE PROPS        
        rci=rhoiT(Tc).*cpiT(Tc);
        ki=kiT(Tc);
% MATRIX PROPS WITH LITHOSTATIC PRESSURE 
        rcm=rcmT(rm,cpm,Tc,Pcl);
        km=kmT(k,Tc,kA,kB,Pcl);

        %  permafrost
        if freeze==1,
            [gf dgf]=ftheta(Tc,Tf,w);
        else
            gf=one;dgf=zero;
        end

        porm=one-por;
        porf=por.*gf;
        pori=por-porf;

%         switch lower(mean)
%             case {'a' ,'ari', 'arithmetic'}
%                 keff  = porf.*kf +  pori.*ki + porm.*km;
%             case {'g','geo','geometric'}
                keff= exp(log(kf).*porf+log(ki).*pori+log(km).*porm);
%             case {'h','har','harmonic'}
%                 keff= 1./ (porf./kf +pori./ki + porm./km);
%             case {'s','sqr','sqrmean'}
%                 keff=(porf.*sqrt(kf)+pori.*sqrt(ki)+porm.*sqrt(km)).^2;
%             otherwise
%                 disp(['WMEAN: mode set to arithmetic, >', mean, '<  not defined'])
%                 keff  = porf.*kf +  pori.*ki +  porm.*km;
%         end
        rceff = porm.*rcm + ...
            pori.*rci +  ...
            porf.*rcf + por.*rhof.*Lh.*dgf;

        %         if nargout > 4 , k_eff(:,itime)=keff'; end
        %         if nargout > 5 , rc_eff(:,itime)=rceff'; end
        %         if nargout > 6 , ipor(:,itime)=(pori./por)'; end
        %         if nargout > 7 , lheat(:,itime)=(por.*rhof.*Lh.*dgf)'; end
        %         if nargout > 8 , rci(:,itime)=(rhoi.*cpi)'; end
        %         ipor(:,itime)=(pori./por)';

        dl = keff./dz;
        N(itime)=max((keff*dt(itime))./(rceff.*dz.^2));

        %        define matrix coefficienGST for interior poinGST
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
        rhs(nz)= rhs(nz) - 2*acn*dz(nc)*qb/keff(nc);
        
     
        F= spdiags( 1./rc', 0, nz, nz);
        A=F*A; rhs=rhs./rc';
        L = I-dt(itime)*thetstep(itime)*A;
        %       right hand side
        r =(I+dt(itime)*(1-thetstep(itime))*A)*Tlast+dt(itime)*rhs;
            
        %       modify for top (dirichlet) bc:
        r(1) = GST(it(itime)); 
        L(1,1)=1.;
        %       solve by LU
%         [row,col]=find(L==min(nonzeros(L)));
%         A=L(row,col);
%         [row,col]=find(L==max(nonzeros(L)));
%         B=L(row,col);
%         disp([num2str(A),'    ',num2str(B)])
%         [row,col]=find(r==min(nonzeros(r))) 
%         r(row,col)
%         [row,col]=find(r==max(nonzeros(r)))
%         r(row,col)
        Tnew = (L\r);
        
        checktol=norm(Tnew-Titer,inf);
                         
        Titer=Tnew;
        
        if debug > 1, disp([' maximal deviation at timestep ',...
                num2str(itime),' is ', num2str(checktol), ' K  (',...
                num2str(iter),')']);
        end
        
        if checktol <= tol || iter >= maxiter,
            break
        end



        %    end NONLINEAR ITERATION
    end
    
    Tlast=Tnew;
    
    % end TIME STEPPING
    if out ~=0,
        if mod(itime,out) ==0
           
            T(:,itime+1)=Tnew;
        end
    end
    % end TIME STEPPING
   
end

if out==0;T=Tlast;end

if nargout > 1
    dT=diff(T0);
    dT=dT(:)./dz(:);
    Q = keff(:).*dT(:);
    ipor=pori;
    kbulk=keff;
end
    






