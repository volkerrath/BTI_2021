function [Theta,dTheta]=ftheta(T,freeze)
% [theta,dtheta]=ftheta(T,freeze) calculates the fluid/ice
% partition function used for the apparent heat capacity approach
% to phase change.
%
% On Input:
% T                 =     temperature (Celsius)
% freeze.fun        =     function type used (lun/gal/nic/hin)
% freeze par        =     paramters passsed to freeze.fun
%
% If freeze.fun  = 'lun':
% Tf                =     phase boundary temperature
% w                 =     scaling factor for "mushy" region
%
%
% If freeze.fun  = 'gal':
% Tf                =     phase boundary temperature
%
% If freeze.fun  = 'nic':
% Tf                =     phase boundary temperature
% p                 =     characterisitic exponent
%
% If freeze.fun  = 'hin':
% Tb          =     freezing interval
%                         Tb(1) = T_liquidus (-0.001)
%                         Tb(2) = T_solidus  (-2)
% b                 =      scaling coefficient
%
% On Ouput:
% Theta            =    value of partition function (0 <= Theta <=1)
% dTheta            =   derivative of partition function
%                       with respect to temperature
%
% REFERENCES:
% Lunardini, V. J.
%       Heat Transfer in Cold Climates
%       Litton Educational Publishing, 1981
% Galushkin, Y.
%       Numerical simulation of permafrost evolution as a part of
%       sedimentary basin modeling: permafrost in the
%       Pliocene-Holocene climate history of the Urengoy field
%       in the West Siberian basin
%       Canadian Journal of Earth Sciences, 1997, 34, 935-948
% Nicolsky, D. J.; Romanovsky, V. E. & Tipenko, G.
%       Estimation of thermal properties of saturated soils using
%       in-situ temperature measurements
%       The Cryosphere, 2007, 1, 41-58
% Hinzman, L.; Goering, D. J. & Kane, D.
%       A distributed thermal model for calculating soil temperature
%       profiles and depth of thaw in permafrost regions
%       Jour. Geophys. Res. 1998, 103, D22, 28975-28991
%
% vr  Nov 2010



Theta=ones(size(T));if nargout>1, dTheta=zeros(size(T)); end



switch lower(freeze.fun)

    case {'lun' 'l'}
        if isempty(freeze.par),Tf=0;w=.5;
        else
            Tf=freeze.par(1);w=freeze.par(2);
        end
        
        
        Theta=exp(-((T-Tf)./w).^2);
        b=find(T>Tf);Theta(b)=1.;

        if nargout>1,
            a=-2/(w*w);
            dTheta=a*(T-Tf).*Theta;
            dTheta(b)=0.;
        end

    case {'gal' 'g'}
        if isempty(freeze.par),Tf=0;
        else
            Tf=freeze.par(1);
        end

        a0=   1.0;
        a1=   0.60152823763179;
        a2=   0.23218232347212;
        a3=   0.04669792788297;
        a4=   0.00535597924776;
        a5=   0.00036415588418;
        a6=   0.00001450956751;
        a7=   0.00000031279149;
        a8=   0.00000000281502;
%
        Theta=a0+a1*T+a2*T.^2+a3*T.^3+a4*T.^4+a5*T.^5+a6*T.^6 ...
            +  a7*T.^7 + a8*T.^8;
        b=find(T>Tf);Theta(b)=1.;

        if nargout>1,
            dTheta=a1+2*a2*T+3*a3*T.^2+4*a4*T.^3+5*a5*T.^4 ...
                +6*a6*T.^5 +  7*a7*T.^6 + 8*a8*T.^7;
            dTheta(b)=0.;
        end

    case {'nic' 'n'}
        if isempty(freeze.par),Tf=-0.001;p=0.5;
        else
            Tf=freeze.par(1);p=freeze.par(2);
        end
        if p >= 0.5 && p <= 0.8 && Tf < 0,
            b=find(T<Tf);Theta=ones(size(T));
            Theta(b)=abs(Tf)^p * abs(T(b)).^(-p) ;

            if nargout>1,
                dTheta=zeros(size(T));
                dTheta(b)=abs(Tf)^p*(-p)*abs(T(b)).^(-p-1).*sign(T(b));
            end
        else
            error ('ftheta_nicolsky: incorrect values for parameters p/Tf ')
        end


     case {'hin' 'h'}
%
         if isempty(freeze.par),
             Tb=[-0.001 -2];p=10;
         else
             Tb=freeze.par(1:2);p=freeze.par(3);
         end
         b1=find(T>=Tb(1));
         b2=find(T<=Tb(2));
         b3=find(T>Tb(2) & T<Tb(1));

         f=p./(Tb(1)-Tb(2));
         Theta(b3)=(exp(-f*T(b3))-exp(-f*Tb(2)))./(exp(-f*Tb(1))-exp(-f*Tb(2)));
         Theta(b1)=1.;Theta(b2)=0.;


         if nargout>1,
             fp=f*p;
             dTheta=zeros(size(T));

         end

     end

