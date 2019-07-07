function qz=hflux(A)

%====================================================================== 
%    calculate z heat flux at cell center 
%
%     input:
%       i,j,k                               grid indices    
%     output:
%       qzc                                 (W/m^2)
%
%     last change:  vr    nov 26,2004
%======================================================================

% hier: input nur k (z-Richtung)
% A ist Matrix mit Spalten: k-Index, temp, wlz, delz

    k=A(:,1);
    temp=A(:,2);
    wlz=A(:,3);
    delz=A(:,4);
    qz=0.d0;
    k0=length(k);
    
    for k=1:k0-1
%       lk=0.d0
       if (k0>1) & (k<=k0-1); 
           f1=wlz(k);  
           f2=wlz(k+1);
           prod  = f1*f2;
           summ  = f1*delz(k+1)+f2*delz(k);
           if (summ>0.d0) lk(k) = 2.d0*prod/summ;
           end 
       
       end

           
    end

    for k=1:k0-1

       if (k>1) & (k<=k0) 
          d1 = temp(k+1) - temp(k);
          d2 = temp(k) - temp(k-1);
      	  qz(k) = 0.5d0*(-lk(k)*d1-lk(k-1)*d2);
       end
       if (k==1)  	     
          qz(k)=-lk(1)*(temp(2)-temp(1));
       end
       if (k==k0)  	  
          qz(k)=-lk(k0-1)*(temp(k0)-temp(k0-1));
       end

    end 