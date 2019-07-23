function [unit]=get_props(file)
% GET_PROPS reads geological units and associated physical parameter 
% values from file
%
% Arguments:
% file                  filename for temperature log
% 
% unit                  structure associated with units                  
%
% unit.name             unit name
% unit.number           unit number
% unit.k                thermal conductivity
% unit.kA,kB            coefficients of temperature dependence of condictivity
% unit.h                heat production
% unit.r                density
% unit.c                heat capacity
% unit.d                diffusivity
% unit.p                porosity
% unit.dk, .dh, .dr,   errors
%     .dc, .dd, .dp

%% V. R., May. 6, 2003 
header=0;

% read from file
[number,k,dk,kA,kB,h,dh,r,dr,c,dc,d,dd,p,dp,name]=...
        textread(file,'%d   %f %f %f %f   %f %f %f %f %f %f %f %f  %f %f %q',...
        'commentstyle','matlab','headerlines',header);

% store into structure unit
units=length(number);
for i=1:units
 unit(i).name=name{i};
 unit(i).number=number(i);
 unit(i).k=k(i); unit(i).dk=dk(i);unit(i).kA=kA(i);unit(i).kB=kB(i);
 unit(i).h=h(i); unit(i).dh=dh(i);
 unit(i).r=r(i); unit(i).dr=dr(i);
 unit(i).c=c(i); unit(i).dc=dc(i);
 unit(i).d=d(i); unit(i).dd=dd(i);
 unit(i).p=p(i); unit(i).dp=dp(i);

end

