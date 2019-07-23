function [err]=put_props(file,unit,text)
% PUT_PROPS writes geological units and associated physical parameter 
% values to file
%
% Arguments:
% file                  filename for temperature log
% text                  header             
% unit                  structure associated with units                  
%
% unit.number           unit number
% unit.name             unit name
% unit.k                thermal conductivity
% unit.kA,kB            coefficients of temperature dependence of condictivity
% unit.h                heat production
% unit.r                density
% unit.c                heat capacity
% unit.dk,dh,dr,dc         errors
%
%% V. R., Sept. 29, 2002 

if nargin<3, text=strcat([filename ' written on ' date]); end

err=0;
fid=fopen(file,'w');

fprintf(fid,'%% %s\n',text);

units=length(unit);
for i=1:units 
 number=unit(i).number;
 k=unit(i).k;dk=unit(i).dk; kA=unit(i).kA; kB=unit(i).kB;
 h=unit(i).h;dh=unit(i).dh;
 r=unit(i).r;dr=unit(i).dr;
 c=unit(i).c;dc=unit(i).dc;
 d=unit(i).d;dd=unit(i).dd;
 p=unit(i).p;dp=unit(i).dp;
name=unit(i).name;
 fprintf(fid,'%5i %10g %10g %10g %10g  %10g %10g  %10g %10g  %10g %10g  %10g %10g  %10g %10g \t "%s" \n',...     
         number,k,dk,kA,kB,h,dh,r,dr,c,dc,d,dd,p,dp,name);

end                                                                         
fclose(fid);
