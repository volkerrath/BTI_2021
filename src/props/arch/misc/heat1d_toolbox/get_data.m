function [data]=get_data(file,z,header)
% GET_DATA constructs inversion data from temperature logs
% by interpolation to model depths fiven in z
%
% Arguments:
% file        filename for temperature log
% z           depths for model nodes
%
% data                   structure associated with observed data          
% 
% data.t                 temperature
% data.id                points to mesh nodes associated with observation
% data.err               mean error of observation
%
% V. R., Oct. 30, 2002 

if nargin <3 header=0; end

nz=length(z);

[zd,td,err,flag]  = textread(file,'%f%f%f%f%*[^\n]',...
          'commentstyle','matlab','headerlines',header);

good=find(flag==0);zdf=zd(good);tdf=td(good);errf=err(good);
t=interp1(zdf,tdf,z,'linear');
e=interp1(zdf,errf,z,'linear');
id=find(z > min(zdf) & z < max(zdf));nd=length(id);
% store into structure data 
data.T=t;data.id=id;data.nd=nd;data.z=z;data.Err=e;
data.Cov=e.^2;
 

