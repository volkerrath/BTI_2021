function [model]=get_model(file,zin)
% GET_MODEL constructs model parameter index ip
% from file with layer boundaries
%
% Arguments:
% file        filename for temperature log
% z           depths for model nodes
%
% model                  structure associated with model            
%
% model.ip               unit numbers associated to mesh cells
% model.z                depth of mesh points
% model.is               unit number associated with depths zs 
% model.zs               depth of unit boundaries
% model.gt, .dgt         surface temperature and error
% model.qb, .dqb         basal heat flow density and error
% 
% V. R., Oct. 29, 2002 

header=1;
[gt dgt qb dqb] = textread(file,'%f %f %f %f ',header);
[zs,is]         = textread(file,'%f %f %*[^\n]','commentstyle','matlab','headerlines',header);


if nargin <2,
    zstart=10;zend=zs(length(zs));
    z=[0,logspace(log10(zstart),log10(zend),201)];
else
    z=zin;
end

nz=length(z);
lower=z(1);
for i=1:length(zs)
    upper=min([zs(i),z(nz)]);    
    ip(find(z<upper&z>=lower))=is(i);
    lower=upper;
end
ip(length(ip):nz-1)=is(length(zs));


% store into structure model 
model.ip=ip;model.z=z;model.is=is;model.zs=zs;model.gt=gt;model.dgt=dgt;model.qb=qb;model.dqb=dqb;