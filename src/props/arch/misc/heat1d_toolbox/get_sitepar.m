function [k,kA,kB,h,p,c,r,rc,ip,z,qb,gt,site]= get_sitepar(model)
% extracts model parameters from model structure
k=model.k;
h=model.h;
kA=model.kA;kB=model.kB;
p=model.p;
c=model.c;r=model.r;rc=model.rc;
qb=model.qb;
gt=model.gt;
ip=model.ip;
z=model.z;
site=model.site;