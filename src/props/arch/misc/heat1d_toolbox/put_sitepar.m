function [model]= put_sitepar(k,kA,kB,h,p,c,r,rc,z,ip,qb,gt,site)
% stores model parameters from model structure
model.k=k;
model.kA=kA;
model.kB=kB;
model.h=h;
model.p=p;
model.c=c;
model.r=r;
model.rc=rc;
model.qb=qb;
model.gt=gt;
model.ip=ip;
model.z=z;
model.site=site;