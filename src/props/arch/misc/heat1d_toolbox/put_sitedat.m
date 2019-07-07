function [data]= put_sitedat(Tobs,id,zd,cov,err,site)
% extracts site data from data structure
data.Tobs=Tobs;
data.id=id;
data.cov=cov;
data.err=err;
data.z=zd;
data.n=length(zd);
data.site=site;