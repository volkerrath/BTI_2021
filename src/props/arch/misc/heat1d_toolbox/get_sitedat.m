function [Tobs,id,zd,cov,err,site]= get_sitedat(data)
% extracts site data from data structure
Tobs=data.Tobs';
err=data.err';
id=data.id;
cov=data.cov';
zd=data.z';
site=data.site;