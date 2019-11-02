function [T,Q]=heat1dat(parvec,time,kappa,k,z,tlog,refyr,out)  
% [T]=heat1dat(parvec,gtime,tlog,kappa,k,z,refyr,nz,out); 
% time:     times at which the GTH is evaluated
% k:   average thermal conductivity of the log
% kappa:    diffusivity
% tlog:    year the log was recorded
% refyr:   Reference year for time scale 
% nstep:     Number of time steps in the model
% Method  is described in:
% Beltrami, H. & Mareschal, J.-C. 
% Recent warming in eastern Canada inferred from Geothermal Measurements 
% Geophys. Res. Lett., 1991, 18(4), 605-608

nprof=1;nz=length(z);nt=length(time);
% setup matrix for forward modeling
[M,nmod] = msetup(time,nprof,kappa,k,z,nz,tlog,refyr,nt,out);
% solve forward problem
Temp=M*parvec;
T.val=Temp;T.z=z;T.grd=diff(Temp)./diff(z);
T.zm=0.5*(z(1:nz-1)+z(2:nz));

if nargout>1,
   Q.val=k*T.grd;
   Q.z=T.zm;
end

