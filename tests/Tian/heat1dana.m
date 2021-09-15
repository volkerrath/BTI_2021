function [T,Q]=heat1dana(parvec,time,kappa,k,z,tlog,refyr,out)  
%            [T]=Heat1Dana(parvec,gtime,tlog,kappa,k,z,refyr,nz,out); 
% time:     times at which the GTH is evaluated
% k:   average thermal conductivity of the log
% kappa:    diffusivity
% tlog:    year the log was recorded
% refyr:   Reference year for time scale 
% nstep:     Number of time steps in the model
%
nprof=1;nz=length(z);nt=length(time);
% setup matrix for forward modeling
[M,nmod] = msetup(time,nprof,kappa,k,z,nz,tlog,refyr,nt,out);
% solve forward problem
Temp=M*parvec;
T.val=Temp;T.z=z;
if nargout>1,
   T.grd=diff(Temp)./diff(z);
   Q.z=0.5*(z(1:nz-1)+z(2:nz));
end
if nargout>2,
   Q.val=k*T.grd;
end

