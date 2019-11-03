function [T,Q]=heat1dat(k,r,c,Qb,z,t,gst,T0,tlog,refyr,out) 
% time:     times at which the GTH is evaluated
% k:   average thermal conductivity of the log
% kappa:    diffusivity
% tlog:    year the log was recorded
% refyr:   Reference year for time scale 
% nstep:     Number of time steps in the model
% Method  is described in:
% 
% Beltrami, H. & Mareschal, J.-C. 
%   Recent warming in eastern Canada inferred from Geothermal Measurements 
%   Geophys. Res. Lett., 1991, 18(4), 605-608
% Beltrami, H.; Jessop, A. M. & Mareschal, J.-C. 
%   Ground temperature histories in eastern and central Canada from 
%   geothermal measurements: evidence of climate change 
%   Palaeogeography, Palaeoclimatology, Palaeoecology, 1992, 98, 167-184
% Mareschal, J.-C. & Beltrami, H. 
%   Evidence for recent warming from perturbed thermal gradients: examples 
%   from eastern Canada Clim. Dyn., 1992, 6, 135-143

nz=length(z);nt=length(time);

gtime=-flipud(t);
gsth=flipud(gst);

kappa = k/(r*c);
parvec=[gsth;T0;Qb];

% setup matrix for forward modeling
[M,nmod] = msetup(t,kappa,k,z,nz,tlog,refyr,nt,out);
% solve forward problem
Temp=M*parvec;
T.val=Temp;T.z=z;T.grd=diff(Temp)./diff(z);
T.zm=0.5*(z(1:nz-1)+z(2:nz));

if nargout>1
   Q.val=k*T.grd;
   Q.z=T.zm;
end

