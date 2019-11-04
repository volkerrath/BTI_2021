function [M] = msetup(time,diffu,conduc,zdat,tlog,refyr,out)
% setting up matrix for inversion.
% This function will be called from the heat1dat
% and sets up the model matrix for the inversion.
% 
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
%
% Original code by J.C. Mareschal
% Last modified:
% A. Hartmann, 3-May-01, created a separate function
% V. Rath, 8-February-04, vectorized
% V. Rath, 21-June-07
% V. Rath,  1-Nov-19

% t0=cputime;

nt=length(time);npar=nt+2;
nz = length(zdat);
M=zeros(nz,npar);

% disp(strcat(['Mstetup: ',num2str([nt ndat alldat allpar])]));
% begin loop on temperature profiles
mindat=maxdat+1;
maxdat=mindat+nz-1;
shift=refyr-tlog;
active=find(time > shift);na=length(active);
tinv(active)= .5./sqrt(diffu*(time(active)-shift));
for idata=mindat:maxdat
    % ...................number of steps after logging date
    M(idata,active)=erfc(zdat(idata)*tinv(active));
    M(idata,active(2:na))=M(idata,active(2:na))-M(idata,active(1:na-1));
    M(idata,npar-1)=1.;
    M(idata,npar)=zdat(idata)/conduc;
end
%     if out>1, disp(['... row ',num2str(iprof),' finished, ', ...
%             num2str(cputime-t0,'%0.5g'),' seconds']); end
% finish setting up matrix: end loop on temperature profiles
% if out>0, disp([' finished, ', num2str(cputime-t0,'%0.5g'),' seconds']); end