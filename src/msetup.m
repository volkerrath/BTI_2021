function [a,nmod] = ...
    msetup(time,nprof,diffu,conduc,zdat,ndat,tlog,refyr,nu,out)
% setting up matrix for inversion.
% This function will be called from the heat1dat
% and sets up the model matrix for the inversion.
% Method  is described in:
% Beltrami, H. & Mareschal, J.-C. 
% Recent warming in eastern Canada inferred from Geothermal Measurements 
% Geophys. Res. Lett., 1991, 18(4), 605-608

% Original code by J.C. Mareschal
% Last modified:
% A. Hartmann, 3-May-01, created a separate function
% V. Rath, 8-February-04, vectorized
% V. Rath, 21-June-07
% V. Rath,  1-Nov-19

% t0=cputime;

nmod=length(time);
alldat=sum(ndat);allpar=nu+2*nprof;a=zeros(alldat,allpar);
maxdat=0;
% disp(strcat(['Mstetup: ',num2str([nmod ndat alldat allpar])]));
for iprof=1:nprof
    % begin loop on temperature profiles
    mindat=maxdat+1;
    maxdat=mindat+ndat(iprof)-1;
    nmod=nmod+2;
    shift=refyr-tlog(iprof);
    active=find(time > shift);na=length(active);
    tinv(active)= .5./sqrt(diffu(iprof)*(time(active)-shift));
    for idata=mindat:maxdat
        % ...................number of steps after logging date
        a(idata,active)=erfc(zdat(idata)*tinv(active));
        a(idata,active(2:na))=a(idata,active(2:na))-a(idata,active(1:na-1));
        a(idata,nmod-1)=1.;
        a(idata,nmod)=zdat(idata)/conduc(iprof);
    end
%     if out>1, disp(['... row ',num2str(iprof),' finished, ', ...
%             num2str(cputime-t0,'%0.5g'),' seconds']); end
    % finish setting up matrix: end loop on temperature profiles
end
% if out>0, disp([' finished, ', num2str(cputime-t0,'%0.5g'),' seconds']); end