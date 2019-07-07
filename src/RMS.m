function [rms]=RMS(res,err)
res=res(isfinite(res));res=res(:);nr=length(res);

if nargin==1, err=ones(size(res));end

if length(err)==1, 
    err=err*ones(size(res));
else
    err=err(isfinite(res));err=err(:);
end

res=res./err;
rms=sqrt(sum(res.*res)/(nr-1));

end