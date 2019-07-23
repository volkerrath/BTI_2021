function[L,W]=L_1D(k,flag)
% REG1D generates the regularization matrix L
%
% [L,W]=reg1d(k,flag)
% Arguments:
%     k    = dimension ofd parameter space
%     flag =
%           'l0' or 'marq'   norm(m) = smallest m regularization
%           'l1' or 'grad'   norm(grad)= flattest model regularization
%           'l2' or 'lap'    norm(lap) = smoothest model regularization
%
% V. R., June 2000

if nargin < 2 flag=0; end


e=ones(k,1);

switch lower(flag)
    case {'l0','marq'}
%        disp(' ... norm(m) = smallest m regularization used')
        L = spdiags([e],0,k,k);
        if nargout > 1, W=zeros(k,0); end
    case {'l1','grad'}
%        disp(' ... norm(grad)= flattest model regularization used')
        L = spdiags([-1*e,e],0:1,k-1,k);
        if nargout > 1, 
             W(:,1) = ones(k,1);
             W(:,1) = W(:,1)/norm(W(:,1));
        end 
    case {'l1x','gradx'}
%        disp(' ... norm(grad)= flattest model regularization used')
%        disp(' ... quadratic form incl. Neumann BC')
        L = spdiags([-1*e,e],0:1,k,k);
        if nargout > 1, 
             W(:,1) = ones(k,1);
             W(:,1) = W(:,1)/norm(W(:,1));
        end 
    case {'l2','lap'}
%        disp(' ... norm(lap)= smoothest model regularization used')
        L = spdiags([e,-2*e,e],0:2,k-2,k);
        if nargout > 1, 
          W(:,1) = ones(k,1); 
          W(:,2) = W(:,1).*[1:k]'; W(:,1) = W(:,1)/norm(W(:,1)); 
          W(:,2) = W(:,2) - W(:,1)*(W(:,1)'*W(:,2)); W(:,2) = W(:,2)/norm(W(:,2)); 
        end
    case {'l2x','lapx'}
%        disp(' ... norm(lap)= smoothest model regularization used')
%        disp(' ... quadratic form incl. Neumann BC')
        L = spdiags([e,-2*e,e],0:2,k,k);
        if nargout > 1, 
          W(:,1) = ones(k,1); 
          W(:,2) = W(:,1).*[1:k]'; W(:,1) = W(:,1)/norm(W(:,1)); 
          W(:,2) = W(:,2) - W(:,1)*(W(:,1)'*W(:,2)); W(:,2) = W(:,2)/norm(W(:,2)); 
        end
   otherwise
        disp(' Unknown method.')
        L = zeros(k,k);
end
