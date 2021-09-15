function[L]=reg1d(k,flag)
% REG1D generates the regularization matrix L 
%
% [L]=reg1d(k,flag)
% Arguments:
%     k    = dimension ofd parameter space
%     flag =
%           'l0' or 'marq'   norm(m) = smallest m regularization
%           'l1' or 'grad'   norm(grad)= flattest model regularization
%           'l2' or 'lap'    norm(lap) = smoothest model regularization 
%
% V. R., June 2000 

if nargin < 2 flag=0; end

  N=k;
  
  e=ones(k,1);
  
  switch flag
        case {'l0','marq'}
           disp(' ... norm(m) = smallest m regularization used')
           Op = spdiags([e],0,k,k); 
        case {'l1','grad'}
             disp(' ... norm(grad)= flattest model regularization used')
             Op = spdiags([-1*e,e],0:1,k,k);
             Op(k,k)=0;   
        case {'l2','lap'}
           disp(' ... norm(lap)= smoothest model regularization used')
           Op = spdiags([e,-2*e,e],-1:1,k,k);
           Op(1,1)=0;Op(1,2)=0;
           Op(k,k-1)=0;Op(k,k)=0;
        otherwise
             disp(' Unknown method.')
           Op = zeros(k,k)
        end
L = Op;       
%  L = sparse(Op);       
%