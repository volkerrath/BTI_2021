function[W,Wi]=we1d(m,x,mode,eps,bet)
% [W,Wi]=function we1d(m,x,mode,eps,bet)
% caclculates the diagonal parameter weighting matrix
% W and its inverse, Wi.
% =>     m           parameter vector
%        x           patial coordinate vector
%        mode        Penalty type:
%                       'SM'    =  maxsmoothness
%                       'TV'    =  min total deviation 
%                       'ME'    =  min entropy
%                       'MS'    =  min support
%                       'MGS'   = min gradient support
% <=     W          sparse weight diagonal matrix
%        Wi         inv(W) diagonal matrix
% see:
% Zhdanov, M. S. (2002): Geophysical inverse theory and 
%       regularization problems, p. 155ff 

% defaults
if nargin<5, bet=1.e-6;end 
if nargin<4, eps=1.e-6;end 
if nargin<3, mode='sm';end 

M=length(m);
e=ones(M,1);

[Mx Nx]=size(x);
if       Mx==M, 
    dx=[diff(x)' 1]'; e=e./dx;
elseif  Nx==M,
    dx=[diff(x) 1]'; e=e./dx;
else
    disp([' dimension of m (' num2str(M) ...
            ') not konsistent with size of x (' num2str(Mx) num2str(Nx) ')'])
end

disp(['... calculating <' lower(mode) '> model weights'])
    
    switch lower(mode)
        case {'sm'}
            % maximum smoothness 
            eps2=eps^2;
            g = spdiags([-1*e,e],0:1,M,M);g(M,M)=0;
            gm=g*m;
            wgt=gm./sqrt(m.*m+eps2);
        case {'me'}
            % min entropy 
            eps2=eps^2;q=sum(abs(m));
            wgt=sqrt(abs(m).*log(q./abs(m))/q)...
                ./sqrt(m.*m+eps2);
        case {'tv'}
            % min total variation 
            eps2=eps^2;bet2=bet^2;
            g = spdiags([-1*e,e],0:1,M,M);g(M,M)=0;
            gm=g*m;
            wgt=sqrt(sqrt(abs(gm).^2+bet2))./sqrt(m.*m+eps2);
        case {'ms'}
            % min support
            eps2=eps^2;
            wgt=1./sqrt(m.*m+eps2);
        case {'mgs'}
            % min gradient support
            eps2=eps^2;bet2=bet^2;
            g = spdiags([-1*e,e],0:1,M,M);g(M,M)=0;
            gm=g*m;
            wgt=gm./(sqrt(gm.*gm+eps2).*sqrt(m.*m+bet2));
        otherwise
            disp(' Unknown method. Unit matrix used')
            wgt = e;
    end
    
    W=spdiags(wgt,[0],M,M);
    if nargout > 1, Wi=spdiags(1./wgt,[0],M,M); end
    
    