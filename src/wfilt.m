function [y]= wfilt(x,M,N,fmode,debug)
% WFILT(X,M,N,MODE,DEBUG)  runs a weighted average filter of length 2*M+1
% N times on the input time series X. MODE can currently be either
% 'b' (boxcar, default) or 't' (triangular).
% V. R., modified June 5, 2016
if nargin < 5, debug =   0; end
if nargin < 4, fmode  = {'b' 'mir'}; end
if nargin < 3, N=1; end

x=x(:);nx=length(x);

fil  =char(fmode(1));
pad  =char(fmode(2));

if rem(M,2)==0, 
M=M+1; 
disp([' filter length M even, M+1 assumed']);
end
L=(M-1)/2;

%DEFINE FILTER WEIGHTS
switch lower(fil)
    case{'t' 'tri' 'triangular'}
        w = tri(2*L+1);
    case{'b' 'box' 'boxcar'}
        w = box(2*L+1);
    otherwise
     disp(['filter ',lower(fil),' not implemented, box assumed']);
      w = box(2*L+1);
end


% PAD
switch lower(pad)
    case{'a' 'avg' 'average'}
        x1=mean(x(1:L))    *ones(L,1);
        x2=mean(x(nx-L+1:nx))*ones(L,1);
    case{'m' 'mir' 'mirror'}
        x1=flipud(x(1:L));x2=flipud(x(nx-L+1:nx));
    otherwise
     disp(['padding ',lower(pad),' not implemented, avg assumed']);
        x1=mean(x(1:L))    *ones(L,1);
        x2=mean(x(nx-L+1:nx))*ones(L,1);
end
xa=[x1;x;x2];
nxa=length(xa);
xb=xa;


% RUN FILTER
for n=1:N
    
    for k = L+1:nxa-L
        aux = xa(k-L:k+L);
        xb(k) = sum(w.*aux);
    end;
    xa=xb;
end

y=xa(L+1:nxa-L);

if debug,
    figure;
    plot(x, 'LineWidth',2,'Color','r');hold on;
    if L ~= 0, plot(y, 'LineWidth',2,'Color','b','LineStyle','--');end
    grid on;
    xlabel('samples','FontSize',14);ylabel('x','FontSize',14);
    if L ~= 0, grid on;legend('original','filtered','Location','NorthWest');end
    title(['filter:',smode,'/',num2str(M)],'FontSize',14)
    %xlim([0 20])
end






