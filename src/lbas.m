function [X,rho,eta] = lbas(A,b,W,k,reorth)
%LBAS Lanczos bidiagonalization solver with augmented subspace
%
% [X,rho,eta] = lbas(A,b,W,k,reorth)
%
% A - coeffient matrix
% b - right-hand side
% W - a matrix with orthonormal columns (or a positive integer)
% k - the number of iterations
% reorth - 0: no reorthogonalization, 1: MGS ditto (default)
%
% X - matrix with k columns, each column holds an iteration vector
% rho - residual norms
% eta - solution norms
%
% If W is a positive integer p, then W has p columns corresponding
% to the polynomials of degree 0,1,...,p-1.
% Per Christian Hansen, DTU Compute and Kuniyoshi Abe, Gifu Shotoku
% Gakuen University, Sept. 4, 2015.
if nargin < 4, error('Too few input arguments'), end
if nargin < 5, reorth = 1; end
[m,n] = size(A);
[nW,p] = size(W);
if nW == 1 && p == 1 && W > 0 && ~rem(W,1); % Set W.
    p = W;
    W = ones(n,p);
    for i = 1:p-1
        W(:,i+1) = (1:n).ˆi;
    end
    W = orth(W);
elseif nW ~= size(A,2);
    error('No. rows in W must equal no. columns in A');
end
% Initialize.
X = zeros(n,k);
U = zeros(m,k+1);
V = zeros(n,k);
Bk = zeros(k+1,k);
Gk = zeros(k+1,p);
g = zeros(k+p+1,1);
% Prepare for iterations.
beta = norm(b);
u = b/beta;
U(:,1) = u;
v = zeros(n,1);
AW = A*W;
Gk(1,:) = u'*AW;
normb = beta;
g(1) = normb;
% Commence iterations.
for j = 1:k
    % Next (rightmost) column in lower bidiagonal part of matrix
    % via Lanczs bidiagonalization process.
    r = A'*u - beta*v;
    if reorth==1
        for i=1:j-1, r = r - (V(:,i)'*r)*V(:,i); end
    end
    alpha = norm(r); v = r/alpha;
    V(:,j) = v;
    Bk(j,j) = alpha;
    pp = A*v - alpha*u;
    if reorth==1
        for i=1:j, pp = pp - (U(:,i)'*pp)*U(:,i); end
    end
    beta = norm(pp); u = pp/beta;
    U(:,j+1) = u;
    Bk(j+1,j) = beta;
    % Apply stored orthog. transf. to new column of Bk.
    if j>1
        Bk(j-1:j,j) = [-si;conj(co)]*Bk(j,j);
    end
    % Determine new orthog. transf. to make Bk upper triangular.
    nu = norm(Bk(j:j+1,j));
    if nu==0 , error('Breakdown'), end
    co = Bk(j,j)/nu;
    si = -Bk(j+1,j)/nu;
    Bk(j,j) = co*Bk(j,j) - si*Bk(j+1,j);
    Bk(j+1,j) = 0;
    % Apply the orthog. transf. to updated G and to rhs.
    if j>1
        Gk(j,:) = saveG;
        g(j) = saveg;
        g(j+1:j+p) = 0;
    end
    Gk(j+1,:) = u'*AW;
    Gk(j:j+1,:) = [co,-si;si,conj(co)]*Gk(j:j+1,:);
    saveG = Gk(j+1,:); % To be used in next iteration.
    g(j:j+1) = [co,-si;si,conj(co)]*g(j:j+1);
    saveg = g(j+1);
    % Ditto.
    % Needed for bottom right block Fk.
    if j==1
        Utilde = ort(U(:,1:j+1),AW);
    else
        %Utilde = ort(U(:,j+1),Utilde);
        Utilde = ort(U(:,1:j+1),Utilde);
    end
    6
    % QR factorization of bottom right block.
    [qq,rr] = qr([Gk(j+1,:);Utilde'*AW]);
    Gk(j+1,:) = rr(1,:);
    Fk = rr(2:end,:);
    g(j+1:j+1+p) = qq'*g(j+1:j+1+p);
    % Compute solution.
    y = [Bk(1:j+1,1:j) Gk(1:j+1,:); zeros(p,j) Fk] \ g(1:j+p+1);
    X(:,j) = [V(:,1:j),W]*y;
    % svd([Bk(1:j+1,1:j) Gk(1:j+1,:); zeros(p,j) Fk])'
    disp(y')
end
if nargout > 1
    rho = sqrt(sum(abs(A*X-repmat(b,1,k)).ˆ2));
end
if nargout > 2
    eta = sqrt(sum(abs(X).ˆ2));
end
% Subfunction ==================================================
function Vw = ort(V,W)
% Orthonormalize W with respect to V (remove components along V).
% Henrik Garden & Per Chr. Hansen, DTU Compute, July 30, 2013.
k = size(V,2)-1;
p = size(W,2);
for s = 1:p
    w = W(:,s);
    for i = 1:k+s
        vi = V(:,i);
        w = w-vi'*w*vi;
    end
    w = w/norm(w);
    V(:,k+s+1) = w;
end
Vw = V(:,end-p+1:end);
