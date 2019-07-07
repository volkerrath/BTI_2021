function varargout = mcgls(varargin);
%MCGLS(A,B,TAU) 
%  solves the systems of linear equations 
%            (A'*A+TAU(I)*EYE(N))*X(:,I)=A'*B    (I=1:L) 
%  for X with the multishift CGLS algorithm. The right hand side (column)
%  vector B must have length M, where the coefficient matrix A is M-by-N.
%  TAU is a vector of reals (the shifts) of length L.
%
%  [X,NORM_R] = MCGLS(A,B,TAU) returns the solution X and the convergence
%  history, that is, an array of the ESTIMATED relative residual
%  norms NORM_R at each iteration step. Here
%    NORM_R(:,I) = NORM(A'*B-(A'*A+TAU(I)*EYE(N))*X(:,I))/NORM(A'*B) 
%
%  MCGLS(A,B,TAU)
%  MCGLS('Afun',B,TAU)
%  The first input argument is either an M-by-N matrix (which can be full
%  or sparse, symmetric or nonsymmetric, real or complex), or a string
%  containing the name of an M-file which applies a linear operator to a
%  N-vector. In the latter case, the M-file, say Afun.m, must return an
%  N-vector with Y = Afun(X,'')  (for Y=A*X)  and with 
%  Y = Afun(X,'transpose')  (for Y=A'*X). If Y = Afun(X,'') results in an
%  error then the operator in Afun is supposed to be symmetric and 
%  Y = Afun(X) is taken for Y = A*X as well as for Y = A'*X.
%
%  X = MCGLS(A,B,TAU,TOL,MAX_IT,DISP) specifies parameters as tolerance TOL
%  and maximum number of iterations MAX_IT of the method. If DISP == 1
%  then the options are printed and convergence history is plotted. The
%  iteration is started from the trivial initial guess. Iterates for each
%  system are produced until the respective system either converges, or
%  the maximum number MAX_IT is reached. Convergence for the I-th vector
%  X(:,I) is achieved when the iterate X(:,I) has relative residual norm
%  NORM_R(:,I) less than or equal the tolerance TOL(I). TOL is a vector 
%  of positive reals of length <= L.
%    Defaults:
%                TOL = 1e-9*ONES(SIZE(TAU))
%             MAX_IT = 200
%                TAU = 0
%               DISP = 0
%    If  LT=LENGTH(TOL)<L then TOL(I)=TOL(LT) for I=LT:L
%  
%
%  The parameters can also be defined in a struct:
%  X = MCGLS(A,B,TAU,OPTIONS) or X = MCGLS(A,B,OPTIONS0), where
%   OPTIONS = STRUCT('residual_acc',TOL,'max_it',MAX_IT,'disp',DISP)
%  and  
%   OPTIONS0 = ...
%     STRUCT('Shifts',TAU,'residual_acc',TOL,'max_it',MAX_IT,'disp',DISP)
%
%  X = MCGLS(A,B,TAU,OPTIONS,PAR1,PAR2,..)
%  Additional parameters are passed to the M-file for the operator
%  Y = Afun(Y,'',PAR1,PAR2,...) and Y = Afun(Y,'transpose',PAR1,PAR2,..)
%
%  MCGLS returns the default values.
%
%
% Reference:
%   Jasper van den Eshof and Gerard Sleijpen,
%     "Accurate conjugate gradient methods for families of shifted systems",
%     Applied Numerical Mathematics, Volume 49, Issue 1, 2004. Pages 17-37. 


%  [X,NORM_R,J] = MCGLS(A,B,TAU) returns also the indices J of unconverged 
%  X(:,I) (for I=J). 
%
%  The fieldnames in the OPTIONS struct are identified by the first two
%  letters. The identification is case insensitive

% Jasper van den Eshof, Gerard Sleijpen
% Copyright (c), January 13, 2003


%-----------------------------------------------------------------------------
%%% Set default values and names of default values
defaultoptions=struct('Shifts',0,'residual_acc',1.e-9,'max_it',200,'disp',0);
%-----------------------------------------------------------------------------
%%% Read input parameters
[n,b,tau,tol,max_it,show]=setparameters(defaultoptions,varargin{:});
%-----------------------------------------------------------------------------
%%% Display options
options=struct('Shifts',tau,'residual_acc',tol,'max_it',max_it,'disp',show);
if nargin==0, default_options=options, return, end
if show,  options=options, end
if n<1, [varargout{1:nargout}]=output([],[],[],show); return, end
%-----------------------------------------------------------------------------

if size(tau,1)>1, tau=tau'; end
if size(tol,1)>1, tol=tol'; end
l=length(tau);  J=1:l;
lt=length(tol); tol(1,lt:l)=tol(1,lt);


z=b;
r=mvt0(z); rho = norm(r); 
u = r; 
x=zeros(n,l); p=r*ones(1,l);
gamma=ones(1,l); tt=tau;
k=1; tol=tol*rho; rho=rho*rho; 

NORM_r(k,:)=sqrt(rho)./abs(gamma);

while  k <= max_it 
             
  y=mv(u); 
  sigma=y'*y; alpha=rho/sigma;  
  z=z-alpha*y; 
  r=mvt(z); sigma=rho; rho=r'*r;  
  beta=rho/sigma;   
  u=r+beta*u; 

  % Hist(x)

  for j=J
    %%% Differential QD step
    el=1+alpha*tt(j);
    tt(j)=tau(j)+(beta/el)*tt(j);  gamma(j)=gamma(j)*el;
    %%% update additional system
    x(:,j)=x(:,j)+p(:,j)*(alpha/gamma(j));
    p(:,j)=r+p(:,j)*(beta/el); 
  end

  normr=sqrt(rho)./abs(gamma(J));
  NORM_r(k+1,:)=NORM_r(k,:);
  k=k+1; NORM_r(k,J)=normr; 

  J=J(find(normr>tol(J)));      % find unconverged
  if isempty(J), break, end     % stop if all are converged

end % while

% Hist(x), Hist


if nargout == 0
   solution = x
else
  [varargout{1:nargout}]=output(x,NORM_r,J,show);
end

return
%====================================================================
%  ACTIONS OPERATORS
%====================================================================
function y=mv(y)
% Action of the matrix/operator A

  global A_operator Operator_params

  if ischar(A_operator)
    y=feval(A_operator,y,Operator_params{:});
  else 
    y=A_operator*y;
  end

return
%-------------------------------------------------------------------
function y=mvt(y)
% Action of the transpose of the matrix/operator A

  global A_operator OperatorTr_params

  if ischar(A_operator)
    y=feval(A_operator,y,OperatorTr_params{:});
  else 
    y=(A_operator)'*y;
  end

return
%-------------------------------------------------------------------
function y=mvt0(y)
% Check definition of the transpose of the matrix/operator A.
%
% If A an operator (function) then y=A(y,'transpose',...)
% should return a vector. If not, then A is supposed to be
% symmetric, i.e. A'*y=A*y.

  global A_operator Operator_params OperatorTr_params

  if ischar(A_operator) 
    try %%% is A(y,'transpose',...) defined?
      OperatorTr_params={'transpose',Operator_params{:}};
      y=feval(A_operator,y,OperatorTr_params{:});
      Operator_params={'',Operator_params{:}};
    catch %%% if not A is symmetric; 
      OperatorTr_params=Operator_params;
      y=feval(A_operator,y,OperatorTr_params{:});
    end
  else 
    y=(A_operator)'*y;
  end

return
%====================================================================
%      OUTPUT
%====================================================================
function varargout=output(x,hist,UnConv,show)

  if show
    J=0:size(hist,1)-1; plot(J,log10(hist),'.-')
    title('Convergence history of Multishift CGLS')
    xlabel('iteration step k')
    ylabel('log_{10}(|r_k|_2), where r_k\equiv A''*b-(A''*A+tau(i)*I)*x_k(:,i), i=1,..')
  end

  varargout{1}=x; 
  if nargout>1, varargout{2}=hist; end
  if nargout>2, varargout{3}=UnConv; end

return
%====================================================================
%        INPUT
%====================================================================
function varargout=setparameters(dopt,varargin);
% [n,b,tau,tol,max_it,show]=setparameters(dopt,varargin);

global A_operator Operator_params

  na=length(varargin); 
  if na == 0
     n=0; b=[];
  else %%% find the matrix
    A_operator=varargin{1};
    if ischar(A_operator)
      n=-1;
      if exist(A_operator) ~=2
        msg=sprintf('  Can not find the M-file ''%s.m''  ',A_operator);
        errordlg(msg,'MATRIX'),n=-2; return
      end
    else
      [m,n] = size(A_operator);
    end

    %%% find the right hand side vector
	 b=varargin{2}; 
    [nb,mb]=size(b);
    if n > 0 
      if nb ~= m | mb ~=1
        msg=sprintf(' The right hand side vector b')
        msg=sprintf('%s\n is %i by %i, b must be %i by 1.',msg,nb,mb,m);
        errordlg(msg,'INPUT VECTOR'),n=-3;  return
      end
    else
      n=nb; 
    end
  end

  varargout{1}=n;
  varargout{2}=b;

  %%% find the other parameters
  [varargout{3:nargout+1}]=findparameters(1,dopt,varargin{3:na});
  Operator_params=varargout{nargout+1};


return
%======================================================================
%======================================================================
%======================================================================
function varargout=findparameters(m,dopt,varargin);
%[PAR1,PAR2,...,OTHERPARS]=...
%      findparameters(M,DEFAULT_OPTIONS,INPAR_1,INPAR_2,INPAR_3,...);
% Let K be the largest index <= M such that INPAR_1,...,INPAR_K are not structures.
% Then these values are assigned in this order to the first K output parameters;
% they correspond to the first K values in the DEFAULT_OPTIONS.
% Th other output arguments are determined according to FINDPARAMETER

% Jasper van den Eshof, Gerard Sleijpen
% Copyright (c), January 13, 2003

  na=length(varargin); no=nargout;
  fopts = fieldnames(dopt);
  j=1;
  while j<=min(m,na)
    if isstruct(varargin{j}), break, end
    varargout{j}=varargin{j}; dopt=rmfield(dopt,char(fopts(j,:))); j=j+1; 
  end
  [varargout{j:no}]=findparameter(dopt,varargin{j:na});

return
%======================================================================
function varargout=findparameter(dopt,varargin);
%[PAR1,PAR2,...,OTHERPARS]=...
%      mergeopt(DEFAULT_OPTIONS,INPAR1,INPAR2,INPAR3,...);
%
% Replaces DEFAULT_OPTIONS values by new values from the input parameter list
% INPAR1,INPAR2,INPAR3,...
% and list the values in order of occurence in DEFAULT_OPTIONS list.
%
% INPAR1 is not a structure then:
% Input parameters are classified as strings, booleans, pos. integers, doubles
% (in this order) and within a class they are assigned in order of occurence.
% The parameters that are excessive with respect to DEFAULT_OPTIONS
% and the ones that have not been classified are collected
% in OTHERPARS.
%
% If INPAR1 is a structure then:
% INPAR1 is taken as an OPTIONS list and DEFAULT_OPTIONS values are 
% replaced by the new values from OPTIONS and the values from the 
% new "merged" option list are listed in PAR1,...  in order of 
% occurence in DEFAULT_OPTIONS list.
% Field names are case insensative and only the first 2 characters 
% are checked.
%

% Matlab does not distinguish between 1.0 and 1. As a consequence
% the double 1.0 will be interpreted as boolean. 
% The same remark aplies to integers. 

% N is the number of first characters to be compared in field names
  N=2;

  if isempty(varargin) | ~isstruct(varargin{1})
     [varargout{1:nargout}]=findopt(dopt,varargin{1:nargin-1});
  else
     [varargout{1:nargout-1}]=mergeopt(dopt,varargin{1},N);
     if nargin > 2,  %%% collect the additional parameters in one cell
        varargout{nargout}=varargin(2:nargin-1); 
     else, varargout{nargout}={};
     end
  end

return
%====================================================================
function varargout=findopt(dopt,varargin);
%[PAR1,PAR2,...,OTHERPARS]=...
%      mergeopt(DEFAULT_OPTIONS,INPAR1,INPAR2,INPAR3,...);
%
% Replaces DEFAULT_OPTIONS values by new values from the input parameter list
% INPAR1,INPAR2,INPAR3,...
% and list the values in order of occurence in DEFAULT_OPTIONS list.
% Input parameters are classified as strings, booleans, pos. integers, vectors 
% (in this order) and within a class they are assigned in order of occurence.
% The parameters that are excessive with respect to DEFAULT_OPTIONS
% and the ones that have not been classified are collected
% in OTHERPARS.

  fopts = fieldnames(dopt);
  ndo=length(fopts);
  jds=[]; jdb=[]; jdi=[]; jdv=[]; jdm=[];
  jde=[]; % To check for emtpty default options
  for j=1:ndo
    ar=getfield(dopt,char(fopts(j,:)));
    if isempty(ar),       jde=[jde,j]; % empty
    elseif ischar(ar),    jds=[jds,j]; % strings
    elseif isboolean(ar), jdb=[jdb,j]; % booleans
    elseif isinteger(ar), jdi=[jdi,j]; % positive integers
    elseif isvector(ar),  jdv=[jdv,j]; % vectors
    else,                 jdm=[jdm,j]; % matrices
    end
    varargout{j}=ar;
  end

  if ~isempty(jde), 
    str=char(fopts(jde(1),:));
    str=sprintf('The option ''%s'' is empty.',str);
    str=[str,...
  sprintf('\nEmpty strings and empty non-strings are indistinguishable!')];
    error(str); 
  end

  % defaults=varargout
  % jde,jds,jdb,jdi,jdv,jdm

  js=0; jb=0; ji=0; jv=0; jm=0;
  no=length(varargin); J=1:no;
  for j=1:no
    ar=varargin{j};
    if isempty(ar)
    elseif (ischar(ar) & js<length(jds))         % strings
      js=js+1; varargout{jds(js)}=ar; J=J(find(J~=j));
    elseif ~isnumeric(ar)        % rule out cells, structs, strings
    elseif (isboolean(ar) & jb<length(jdb))      % boolean 
      jb=jb+1; varargout{jdb(jb)}=ar; J=J(find(J~=j));
    elseif (isinteger(ar) & ji<length(jdi))      % positive integers
      ji=ji+1; varargout{jdi(ji)}=ar; J=J(find(J~=j));
    elseif (isvector(ar) & jv<length(jdv))       % vectors
      jv=jv+1; varargout{jdv(jv)}=ar; J=J(find(J~=j));
    elseif jm<length(jdm),                       % matrices
      jm=jm+1; varargout{jdm(jm)}=ar; J=J(find(J~=j));
    end
  end

  %%% collect the additional parameters in one cell
  varargout{ndo+1}=varargin(J);

return
%====================================================================
function varargout=mergeopt(dopt,opt,N);
%[PAR1,PAR2,...]=mergeopt(DEFAULT_OPTIONS,OPTIONS,N);
%
% Replaces DEFAULT_OPTIONS values by new values from OPTIONS
% and list the values in order of occurence in DEFAULT_OPTIONS list.
% Field names are case insensative and only the first N
% characters are checked.

  fdopts = fieldnames(dopt);
  fopts=fieldnames(opt);
  ndo=length(fdopts);
  for j=1:ndo
    name=char(fdopts(j,:)); 
    ar=getfield(dopt,name);
    name=lower(name(1:N));
    [varargout{j},ok] = findfield(opt,fopts,name,ar);
  end

return
%======================================================================
function [a,ok]=findfield(options,fopts,str,default)
% Searches the fieldnames in FOPTS for the string STR.
% The field is detected if only the first part of the fieldname
% matches the string STR. The search is case insensitive.
% If the field is detected, then OK=1 and A is the fieldvalue.
% Otherwise OK=0 and A=DEFAULT

   l=size(str,2); j=min(find(strncmpi(str,fopts,l)));
   if ~isempty(j)
      a=getfield(options,char(fopts(j,:))); ok=1;
   elseif nargin>3
      a=default; ok=0;
   else
      a=[]; ok=0;
   end

return
%======================================================================
function ok=isboolean(ar)
  ok = (length(ar)==1 & (ar == 0 | ar == 1));
return
%======================================================================
function ok=isinteger(ar)
  ok = (length(ar)==1 & round(ar) == ar & ar > 0);
return
%======================================================================
function ok=isvector(ar)
  ok = (min(size(ar))==1);
return
%======================================================================


