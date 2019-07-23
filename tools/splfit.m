function  [fail, S] = splfit(tyw,x,BC)
%SPLFIT  Weighted least squares fit with cubic spline
%    [fail, S] = splfit(tyw,x{,BC})
%
% Input
% tyw :  Data points and weights.  Array with 2 or 3 columns,
%        tyw(:,1) :  Abscissae, must be in increasing order.
%        tyw(:,2) :  Ordinatates.
%        tyw(:,3) :  Weights.  If  tyw  holds less than 3 columns
%                    all weights are set to 1.
% x   :  Knots.  Must satisfy  
%            x(1) < x(2) <= ... <= x(end-1) < x(end) .
% BC  :  If present, then the fit is computed with respect to
%        boundary conditions
%          BC(k,1:3)*[s(zk); s'(zk); s"(zk)] = BC(k,4)  for  k=1,2
%        with  z1=x(1), z2 = x(end).
%
% Output
% fail :  Performance indicator
%         fail = 0 :  No problems (otherwise, the contents of
%                     S of no significance).
%                1 :  length(x) < 2.
%                2 :  Knots are not in increasing order.
%                3 :  Some data abscissa is outside [x(1), x(end)].
%                4 :  Too few active data points (points with
%                     strictly positive weights).
%                5 :  Schoenberg-Whitney condition is not satisfied.
% S    :  Struct representing the spline.  Fields
%         x   :  Knots.
%         c   :  Coefficients in B-spline representation of  s .
%         pp  :  Piecewise polynomial representation of  s .
%         sdy :  Estimated standard deviation of data.
%         sdc :  Estimated standard deviation of c.
%
% Hans Bruun Nielsen, IMM, DTU.  00.08.25 / 08.31

  % Initialize
	fail = 0;  S = [];
	%  Check knots
  n = length(x) - 1;  [m nw] = size(tyw);
  if  n < 1,  fail = 1;  return, end
  if  any(diff(x) < 0) | x(2) == x(1) | x(end-1) == x(end)
    fail = 2;  return, end
  % Check (and order) data
  [m q] = size(tyw);
  if  q > 2
    j = find(tyw(:,3) > 0);   ma = length(j);
    if  ma < m,  tyw = tyw(j,1:3); end
  else,  ma = m;   tyw = [tyw ones(m,1)]; end
  [tyw(:,1) js] = sort(tyw(:,1));  tyw(:,2:3) = tyw(js,2:3);
  if  tyw(1,1) < x(1) | tyw(end,1) > x(end)
    fail = 3;  return, end  
  % Extended knot set 
  x = x(:)';   xx = [x(1)*[1 1] x x(end)*[1 1]];
  
  %  Boundary conditions
  con = zeros(1,2);    alfa = zeros(1,6);
  if  nargin > 2
    B = splbsd(x(1), xx(1:6));
    alfal = BC(1,1:3)*B(:,1:3);
    if  alfal(1)    % Active constraint
      alfa(1:3) = [-alfal(2:3) BC(1,4)]/alfal(1);   con(1) = 1;  
    end
    B = splbsd(x(n+1), xx(n:n+5));
    alfar = BC(2,1:3)*B(:,2:4);
    if  alfar(3)    % Active constraint
      alfa(4:6) = [-alfar(1:2) BC(2,4)]/alfar(3);   con(2) = 1;
    end
  end  %  Boundary conditions
  
  % Find active points
  if  nw < 3,  tyw = [tyw  ones(m,1)]; end
  j = find(tyw(:,3) > 0);   ma = length(j);
  ns = n+3 - sum(con);
  if  ma < ns,  fail = 4;  return, end
  if  ma < m,  tyw = tyw(j,1:3); end    % Compress data

  if  n == 1    % One spline interval 
    [c  sdy  sdc] = splfit1(tyw,x,alfa,con); 
    
  else   % General case.  Successive orthogonal transformation
    R = zeros(ns,5);  kr = 0;    % Rhs in last column
    i2 = 0;   G = zeros(5,5);   lg = 0;  sumr = 0;
    for  j = 1 : n
      i = find(tyw(i2+1:ma,1) <= x(j+1));   li = length(i);
      if  li    %  Contributions
        i1 = i2+1;  i2 = i2+li;
        A = [repmat(tyw(i1:i2,3),1,5) .* ...
            [splbsv(tyw(i1:i2,1),xx(j:j+5)) tyw(i1:i2,2)]
             zeros(lg,5)];
        if  (j == 1) & con(1)    % Apply boundary constraints
          A(1:li,:) = A(1:li,:) * [alfa(1:2) 0 0 -alfa(3)
            1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 0 1];
        end
        if  (j == n) & con(2)
          A(1:li,:) = A(1:li,:) * [1 0 0 0 0; 0 1 0 0 0
            0 0 1 0 0; 0 alfa(4:5) 0 -alfa(6); 0 0 0 0 1];
        end
        if  j > 1    % Add  G
          A(li+[1:4],:) = G(1:4,:);
        else,  lg = 4; end
        %  Orthogonal transformation
        A = qr(A);   m = min(5,size(A,1));
        G(1:m,:) = triu(A(1:m,:));
      end  % Contribution
      if  (j > 1) | ~con(1)    % Transfer to R and shift G
        kr = kr + 1;  R(kr,:) = G(1,:);
        if  R(kr,1) == 0,  fail = 5;  return, end
        G(1:4,:) = [G(2:5,2:4) zeros(4,1) G(2:5,5)];
      else    % Update residual info from 1st interval
        if  li > 4,  G(4,5) = norm(G(4:5,5)); end
      end    
      G(5,:) = zeros(1,5);
    end % Knot intervals

    %  Transfer last rows to R
    for  k = 1 : ns-kr
      R(kr+k,:) = [G(k,k:4) zeros(1,k-1) G(k,5)];
      if  R(kr+k,1) == 0,  fail = 5;  return, end
    end
  
    % Get solution
    if  ma > ns,   sdy = norm(G(ns+1-kr:4,5))/sqrt(ma - ns);
    else,  sdy = 0; end
    d = zeros(ns+3,1);
    for  i = ns : -1 : 1
      d(i) = (R(i,5) - R(i,2:4)*d(i+[1:3]))/R(i,1);
    end
    if  con(2),  c = [d(1:ns); alfa(4:6)*[d(ns-1:ns); 1]];
    else,  c = d(1:ns); end
    if  con(1),  c = [alfa(1:3)*[c(1:2); 1]; c]; end
    
    %  Estimate standard deviation of the cj
    if  ma > ns
      sdc = splvbc(R,sdy);   lc = length(sdc);
      if  con(2)
        sdc = [sdc; norm(alfa(4:5).*sdc(lc-1:lc)')];
      end
      if  con(1)
        sdc = [norm(alfa(1:2).*sdc(1:2)'); sdc];
      end
    else, sdc = zeros(n+3,1); end
  end  % General case
  
  % Return results
  S = struct('x',x, 'c',c(:).', 'pp',splmpp(x,c), ...
             'sdy',sdy, 'sdc',sdc(:)');
  
% ------------------   Auxiliary functions

function  N = splbsd(t,z)
% Values and first two derivatives of the four nonzero B-splines
% at (scalar)  t,   z3 <= t <= z4
  M = zeros(3,5);   N = zeros(3,4);
  M(1,2) = 1/(z(4) - z(3));
  for  r = 1 : 2
    s = 0 : r;   jj = 3 - r;
    M(r+1,2+s) = ((t - z(jj+s)).*M(r,1+s) + ...
          (z(4+s) - t).*M(r,2+s)) ./ (z(4+s) - z(jj+s));
  end
  N(1,1:3) = (z(4:6) - t) .* M(3,2:4);
  N(1,2:4) = N(1,2:4) + (t - z(1:3)) .* M(3,2:4);
  %  Differentiate
  N(2,:) = -3*diff(M(3,:));
  M(2,2:4) = diff(M(2,1:4)) ./ (z(4:6) - z(1:3));
  N(3,:) = 6*diff(M(2,:));

function  N = splbsv(t,z)
% Values of the four nonzero B-splines at (vector)  t,  
% z3 <= t <= z4
  q = length(t);   N = zeros(q,4);
  z3 = z(3);   z4 = z(4);   hj = z4 - z3;
  u = t - z3;   v = z4 - t;
  N(:,1) = v/(hj*(z4 - z(2)));
  N(:,2) = u/(hj*(z(5) - z3));    % M1
  N(:,3) = u.*N(:,2)/(z(6) - z3);
  N(:,2) = ((t - z(2)).*N(:,1) + (z(5)-t).*N(:,2))/(z(5)-z(2));
  N(:,1) = v.*N(:,1)/(z(4)-z(1));    % M2
  N(:,4) = u.*N(:,3);
  N(:,3) = (t - z(2)).*N(:,2) + (z(6) - t).*N(:,3);
  N(:,2) = (t - z(1)).*N(:,1) + (z(5) - t).*N(:,2);
  N(:,1) = v.*N(:,1);
  
function [c,sdy,sdc] = splfit1(tyw,x,alfa,con)
% Fit with cubic spline over one knot interval.
  m = size(tyw,1);   h = diff(x); 
  u = (tyw(:,1) - x(1))/h;   v = 1 - u;
  uu = u.*u;   vv = v.*v;   wu = tyw(:,3).*u;   wv = tyw(:,3).*v;
  A = [wv.*vv 3*wu.*vv 3*wv.*uu wu.*uu tyw(:,3).*tyw(:,2)];
  if  any(con)  % Constraints
    I = eye(4);
    if  con(1) & con(2)
      B = [alfa(1:3).'  I(1:3,1:2)  alfa(4:6).'  I(1:3,3)];
    elseif  con(1)
      B = [[alfa(1:2).'; 0; alfa(3)] I];
    else
      B = [I(:,1:3)  [0; alfa(4:6).']  I(:,4)];
    end
    A = A*B';  
  end
  ns = 4 - sum(con);   ns1 = ns + 1;
  A = qr(A);   R = triu(A(1:ns,1:ns));
  if  m > ns,   sdy = abs(A(ns1,end))/sqrt(m-ns);
  else,         sdy = 0; end  
  d = R\A(1:ns,end);
  if  con(2),  c = [d; alfa(4:6)*[d(ns-1:ns); 1]];
  else,        c = d; end
  if  con(1),  c = [alfa(1:3)*[c(1:2); 1]; c]; end
  if  sdy    % Estimate deviation of coefficients
    V = inv(R);  dc2 = diag(V'*V);
    a2 = abs(alfa).^2;
    if  con(2),  sdc = [dc2; a2(4:5)*dc2(end-1:end)];
    else,        sdc = dc2; end
    if  con(1),  sdc = [a2(1:2)*sdc(1:2); sdc]; end
    sdc = sdy * sqrt(sdc);
  else,  sdc = zeros(size(c));  end
  
function  pp = splmpp(x,c)
% Given knots  x(1:n+1) and B-spline coefficients c(1:n+3).
% Compute piecewise polynomial representation.
  n = length(x) - 1;   pp = zeros(5,n+1);   p = 0;
  xx = [x(1) x(1) x(:)' x(n+1) x(n+1)];
  for  j = 1 : n
    z = xx(j:j+5);
    if  x(j+1) > x(j)  % Non empty interval
      p = p+1;   pp(1,p) = x(j);
      if  z(2) == x(j)  % After empty interval
        B = splbsd(x(j),z);   q = 1;
      else,  q = 2;  end
      pp(2:4,p) = B(:,q:q+2)*c(j:j+2);   pp(4,p) = .5*pp(4,p);
      B = splbsd(x(j+1),z); 
      d2s = .5*(B(3,2:4)*c(j+1:j+3));
      pp(5,p) = (d2s - pp(4,p))/(3*(x(j+1) - x(j)));
    end  % Non empty interval
  end  % j
  p = p+1;   pp(1,p) = x(n+1);    pp(2,p) = c(n+3);
  pp(3,p) = B(2,3:4)*c(n+2:n+3);  pp(4,p) = d2s;
  pp = pp(:,1:p); 
  
function  sdc = splvbc(R,sdy)
% Given upper triangular, banded R and standard deviation
% sdy  of data points.  Compute standard deviation of 
% B-spline coefficients.
  n = size(R,1);    k = n-1;   
  V = [zeros(n,2) 1./R(:,1)];   sdc = V(:,3).^2;
  V(1:k,2) = -(R(1:k,2).*V(2:n,3))./R(1:k,1);
  sdc(1:k) = sdc(1:k) + V(1:k,2).^2;
  k = n-2;
  V(1:k,1) = -(R(1:k,2).*V(2:k+1,2) + R(1:k,3).*V(3:n,3))./R(1:k,1);
  sdc(1:k) = sdc(1:k) + V(1:k,1).^2;   jj = 1:3;   
  for  k = n-3 : -1 : 1
    j1 = jj(1);  j2 = jj(2);  j3 = jj(3);   jj = jj([3 1 2]);
    V(1:k,j3) = -(R(1:k,2).*V(2:k+1,j1) + R(1:k,3).*V(3:k+2,j2) ...
                  + R(1:k,4).*V(4:k+3,j3))./R(1:k,1);
    sdc(1:k) = sdc(1:k) + V(1:k,j3).^2;
  end
  sdc = sdy * sqrt(sdc);
  
  