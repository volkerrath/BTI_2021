function  f = splval(S,t,sp)
%SPLVAL  Values of cubic spline.    f = splval(S,t {,sp})
% The spline  s  is given by the struct  S (as descibed in the 
% output from  SPLFIT).
% If  sp  is present, then  
% sp = 1 :  t  must be a scalar in the range  [S.x(1), S.x(end)]
%           f  is a 2*3 matrix with
%               f(1,:) = [s(t), s'(t), s"(t)]
%           and  f(2,:)  holds their standard deviation.
% sp = 2 :  t  is not used.
%           f  is a vector with zeros of  s' .
% Otherwise,  t  can be a vector, and  f  is a vector of
%             the same type with  f(i) = s(t(i)) .
%
% Hans Bruun Nielsen, IMM, DTU.  00.08.24 / 08.29

  if  nargin > 2    % Special cases
    switch  sp
    case 1  % s, s', s''  and standard deviation
      if  any(size(t) ~= 1)
        error('Argument  t  must be a scalar'),  end
      if  (t < S.x(1)) | (t > S.x(end))
        error('Argument  t  is outside knot range'), end
      f = splvdv(S,t);  
  
    case  2 % zeros of  s'
      pp = S.pp;
      if  any(imag(pp))
        error('Can only be used with real-valued S.pp'), end
      p = size(pp,2);   f = zeros(1,2*p);   kf = 0;
      for  j = 1 : p-1
        dc = fliplr([1:3].*pp(3:5,j).');
        r = roots(dc);
        if  ~any(imag(r))  % all real
          k = find(0 <= r & r <= diff(pp(1,j:j+1)));
          for  l = 1 : length(k)
            kf = kf +1;  f(kf) = pp(1,j) + r(k(l));
          end
        end
      end
      % Remove multiple entries
      thr = max(10*eps*(pp(1,end) - pp(1,1)), .001*min(diff(pp(1,:))));
      f = sort(f(1:kf));   k = 1;
      for  j = 2 : kf
        if  f(j) <= f(k)+thr
          f(k) = (f(k) + f(j))/2;
        else,  k = k+1;   f(k) = f(j); end
      end
      f = f(1:k);
        
    otherwise
      error('sp  should have the value 1 or 2')
    end % switch
    return
  end % special
  
  %  Simple case
  if  min(size(t)) ~= 1
    error('Argument  t  must be a vector'),  end
  % Sort arguments to compute from left to rights
  [t,js] = sort(t);   f = zeros(size(t));
  pp = S.pp;   p = size(pp,2);   q = length(t);
  %  Points to the left
  i = find(t <= pp(1,1));   usd = length(i);
  if  usd
    v = t(i) - pp(1,1);
    f(js(i)) = pp(2,1) + v.*(pp(3,1) + pp(4,1)*v);
  end
  for  j = 1 : p-1
    i = find(t(usd+1:q) <= pp(1,j+1));   li = length(i);
    if  li    % Unused arguments
      v = t(usd+i) - pp(1,j);
      f(js(usd+i)) = pp(2,j) + v.*(pp(3,j) + v.*(pp(4,j) + pp(5,j)*v));
      usd = usd + li;
    end
  end
  if  usd < q
    v = t(usd+1:q) - pp(1,p);
    f(js(usd+1:q)) = pp(2,p) + v.*(pp(3,p) + pp(4,p)*v);
  end
  
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

function  vdv = splvdv(S,t)
% s, s', s" and their standard deviation, evaluated at  t
  j = max(find(t > S.x));  
  if  length(j) == 0,  j = 1; end
  xx = [S.x(1)*[1 1] S.x S.x(end)*[1 1]];
  B = splbsd(t,xx(j:j+5));   C = B * diag(S.sdc(j:j+3));
  vdv = [B*S.c(j:j+3).' sqrt(diag(C*C'))].';
