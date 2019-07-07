function w = tri(n)
% W = TRI(N) returns the N-point triangular window.
%  needs an odd length sequence
w = 2*(1:(n+1)/2)/(n+1);
w = [w w((n-1)/2:-1:1)]';
w = w/sum(w);
