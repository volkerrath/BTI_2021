function w = box(n)
% W = (N) returns the N-point boxcar window.
w =ones(n,1);w = w/sum(w);
