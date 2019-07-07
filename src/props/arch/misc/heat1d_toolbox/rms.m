  function [rms]=rms(r)
  % calulates the rms of vector r
  rms=sqrt(sum(r.*r)/(length(r)-1));