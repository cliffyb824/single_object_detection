function [i,j] = cartesian_to_pixel(u,v,nrow,scalefac)

i=nrow-v/scalefac;
j=u+1;
