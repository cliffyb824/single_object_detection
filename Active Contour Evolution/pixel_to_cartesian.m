function [u,v] = pixel_to_cartesian(i,j,nrow,scalefac)

u=j-1;
v=scalefac*(nrow-i);