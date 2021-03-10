function [uu,vv] = generate_coordinates(nrow,ncol,scalefac)

[jj,ii]=meshgrid(1:ncol,1:nrow);
[uu,vv]=pixel_to_cartesian(ii,jj,nrow,scalefac);