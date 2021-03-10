function C = curve_to_regions(beta,scalefac,nrow,ncol)

[betapix(:,2),betapix(:,1)]=cartesian_to_pixel(beta(:,1),beta(:,2),nrow,scalefac);
C=poly2mask(betapix(:,1)',betapix(:,2)',nrow,ncol);


