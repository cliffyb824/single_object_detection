function [dshape,dbin] = compare_to_ground_truth(beta,betatrue,nrow,ncol,scalefac)

[betapix(:,2),betapix(:,1)]=cartesian_to_pixel(beta(:,1),beta(:,2),nrow,scalefac);
C=poly2mask(betapix(:,1)',betapix(:,2)',nrow,ncol);

[betatruepix(:,2),betatruepix(:,1)]=cartesian_to_pixel(betatrue(:,1),betatrue(:,2),nrow,scalefac);
Ctrue=poly2mask(betatruepix(:,1)',betatruepix(:,2)',nrow,ncol);

union=double(C | Ctrue);
intersection=double(C & Ctrue);

% Shape distance and binary image distance measures
[~,dshape,~]=inverseExp_Coord(beta,betatrue,true);
dbin=sum(sum(union-intersection))/sum(sum(union));




