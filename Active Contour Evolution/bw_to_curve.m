function beta = bw_to_curve(BW,T,nrow,scalefac)

dx=bwboundaries(BW);
betapix=dx{1};
[beta(:,1),beta(:,2)]=pixel_to_cartesian(betapix(:,1),betapix(:,2),nrow,scalefac);
beta=ReSampleCurve(beta,T);
% beta=ReSampleCurve(beta,T2);
% beta=center_curve(beta);

