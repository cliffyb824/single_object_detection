function plot_curves_on_image(curves,I,nrow,scalefac,fig,linew,cmax,color)
  
N=length(curves);
[T,~]=size(curves{1});
betapix=zeros(T,2);
figure(fig); clf; imagesc(uint8(I)); colormap gray; axis off tight; caxis([0,cmax]); daspect([1 1/scalefac 1])
hold on;
for i=1:N
    [betapix(:,2),betapix(:,1)]=cartesian_to_pixel(curves{i}(:,1),curves{i}(:,2),nrow,scalefac);
    plot(betapix(:,1),betapix(:,2),color,'LineWidth',linew);
end

    