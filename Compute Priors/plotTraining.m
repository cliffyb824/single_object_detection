function plotTraining(Sprior)

close all;

fig=0;
[nrow,~,ntrain]=size(Sprior.snips);
scalefac=Sprior.scalefac;

fig=fig+1; figure(fig);
h1=tight_subplot(2,5,0.01,0.01,0.01); % Hard coded for 10 training images
for i=1:ntrain
    axes(h1(i));
    imagesc(Sprior.snips(:,:,i)); axis off tight; colormap gray; caxis([0,5]); % Hard coded caxis limits
    daspect([1 1/Sprior.scalefac 1]);
end

fig=fig+1; figure(fig);
h2=tight_subplot(2,5,0.01,0.01,0.01);
for i=1:ntrain
    axes(h2(i));
    imagesc(Sprior.snips(:,:,i)); axis off tight; colormap gray; caxis([0,5]); % Hard coded caxis limits
    daspect([1 1/Sprior.scalefac 1]);
    hold on; 
    [betapix(:,2),betapix(:,1)]=cartesian_to_pixel(Sprior.Sbetatrain(i).beta(:,1),Sprior.Sbetatrain(i).beta(:,2),nrow,scalefac);
    plot(betapix(:,1),betapix(:,2),'r','LineWidth',2);
end
    
fig=fig+1; figure(fig);
h3=tight_subplot(2,5,0.01,0.01,0.01);
for i=1:ntrain
    axes(h3(i));
    imagesc(Sprior.Ctrue(:,:,i)); axis off tight; daspect([1 1/Sprior.scalefac 1])
end

pix=linspace(-4,4,400); % Hard coded log-pixel distribution domain
f0=skewLogisticpdf(pix,Sprior.theta0(1),Sprior.theta0(2),Sprior.theta0(3));
f1=skewLogisticpdf(pix,Sprior.theta1(1),Sprior.theta1(2),Sprior.theta1(3));

fig=fig+1; figure(fig);
plot(pix,f0,pix,f1,'LineWidth',2);
% title('Pixel Densities per Region (skew-logistic)')
xlabel('Log of Pixel Value')
legend('Background','Target')

lengths=linspace(0,600,200); % Hard coded length distribution domain
g=normpdf(lengths,Sprior.muL,Sprior.sigL);

fig=fig+1; figure(fig);
plot(lengths,g,'LineWidth',2);
title('Target Boundary Length Density (Gaussian)')
xlabel('Boundary Length')

% Plot training target shapes 
% n=ceil(sqrt(ntrain));
fig=fig+1; figure(fig);
h4=tight_subplot(2,5,0.01,0.01,0.01); % Hard coded for 10 training shapes
for i=1:ntrain
    axes(h4(i));
    plotCurve(Sprior.Sbetatrain(i).beta,fig,2,0,[17 44]);
end

% Plot Karcher mean
fig=fig+1;
plotCurve(Sprior.Sbetamean.beta,fig,2,0,[17 44]);
title('Target Mean Shape')
