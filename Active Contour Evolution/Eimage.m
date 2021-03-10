function [E,gradEbeta] = Eimage(Sbeta,Sopt,Sprior)

I=Sopt.I; % log of image
uu=Sopt.uu;
vv=Sopt.vv;
scalefac=Sopt.scalefac;
theta0=Sprior.theta0;
theta1=Sprior.theta1;

[nrow,ncol]=size(I);
C=curve_to_regions(Sbeta.beta,scalefac,nrow,ncol); 
x0=I(C==0);
x1=I(C==1);
E=-skewLogisticLikelihood(theta0,x0)-skewLogisticLikelihood(theta1,x1);
E=E*scalefac;

val=interp2(uu,vv,I,Sbeta.beta(:,1),Sbeta.beta(:,2));
ratio=log(skewLogisticpdf(val,theta1(1),theta1(2),theta1(3))./...
    skewLogisticpdf(val,theta0(1),theta0(2),theta0(3)));
gradEbeta=-repmat(ratio,1,2).*Sbeta.n;

% For debugging:

% keyboard;
% figure(4); hold on; plot(Sbeta.beta(:,1),Sbeta.beta(:,2),'LineWidth',2); axis equal off;
% quiver(Sbeta.beta(:,1),Sbeta.beta(:,2),-gradEbeta(:,1),-gradEbeta(:,2),'b');
% curves{1}=Sbeta.beta;
% plot_curves_on_image(curves,exp(I),nrow,scalefac,5,2,Sopt.cmax,'r')
% plot_curves_on_image(curves,C,nrow,scalefac,6,2,2,'r')
