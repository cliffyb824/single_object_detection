function [E,gradE] = Elength(Sbeta,Sprior)

E=log(sqrt(2*pi)*Sprior.sigL)+(Sbeta.length-Sprior.muL)^2/(2*Sprior.sigL^2);
gradE=(Sbeta.length-Sprior.muL)/(Sprior.sigL^2)*diag(Sbeta.kappa)*Sbeta.n;

% % For debugging
% keyboard;
% beta=Sbeta.beta;
% figure(1); clf; hold on;
% plot(beta(:,1),beta(:,2),'b','LineWidth',2); 
% quiver(beta(:,1),beta(:,2),-gradE(:,1),-gradE(:,2),'b'); axis equal off;

