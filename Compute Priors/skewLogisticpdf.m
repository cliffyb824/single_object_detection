function fx = skewLogisticpdf(x,mu,sig,alpha)

z=(x-mu)/sig;
phi=1/sig*exp(-z)./(1+exp(-z)).^2;
Phi=1./(1+exp(-alpha*z));
fx=2*phi.*Phi;