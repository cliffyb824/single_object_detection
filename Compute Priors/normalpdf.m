function fx = normalpdf(x,mu,sig)

z=(x-mu)/sig;
fx=(2*pi*sig^2)^(-1/2)*exp(-1/2*z^2);