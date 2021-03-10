function L = skewLogisticLikelihood(theta,x)

n=length(x);
mu=theta(1);
sig=theta(2);
alpha=theta(3);
z=(x-mu)/sig;
L=n*log(2)-n*log(sig)-sum(z)-2*sum(log(1+exp(-z)))-sum(log(1+exp(-alpha*z)));
