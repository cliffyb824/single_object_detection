function L = normallikelihood(theta,x)

n=length(x);
mu=theta(1);
sig=theta(2);
z=(x-mu)/sig;
L=-n/2*log(2*pi)-n*log(sig)-1/(2*sig^2)*sum(z.^2);