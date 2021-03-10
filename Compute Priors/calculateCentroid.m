function c = calculateCentroid(beta)

% Compute the centroid (center of mass) of a continuous curve beta. The
% formula for center of mass is c = (1/L) * \int_0^1 \beta(t)|\dot{\beta}(t)|dt,
% where L = \int_0^1 |\dot{\beta}(t)|dt is the length of the curve.
%
% Inputs: 
%   beta - n x p matrix, i.e., an n-tuple in R^p. 
% Outputs: 
%   c - 1 x p vector.

beta=beta';
n=length(beta);
betadot=gradient(beta,1/(n-1));
normbetadot=sqrt(sum(betadot.*betadot));
integrand=beta*diag(normbetadot);
L=trapz(linspace(0,1,n),normbetadot);
c=trapz(linspace(0,1,n),integrand,2)/L;
c=c';