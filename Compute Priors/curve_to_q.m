function q = curve_to_q(beta)

% Compute the square root velocity function (srvf) q(t) of the continuous 
% curve \beta(t). The formula for q is q(t) = \dot{\beta}(t) / \sqrt(|\dot{\beta}(t)|}.
%
% Inputs: 
%   beta - An n x p matrix, i.e., an n-tuple in R^p. 
% Outputs: 
%   q - An n x p matrix representing the srvf of beta. q is normalized to
%       have unit norm in L2 function space, which is equivalent to beta 
%       having unit length.

beta=beta';
n=length(beta);
betadot=gradient(beta,1/(n-1));
normbetadot=sqrt(sum(betadot.*betadot));
q=betadot*diag(1./sqrt(normbetadot)); 
q=q';
q=q/sqrt(InnerProd_Q(q,q));



