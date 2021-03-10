function beta = q_to_curve(q)

% Recover the continuous coordinate function \beta(t) from the square root 
% velocity function (srvf) q(t). The formula for q is 
% q(t) = \dot{\beta}(t) / \sqrt(|\dot{\beta}(t)|}.
%
% Inputs: 
%   q - An n x p matrix representing the srvf of beta.
%   
% Outputs: 
%   beta - An n x p matrix, i.e., an n-tuple in R^p. 

q=q';
n=length(q);
qnorm=sqrt(sum(q.*q));
integrand=q*diag(qnorm);
beta=cumtrapz(integrand,2)/n;
beta=beta';