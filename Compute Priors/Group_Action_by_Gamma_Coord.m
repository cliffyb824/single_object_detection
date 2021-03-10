function beta_new = Group_Action_by_Gamma_Coord(beta,gamma)

% Compute the group action of \gamma, an element of the diffeomorphism 
% (re-parameterization) group \Gamma, on the curve \beta. The group action 
% in this case is simply the function composition \beta \mapsto
% \beta(\gamma(t)). 
%
% Inputs: 
%   beta - An n x p matrix representing a continuous curve in R^p. 
%   gamma - An n x 1 array of strictly increasing values of t with
%       gamma(1)=0 and gamma(n)=1. 
%
% Outputs:
%   beta_new - An n x p matrix representing the re-parameterized curve 
%       \beta(\gamma(t)).

beta=beta';
gamma=gamma';
n=length(beta);
beta_new=spline(linspace(0,1,n),beta,gamma);
beta_new=beta_new';


