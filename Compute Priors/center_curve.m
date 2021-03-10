function beta_new = center_curve(beta)

% Translate a continuous curve beta to have centroid (center of mass) at 
% the origin.
%
% Inputs: 
%   beta - n x p matrix, i.e., an n-tuple in R^p. 
% Outputs: 
%   beta_new - n x p matrix. A centered version of beta.

[n,~]=size(beta);
beta_new=beta-repmat(calculateCentroid(beta),n,1);
