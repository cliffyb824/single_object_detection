function beta_std = standardizeCurve(beta)

% Standardize a curve beta, i.e. rescale to have unit length and translate 
% so that the centroid is at the origin
%
% Inputs: 
%   beta - n x p matrix, i.e., an n-tuple in R^p. 
% Outputs: 
%   beta_std - n x p matrix. A standardized version of beta.

beta_std=scale_curve(center_curve(beta));