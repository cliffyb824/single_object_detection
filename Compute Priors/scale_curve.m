function [beta0,L] = scale_curve(beta)

% Computes length of a continuous curve beta, and rescales it to have 
% length 1. The formula for the length of the curve is 
% L = \int_0^1 |\dot{\beta}(t)|dt.
%
% Inputs: 
%   beta - n x p matrix, i.e., an n-tuple in R^p. 
% Outputs: 
%   beta0 - beta rescaled to have unit length.
%   L - The original length of beta. 

beta=beta';
[~,n]=size(beta);
betadot=gradient(beta,1/(n-1));
normbetadot=sqrt(sum(betadot.*betadot));
L=trapz(linspace(0,1,n),normbetadot);
beta0=beta'/L;
