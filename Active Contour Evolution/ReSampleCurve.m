function betanew = ReSampleCurve(beta,N)

% Re-sample a continuous curve beta so that it has N uniformly spaced 
% sample points. 
%
% Inputs: 
%   beta - An n x p matrix representing the n sample points on the curve. 
%   N - A number of sample points for the output curve.
% Outputs: 
%   betanew - An N x p matrix representing a uniformly sampled version of
%       beta. 

beta=beta';
n=length(beta);
d=beta(:,2:n)-beta(:,1:end-1);
w=[0 sqrt(sum(d.*d))];
gamma=cumsum(w)/sum(w);   
gammanew=linspace(0,1,N);
betanew=spline(gamma,beta,gammanew);
betanew=betanew';