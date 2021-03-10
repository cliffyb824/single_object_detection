function gammaI = invertGamma(gamma)

% Compute \gamma^{-1}(t) for \gamma(t)\in\Gamma, the diffeomorphism group
% on [0,1]. 
%
% Inputs: 
%   gamma - An n x 1 array of strictly increasing values of t with
%       gamma(1)=0 and gamma(n)=1. 
% 
% Outputs:
%   gammaI - The inverse of gamma is also an n x 1 array of strictly 
%       increasing values of t with gammaI(1)=0 and gammaI(n)=1. 

gamma=gamma';
n=length(gamma);
x=linspace(1/n,1,n);
gammaI=interp1(gamma,x,x,'linear');
gammaI=(gammaI-gammaI(1))/(gammaI(n)-gammaI(1));
gammaI=gammaI';
