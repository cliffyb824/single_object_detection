function [E,gradE,normgrad] = Energy_KarcherMean(betamean,Data,mu,basis,register)

% Compute the objective function to be minimized, and its gradient, in the 
% gradient descent Karcher mean algorithm implemented in findKarcherMean.m. 
% The objective function is E(q) = \sum_{i=1}^N d(q,q_i)^2, where q_i is 
% the srvf of \beta_i, and d(.,.) is the geodesic distance on closed curve, 
% elastic shape space. The gradient of E at q is given by 
% \nabla E(q) = -(2/N) \sum_{i=1}^N \exp_{q}^{-1}(q_i), negative 2 times 
% the average shooting vector from q to all other data points q_i. 
%
% Inputs: 
%   betamean - A T x 2 matrix representing the current estimate of the 
%       mean curve.
%   Data - A 1 x N cell array where Data{i} is \beta_i. 
%   mu - A T x 2 matrix representing the SRVF of betamean.
%   basis - A cell array consisting of the two basis vectors of the normal
%       space of the closed curve manifold at mu.
%   register - A logical value indicating whether or not to perform elastic
%       registration. If not, then the analysis is equivalent to a
%       landmark-based analysis. 
%
% Outputs: 
%   E - The scalar value of the energy function evaluated at the srvf of 
%       betamean.
%   gradE - The n x p negative gradient vector in the tangent space of the 
%       shape manifold.
%   normgrad - The norm of gradE using the L2 inner product. 

N=length(Data);
[n,p]=size(Data{1});
basis=gramSchmidt2(basis);

% Compute shooting vector from betamean to each data point. Energy is sum
% of squared geodesic distances. 
sumv=zeros(n,p);
E=0;
for i=1:N
    beta=Data{i};
    [v,d,~]=inverseExp_Coord(betamean,beta,register); % Compute shooting vector in open curve space
    v=projectTangent(v,mu,basis); % Project to closed curve space
    sumv=sumv+v;
    E=E+d^2;
end

% Negative gradient of energy is 2 times the average shooting vector.
gradE=2*sumv/N;

% Project to horizontal space of O(p) orbit of q. In other words, remove
% the rotational component of the negative gradient vector. (Optional)
mu=curve_to_q(betamean);
b=mu*[0 -1; 1 0];
b0=b/sqrt(InnerProd_Q(b,b));
gradE=gradE-InnerProd_Q(gradE,b0)*b0;

% Norm of gradient
normgrad=sqrt(InnerProd_Q(gradE,gradE));