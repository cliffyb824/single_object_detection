function [v,dist,beta2n] = inverseExp_Coord(beta1,beta2,register)

% Calculate the inverse exponential mapping between two elements of open 
% curve elastic shape space. The shape manifold is described in Chapter 5 
% of the text by Srivastava and Klassen (2016). Here, we fix the first 
% element and optimally rotate and re-parameterize the second to the first.
%
% Inputs: 
%   beta1,beta2 - Two T x 2 matrices representing two parameterized curves 
%       in R^2. Their srvf's in L^2([0,1],R^2) function space are given by 
%       q1 and q2. 
%   register - A logical value indicating whether or not to perform elastic
%       registration of q2 to q1. 
% 
% Outputs: 
%   v - A T x 2 matrix representing the inverse exponential mapping (the 
%       shooting vector) from q1 to an optimally aligned q2. v is in the 
%       tangent space of the unit hypersphere in L^2([0,1],R^2) function 
%       space. 
%   dist - The geodesic distance (arc length) between q1 and the optimally 
%       aligned q2. This value is also equal to the norm of v. 
%   beta2n - The version of beta2 that is optimally aligned to beta1. 

[n,p]=size(beta1);
beta1=standardizeCurve(beta1);
beta2=standardizeCurve(beta2);
q1=curve_to_q(beta1);

% Optimize over O(n) 
[beta2,~]=ProcrustesOpt(beta1,beta2);
q2=curve_to_q(beta2);

if register 
    % Optimize over Gamma
    gam=DynamicProgrammingQ(q1',q2',0,1)';
    gamI=invertGamma(gam);
    beta2=Group_Action_by_Gamma_Coord(beta2,gamI);

    % Optimize over O(n) again
    [beta2n,~]=ProcrustesOpt(beta1,beta2);
    q2n=curve_to_q(beta2);
else
    beta2n=beta2;
    q2n=q2;
end

% Compute geodesic distance
q1dotq2=min([1 InnerProd_Q(q1,q2n)]);
dist=real(acos(q1dotq2)); 

% Compute shooting vector on sphere
u=q2n-q1dotq2*q1;
normu=sqrt(InnerProd_Q(u,u));
if normu>10^-10
    v=u*real(acos(q1dotq2))/normu;
else
    v=zeros(n,p);
end
