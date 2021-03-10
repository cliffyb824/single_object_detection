function [E,beta2new] = Eshape(Sbeta,Sbetamean,register)

beta1=Sbeta.betaStd;
beta2=Sbetamean.betaStd;
q1=Sbeta.q;
q2=Sbetamean.q;
t=Sbeta.t;

% Optimize over O(n) 
[beta2,O]=ProcrustesOpt(beta1,beta2);
q2=q2*O;

if register 
    % Optimize over Gamma
    gam=DynamicProgrammingQ(q1',q2',0,1)';
    gamI=invertGamma(gam);
    beta2=spline(t,beta2',gamI')';

    % Optimize over O(n) again
    [beta2new,O]=ProcrustesOpt(beta1,beta2);
    q2=q2*O;
else
    beta2new=beta2;
end

% Compute geodesic distance
q1dotq2=min([1 inner_prod_q(q1',q2',t)]);
E=real(acos(q1dotq2))^2; 
% keyboard;