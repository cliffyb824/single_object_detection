function [Sbetamean,E,normgrad] = findKarcherMean(Data,Strain,beta0)

% Compute the Karcher mean on closed curve elastic shape space. The shape 
% manifold is described in Chapter 6 of the text by Srivastava and Klassen 
% (2016), and the Karcher mean algorithm is given in Chapter 7. Here, we 
% implement a gradient descent method with Armijo backtracking to minimize 
% the objective function provided in Energy_KarcherMean.m. 
%
% Inputs:
%   Data - A 1 x N cell array where Data{i} is the ith curve \beta_i 
%       represented as an T x 2 matrix. 
%   Strain - A structure array consisting of algorithm options and
%       parameters
%   beta0 - The initializing curve, a T x 2 matrix. 
%
% Outputs:
%   Sbetamean - The final estimate of the mean curve (T x 2 matrix) in its 
%       curve_properties structure array. 
%   E - The evolution of the energy functional over the iterations of the
%       algorithm (1 x numIters array). 
%   normgrad - The evolution of the norm of the gradient over the
%       iterations of the algorithm (1 x numIters array). 

backtracking=Strain.backtracking;
register=Strain.register;
delta0=Strain.delta0;
deltamin=Strain.deltamin;
tol=Strain.tol;
maxit=Strain.maxit;
plotevol=Strain.plotevol;
verbose=Strain.verbose;

% Hard coded variables for Armijo backtracking
c=0.1;
tau=0.5;
% v = VideoWriter('123.avi');
% v.FrameRate = 1;
% open(v);
% Initialize
i=1;
E=zeros(maxit,1);
normgrad=zeros(maxit,1);
T=length(beta0');
t=linspace(0,1,T);
betamean=standardizeCurve(beta0);
mu=curve_to_q(betamean);
mu=mu/sqrt(InnerProd_Q(mu,mu)); 
basis=findBasisNormal(mu);
if plotevol
    fig=10;
    plotCurve(betamean,fig,2,0,[17,44]);
    title(['Iteration: ' num2str(i)]);
end
[E(i),gradE,normgrad(i)]=Energy_KarcherMean(betamean,Data,mu,basis,register);
i=i+1;
delta=delta0;
maxitj=ceil(log(deltamin/delta0)/log(tau));

% Compute the Karcher mean
while normgrad(i-1)>tol && i<=maxit && delta>deltamin
    
    if verbose
        disp([num2str(i-1) ': ' num2str(E(i-1))]);
    end
    
%     % Re-Initialize stepsize. May not be necessary for this manifold.
%     delta=delta0;
    
    % Advance mu and betamean
    mutry=expMapping(gradE,mu,delta);
    mutry=projectCurve(mutry);
    basis=findBasisNormal(mutry);
    betameantry=center_curve(q_to_curve(mutry));
    
    % Compute E and gradE
    [Etry,gradEtry,normgradtry]=Energy_KarcherMean(betameantry,Data,mutry,basis,register);
    
    % Armijo-Goldstein Condition
    ag_cond=E(i-1)-Etry;
    
    % Armijo backtracking
    j=0;
    while backtracking && ag_cond<delta*c*normgrad(i-1)^2 && j<maxitj 
        j=j+1;

        % Minimize step size
        delta=tau*delta; 
        
        % Advance mu and betamean
        mutry=expMapping(gradE,mu,delta);
        mutry=projectCurve(mutry);
        basis=findBasisNormal(mutry);
        betameantry=center_curve(q_to_curve(mutry));

        % Compute E and gradE
        [Etry,gradEtry,normgradtry]=Energy_KarcherMean(betameantry,Data,mutry,basis,register);
        
        % Armijo-Goldstein Condition
        ag_cond=E(i-1)-Etry;
    end
           
    if verbose
        disp(['delta = ' num2str(delta)])
    end
    
    % Update
    mu=mutry;
    betamean=betameantry;
    gradE=gradEtry;
    E(i)=Etry;
    normgrad(i)=normgradtry;
    i=i+1;
    
    if plotevol
        plotCurve(betamean,fig,2,0,[17,44]);
        title(['Iteration: ' num2str(i)]);
%         frame = getframe(gcf);
%         writeVideo(v,frame);
    end
    
end

if verbose
    disp([num2str(i-1) ': ' num2str(E(i-1))]);
end
E=E(1:i-1);
normgrad=normgrad(1:i-1);

if plotevol
    plotCurve(betamean,fig,2,0,[17,44]);
    close;
end

% Return the curve properties structure of a uniformly sampled mean curve
Sbetamean=curve_properties(ReSampleCurve(betamean,T),t);
% close(v);