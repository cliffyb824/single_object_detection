function [betahat,Eevol,error] = active_contour_gradient_descent(Sopt,Sprior)

% Optimize total energy functional via active contour gradient descent
% according to the settings contained in the Sopt structure and the prior
% parameter values contained in the Sprior structure. The method is hard
% coded for one active contour in the image.
% 
% Inputs: 
%   Sopt - Structure array containing algorithm options and settings
%   Sprior - Structure array containing parameter values of prior
%       distributions
%
% Outputs: 
%   betahat - Optimal boundary curve (T x 2 matrix)
%   Eevol - Etotal computed at each iteration of the contour evolution
%   error - Error term computed at each iteration of the contour evolution

% Unpack structure array
beta0=Sopt.beta0;
delta=Sopt.delta;
maxit=Sopt.maxit;
tol=Sopt.tol;
toggleplot=Sopt.toggleplot;
display_iter=Sopt.display_iter;
resample_iter=Sopt.resample_iter;
T=Sopt.T;
t=Sopt.t;
gradtype=Sopt.gradtype;
lambda=Sopt.lambda;

% Initialize active contour algorithm
Eevol=zeros(1,maxit);
Sbeta=curve_properties(beta0,t);
curves{1}=Sbeta.beta;
i=1;
[Eevol(i),Sprior,gradEanalytical]=Etotal(Sbeta,Sopt,Sprior,true,[]);
error=zeros(maxit,1);
error(1)=tol+1;
need_numerical_gradient=false;
if sum(lambda(gradtype==1))>0
    need_numerical_gradient=true;
end

% Active contour gradient descent loop (target and shadow curves, alternating)
while i<maxit %&& error(i)>tol 
    
    % Plot the active contour evolution on top of the image
    if toggleplot
        plot_curves_on_image(curves,exp(Sopt.I),Sopt.nrow,Sopt.scalefac,1,2,Sopt.cmax,'r');
        title(['Iter: ' num2str(i)]);
    end
    frame = getframe(gcf);
    writeVideo(u,frame);
    
    % Display the iteration count and current energy in the command window
    if display_iter
        disp([num2str(i) ': ' num2str(Eevol(i))]);
    end
 
    % Approximate the gradient of target curve beta1 with respect to the
    % terms that require numerical gradient computation
    if need_numerical_gradient
        gradEnumerical=compute_grad_Etotal_numerical(Sbeta,Sopt,Sprior);
    else
        gradEnumerical=zeros(T,2);
    end
    
    % Combine numerical gradient with analtyical gradient
    gradEtotal=gradEnumerical+gradEanalytical;
    
    % Project to the normal direction of the curve (optional)
    gradEtotal=diag(sum(gradEtotal.*Sbeta.n,2))*Sbeta.n;

    % Advance beta according to approx gradient descent, and compute new
    % curve properties
    Sbeta=advance_curve(Sbeta,delta,gradEtotal,resample_iter);
    
    % Update iteration count, energy, and error
    i=i+1;
    [Eevol(i),Sprior,gradEanalytical]=Etotal(Sbeta,Sopt,Sprior,true,[]);
    error(i-1)=abs(Eevol(i)-Eevol(i-1))/Eevol(i-1); % Maybe this needs to be the variance of the relative error within some window
    curves{1}=Sbeta.beta;
    
end

% Plot the final solution on top of the image
if toggleplot
    plot_curves_on_image(curves,exp(Sopt.I),Sopt.nrow,Sopt.scalefac,1,2,Sopt.cmax,'r');
    title(['Iter: ' num2str(i)]);
end

% Display the iteration count and final energy in the command window
if display_iter
    disp([num2str(i) ': ' num2str(Eevol(i))]);
end

error=error(1:i-1);
betahat=Sbeta.beta;


