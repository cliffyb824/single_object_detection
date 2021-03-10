function [betahat,Eevol,error] = active_contour_agm(Sopt,Sprior)

% Experimental optimization algorithm: Accelerated Gradient Method. This is
% the theoretically fastest converging first-order optimization method.
% Code needs to be modified to avoid causing self-intersecting curves,
% (loops) which cause the solution to blow up. 

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
SbetaY_prev=Sbeta;
curves{1}=Sbeta.beta;
i=1;
[Eevol(i),Sprior,gradEanalytical]=Etotal(Sbeta,Sopt,Sprior,true,[]);
error=zeros(maxit,1);
error(1)=tol+1;
lam=0;
need_numerical_gradient=false;
if sum(lambda(gradtype==1))>0
    need_numerical_gradient=true;
end

% Active contour accelerated gradient method loop (target and shadow curves, alternating)
while i<maxit %&& error(i)>tol 
    
    % Plot the active contour evolution
    if toggleplot
        plot_curves_on_image(curves,exp(Sopt.I),Sopt.nrow,Sopt.scalefac,1,2,Sopt.cmax,'r');
        title(['Iter: ' num2str(i)]);
    end
    
    if display_iter
        disp([num2str(i) ': ' num2str(Eevol(i))]);
    end
 
    % Approximate the gradient of curve beta1 with respect to the
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

    % Advance beta according to approx gradient descent with momentum
    SbetaY_cur=advance_curve(Sbeta,delta,gradEtotal,resample_iter);
    lam=(1+sqrt(1+4*lam^2))/2;
    gam=(1-lam)/(lam+1);
    beta=(1-gam)*SbetaY_cur.beta+gam*SbetaY_prev.beta;
    if resample_iter
        beta=ReSampleCurve(beta,Sopt.T);
    end
    Sbeta=curve_properties(beta,t);
        
    % Update iteration count, curves, energy, and error
    i=i+1;
    [Eevol(i),Sprior,gradEanalytical]=Etotal(Sbeta,Sopt,Sprior,true,[]);
    error(i-1)=abs(Eevol(i)-Eevol(i-1))/Eevol(i-1);
    curves{1}=Sbeta.beta;
    SbetaY_prev=SbetaY_cur;
%     pause;
    
end

% Plot the active contour evolution
if toggleplot
    plot_curves_on_image(curves,exp(Sopt.I),Sopt.nrow,Sopt.scalefac,1,2,Sopt.cmax,'r');
    title(['Iter: ' num2str(i)]);
end

if display_iter
    disp([num2str(i) ': ' num2str(Eevol(i))]);
end

error=error(1:i-1);
betahat=Sbeta.beta;
