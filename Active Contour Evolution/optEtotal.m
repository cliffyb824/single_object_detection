function [betahat,Eevol,error] = optEtotal(Sopt,Sprior)

% Optimize Etotal
if strcmp(Sopt.algorithm,'grad') % Gradient descent
    
    [betahat,Eevol,error]=active_contour_gradient_descent(Sopt,Sprior); 
    
elseif strcmp(Sopt.algorithm,'agm') % Accelerated gradient method
    
    [betahat,Eevol,error]=active_contour_agm(Sopt,Sprior);
    
elseif strcmp(Sopt.algorithm,'none') % Just use the initialization
    
    betahat=Sopt.beta0;
    Sbeta=curve_properties(betahat,Sopt.t);
    Eevol=Etotal(Sbeta,Sopt,Sprior,true,[]);
    
end