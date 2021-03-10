function [E,Sprior,gradEbeta] = Etotal(Sbeta,Sopt,Sprior,register,filter_by_gradtype)

% Etotal is the weighted sum of the following energy functionals: 
% 
% E1: image energy 
% E2: shape prior
% E3: length prior
% E4: smoothness prior
%
% where lambda represents the vector of weights.

lambda=Sopt.lambda;
t=Sopt.t;
T=Sopt.T;
E=0;
gradEbeta=zeros(T,2);
if strcmp(filter_by_gradtype,'analytical')
    lambda(Sopt.gradtype==1)=0;
elseif strcmp(filter_by_gradtype,'numerical')
    lambda(Sopt.gradtype==0)=0;
end
        
% Compute the different terms in the total energy

if lambda(1)~=0 
    [E1,gradEpix]=Eimage(Sbeta,Sopt,Sprior);
    E=E+lambda(1)*E1;
    gradEbeta=gradEbeta+lambda(1)*gradEpix;
end

% The shape energy function returns an updated mean shape that has been
% rotated and registered (if register==true) to the current active contour
if lambda(2)~=0
    [E2,betamean_new]=Eshape(Sbeta,Sprior.Sbetamean,register);
    E=E+lambda(2)*E2;
    Sprior.Sbetamean=curve_properties(betamean_new,t);
end

if lambda(3)~=0
    [E3,gradElen]=Elength(Sbeta,Sprior);
    E=E+lambda(3)*E3;
    gradEbeta=gradEbeta+lambda(3)*gradElen;
end

if lambda(4)~=0 
    E4=Esmooth(Sbeta);
    E=E+lambda(4)*E4;
end


    
    
    
    
    
    