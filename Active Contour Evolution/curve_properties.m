function Sbeta = curve_properties(beta,t)

% Computed several relevant properties of a parameterized curve beta and
% store in a structure array Sbeta. 

beta=beta'; 
t=t'; 
T=length(t); 
betadot=grad_closed_curve(beta); 
normbetadot=sqrt(sum(betadot.*betadot)); 
% betadotdot=grad2_closed_curve(beta);  
L=trapz(t,normbetadot); 
c=trapz(t,beta*diag(normbetadot),2)/L; 
betaStd=(beta-repmat(c,1,T))/L; 
q=betadot*diag(1./sqrt(normbetadot)); 
q=q/sqrt(inner_prod_q(q,q,t)); 
u=betadot./repmat(normbetadot,2,1); 
n(1,:)=u(2,:); 
n(2,:)=-u(1,:);

% Calculate signed curvature
udot=grad_closed_curve(u); 
% udot=betadotdot./repmat(normbetadot,2,1); % Check to see if this line yields same udot as above line
s=sum(udot.*n);
normudot=sqrt(sum(udot.*udot)); 
kappa=-sign(s).*normudot;

% % Alternative numerical approximation of curvature. Still need to further 
% % investigate this method versus using normudot for curvature.
% curvature=zeros(1,T);
% x1=beta(1,T-1);
% x2=beta(1,1);
% x3=beta(1,2);
% y1=beta(2,T-1);
% y2=beta(2,1);
% y3=beta(2,2);
% a=sqrt((x1-x2)^2+(y1-y2)^2); % The three sides
% b=sqrt((x2-x3)^2+(y2-y3)^2);
% c=sqrt((x3-x1)^2+(y3-y1)^2);
% A=1/2*abs((x1-x2)*(y3-y2)-(y1-y2)*(x3-x2)); % Area of triangle
% curvature(1)=4*A/(a*b*c); % Curvature of circumscribing circle
% curvature(T)=curvature(1);
% for i=2:T-1
%     x1=beta(1,i-1);
%     x2=beta(1,i);
%     x3=beta(1,i+1);
%     y1=beta(2,i-1);
%     y2=beta(2,i);
%     y3=beta(2,i+1);
%     a=sqrt((x1-x2)^2+(y1-y2)^2); % The three sides
%     b=sqrt((x2-x3)^2+(y2-y3)^2);
%     c=sqrt((x3-x1)^2+(y3-y1)^2);
%     A=1/2*abs((x1-x2)*(y3-y2)-(y1-y2)*(x3-x2)); % Area of triangle
%     curvature(i)=4*A/(a*b*c); % Curvature of circumscribing circle
% end

% Form structure array of curve properties
Sbeta.beta=beta';                       % The curve
Sbeta.betadot=betadot';                 % Tangent vector (first deriv) at each sample point
% Sbeta.betadotdot=betadotdot';           % Second derivative vector at each sample point
Sbeta.normbetadot=normbetadot';         % Norm of tangent vector at each sample point
Sbeta.length=L;                         % Length of the curve
Sbeta.centroid=c;                       % Centroid of the curve
Sbeta.T=T;                              % Number of sample points on the curve
Sbeta.t=t';                             % Uniformly sampled time points from 0 to 1
Sbeta.betaStd=betaStd';                 % Standardized version of the curve (centroid at 0 and unit length)
Sbeta.q=q';                             % Normalized square-root velocity function (SRVF) of the curve
Sbeta.n=n';                             % Unit normal vector at each sample point
Sbeta.kappa=kappa';                     % Curvature (signed)
Sbeta.curvature=normudot';              % Curvature
% Sbeta.curvature=curvature';


    