function betadotdot = grad2_closed_curve(beta)

% Second derivative. Assumed closed curve with beta(:,1)=beta(:,end)

T=length(beta);
dt=1/(T-1);
betadotdot=zeros(2,T);
betadotdot(:,1)=(beta(:,2)-2*beta(:,1)+beta(:,end-1))/dt^2;
betadotdot(:,2:T-1)=(beta(:,3:T)-2*beta(:,2:T-1)-beta(:,1:T-2))/dt^2;
betadotdot(:,T)=betadotdot(:,1);
