function betadot = grad_closed_curve(beta)

% Centered difference. Assumed closed curve with beta(:,1)=beta(:,end)

T=length(beta);
dt=1/(T-1);
betadot=zeros(2,T);
betadot(:,1)=(beta(:,2)-beta(:,end-1))/(2*dt);
betadot(:,2:T-1)=(beta(:,3:T)-beta(:,1:T-2))/(2*dt);
betadot(:,T)=betadot(:,1);
