function qProj = projectCurve(q)

plotstuff=0;
tol=10^-5;
maxit=200;
iter=0;
delta=0.5;
x=q_to_curve(q);

% Algorithm for projection onto closed curve space
r=[x(end,1);x(end,2)];
rnorm=zeros(1,maxit);
rnorm(1)=norm(r);
while (norm(r)>tol) && (iter<maxit)
    
    % Find basis for normal space at q
    basis=findBasisNormal(q);
    
    % Calculate Jacobian
    J=calculateJ(basis);
    
    % Newton-Raphson step to update q
    y=J\(-r);
    dq=delta*(y(1)*basis{1}+y(2)*basis{2});
    normdq=sqrt(InnerProd_Q(dq,dq));
    q=cos(normdq)*q+sin(normdq)*dq/normdq;
    q=q/sqrt(InnerProd_Q(q,q));
    
    % Update x and r from the new q. 
    x=q_to_curve(q);
    r=[x(end,1);x(end,2)];
    rnorm(iter+1)=norm(r);
    iter=iter+1;
    
    % Plot evolution
    if plotstuff
        figure(10);
        plot(x(:,1),x(:,2)); axis equal off;
    end

end

if plotstuff
    rnorm=rnorm(1:iter-1);
    figure; plot(rnorm)
end

qProj=q;