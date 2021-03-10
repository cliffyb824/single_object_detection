function gradEtotal = compute_grad_Etotal_numerical(Sbeta,Sopt,Sprior)

fbasis=Sopt.fbasis;
epsilon=Sopt.epsilon;
J=length(fbasis);
T=length(Sbeta.beta');
gradEtotal=zeros(T,2);
register=false; % No need to call dynamic programming for each small curve perturbation used in gradient approx 
filter_by_gradtype='numerical'; 
Etotal_cur=Etotal(Sbeta,Sopt,Sprior,register,filter_by_gradtype);

for j=1:J
    beta_pert=Sbeta.beta+epsilon*fbasis{j};
    Sbeta_pert=curve_properties(beta_pert,Sbeta.t);
    Etotal_pert=Etotal(Sbeta_pert,Sopt,Sprior,register,filter_by_gradtype);
    partial=(Etotal_pert-Etotal_cur)/epsilon;
    gradEtotal=gradEtotal+partial*fbasis{j};
end

% % For debugging:
% keyboard;
% beta=Sbeta.beta;
% figure(1); clf; hold on;
% plot(beta(:,1),beta(:,2),'b','LineWidth',2); 
% quiver(beta(:,1),beta(:,2),-gradEtotal(:,1),-gradEtotal(:,2),'b'); axis equal off;