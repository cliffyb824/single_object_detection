function Sbeta_new = advance_curve(Sbeta,delta,gradEtotal,resample_iter)

% Advance curve via gradient descent

beta_new=Sbeta.beta-delta*gradEtotal;
if resample_iter
    beta_new=ReSampleCurve(beta_new,Sbeta.T);
end
Sbeta_new=curve_properties(beta_new,Sbeta.t);