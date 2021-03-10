function E = Esmooth(Sbeta)

Sbeta=curve_properties(Sbeta.betaStd,Sbeta.t); % Use unit length curve to compute smoothness penalty
E=trapz(Sbeta.t,Sbeta.curvature.^2);

