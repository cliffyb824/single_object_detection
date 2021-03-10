function beta0 = find_best_rotation_for_init(Sprior,Sinit)

ntheta=Sinit.ntheta;

% Shape prior mean with average length, centered at origin
beta0c=Sprior.muL*Sprior.Sbetamean.betaStd;

% Optimize over rotation when located at average centroid
Sinit.lambda(2:4)=0; % Hard coded lambda indices
theta=linspace(0,2*pi,ntheta);
E=zeros(ntheta,1);
for i=1:ntheta
    R=rotMat(theta(i));
    beta0try=beta0c*R+repmat(Sprior.muc,Sinit.T,1);
    Sbeta0try=curve_properties(beta0try,Sinit.t);
    E(i)=Etotal(Sbeta0try,Sinit,Sprior,false,[]);
end
[~,idx]=min(E);
beta0=beta0c*rotMat(theta(idx))+repmat(Sprior.muc,Sinit.T,1);



