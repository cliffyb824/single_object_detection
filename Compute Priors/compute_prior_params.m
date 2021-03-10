function Sprior = compute_prior_params(Strain)

% Unpack structure array
filename=Strain.filename;
options=Strain.options;
idxtrain=Strain.idxtrain;

% Load training data and downselect according to training indices idxtrain
load(filename); 
Sprior.snips=snips(:,:,idxtrain);
Sprior.Ctrue=Ctrue(:,:,idxtrain);
betatrain=betatrain(idxtrain);
Sprior.scalefac=scalefac;

% Initialize variables for training algorithms
ntrain=length(idxtrain);
x0=[];
x1=[];
L=zeros(ntrain,1);
T=length(betatrain{1}');
t=linspace(0,1,T);
Data{ntrain}=[];
centroid=zeros(ntrain,2);

disp('Collecting Data...');
for i=1:10%ntrain
    
    % Image
    I=snips(:,:,i);
    
    % Compute centroid
    centroid(i,:)=calculateCentroid(betatrain{i});
    
    % Collect log-transformed pixels for each region 
    x0=[x0; log(double(I(Ctrue(:,:,i)==0)))];
    x1=[x1; log(double(I(Ctrue(:,:,i)==1)))];
    
    % Compute curve properties of training curves
    Sprior.Sbetatrain(i)=curve_properties(betatrain{i},t);
    Data{i}=Sprior.Sbetatrain(i).betaStd;
    
    % Assign curve length vectors
    L(i)=Sprior.Sbetatrain(i).length; 
    
end

Sprior.muc=mean(centroid);

% Pixel distribution MLEs (skew-logistic) for background and target regions
disp('Computing Pixel Distribution MLEs...');
thetainit=[0,1,0];
% Sprior.theta0=fminunc(@(theta) -skewLogisticLikelihood(theta,x0),thetainit,options);
Sprior.theta0=fminunc(@(theta) -skewLogisticLikelihood(theta,x0),thetainit,options);
thetainit=[1,1,1];
% Sprior.theta1=fminunc(@(theta) -skewLogisticLikelihood(theta,x1),thetainit,options);
Sprior.theta1=fminunc(@(theta) -skewLogisticLikelihood(theta,x1),thetainit,options);

% Boundary length distribution MLEs (gaussian)
disp('Computing boundary length MLEs...');
Sprior.muL=mean(L);
Sprior.sigL=std(L);

% Target and shadow boundary curve means
disp('Computing Target Karcher Mean...');
beta0=Data{1};
[Sprior.Sbetamean,~,~]=findKarcherMean(Data,Strain,beta0);

disp('Done');




