function basis = findBasisNormal(q)

% Return basis vectors for normal space of closed curve manifold at q. 
% "basis" is a cell array of size 2 with each entry being a Tx2 matrix.  

q=q';
T=length(q);
t=linspace(0,1,T);
normq=sqrt(sum(q.*q));
qq=q./[normq;normq];
f1=qq*diag(q(1,:))+[normq;zeros(1,T)];
f2=qq*diag(q(2,:))+[zeros(1,T);normq];
I1=sum(q.*f1);
I2=sum(q.*f2);
b1=f1-q*trapz(t,I1);
b2=f2-q*trapz(t,I2);
basis{1}=b1';
basis{2}=b2';

