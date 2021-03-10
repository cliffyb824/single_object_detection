function J = calculateJ(basis)

b1=basis{1};
b2=basis{2};
T=length(b1');
t=linspace(0,1,T);
integrand11=zeros(1,T);
integrand12=zeros(1,T);
integrand22=zeros(1,T);
for i=1:T
    integrand11(i)=b1(i,:)*b1(i,:)';
    integrand12(i)=b1(i,:)*b2(i,:)';
    integrand22(i)=b2(i,:)*b2(i,:)';
end
J(1,1)=trapz(t,integrand11);
J(1,2)=trapz(t,integrand12);
J(2,2)=trapz(t,integrand22);
J(2,1)=J(1,2);

