function b0 = compute_similarity_basis(Sbeta)

T=Sbeta.T;
t=Sbeta.t;
beta=Sbeta.beta;
b{1}=[ones(T,1) zeros(T,1)];
b{2}=[zeros(T,1) ones(T,1)];
b{3}=beta*[0 -1; 1 0];
b{4}=beta; 
b0=gramSchmidt(b,t);


