function b0 = gramSchmidt(b,t)

n=length(b);
b0{n}=[];
b0{1}=b{1}/sqrt(inner_prod_q(b{1}',b{1}',t));
if n>1
    for i=2:n
        s=0;
        for j=1:i-1
            s=s+inner_prod_q(b{j}',b{i}',t)*b{j};
        end
        b0{i}=b{i}-s;
        b0{i}=b0{i}/sqrt(inner_prod_q(b0{i}',b0{i}',t));
    end
end
