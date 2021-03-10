function fbasis = create_basis(d,T)

t=linspace(0,1,T)';
fbasis{1,4*d+2}=[];
fbasis{1}=[ones(T,1) zeros(T,1)];
fbasis{2}=[zeros(T,1) ones(T,1)];
i=3;
for n=1:d
    fbasis{i}=[sqrt(2)*sin(2*n*pi*t) zeros(T,1)];
    i=i+1;
    fbasis{i}=[zeros(T,1) sqrt(2)*sin(2*n*pi*t)];
    i=i+1;
    fbasis{i}=[sqrt(2)*cos(2*n*pi*t) zeros(T,1)];
    i=i+1;
    fbasis{i}=[zeros(T,1) sqrt(2)*cos(2*n*pi*t)];
    i=i+1;
end