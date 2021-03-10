function basis_o = gramSchmidt2(basis)

b1=basis{1};
b2=basis{2};
basis1=b1/sqrt(InnerProd_Q(b1,b1));
b2=b2-InnerProd_Q(basis1,b2)*basis1;
basis2=b2/sqrt(InnerProd_Q(b2,b2));
basis_o{1}=basis1;
basis_o{2}=basis2;
 