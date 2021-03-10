function val = inner_prod_q(q1,q2,t)

val=trapz(t,sum(q1.*q2));