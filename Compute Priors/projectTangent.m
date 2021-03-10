function wproj = projectTangent(w,q,basis)

% w is any vector, q is srvf of closed curve, basis is orthonormal basis of
% normal space to closed curve manifold at q

w=w-InnerProd_Q(w,q)*q; % First project to tangent space of open curve space
wproj=w-InnerProd_Q(w,basis{1})*basis{1}-InnerProd_Q(w,basis{2})*basis{2};

