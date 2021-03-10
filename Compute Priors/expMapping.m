function q2 = expMapping(v,q1,delta)

% Compute the exponential mapping on the unit hypersphere in L2([0,1],R^2) 
% space from q1 in the direction of the scaled shooting vector \delta*v. 
%
% Inputs: 
%   v - A T x 2 matrix representing a tangent vector to the hypersphere.
%   q1 - A T x 2 matrix representing an element of the hypersphere.
%   delta - A scalar in (0,1].
%
% Outputs: 
%   q2 - A T x 2 matrix representing the above exponential mapping.

normv=sqrt(InnerProd_Q(v,v));
if normv<10^-10 % Prevent from dividing by zero (hard-coded tolerance). 
    q2=q1;
else
    q2=cos(delta*normv)*q1+sin(delta*normv)*v/normv;
end

    
