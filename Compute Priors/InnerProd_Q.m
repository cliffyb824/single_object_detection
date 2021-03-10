function val = InnerProd_Q(q1,q2)

% Compute the L^2([0,1],R^2) inner product given by 
% <q_1,q_2> = \int_0^1 q_1(t)q_2(t)dt using trapezoidal integration. 
% 
% Inputs: 
%   q1,q2 - Two T x 2 matrices representing functions (srvf's) in 
%       L^2([0,1],R^2)
%
% Outputs: 
%   val - The inner product <q_1,q_2>    

q1=q1';
q2=q2';
T=length(q1);
val=trapz(linspace(0,1,T),sum(q1.*q2));