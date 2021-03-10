function [X2new,O] = ProcrustesOpt(X1,X2)

% Optimally rotate/reflect the point set X2 to another point set X1 using 
% Procrustes alignment. 
%
% Inputs:
%   X1 and X2 - Two n x p matrices. 
% Outputs: 
%   X2new - The optimally rotated/reflected version of X2.
%   O - The optimal member of O(p), the orthogonal group, that aligns X2 to
%       X1. 

[O1,~,O2] = svd(X2'*X1);
O=O1*O2';
X2new=X2*O;

