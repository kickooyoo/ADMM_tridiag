function cost = ADMM_tridiag_cost(y, F, S, CH, CV, beta, x)
%function cost = ADMM_tridiag_cost(y, F, S, CH, CV, beta, x)
% 
% inputs:
%       y [Ns*Nc 1] undersampled data vector
%       F (fatrix2) undersampled Fourier encoding matrix
%       S (fatrix2) sensitivity map diagonal matrix
%       CH (fatrix2) horizontal finite differences
%       CV (fatrix2) vertical finite differences
%       beta (real scalar) spatial regularization parameter
%       x [Nx Ny] initial guess for x
% outputs:
%       cost [1] objective value of original cost func
%
% 01/26/2016 Mai Le
% University of Michigan

cost = norm(col(y) - col(F * (S * x)),2)^2/2 + beta * norm(col(CH * x),1) + ...
        beta * norm(col(CV * x),1);

end

