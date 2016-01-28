function [x, xsaved, err, cost, time] = tridiag_ADMM(y, F, S, C1, C2, ...
        alph, beta, xinit, xtrue, niters, varargin)
% function [x, xsaved, err, cost, time] = tridiag_ADMM(y, F, S, C1, C2, ...
%         alph, beta, mu, Nx, Ny, xinit, xtrue, niters, varargin)
% 
% implements tridiag_ADMM algo
% inputs:
%       y [Ns*Nc 1] undersampled data vector
%       F (fatrix2) undersampled Fourier encoding matrix
%       S (fatrix2) sensitivity map diagonal matrix
%       C1 (fatrix2) horizontal finite differences
%       C2 (fatrix2) vertical finite differences
%       alph (scalar, \in [0, 1]) parameter balancing between u3 and x
%       beta (real scalar) spatial regularization parameter
%       xinit [Nx Ny] initial guess for x
%       xtrue [Nx Ny] true solution, necessary for NRMSE comparison
%       niters (integer) number of iterations
% varargin:
%       mu [5 1] AL tuning parameters, default all 1
%       mask (logical) [Nx Ny] mask for NRMSE calculation
% outputs:
%       x [Nx Ny] reconstructed image
%       xsaved [Nx Ny niters] estimated image at each iteration
%       err [niters 1] NRMSE at each iteration
%       cost [niters 1] objective value of original cost func at each iter
%       time [niters 1] wall time per each iteration
%
% 01/26/2016 Mai Le
% University of Michigan

Nx = size(xinit,1);
Ny = size(xinit,2);

Nc = S.arg.Nc;
arg.compile_mex = false;
arg.mask = true(Nx, Ny); %generate_mask('slice67',1,Nx,Ny);
arg.mu = ones(1,5);
arg = vararg_pair(arg, varargin);

iter = 0;
x = xinit(:);
u0 = C1*x;
u3 = x;
u1 = C2*u3;
u2 = (1-alph)*S*u3+alph*S*x;
v3 = zeros(size(u3));
v4 = -v3;
eta0 = zeros(size(u0));
eta1 = zeros(size(u1));
eta2 = zeros(size(u2));
eta3 = zeros(size(u3));
eta4 = zeros(size(x));
% renaming to match paper indeces
mu0 = arg.mu(1);
mu1 = arg.mu(2);
mu2 = arg.mu(3);
mu3 = arg.mu(4);
mu4 = arg.mu(5);

[eig_FF,Qbig] = get_eigs(F,Nc);
SS = S'*S;
eig_SS = reshape(SS*ones(prod(S.idim),1), Nx, Ny); % diagonal, not BCCB

% pass tridiag of CC into mex
subCCT = single(-mu1 * ones(Ny - 1, Nx));
diagCCT = single(mu1 * cat(1, ones(1, Nx), 2*ones(Ny-2, Nx), ones(1, Nx)));
subCC = single(-mu0 * ones(Nx - 1, Ny));
diagCC = single(mu0 * cat(1, ones(1, Ny), 2*ones(Nx-2, Ny), ones(1, Ny)));

% to do: make sure mex compiled
nthread = int32(jf('ncore'));
if arg.compile_mex
        curr = cd;
        if ~strcmp(curr(end-11:end), 'ADMM_tridiag')
                display('must be in ADMM_tridiag dir to compile tridiag_inv_mex_noni.c');
                keyboard;
        end
        mex -O CFLAGS="\$CFLAGS -std=c99 -DMmex" -I./pthread_tutor/def/ ./pthread_tutor/tridiag_inv_mex_noni.c
end

time = zeros(niters,1);
err(1) = calc_NRMSE_over_mask(x, xtrue, arg.mask);
cost(1) = calc_cost(beta,C1,C2,F,S,y,x);

while(iter < niters)
        iter_start = tic;
        u0 = soft(C1 * double(x) - eta0, beta/mu0);
        u1 = soft(C2 * double(u3) - eta1, beta/mu1);
        u2 = u2_update(mu2, alph, eig_FF, Qbig, F, S, y, u3, x, eta2, Nx, Ny);
        u3 = u3_update_mex(mu1, mu2, mu3, alph, eig_SS, C2, S, u1, u2, x, v3, eta1, eta2, eta3, Nx, Ny, subCCT, diagCCT, nthread);
        x = x_update(mu0, mu2, mu4, alph, eig_SS, C1, S, u0, u2, u3, v4, eta0, eta2, eta4, Nx, Ny, subCC, diagCC, nthread);
        
        % skip v0, v1, v2 because they are constrained to be zero
        v3 = (mu3*(-u3 - eta3) + mu4*(-x + eta4)) ./ (mu3+mu4);
        v4 = -v3;
        
        % eta updates
        eta0 = eta0 - (-u0 + C1 * x);
        eta1 = eta1 - (-u1 + C2 * u3);
        eta2 = eta2 - (-u2 + (1 - alph) * S * u3 + alph * S * x);
        eta3 = eta3 - (-u3 - v3);
        eta4 = eta4 - (x - v4);
        
        iter = iter+1;
        iter_time = toc(iter_start);
        time(iter) = iter_time;
        err(iter) = calc_NRMSE_over_mask(x, xtrue, arg.mask);
        
        if mod(iter,10) == 0
                printf('%d/%d iterations',iter,niters)
        end
        xsaved(:,:,iter) = reshape(x,Nx,Ny);
        cost(iter+1) = calc_cost(beta, C1, C2, F, S, y, x);
end

end

function cost = calc_cost(beta, C1, C2, F, S, y, x)
cost = norm(col(y) - col(F * (S * x)),2)^2/2 + beta * norm(C1 * x,1) + ...
        beta * norm(C2 * x,1);
end

function out = soft(in,thresh)
out = (in - thresh * sign(in)) .* (abs(in) > thresh);
end

function u2 = u2_update(mu2, alph, eig_FF, Q, F, S, y, u3, x, eta2, Nx, Ny)
arg_u2 = F'*y + mu2 * ((1-alph) * S * u3 + alph * S * x - eta2);
u2 = Q'*((Q*arg_u2)./(eig_FF+mu2))/(Nx*Ny);
end

function u3 = u3_update_mex(mu1, mu2, mu3, alph, eig_SS, C2, S, u1, u2, x, v3, eta1, eta2, eta3, Nx, Ny, subCCT, diagCCT, nthread)
u3arg = mu1 * C2' * (u1 + eta1) + mu2 * (1-alph) * S' * (u2 - alph * S * x + eta2) + ...
        mu3 * (-v3 - eta3);

% transpose to make Hessian tridiagonal, size now Ny Nx
flipu3arg = single(reshape(u3arg, Nx, Ny));
flipu3arg = flipu3arg.';
flipSS = eig_SS.';

diagvals = single(diagCCT + mu2 * (1-alph)^2 * flipSS + mu3);
u3out = tridiag_inv_mex_noni(subCCT, diagvals, subCCT, flipu3arg, nthread);

% transpose solution back
flipu3 = reshape(u3out, Ny, Nx); 
u3 = flipu3.';
u3 = u3(:);
end

function x = x_update(mu0, mu2, mu4, alph, eig_SS, C1, S, u0, u2, u3, v4, eta0, eta2, eta4, Nx, Ny, subCC, diagCC, nthread)
xarg = mu0 * C1' * (u0 + eta0) + mu2 * alph * S' * (u2 - (1-alph) * S * u3 + eta2) + ...
        mu4 * (v4 + eta4);
xarg = single(reshape(xarg, Nx, Ny));
diagvals = diagCC + mu2 * alph^2 .* eig_SS + mu4;
x = tridiag_inv_mex_noni(subCC, diagvals, subCC, xarg, nthread);
end