function [x, xsaved, err, cost, time] = tridiag_ADMM_W(y, F, S, C1, C2, ...
        alph, beta, xinit, xtrue, niters, varargin)
% function [x, xsaved, err, cost, time] = tridiag_ADMM_W(y, F, S, C1, C2, ...
%         alph, beta, mu, Nx, Ny, xinit, xtrue, niters, varargin)
%
% implements tridiag ADMM algo with wavelets
% inputs:
%       y [Ns*Nc 1] undersampled data vector
%       F (fatrix2) undersampled Fourier eNcoding matrix
%       S (fatrix2) sensitivity map diagonal matrix
%       C1 (fatrix2) horizontal finite differeNces
%       C2 (fatrix2) vertical finite differeNces
%       alph (scalar, \in [0, 1]) parameter balaNcing between u3 and x
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
%       cost [niters 1] objective value of original cost fuNc at each iter
%       time [niters 1] wall time per each iteration
%
% 01/28/2016 Mai Le
% University of Michigan

Nx = size(xinit,1);
Ny = size(xinit,2);

Nc = S.arg.Nc;
arg.compile_mex = false;
arg.mask = true(Nx, Ny);
arg.mu = ones(1,5);
arg.alph2 = 0.5;
arg.beta2 = beta;
arg = vararg_pair(arg, varargin);

% construct wavelets
%W = Gwave1('mask',true(Nx,Ny),'wname','haar'); % how do i make sure it is undecimated?
W = Godwt1(true(Nx,Ny));
CHW = [C1; arg.beta2*arg.alph2/beta*W];
CVW = [C2; arg.beta2*arg.alph2/beta* W];

iter = 0;
x = single(xinit(:));
y = single(y);
u0 = CHW * x;
u3 = x;
u1 = CVW * u3;
u2 = (1-alph) * S * u3 + alph * S * x;
v3 = zeros(size(u3));
v4 = -v3;
eta0 = zeros(size(u0));
eta1 = zeros(size(u1));
eta2 = zeros(size(u2));
eta3 = zeros(size(u3));
eta4 = zeros(size(x));
% rename to match paper indeces
mu0 = arg.mu(1);
mu1 = arg.mu(2);
mu2 = arg.mu(3);
mu3 = arg.mu(4);
mu4 = arg.mu(5);

[eig_FF,Qbig] = get_eigs(F,Nc);
SS = S'*S;
eig_SS = SS*ones(S.idim,1);

% pass tridiag of CHW'CHW into mex
Wconst = arg.beta2*arg.alph2/beta;
subCCT = single(-mu1 * ones(Ny - 1, Nx));
diagCCT = single(mu1 * (cat(1, ones(1, Nx), 2*ones(Ny-2, Nx), ones(1, Nx)) + Wconst.^2));
subCC = single(-mu0 * ones(Nx - 1, Ny));
diagCC = single(mu0 * (cat(1, ones(1, Ny), 2*ones(Nx-2, Ny), ones(1, Ny)) + Wconst.^2));

% to do: make sure mex compiled
nthread = int32(jf('Ncore'));
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
cost(1) = calc_cost(beta, arg.beta2, C1, C2, F, S, W, y, x);
while(iter < niters)
        iter_start = tic;
        u0 = soft(CHW * double(x) - eta0, beta/mu0);
        u1 = soft(CVW * double(u3) - eta1, beta/mu1);
        u2 = u2_update(mu2, alph, eig_FF, Qbig, F, S, y, u3, x, eta2, Nx, Ny, Nc);
        u3 = u3_update_mex(mu, alph, arg.alph2, beta, arg.beta2, eig_SS, ...
                C1, CVW, S, u1, u2, x, v3, eta1, eta2, eta3, Nx, Ny);
        x = x_update(mu, alph, arg.alph2, beta, arg.beta2, eig_SS, CHW, S, ...
                u0,u2,u3,v4,eta0,eta2,eta4,Nx,Ny);
        
        % skip v0, v1, v2 because they are constrained to be zero
        v3 = (mu3 * (-u3 - eta3) + mu4 * (-x + eta4)) ./ (mu3 + mu4);
        v4 = -v3;
        
        % eta updates
        eta0 = eta0 - (-u0 + CHW * x);
        eta1 = eta1 - (-u1 + CVW * u3);
        eta2 = eta2 - (-u2 + (1-alph) * S * u3 + alph * S * x);
        eta3 = eta3 - (-u3 - v3); 
        eta4 = eta4 - (x - v4);

        iter = iter+1;
        time(iter) = toc(iter_start);
        err(iter) = calc_NRMSE_over_mask(x, xtrue, arg.mask);
        
        if mod(iter,10) == 0
                printf('%d/%d iterations',iter,niters)
        end
        xsaved(:,:,iter) = reshape(x, Nx, Ny);
        cost(iter + 1) = calc_cost(beta, arg.beta2, C1, C2, F, S, W, y, x);
end

end

function cost = calc_cost(beta, arg.beta2, C1, C2, F, S, W, y, x)
cost = norm(col(y) - col(F * (S * x)),2)^2/2 + beta * norm(C1 * x,1) + ...
        beta * norm(C2 * x,1) + arg.beta2 * norm(W * x,1);
end

function out = soft(in,thresh)
out = (in - thresh * sign(in)) .* (abs(in) > thresh);
end

function u2 = u2_update(mu2, alph, eig_FF, Q, F, S, y, u3, x, eta2, Nx, Ny)
arg_u2 = F' * y + mu2 * ((1-alph) * S * u3 + alph * S * x - eta2);
u2 = Q' * ((Q * arg_u2) ./ (eig_FF + mu2)) / (Nx * Ny);
end

function u3 = u3_update_mex(mu1, mu2, mu3, alph, alph2, beta1, beta2, ...
        eig_SS, C1, CVW, S, u1, u2, x, v3, eta1, eta2, eta3, Nx, Ny, ...
        subCCT, diagCCT, nthread)
u3arg = mu1 * CVW' * (u1 + eta1) + mu2 * (1-alph) * S' * (u2 - alph * S * x + eta2) + ...
        mu3 * (-v3 - eta3);

% transpose to make Hessian tridiagonal, size now Ny Nx
flipu3arg = single(reshape(u3arg, Nx, Ny));
flipu3arg = flipu3arg.';
flipSS = eig_SS.';

diagvals = diagCCT + mu2 * (1-alph)^2 * flipSS + mu3 + mu1 * (beta2/beta1)^2 * (1-alph2)^2;
u3out = tridiag_inv_mex_noni(subCCT, diagvals, subCCT, flipu3arg, nthread);

% transpose solution back
flipu3 = reshape(u3out, Ny, Nx); 
u3 = col(flipu3.');
end

function x = x_update(mu0, mu2, mu4, alph, alph2, beta1, beta2, eig_SS, ...
        CHW, S, u0, u2, u3, v4, eta0, eta2, eta4, Nx, Ny, subCC, diagCC, nthread)
xarg = mu0 * CHW' * (u0 + eta0) + mu2 * alph * S' * (u2 - (1-alph) * S * u3 + eta2) + ...
        mu4 * (v4 + eta4);
xarg = reshape(xarg, Nx, Ny);
diagvals = diagCC + mu2 * alph^2 .* eig_SS + mu4 + mu0 * (beta2/beta1)^2 * alph2^2;
x = tridiag_inv_mex_noni(subCC, diagvals, subCC, xarg, nthread);
end





