function [x, xsaved, err, cost, time] = ADMM_tridiag(y, F, S, CH, CV, ...
        beta, xinit, xtrue, niters, varargin)
% function [x, xsaved, err, cost, time] = ADMM_tridiag(y, F, S, CH, CV, ...
%        beta, xinit, xtrue, niters, varargin)
% 
% implements tridiag ADMM algo
% inputs11:
%       y [Ns*Nc 1] undersampled data vector
%       F (fatrix2) undersampled Fourier encoding matrix
%       S (fatrix2) sensitivity map diagonal matrix
%       CH (fatrix2) horizontal finite differences
%       CV (fatrix2) vertical finite differences
%       beta (real scalar) spatial regularization parameter
%       xinit [Nx Ny] initial guess for x
%       xtrue [Nx Ny] true solution, necessary for NRMSE comparison
%       niters (integer) number of iterations
% varargin:
%       nthread (int32) number of threads
%	mu [5 1] AL tuning parameters, default all 1
%       mask (logical) [Nx Ny] mask for NRMSE calculation
%       alph (scalar, \in [0, 1]) parameter balancing between u3 and x,
%               default: 0.5
%       alphw (scalar, \in [0, 1]) [[orthonormal wavelet case]]
%               parameter balancing between u0 and u1, default: 0.5
%       betaw (real scalar) [[orthonormal wavelet case]] 
%               spatial regularization parameter for wavelets
%	pot (potential_fun) for penalty, default l1
%       compile_mex recompiles tridiag_inv_mex_noni, default: false
%	debug
%		plots all aux vars in each iteration
%	timing (string) 'all' or 'tridiag'
%	fancy_mu34 (boolean), uses spatially varying mu3, mu4
%	parFFT (boolean), parfor for u2 update
%	save_progress (string) 
%		if not empty, will save tmp file with string suffix 
%		every 1k iters
% outputs:
%       x [Nx Ny] reconstructed image
%       xsaved [Nx Ny niters] estimated image at each iteration
%       err [niters 1] NRMSE at each iteration
%       cost [niters 1] objective value of original cost func at each iter
%       time [niters 1] wall time per each iteration
%
% 01/26/2016 Mai Le
% 02/15/2016 bells and whistles 
% University of Michigan

Nx = size(xinit,1);
Ny = size(xinit,2);
Nc = S.arg.Nc;

arg.compile_mex = false;
arg.mask = true(Nx, Ny); 
arg.mu = [];
arg.alph = 0.5;
arg.alphw = 0.5;
arg.betaw = 0;
arg.debug = false;
arg.nthread = int32(jf('ncore'));
arg.timing = 'all'; % 'tridiag'
arg.fancy_mu34 = false;
arg.parFFT = false;
arg.save_progress = [];
arg.pot = [];
arg = vararg_pair(arg, varargin);

if arg.parFFT
        tmp = gcp('nocreate');
        if prod(size(tmp)) == 0
                parpool()
        end
end

% eigvals for SS, get mus
SS = S'*S;
eig_SS = reshape(SS * ones(prod(S.idim),1), Nx, Ny); 

if isempty(arg.mu)
	arg.mu = get_mu(eig_SS, [], Nx*Ny, beta, ...
                'split', 'ADMM-tridiag', 'alph', arg.alph, 'fancy_mu34', arg.fancy_mu34);
        % do we really want to mask for mu?
end

tic

if length(arg.mu) ~= 5
	display('wrong size for mu convergence parameters');
	keyboard;
end

if ~strcmp(class(arg.nthread), 'int32')
	display('nthread must be int32');
	arg.nthread = int32(arg.nthread);
end

%iter = 0;
x = single(xinit(:));
y = single(y);
u0 = CH * x;
u3 = x;
u1 = CV * u3;
u2 = (1-arg.alph) * S * u3 + arg.alph * S * x;
v3 = zeros(size(u3));
v4 = -v3;
eta0 = zeros(size(u0));
eta1 = zeros(size(u1));
eta2 = zeros(size(u2));
eta3 = zeros(size(u3));
eta4 = zeros(size(x));
% renaming to match paper indeces
mu0 = arg.mu{1};
mu1 = arg.mu{2};
mu2 = arg.mu{3};
mu3 = col(arg.mu{4});
mu4 = col(arg.mu{5});

if mu0 ~= mu1
	display('using warning: otential function for mu0 everywhere');
	keyboard
end

if isempty(arg.pot)
	arg.pot = potential_fun('l1', beta/mu0);
end
shrink = @(a, t) arg.pot.shrink(a, t);
calc_cost = @(beta, CH, CV, F, S, y, x) norm(col(y) - col(F * (S * x)),2)^2/2 + ...
        beta * sum(col(arg.pot.potk(col(CH * x)))) + beta * sum(col(arg.pot.potk(col(CV * x))));

[eig_FF,Qbig] = get_eigs(F,Nc);

% pass tridiag of C'C into mex
Wconst = arg.betaw * arg.alphw / beta;
WVconst = arg.betaw * (1-arg.alphw) / beta;
subCCT = single(-mu1 * ones(Ny - 1, Nx)); % tranpose dims 
subCC = single(-mu0 * ones(Nx - 1, Ny));
if arg.fancy_mu34
        diagCCT = single(mu1 * (cat(1, ones(1, Nx), 2*ones(Ny-2, Nx), ones(1, Nx)) + WVconst.^2) + ...
                reshape(mu3, Nx, Ny).' + mu1 * WVconst);
        diagCC = single(mu0 * (cat(1, ones(1, Ny), 2*ones(Nx-2, Ny), ones(1, Ny)) + Wconst.^2) + ...
                reshape(mu4, Nx, Ny) + mu0 * Wconst);
else
        diagCCT = single(mu1 * (cat(1, ones(1, Nx), 2*ones(Ny-2, Nx), ones(1, Nx)) + WVconst.^2) + ...
                mu3 + mu1 * WVconst);
        diagCC = single(mu0 * (cat(1, ones(1, Ny), 2*ones(Nx-2, Ny), ones(1, Ny)) + Wconst.^2) + ...
                mu4 + mu0 * Wconst);
end

% to do: make sure mex compiled
if arg.compile_mex
        curr = cd;
        if ~strcmp(curr(end-11:end), 'ADMM_tridiag')
                display('must be in ADMM_tridiag dir to compile tridiag_inv_mex_noni.c');
                keyboard;
        end
        mex -O CFLAGS="\$CFLAGS -std=c99 -DMmex" -I./pthread_tutor/def/ ./pthread_tutor/tridiag_inv_mex_noni.c
end

err(1) = calc_NRMSE_over_mask(x, xtrue, arg.mask);
cost(1) = calc_cost(beta, CH, CV, F, S, y, x);
time(1) = toc;
tridiag_time(1) = 0;
%while(iter < niters)
for iter = 1:niters
        iter_start = tic;
        u0 = shrink(CH * x - eta0, beta/mu0);
        u1 = shrink(CV * u3 - eta1, beta/mu1);
        u2 = u2_update(mu2, arg.alph, eig_FF, Qbig, F, S, y, u3, x, eta2, Nx, Ny, arg.parFFT);
	tridiag_tic = tic;
        u3 = u3_update_mex(mu1, mu2, mu3, arg.alph, eig_SS, CV, S, u1, u2, x, v3, eta1, eta2, eta3, Nx, Ny, subCCT, diagCCT, arg.nthread);
        x = x_update(mu0, mu2, mu4, arg.alph, eig_SS, CH, S, u0, u2, u3, v4, eta0, eta2, eta4, Nx, Ny, subCC, diagCC, arg.nthread);
        tridiag_time(iter + 1) = toc(tridiag_tic);

        % skip v0, v1, v2 because they are constrained to be zero
        v3 = (mu3 .* (-u3 - eta3) + mu4 .* (-x + eta4)) ./ (mu3 + mu4);
        v4 = -v3;

	if arg.debug
		subplot(2,2,1); im(reshape(abs(x), Nx, Ny));
		subplot(2,2,2); im(reshape(abs(u2), Nx, Ny, Nc));
		subplot(2,2,3); im(cat(1, reshape(abs(u0), Nx, Ny), reshape(abs(u1), Nx, Ny)));
		subplot(2,2,4); im(reshape(abs(u3), Nx, Ny));
		drawnow;
		pause(1);
	end
        
        % eta updates
        eta0 = eta0 - (-u0 + CH * x);
        eta1 = eta1 - (-u1 + CV * u3);
        eta2 = eta2 - (-u2 + (1 - arg.alph) * S * u3 + arg.alph * S * x);
        eta3 = eta3 - (-u3 - v3);
        eta4 = eta4 - (x - v4);
        
        time(iter + 1) = toc(iter_start);
	err(iter + 1) = calc_NRMSE_over_mask(x, xtrue, arg.mask);

        if mod(iter,10) == 0
                printf('%d/%d iterations',iter,niters)
        end

	if (mod(iter,1000) == 0) && ~isempty(arg.save_progress)
		save(sprintf('tmp_%s',arg.save_progress), 'x');
	end

        xsaved(:,:,iter) = reshape(x, Nx, Ny);
        cost(iter + 1) = calc_cost(beta, CH, CV, F, S, y, x);
end
x = reshape(x, Nx, Ny);
if strcmp(arg.timing, 'tridiag')
	time = tridiag_time;
end
end

function u2 = u2_update(mu2, alph, eig_FF, Q, F, S, y, u3, x, eta2, Nx, Ny, par)
arg_u2 = F' * y + mu2 * ((1-alph) * S * u3 + alph * S * x - eta2);
if par
        arg_u2 = reshape(arg_u2, F.arg.Nx, F.arg.Ny, F.arg.Nc);
        parfor ii = 1:F.arg.Nc
                U2(:,:,ii) = fft2(arg_u2(:,:,ii));
        end
        U2 = U2./(reshape(eig_FF, F.arg.Nx, F.arg.Ny, F.arg.Nc) + mu2);
        parfor ii = 1:F.arg.Nc
                u2(:,:,ii) = ifft2(U2(:,:,ii));
        end
        u2 = col(u2);
else
        u2 = Q' * ((Q * arg_u2) ./ (eig_FF + mu2)) / (Nx * Ny);
end
end

function u3 = u3_update_mex(mu1, mu2, mu3, alph, eig_SS, CV, S, u1, u2, x, ...
        v3, eta1, eta2, eta3, Nx, Ny, subCCT, diagCCT, nthread)
u3arg = mu1 * CV' * (u1 + eta1) + mu2 * (1-alph) * S' * (u2 - alph * S * x + eta2) + ...
        mu3 .* (-v3 - eta3);

% transpose to make Hessian tridiagonal, size now Ny Nx
flipu3arg = reshape(u3arg, Nx, Ny);
flipu3arg = flipu3arg.';
flipSS = eig_SS.';

diagvals = diagCCT + mu2 * (1-alph)^2 * flipSS;
u3out = tridiag_inv_mex_noni(subCCT, diagvals, subCCT, flipu3arg, nthread);

% transpose solution back
flipu3 = reshape(u3out, Ny, Nx); 
u3 = col(flipu3.');
end

function x = x_update(mu0, mu2, mu4, alph, eig_SS, CH, S, u0, u2, u3, v4, ...
        eta0, eta2, eta4, Nx, Ny, subCC, diagCC, nthread)
xarg = mu0 * CH' * (u0 + eta0) + mu2 * alph * S' * (u2 - (1-alph) * S * u3 + eta2) + ...
        mu4 .* (v4 + eta4);
xarg = reshape(xarg, Nx, Ny);
diagvals = diagCC + mu2 * alph^2 .* eig_SS;
x = col(tridiag_inv_mex_noni(subCC, diagvals, subCC, xarg, nthread));
end
