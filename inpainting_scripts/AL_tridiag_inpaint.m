function [x, xsaved, err, cost, time] = AL_tridiag_inpaint(y, D, CH, CV, ...
        beta, xinit, xtrue, niters, varargin)
% function [x, xsaved, err, cost, time] = AL_tridiag(y, D, CH, CV, ...
%        beta, xinit, xtrue, niters, varargin)
% 
% implements tridiag AL algo for inpainting
% inputs:
%       y [Ns 1] undersampled data vector
%       D (fatrix2) diagonal inpainting matrix
%       CH (fatrix2) horizontal finite differences
%       CV (fatrix2) vertical finite differences
%       beta (real scalar) spatial regularization parameter
%       xinit [Nx Ny] initial guess for x
%       xtrue [Nx Ny] true solution, necessary for NRMSE comparison
%       niters (integer) number of iterations
% varargin:
%       nthread (int32) number of threads
%	mu [3 1] AL tuning parameters, default all 1
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
%	fancy_mu34 (double), if specified, uses spatially varying mu3, mu4
%		sets mu0 param relative to sense maps
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

[Nx, Ny] = size(D.arg.samp);

arg.compile_mex = false;
arg.mask = true(Nx, Ny); 
arg.mu = [];
arg.alph = 0.5;
arg.alphw = 0.5;
arg.betaw = 0;
arg.debug = false;
arg.nthread = int32(jf('ncore'));
arg.timing = 'all'; % 'tridiag'
arg.mu_args = {};
arg.save_progress = [];
arg.potx = [];
arg.poty = [];
arg.kapp = 4;
arg.pos = 0.1;
arg.output_xsaved = 0;
arg.Dfatrix = false;
arg = vararg_pair(arg, varargin);

% eigvals for DD, get mus
DD = D'*D;
eig_DD = reshape(DD * ones(D.idim), Nx, Ny); 

if isempty(arg.mu)
        CCx = sparse(diag(-1*ones(1,Nx-1),-1) + diag(cat(2, 1, 2*ones(1,Nx - 2), 1)) + diag(-1*ones(1,Nx-1),1));
	CCy = sparse(diag(-1*ones(1,Ny-1),-1) + diag(cat(2, 1, 2*ones(1,Ny - 2), 1)) + diag(-1*ones(1,Ny-1),1));
	[eigV, eigD] = eigs(CCx);
	RRmaxx = max(eigD(:));
	[eigV, eigD] = eigs(CCy);
	RRmaxy = max(eigD(:));
	c = max(arg.alph.^2, (1-arg.alph).^2) + arg.pos;
	w = (arg.betaw * arg.alphw / beta ).^2; % hardcoded for alphw = 0.5 :(
	mu1 = c*(arg.kapp - 1)./(RRmaxy - (arg.kapp - 1) * w);
	mu0 = c*(arg.kapp - 1)./(RRmaxx - (arg.kapp - 1) * w);
	mu2 = c - max(arg.alph.^2, (1-arg.alph).^2)*eig_DD;
else
	mu0 = arg.mu{1};
	mu1 = arg.mu{2};
	mu2 = arg.mu{3};
end

tic

if ~strcmp(class(arg.nthread), 'int32')
	display('nthread must be int32');
	arg.nthread = int32(arg.nthread);
end

%iter = 0;
x = single(xinit(:));
y = single(y);
u0 = CH * x;
u2 = x;
u1 = CV * u2;
eta0 = zeros(size(u0));
eta1 = zeros(size(u1));
eta2 = zeros(size(u2));
% renaming to match paper indices

if isempty(arg.potx) || isempty(arg.poty) || (mu0 ~= mu1)
	if ~isscalar(mu0) || ~isscalar(mu1) || ~isscalar(beta)
		arg.potx.shrink = @(x, t) sign(col(x)).*max(abs(col(x)) - col(t) ,0);
		arg.potx.potk = @(x) abs(x);
		arg.poty.shrink = @(x, t) sign(col(x)).*max(abs(col(x)) - col(t) ,0);
		arg.poty.potk = @(x) abs(x);
	else 	
		arg.potx = potential_fun('l1', beta/mu0);
                arg.potx = potential_fun('l1', beta/mu0);
		arg.poty = potential_fun('l1', beta/mu1);
	end
end
shrinkx = @(a, t) arg.potx.shrink(a, t);
shrinky = @(a, t) arg.poty.shrink(a, t);
calc_cost_tridiag_inpaint = @(beta, CH, CV, D, y, x) norm(double(col(y) - col(D * x)),2)^2/2 + ...
        sum(col(beta(:) .* arg.potx.potk(col(CH * x))), 'double') + sum(col(beta(:) .* arg.poty.potk(col(CV * x))),'double');

% pass tridiag of C'C into mex
% subCC = - mu0 * I, subCCT = - mu1 * I
% diagCC = mu0 * Ch'*Ch + mu2 + mu0 * betaw * alphaw /beta
% diagCCT = mu1 * Cv'*Cv + mu2 + mu1 * betaw * (1-alphw) / beta
[subCC, subCCT, diagCC, diagCCT] = construct_Hessian_diags(mu0, mu1, mu2, mu2, Nx, Ny, beta, 'betaw', arg.betaw, 'alphw', arg.alphw); 

if arg.compile_mex
        confirm_compile('tridiag_inv_mex_noni');
end

if ~arg.Dfatrix 
        Dsamp = D.arg.samp;
end

% precompute
flipDD = eig_DD.';
diagvalsCCT = diagCCT + (1-arg.alph)^2 * flipDD; % mu1 Cv'Cv + mu2 I
diagvalsCC = diagCC + arg.alph^2 .* eig_DD; % diag CC = mu0 Ch'Ch + mu2 I

err(1) = calc_NRMSE_over_mask(x, xtrue, true(size(arg.mask)));
cost(1) = calc_cost_tridiag_inpaint(beta, CH, CV, D, y, x);
time(1) = toc;
tridiag_time(1) = 0;
%while(iter < niters)
for iter = 1:niters
        iter_start = tic;
        tridiag_tic = tic;
        if arg.Dfatrix
                u2 = u2_update(mu1, mu2, arg.alph, eig_DD, CV, D, y, u1, x, ...
                        eta1, eta2, subCCT, diagCCT, arg.nthread);
                if any(isnan(u2(:))) || any(u2(:) > 1e5), keyboard; end
                x = x_update(mu0, mu2, arg.alph, eig_DD, CH, D, y, u0, u2, ...
                        eta0, eta2, subCC, diagCC, arg.nthread);
                if any(isnan(x(:))) || any(x(:) > 1e5), keyboard; end
        else
                u2 = u2_update_nf(mu1, mu2, arg.alph, CV, Dsamp, y, u1, x, ...
                        eta1, eta2, subCCT, diagvalsCCT, arg.nthread, Nx, Ny);
                if any(isnan(u2(:))) || any(u2(:) > 1e5), keyboard; end
                x = x_update_nf(mu0, mu2, arg.alph, CH, Dsamp, y, u0, u2, ...
                        eta0, eta2, subCC, diagvalsCC, arg.nthread, Nx, Ny);
                if any(isnan(x(:))) || any(x(:) > 1e5), keyboard; end
        end
        tridiag_time(iter + 1) = toc(tridiag_tic);
        % compute once
        CHx = CH*x;
        CVu2 = CV*u2;
        u0 = shrinkx(CHx - eta0, beta./mu0);
        if any(isnan(u0(:))) || any(u0(:) > 1e5), keyboard; end
        u1 = shrinky(CVu2 - eta1, beta./mu1);
        if any(isnan(u1(:))) || any(u1(:) > 1e5), keyboard; end
        
        % eta updates
        eta0 = eta0 - (-u0 + CHx);
        eta1 = eta1 - (-u1 + CVu2);
        eta2 = eta2 - (-u2 + x);
        
        time(iter + 1) = toc(iter_start);
	err(iter + 1) = calc_NRMSE_over_mask(x, xtrue, true(size(arg.mask)));
        
        if arg.debug
                subplot(2,2,1); im(reshape(x, Nx, Ny));
                subplot(2,2,2); im(reshape(u2, Nx, Ny));
                subplot(2,2,3); im(reshape(u0, Nx, Ny));
                subplot(2,2,4); im(reshape(u1, Nx, Ny));
		drawnow;
		pause(1);
        end
        
        if mod(iter,100) == 0
                printf('%d/%d iterations',iter,niters)
        end

	if (mod(iter,1000) == 0) && ~isempty(arg.save_progress)
		save(sprintf('tmp_%s',arg.save_progress), 'x');
	end
	if arg.output_xsaved
		xsaved(:,:,iter) = reshape(x, Nx, Ny);
	end
        cost(iter + 1) = calc_cost_tridiag_inpaint(beta, CH, CV, D, y, x);
end
x = reshape(x, Nx, Ny);
if strcmp(arg.timing, 'tridiag')
	time = tridiag_time;
end

if ~arg.output_xsaved
	xsaved = 0;
end
end

function u2 = u2_update(mu1, mu2, alph, eig_DD, CV, D, y, u1, x, ...
        eta1, eta2, subCCT, diagCCT, nthread)
u2arg = mu1(:) .* (CV' * (u1 + eta1)) + (1-alph) * D' * (y - alph * D * x) + ...
        mu2(:) .* (x - eta2);

% transpose to make Hessian tridiagonal, size now Ny Nx
flipu2arg = reshape(u2arg, D.arg.Nx, D.arg.Ny);
flipu2arg = flipu2arg.';
flipDD = eig_DD.';

diagvals = diagCCT + (1-alph)^2 * flipDD; % mu1 Cv'Cv + mu2 I
u2out = tridiag_inv_mex_noni(subCCT, diagvals, subCCT, flipu2arg, nthread);

% transpose solution back
flipu2 = reshape(u2out, D.arg.Ny, D.arg.Nx); 
u2 = col(flipu2.');
end

function x = x_update(mu0, mu2, alph, eig_DD, CH, D, y, u0, u2, ...
        eta0, eta2, subCC, diagCC, nthread)
xarg = mu0(:) .* (CH' * (u0 + eta0)) + alph * D' * (y - (1- alph) * D * u2) + ...
        mu2(:) .* (u2 + eta2);
xarg = reshape(xarg, D.arg.Nx, D.arg.Ny);
diagvals = diagCC + alph^2 .* eig_DD; % diag CC = mu0 Ch'Ch + mu2 I

x = col(tridiag_inv_mex_noni(subCC, diagvals, subCC, xarg, nthread));

end

function u2 = u2_update_nf(mu1, mu2, alph, CV, samp, y, u1, x, ...
        eta1, eta2, subCCT, diagvals, nthread, Nx, Ny)
u2arg = mu1(:) .* (CV' * (u1 + eta1)) + (1-alph) * col(embed(y - alph * masker(x, samp), samp)) + ...
        mu2(:) .* (x - eta2);

% transpose to make Hessian tridiagonal, size now Ny Nx
flipu2arg = reshape(u2arg, Nx, Ny);
flipu2arg = flipu2arg.';
u2out = tridiag_inv_mex_noni(subCCT, diagvals, subCCT, flipu2arg, nthread);

% transpose solution back
flipu2 = reshape(u2out, Ny, Nx); 
u2 = col(flipu2.');
end

function x = x_update_nf(mu0, mu2, alph, CH, samp, y, u0, u2, ...
        eta0, eta2, subCC, diagvals, nthread, Nx, Ny)
xarg = mu0(:) .* (CH' * (u0 + eta0)) + alph * col(embed(y - (1- alph) * masker(u2, samp), samp)) + ...
        mu2(:) .* (u2 + eta2);
xarg = reshape(xarg, Nx, Ny);

x = col(tridiag_inv_mex_noni(subCC, diagvals, subCC, xarg, nthread));

end
