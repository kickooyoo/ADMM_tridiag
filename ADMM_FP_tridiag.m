% function [x, xsaved, cost] = SENSE_FP(y,F,S,CH,CV,alph,beta,mu,Nx,Ny,xinit,niters,arg.attempt_par)
function [x, xsaved, err, cost, time] = ADMM_FP_tridiag(y, F, S, CH, CV, ...
        beta, xinit, xtrue, niters, varargin)
% function [x, xsaved, err, cost, time] = ADMM_tridiag(y, F, S, CH, CV, ...
%        beta, xinit, xtrue, niters, varargin)
%
% implements tridiag ADMM algo
% inputs11:
%       y [Ns*Nc 1] undersampled data vector
%       F (fatrix2) undersampled Fourier eNcoding matrix
%       S (fatrix2) sensitivity map diagonal matrix
%       CH (fatrix2) horizontal finite differeNces
%       CV (fatrix2) vertical finite differeNces
%       beta (real scalar) spatial regularization parameter
%       xinit [Nx Ny] initial guess for x
%       xtrue [Nx Ny] true solution, necessary for NRMSE comparison
%       niters (integer) number of iterations
% varargin:
%       nthread (int32) number of threads
%	mu [5 1] AL tuning parameters, default all 1
%       mask (logical) [Nx Ny] mask for NRMSE calculation
%       alph (scalar, \in [0, 1]) parameter balaNcing between u3 and x,
%               default: 0.5
%       alphw (scalar, \in [0, 1]) [[orthonormal wavelet case]]
%               parameter balaNcing between u0 and u1, default: 0.5
%       betaw (real scalar) [[orthonormal wavelet case]]
%               spatial regularization parameter for wavelets
%	pot (potential_fun) for penalty, default l1
%       compile_mex recompiles tridiag_inv_mex_noni, default: false
%	debug
%		plots all aux vars in each iteration
%	timing (string) 'all' or 'tridiag'
%	fancy_mu34 (double), if specified, uses spatially varying mu3, mu4
%		sets mu0 param relative to sense maps
%	parFFT (boolean), parfor for u2 update
%	save_progress (string)
%		if not empty, will save tmp file with string suffix
%		every 1k iters
% outputs:
%       x [Nx Ny] reconstructed image
%       xsaved [Nx Ny niters] estimated image at each iteration
%       err [niters 1] NRMSE at each iteration
%       cost [niters 1] objective value of original cost fuNc at each iter
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
arg.mu_args = {};
arg.parFFT = false;
arg.save_progress = [];
arg.potx = [];
arg.poty = [];
arg.attempt_par = 0;
arg.pmethod = 'feval'; %'pfor'; % or 'spmd' or 'feval'
arg.prof = 0;
arg.timing = false;
arg.warmup = 0; %5; % number of warmup iters
arg = vararg_pair(arg, varargin);

% eigvals for SS, get mus
SS = S'*S;
eig_SS = reshape(SS * ones(prod(S.idim),1), Nx, Ny);

if isempty(arg.mu)
        arg.mu = get_mu(eig_SS, [], Nx*Ny, beta, 'mask', arg.mask, ...
                'split', 'ADMM-FP-tridiag', 'alph', arg.alph, ...
                arg.mu_args{:});
%         arg.mu = num2cell(ones(1,9));
end

tic

if length(arg.mu) ~= 9
        display('wrong size for mu convergence parameters');
        keyboard;
end

if ~strcmp(class(arg.nthread), 'int32')
        display('nthread must be int32');
        arg.nthread = int32(arg.nthread);
end

y = single(y);
[x, u0, u1, u2, u3, v0, v2, v4, v5, v7, eta0, eta1, eta2, eta3, ...
        eta4, eta5, eta6, eta7, eta8] = init_vars(xinit, CH, CV, S, arg.alph);
mu0 = arg.mu{1};
mu1 = arg.mu{2};
mu2 = arg.mu{3};
mu3 = arg.mu{4};
mu4 = arg.mu{5};
mu5 = arg.mu{6};
mu6 = arg.mu{7};
mu7 = col(arg.mu{8});
mu8 = col(arg.mu{9});

if isempty(arg.potx) || isempty(arg.poty) || (mu0 ~= mu1)
        if ~isscalar(mu0)
                arg.potx.shrink = @(x, t) sign(col(x)).*max(abs(col(x)) - col(t) ,0);
                arg.potx.potk = @(x) abs(x);
                arg.poty.shrink = @(x, t) sign(col(x)).*max(abs(col(x)) - col(t) ,0);
                arg.poty.potk = @(x) abs(x);
        else
                arg.potx = potential_fun('l1', beta/mu0);
                arg.poty = potential_fun('l1', beta/mu1);
        end
end
shrinkx = @(a, t) arg.potx.shrink(a, t);
shrinky = @(a, t) arg.poty.shrink(a, t);

[eig_FF, Qbig] = get_eigs(F, Nc);

% pass tridiag of C'C into mex
[subCC, subCCT, diagCC, diagCCT] = construct_Hessian_diags(mu3, mu1, mu7, mu8, Nx, Ny, beta);

if arg.compile_mex
        confirm_compile('tridiag_inv_mex_noni');
end

if(arg.attempt_par)
        %if matlabpool('size') == 0
        %matlabpool('open',Ncore);
        %		pool = parpool(Ncore); % only in 2013b
        %		en
        pool = gcp('nocreate');
        if numel(pool) == 0
                pool = parpool();
        end
end

err(1) = calc_NRMSE_over_mask(x, xtrue, arg.mask);
cost(1) = calc_cost(beta, CH, CV, F, S, y, x, arg);
time(1) = toc;
if arg.prof
        profile on
end
warmup_iter = 0;
iter = 1;
while iter <= niters
        iter_start = tic;
        if arg.timing, u_start = tic; end
        if (arg.attempt_par)
                u_start = tic;
                if strcmp(arg.pmethod,'pfor')
                        parfor ui = 1:5
                                switch ui
                                        case 1
                                                %u0 = soft(-v0-eta0,beta/mu(1)); %mu0
                                                bigu{ui} = u0_update(v0, eta0, beta, mu0, arg.potx);
                                        case 2
                                                %u1 = soft(-v2-eta2,beta/mu(3));
                                                bigu{ui} = u1_update(v2, eta2, beta, mu2, arg.poty);
                                        case 3
                                                bigu{ui} = u2_update(mu4, eig_FF, Qbig, F, y, v4, eta4);
                                        case 4
                                                bigu{ui} = u3_update(mu3, mu5, mu7, arg.alph, eig_SS, CV, S, -v2, ...
                                                        v5, v7, eta3, eta5, eta7, subCCT, diagCCT, arg.nthread);
                                        case 5
                                                bigu{ui} = x_update(mu1, mu6, mu8, arg.alph, eig_SS, CH, S, -v0, ...
                                                        -v4 - v5, -v7, eta1, eta6, eta8, subCC, diagCC, arg.nthread);
                                        otherwise
                                                display('no such case for ui parfor loop');
                                end
                        end
                        u0 = bigu{1};
                        u1 = bigu{2};
                        u2 = bigu{3};
                        u3 = bigu{4};
                        x = bigu{5};
                elseif strcmp(arg.pmethod,'spmd');
                        if pool.NumWorkers < 4
                                display('not enough workers in pool, use another pmethod');
                                keyboard;
                        end
                        spmd;
                                bigu = update_u(labindex, v0, v2, v4, v5, ...
                                        v7, mu0, mu1, mu2, mu3, mu4, mu5, ...
                                        mu6, mu7, mu8, eta0, eta1, eta2, ...
                                        eta3, eta4, eta5, eta6, eta7, eta8, ...
                                        beta, subCC, subCCT, ...
                                        diagCC, diagCCT, ...
                                        CH, CV, S, eig_FF, eig_SS, F, Qbig, y, arg);

                        end
                        % extract vectors from composites
                        tmp = bigu{1};
                        u0 = tmp.H;
                        u1 = tmp.V;
                        u2 = bigu{2};
                        u3 = bigu{3};
                        x = bigu{4};
                elseif strcmp(arg.pmethod,'feval')
                        u0 = parfeval(pool,@u0_update, 1, v0, eta0, beta, mu0, arg.potx);
                        u1 = parfeval(pool,@u1_update, 1, v2, eta2, beta, mu2, arg.poty);
                        u2 = parfeval(pool,@u2_update,1, mu4, eig_FF, Qbig, F, y, v4, eta4);
                        u3 = parfeval(pool,@u3_update,1, mu3, mu5, mu7, arg.alph, eig_SS, CV, S, -v2, ...
                                v5, v7, eta3, eta5, eta7, subCCT, diagCCT, arg.nthread);
                        x = parfeval(pool,@x_update,1,mu1, mu6, mu8, arg.alph, eig_SS, CH, S, -v0, ...
                                -v4 - v5, -v7, eta1, eta6, eta8, subCC, diagCC, arg.nthread);
                        u0 = fetchOutputs(u0);
                        u1 = fetchOutputs(u1);
                        u2 = fetchOutputs(u2);
                        u3 = fetchOutputs(u3);
                        x = fetchOutputs(x);
                end
                u_times(iter) = toc(u_start);
        else 
                if arg.timing, u_start = tic; end
                u0 = u0_update(v0, eta0, beta, mu0, arg.potx);
                if arg.timing
                        u_times(iter, 1) = toc(u_start);
                        u_start = tic;
                end
                u1 = u1_update(v2, eta2, beta, mu2, arg.poty);
                if arg.timing
                        u_times(iter, 2) = toc(u_start);
                        u_start = tic;
                end
                u2 = u2_update(mu4, eig_FF, Qbig, F, y, v4, eta4);
                if arg.timing
                        u_times(iter, 3) = toc(u_start);
                        u_start = tic;
                end
                u3 = u3_update(mu3, mu5, mu7, arg.alph, eig_SS, CV, S, -v2, ...
                        v5, v7, eta3, eta5, eta7, subCCT, diagCCT, arg.nthread);
                if arg.timing
                        u_times(iter, 4) = toc(u_start);
                        u_start = tic;
                end
                x = x_update(mu1, mu6, mu8, arg.alph, eig_SS, CH, S, -v0, ...
                        -v4 - v5, -v7, eta1, eta6, eta8, subCC, diagCC, arg.nthread);
                if arg.timing
                        u_times(iter, 5) = toc(u_start);
                end
        end
        if arg. timing, v_start = tic; end
        Sx = S * x;
        Su3 = S * u3;
        AWy1 = mu4 * (-u2 - eta4) + mu6 *(-arg.alph * Sx + eta6);
        AWy2 = mu5 * ((1-arg.alph) * Su3 - eta5) + mu6 * (-arg.alph * Sx + eta6);
        v45det = (mu4 + mu6) * (mu5 + mu6) - mu6.^2;
        if (arg.attempt_par)
                if strcmp(arg.pmethod, 'pfor')
                        parfor vi = 1:5
                                switch vi
                                        case 1
                                                bigv{vi} = v0_update(mu0, mu1, u0, eta0, CH, x, eta1);
                                        case 2
                                                bigv{vi} = v2_update(mu2, mu3, u1, eta2, CV, u3, eta3);
                                        case 3
                                                bigv{vi} = v4_update(AWy1, AWy2, v45det, mu5, mu6);
                                        case 4
                                                bigv{vi} = v5_update(AWy1, AWy2, v45det, mu4, mu6);
                                        case 5
                                                bigv{vi} = v7_update(mu7, mu8, u3, eta7, x, eta8);
                                        otherwise
                                                display('no such case for vi parfor loop');
                                end
                        end
                        v0 = bigv{1};
                        v2 = bigv{2};
                        v4 = bigv{3};
                        v5 = bigv{4};
                        v7 = bigv{5};
                elseif strcmp(arg.pmethod,'spmd')
                        if pool.NumWorkers < 4
                                display('not enough workers in pool, use another pmethod');
                                keyboard;
                        end
                        spmd;
                                bigv = update_v(labindex, u0, u1, u2, u3, ...
                                        x, mu0, mu1, mu2, mu3, mu4, mu5, ...
                                        mu6, mu7, mu8, eta1, eta3, eta7, ...
                                        eta8, CH, CV, AWy1, AWy2, v45det);
                        end
                        tmp = bigv{1};
                        v0 = tmp.H;
                        v2 = tmp.V;
                        v4 = bigv{2};
                        v5 = bigv{3};
                        v7 = bigv{4};
                elseif strcmp(arg.pmethod,'feval')
                        v0 = parfeval(pool,@v0_update, 1, mu0, mu1, u0, eta0, CH, x, eta1);
                        v2 = parfeval(pool,@v2_update,1, mu2, mu3, u1, eta2, CV, u3, eta3);
                        v4 = parfeval(pool,@v4_update,1, AWy1, AWy2, v45det, mu5, mu6);
                        v5 = parfeval(pool,@v5_update,1, AWy1, AWy2, v45det, mu4, mu6);
                        v7 = parfeval(pool,@v7_update,1, mu7, mu8, u3, eta7, x, eta8);
                        v0 = fetchOutputs(v0);
                        v2 = fetchOutputs(v2);
                        v4 = fetchOutputs(v4);
                        v5 = fetchOutputs(v5);
                        v7 = fetchOutputs(v7);
                        
                end
        else % not parallel
                v0 = v0_update(mu0, mu1, u0, eta0, CH, x, eta1);
                v2 = v2_update(mu2, mu3, u1, eta2, CV, u3, eta3);
                v4 = v4_update(AWy1, AWy2, v45det, mu5, mu6);
                v5 = v5_update(AWy1, AWy2, v45det, mu4, mu6);
                v7 = v7_update(mu7, mu8, u3, eta7, x, eta8);
        end
%         v1 = -v0;
%         v3 = -v2;
%         v6 = -v4 - v5;
%         v8 = -v7;
        if arg.timing
                v_times(iter) = toc(v_start);
                eta_start = tic;
        end
        % 	eta updates
        %	eta  = eta - (Au-v)
        
        eta0 = eta0 - (-u0 - v0);
        eta1 = eta1 - (CH * x - (-v0));
        eta2 = eta2 - (-u1 - v2);
        eta3 = eta3 - (CV * u3 - (-v2));
        eta4 = eta4 - (-u2 - v4);
        eta5 = eta5 - ((1-arg.alph) * Su3 - v5);
        eta6 = eta6 - (arg.alph * Sx - (-v4 -v5));
        eta7 = eta7 - (-u3 - v7);
        eta8 = eta8 - (x - (-v7));
        if arg.timing, eta_times(iter) = toc(eta_start); end
        
        time(iter + 1) = toc(iter_start);
	err(iter + 1) = calc_NRMSE_over_mask(x, xtrue, arg.mask);

        if mod(iter,10) == 0
                printf('%d/%d iterations',iter,niters)
        end

	if (mod(iter,1000) == 0) && ~isempty(arg.save_progress)
		save(sprintf('tmp_%s', arg.save_progress), 'x');
        end
        if (any(isnan(col(x))))m
                display('getting NaNs :(');
                keyboard
        end
        xsaved(:,:,iter) = reshape(x, Nx, Ny);
        cost(iter + 1) = calc_cost(beta, CH, CV, F, S, y, x, arg);
        if warmup_iter < arg.warmup
                [x, u0, u1, u2, u3, v0, v2, v4, v5, v7, eta0, eta1, eta2, ...
                        eta3, eta4, eta5, eta6, eta7, eta8] = init_vars(...
                        xinit, CH, CV, S, arg.alph);
                warmup_iter = warmup_iter + 1;
                display('warm up iter!')
        else
                iter = iter + 1;
        end
end
x = reshape(x, Nx, Ny);
if arg.prof
        profile viewer
end

end

function [x, u0, u1, u2, u3, v0, v2, v4, v5, v7, eta0, eta1, eta2, eta3, ...
        eta4, eta5, eta6, eta7, eta8] = init_vars(xinit, CH, CV, S, alph);
x = single(xinit(:));
u0 = CH*x;
u3 = x;
u1 = CV*u3;
u2 = (1-alph) * S * u3 + alph * S * x;
v0 = zeros(size(u0));
% v1 = v0;
v2 = zeros(size(u1));
% v3 = v2;
v4 = zeros(size(u2));
v5 = v4;
% v6 = v4;
v7 = zeros(size(u3));
% v8 = v7;
eta0 = zeros(size(v0));
eta1 = zeros(size(v0));
eta2 = zeros(size(v2));
eta3 = zeros(size(v2));
eta4 = zeros(size(v4));
eta5 = zeros(size(v5));
eta6 = zeros(size(v4));
eta7 = zeros(size(v7));
eta8 = zeros(size(v7));
end

function u = update_u(labindex, v0, v2, v4, v5, v7, mu0, mu1, mu2, ...
        mu3, mu4, mu5, mu6, mu7, mu8, eta0, eta1, eta2, eta3, eta4, eta5, ...
        eta6, eta7, eta8, beta, subCC, subCCT, diagCC, diagCCT, CH, CV, S, ...
        eig_FF, eig_SS, F, Qbig, y, arg)
switch labindex
        case 1
                u.H = u0_update(v0, eta0, beta, mu0, arg.potx);
                u.V = u1_update(v2, eta2, beta, mu2, arg.poty);
        case 2
                u = u2_update(mu4, eig_FF, Qbig, F, y, v4, eta4);
        case 3
                u = u3_update(mu3, mu5, mu7, arg.alph, eig_SS, CV, S, -v2, ...
                        v5, v7, eta3, eta5, eta7, subCCT, diagCCT, arg.nthread);
        case 4
                u = x_update(mu1, mu6, mu8, arg.alph, eig_SS, CH, S, -v0, ...
                        -v4 - v5, -v7, eta1, eta6, eta8, subCC, diagCC, arg.nthread);
        otherwise
                display('invalid labindex');
                keyboard;
end
end

function v = update_v(labindex, u0, u1, u2, u3, x, mu0, mu1, mu2, ...
        mu3, mu4, mu5, mu6, mu7, mu8, eta1, eta3,  ...
        eta7, eta8, CH, CV, AWy1, AWy2, v45det)
switch labindex
        case 1
                v.H = v0_update(mu0, mu1, u0, eta0, CH, x, eta1);
                v.V = v2_update(mu2, mu3, u1, eta2, CV, u3, eta3);
        case 2
                v = v4_update(AWy1, AWy2, v45det, mu5, mu6);
        case 3
                v = v5_update(AWy1, AWy2, v45det, mu4, mu6);
        case 4
                c = v7_update(mu7, mu8, u3, eta7, x, eta8);
        otherwise
                display('invalid labindex');
                keyboard;
end
end

function cost = calc_cost(beta, CH, CV, F, S, y, x, arg)
cost = norm(col(y) - col(F * (S * x)),2)^2/2 + ...
        beta * sum(col(arg.potx.potk(col(CH * x)))) + ...
        beta * sum(col(arg.poty.potk(col(CV * x))));
end

function out = soft(in,thresh)
out = (in - thresh*sign(in)).*(abs(in) > thresh);
end

function u0 = u0_update(v0, eta0, beta, mu, pot)
% u0 = soft(-v0 - eta0, beta/mu);
u0 = pot.shrink(-v0 - eta0, beta/mu);
end

function u1 = u1_update(v2, eta2, beta, mu, pot)
% u1 = soft(-v2 - eta2, beta/mu);
u1 = pot.shrink(-v2 - eta2, beta/mu);
end

function u2 = u2_update(mu4, eig_FF, Q, F, y, v4, eta4)
arg_u2 = F' * y + mu4 * (-v4 - eta4);
u2 = Q' * ((Q * arg_u2) ./ (eig_FF + mu4)) / (F.arg.Nx*F.arg.Ny);
end

function u3 = u3_update(mu3, mu5, mu7, alph, eig_SS, CV, S, v3, v5, v7, ...
        eta3, eta5, eta7, subCCT, diagCCT, nthread)
u3arg = mu3 * CV' * (v3 + eta3) + mu5 * (1-alph) * S' * (v5 + eta5) + ...
        mu7 .* (-v7 - eta7);

% take argument of CV'CV, reshape to rect, transpose,  stretch, apply CH'CH
flipu3arg = reshape(u3arg, S.arg.Nx, S.arg.Ny);
flipu3arg = flipu3arg.';
flipSS = eig_SS.';

% construct diagonal entries
diagvals = diagCCT + mu5 * (1-alph)^2 .* flipSS;
u3out = tridiag_inv_mex_noni(subCCT, diagvals, subCCT, flipu3arg, nthread);

% transpose solution back
flipu3 = reshape(u3out, S.arg.Ny, S.arg.Nx); 
u3 = col(flipu3.');
end

function x = x_update(mu1, mu6, mu8, alph, eig_SS, CH, S, v1, v6, v8, ...
        eta1, eta6, eta8, subCC, diagCC, nthread)
xarg = mu1 * CH' * (v1 + eta1) + mu6 * alph * S' * (v6 + eta6) + ...
        mu8 .* (v8 + eta8);
xarg = reshape(xarg, S.arg.Nx, S.arg.Ny);
diagvals = diagCC + mu6 * alph^2 .* eig_SS;
x = col(tridiag_inv_mex_noni(subCC, diagvals, subCC, xarg, nthread));
end

function v0 = v0_update(mu0, mu1, u0, eta0, CH, x, eta1)
v0 = (mu0 * (-u0 - eta0) + mu1 * (- (CH * x) + eta1)) / (mu0 + mu1);
end

function v2 = v2_update(mu2, mu3, u1, eta2, CV, u3, eta3)
v2 = (mu2 * (-u1 - eta2) + mu3 * (- (CV * u3) + eta3)) / (mu2 + mu3);
end

function v4 = v4_update(AWy1, AWy2, det, mu5, mu6)
v4 = (mu5 + mu6) * AWy1 - mu6 * AWy2;
v4 = v4/det;
end

function v5 = v5_update(AWy1, AWy2, det, mu4, mu6)
v5 = - mu6 * AWy1 + (mu4 + mu6) * AWy2;
v5 = v5/det;
end

function v7 = v7_update(mu7, mu8, u3, eta7, x, eta8)
v7 = (mu7 .* (-u3 - eta7) + mu8 .* (-x + eta8)) ./ (mu7 + mu8);
end
