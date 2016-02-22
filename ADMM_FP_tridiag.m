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
%	faNcy_mu34 (double), if specified, uses spatially varying mu3, mu4
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

x = single(xinit(:));
y = single(y);
u0 = CH*x;
u3 = x;
u1 = CV*u3;
u2 = (1-arg.alph)*S*u3+arg.alph*S*x;
v0 = zeros(size(u0));
v1 = v0;
v2 = zeros(size(u1));
v3 = v2;
v4 = zeros(size(u2));
v5 = v4;
v6 = v4;
v7 = zeros(size(u3));
v8 = v7;
eta0 = zeros(size(v0));
eta1 = zeros(size(v1));
eta2 = zeros(size(v2));
eta3 = zeros(size(v3));
eta4 = zeros(size(v4));
eta5 = zeros(size(v5));
eta6 = zeros(size(v6));
eta7 = zeros(size(v7));
eta8 = zeros(size(v8));
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
calc_cost = @(beta, CH, CV, F, S, y, x) norm(col(y) - col(F * (S * x)),2)^2/2 + ...
        beta * sum(col(arg.potx.potk(col(CH * x)))) + beta * sum(col(arg.poty.potk(col(CV * x))));

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
        Ncore = 4;
        pool = parpool(Ncore); % only in 2013b
end

err(1) = calc_NRMSE_over_mask(x, xtrue, arg.mask);
cost(1) = calc_cost(beta, CH, CV, F, S, y, x);
time(1) = toc;

for iter = 1:niters
        iter_start = tic;
        if (arg.attempt_par)
                if strcmp(arg.pmethod,'pfor')
                        parfor ui = 0:4
                                switch ui
                                        case 0
                                                %u0 = soft(-v0-eta0,beta/mu(1)); %mu0
                                                u0 = u0_update(v0, eta0, beta, mu(1));
                                        case 1
                                                %u1 = soft(-v2-eta2,beta/mu(3));
                                                u1 = u1_update(v2, eta2, beta, mu(3));
                                        case 2
                                                u2 = u2_update(mu4, eig_FF, Qbig, F, y, v4, eta4);
                                        case 3
                                                u3 = u3_update(mu,arg.alph,eig_SS,CH,CV,S,v3,v5,v7,eta3,eta5,eta7,Nx,Ny);
                                        case 4
                                                x = x_update(mu,arg.alph,eig_SS,CH,S,v1,v6,v8,eta1,eta6,eta8,Nx,Ny);
                                        otherwise
                                                display('no such case for ui loop');
                                end
                        end
                elseif strcmp(arg.pmethod,'spmd');
                        spmd;
                                if labindex == 1
                                        u0 = u0_update(v0,eta0,beta,mu(1));
                                elseif labindex == 2
                                        u1 = u1_update(v2,eta2,beta,mu(3));
                                elseif labindex == 3
                                        u2 = u2_update(mu4, eig_FF, Qbig, F, y, v4, eta4);
                                elseif labindex == 4
                                        u3 = u3_update(mu,arg.alph,eig_SS,CH,CV,S,v3,v5,v7,eta3,eta5,eta7,Nx,Ny);
                                elseif labindex == 5
                                        x = x_update(mu,arg.alph,eig_SS,CH,S,v1,v6,v8,eta1,eta6,eta8,Nx,Ny);
                                end
                        end
                        % extract vectors from composites
                        u0 = u0{1};
                        u1 = u1{2};
                        u2 = u2{3};
                        u3 = u3{4};
                        % x
                elseif strcmp(arg.pmethod,'feval')
                        u0 = parfeval(pool,@u0_update, 1, v0, eta0, beta, mu(1));
                        u1 = parfeval(pool,@u1_update, 1, v2, eta2, beta, mu(3));
                        u2 = parfeval(pool,@u2_update,1,mu4, eig_FF, Qbig, F, y, v4, eta4);
                        u3 = parfeval(pool,@u3_update,1,mu,arg.alph,eig_SS,CH,CV,S,v3,v5,v7,eta3,eta5,eta7,Nx,Ny);
                        x = parfeval(pool,@x_update,1,mu,arg.alph,eig_SS,CH,S,v1,v6,v8,eta1,eta6,eta8,Nx,Ny);
                        u0 = fetchOutputs(u0);
                        u1 = fetchOutputs(u1);
                        u2 = fetchOutputs(u2);
                        u3 = fetchOutputs(u3);
                        x = fetchOutputs(x);
                end
        else %not parallel
                u0 = u0_update(v0, eta0, beta, mu0);
                u1 = u1_update(v2, eta2, beta, mu2);
                u2 = u2_update(mu4, eig_FF, Qbig, F, y, v4, eta4);
                u3 = u3_update(mu3, mu5, mu7, arg.alph, eig_SS, CV, S, v3, ...
                        v5, v7, eta3, eta5, eta7, subCCT, diagCCT, arg.nthread);
                x = x_update(mu1, mu6, mu8, arg.alph, eig_SS, CH, S, v1, ...
                        v6, v8, eta1, eta6, eta8, subCC, diagCC, arg.nthread);
        end
        if (arg.attempt_par)
                if strcmp(arg.pmethod,'pfor')
                        parfor vi = 0:4
                                switch vi
                                        case 0
                                                v0 = v0_update(mu0, mu1, u0, eta0, CH, x, eta1);
                                        case 1
                                                v2 = v2_update(mu2, mu3, u1, eta2, CV, u3, eta3);
                                        case 2
                                                v4 = v4_update(u2,u3,x,arg.alph,S,eta4,eta5,eta6,mu4, mu5, mu6);
                                        case 3
                                                v5 = v5_update(u2,u3,x,arg.alph,S,eta4,eta5,eta6,mu4, mu5, mu6);
                                        case 4
                                                v7 = v7_update(mu7, mu8, u3, eta7, x, eta8);
                                        otherwise
                                end
                        end
                elseif strcmp(arg.pmethod,'spmd')
                        spmd;
                                if labindex == 1
                                        v0 = v0_update(mu0, mu1, u0, eta0, CH, x, eta1);
                                elseif labindex == 2
                                        v2 = v2_update(mu2, mu3, u1, eta2, CV, u3, eta3);
                                elseif labindex == 3
                                        v4 = v4_update(u2,u3,x,arg.alph,S,eta4,eta5,eta6,mu4, mu5, mu6);
                                elseif labindex == 4
                                        v5 = v5_update(u2,u3,x,arg.alph,S,eta4,eta5,eta6,mu4, mu5, mu6);
                                elseif labindex == 5
                                        v7 = v7_update(mu7, mu8, u3, eta7, x, eta8);
                                end
                        end
                        
                        % extract vectors from composite
                        v0 = v0{1};
                        v2 = v2{2};
                        v4 = v4{3};
                        v5 = v5{4};
                        % v7
                elseif strcmp(arg.pmethod,'feval')
                        v0 = parfeval(pool,@v0_update,1,mu(1),mu(2),u0,eta0,CH,x,eta1);
                        v2 = parfeval(pool,@v2_update,1,mu(3),mu(4),u1,eta2,CV,u3,eta3);
                        v4 = parfeval(pool,@v4_update,1,u2,u3,x,arg.alph,S,eta4,eta5,eta6,mu4, mu5, mu6);
                        v5 = parfeval(pool,@v5_update,1,u2,u3,x,arg.alph,S,eta4,eta5,eta6,mu4, mu5, mu6);
                        v7 = parfeval(pool,@v7_update,1,mu(8),mu(9),u3,eta7,x,eta8);
                        v0 = fetchOutputs(v0);
                        v2 = fetchOutputs(v2);
                        v4 = fetchOutputs(v4);
                        v5 = fetchOutputs(v5);
                        v7 = fetchOutputs(v7);
                        
                end
        else % not parallel
                v0 = v0_update(mu0, mu1, u0, eta0, CH, x, eta1);
                v2 = v2_update(mu2, mu3, u1, eta2, CV, u3, eta3);
                v4 = v4_update(u2, u3, x, arg.alph, S, eta4, eta5, eta6, mu4, mu5, mu6);
                v5 = v5_update(u2, u3, x, arg.alph, S, eta4, eta5, eta6, mu4, mu5, mu6);
                v7 = v7_update(mu7, mu8, u3, eta7, x, eta8);
        end
        v1 = -v0;
        v3 = -v2;
        v6 = -v4 - v5;
        v8 = -v7;
        
        % 	eta updates
        %	eta  = eta - (Au-v)
        eta0 = eta0 - (-u0 - v0);
        eta1 = eta1 - (CH * x - v1);
        eta2 = eta2 - (-u1 - v2);
        eta3 = eta3 - (CV * u3 - v3);
        eta4 = eta4 - (-u2 - v4);
        eta5 = eta5 - ((1-arg.alph) * S * u3 - v5);
        eta6 = eta6 - (arg.alph * S * x - v6);
        eta7 = eta7 - (-u3 - v7);
        eta8 = eta8 - (x - v8);
        
        time(iter + 1) = toc(iter_start);
	err(iter + 1) = calc_NRMSE_over_mask(x, xtrue, arg.mask);

        if mod(iter,10) == 0
                printf('%d/%d iterations',iter,niters)
        end

	if (mod(iter,1000) == 0) && ~isempty(arg.save_progress)
		save(sprintf('tmp_%s', arg.save_progress), 'x');
	end

        xsaved(:,:,iter) = reshape(x, Nx, Ny);
        cost(iter + 1) = calc_cost(beta, CH, CV, F, S, y, x);
end
x = reshape(x, Nx, Ny);
keyboard

if(arg.attempt_par)
        matlabpool close
end

end

function cost = calc_cost(beta, CH, CV, F, S, y, x)
cost = norm(col(y-F*(S*x)),2)^2/2 + beta*norm(col(CH*x),1) + beta*norm(col(CV*x),1);
end

function out = soft(in,thresh)
out = (in - thresh*sign(in)).*(abs(in) > thresh);
end

function u0 = u0_update(v0, eta0, beta, mu)
u0 = soft(-v0 - eta0, beta/mu);
end

function u1 = u1_update(v2, eta2, beta, mu)
u1 = soft(-v2 - eta2, beta/mu);
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

function v4 = v4_update(u2,u3,x,alph,S,eta4,eta5,eta6,mu4, mu5, mu6)
v4 = (mu5 + mu6) * (mu4 * (-u2 - eta4) + mu6 * (-alph * S * x + eta6)) + ...
        - mu6 * (mu5 * ((1-alph) * S * u3 - eta5) + mu6 * (-alph * S * x + eta6));
det = (mu4 + mu6) * (mu5 + mu6) - mu6.^2;
v4 = v4/det;
end

function v5 = v5_update(u2,u3,x,alph,S,eta4,eta5,eta6,mu4, mu5, mu6,v4);
v5 = - mu6 * (mu4 * (-u2 - eta4) + mu6 * (-alph * S * x + eta6)) + ...
         (mu4 + mu6) * (mu5 * ((1-alph) * S * u3 - eta5) + mu6 * (-alph * S * x + eta6));
det = (mu4 + mu6) * (mu5 + mu6) - mu6.^2;
v5 = v5/det;
end

function v7 = v7_update(mu7, mu8, u3, eta7, x, eta8)
v7 = (mu7 .* (-u3 - eta7) + mu8 .* (-x + eta8)) ./ (mu7 + mu8);
end
