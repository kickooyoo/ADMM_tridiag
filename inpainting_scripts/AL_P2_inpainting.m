function [x, xsave, err, costOrig, time] = AL_P2_inpainting(y, D, R, ...
        xinit, niters, lambda, xtrue, varargin)
% function [x, xsave, err, costOrig, time] = AL_P2_gen(y, D, R, ...
%         xinit, niters, lambda, xtrue, varargin)
% 
% implements AL-P2 with option for CG on x-update for noncirc R
% inputs:
%       y [Ns 1] undersampled data vector
%       D (fatrix2) undersampled Fourier encoding matrix
%       S (fatrix2) sensitivity map diagonal matrix
%       R (fatrix2) stack of finite difference matrices
%       xinit [Nx Ny] initial guess for x
%       niters (integer) number of iterations
%       lambda (real scalar) spatial regularization parameter
%       xtrue [Nx Ny] true solution, necessary for NRMSE comparison
% varargin:
%       mu0 [3 1] AL tuning parameters, default all 1
%       mask (logical) [Nx Ny] mask for NRMSE calculation
%       zmethod (string) 'CG' for non circulant R or 'FFT' for circulant R
%	pot (potential_fun), default l1
%       inner iter (int) inner CG steps for z-update
% 	debug (boolean)
%		plots all aux vars every iteration
%	precon (boolean) uses preconditioner for inner CG
%	save_x (boolean) if false, saves memory
%	parFFT (boolean) uses parfor for u0 update
%		
% outputs:
%       x [Nx Ny] reconstructed image
%       xsave [Nx Ny niters] estimated image at each iteration
%       err [niters 1] NRMSE at each iteratiop

%       costOrig [niters 1] objective value of original cost func at each iter
%       time [niters 1] wall time per each iteration
%
% 01/26/2016 Mai Le
% University of Michigan
Nx = size(xinit,1);
Ny = size(xinit,2);

arg.zmethod = 'CG';
arg.mu = [];
arg.mask = true(Nx, Ny);
arg.inner_iter = 1;
arg.debug = false;
arg.precon = true;
arg.save_x = false;
arg.pot = [];
arg.kapp = 12;
arg.betaw = 0; % only used for tuning
arg.alphw = 0.5; % only used for tuning
arg = vararg_pair(arg, varargin);

tic;
RR = R'*R;
Q = Gdft('mask', true(Nx, Ny));
e0 = zeros(Nx, Ny);
e0(1, 1) = 1; %end/2,end/2) = 1;
eigvalsrr = reshape(Q * RR * e0(:), Nx, Ny); % approximate, for use in precon
% eigvals for DD, get mus
DD = D'*D;
eigvalsdd = reshape(DD * ones(D.idim), Nx, Ny); 

if isempty(arg.mu)
	mu1 = 1/(arg.kapp - 1);
	w = (arg.betaw / lambda).^2 * (arg.alphw .^2 + (1 - arg.alphw).^2);
	mu0 = (arg.kapp - 1) * mu1 / (8 - (arg.kapp - 1) * w);       
else
	mu0 = arg.mu{1};
	mu1 = arg.mu{2};
end
Hu1 = eigvalsdd + mu1;

y = single(y(:));
x = single(xinit(:));
u0 = R*x;
u1 = x;
eta_0 = zeros(size(u0));
eta_1 = zeros(size(u1));

% shrinkage (e.g. soft thresholding) function, t is input, a is threshold
if isempty(arg.pot)
	arg.pot = potential_fun('l1', lambda/mu1);
end
shrink = @(a, t) arg.pot.shrink(a, t);

if strcmpi(arg.zmethod, 'FFT')
	A_CG = [];
	W_CG = [];
	P_CG = [];
else
% 	Qn = (1/sqrt(Nx * Ny)) * Gdft('mask', true(Nx, Ny));
	A_CG = [R; Gdiag(ones(Nx*Ny,1),'mask', true(Nx, Ny))];
	if isfield(R.arg, 'dim')
		Rodim = R.dim(1);
	elseif isfield(R, 'odim')
		Rodim = prod(R.odim);
	else
		display('can not figure out odim of R');
		keyboard;
	end
        W_CG = Gdiag([mu0 * ones(1, Rodim) mu1 * ones(1, Nx*Ny)]); %???
	if arg.precon
		P_CG = qpwls_precon('circ0', {A_CG, W_CG}, Gdiag(zeros(Nx*Ny,1)), true(Nx, Ny)); 
	else
		P_CG = 1;
	end
end

calc_errcost = 1;
if (calc_errcost)
        calc_orig_cost = @(y, D, R, x, lambda) norm(col(y - D*x),2)^2/2 + ...
                lambda* sum(col(arg.pot.potk(col(R*x))));
        err = calc_NRMSE_over_mask(x, xtrue, true(size(arg.mask)));
        costOrig = calc_orig_cost(y, D, R, x, lambda);
end

if arg.save_x
	xsave = zeros(Nx, Ny, niters);
else
	xsave = 0;
end
time(1) = toc;
for ii = 1:niters
        iter_start = tic;
        
        u1 = u1_update(D, y, x, eta_1, mu1, Nx*Ny, Hu1);
        x = x_update(Q, A_CG, W_CG, P_CG, u0, u1, x, mu0, mu1, R, ...
                eta_0, eta_1, eigvalsrr, Nx, Ny, arg);
        Rx = R*x;
        u0 = shrink(Rx + eta_0, lambda/mu0); 
        
        if any(isnan(x)) || any(isinf(x))
                keyboard
        end

        eta_0 = eta_0 - (u0 - Rx);
        eta_1 = eta_1 - (u1 - x);
        time = [time toc(iter_start)];
        if arg.debug
                subplot(2,2,1); im(reshape(x, Nx, Ny));
                title(sprintf('cost: %2.2f', calc_orig_cost(y, D, R, x, lambda)));
                subplot(2,2,2); im(reshape(u0(1:numel(x)), Nx, Ny));
                subplot(2,2,3); im(reshape(u0(numel(x)+1:end), Nx, Ny));
                subplot(2,2,4); im(reshape(u1, Nx, Ny));
                drawnow;
                pause(0.2);
        end
        
        if (calc_errcost)
                err(ii + 1) = calc_NRMSE_over_mask(x, xtrue, true(size(arg.mask)));
                costOrig(ii + 1) = calc_orig_cost(y, D, R, x, lambda);
        end
        if mod(ii,100) == 0
                printf('%d/%d iterations', ii, niters)
        end
   	if arg.save_x
	   	xsave(:,:,ii) = reshape(x, Nx, Ny);
	end
end
x = reshape(x, Nx, Ny);
if (~calc_errcost)
        err = zeros(niters + 1, 1);
        costOrig = zeros(niters + 1, 1);
end
end

function u1 = u1_update(D, y, x, eta_1, mu1, N, Hu1)
rhs = D' * y + mu1 * (x + eta_1);
u1 = rhs ./ Hu1(:);
end

% x-update
% problematic for non-circulant R
% x = argmin (mu0 R'R + u1 I)^(-1)(mu0 R'(u0 - eta_0) + mu1 (x + eta_1))
% x = argmin mu0/2 ||u0 - Rx + eta_0||^2 + mu1/2 ||u1 - x - eta_1||^2
function x = x_update(Q, D, W, P, u0, u1, x, mu0, mu1, R, eta_0, eta_1, ...
        eigvalsrr, Nx, Ny, arg)
% (Q, D, W, P, v, x, z, mu1, u_z, R, eta_1, eta_z, eigvalsrr, Nx, Ny, arg)
switch upper(arg.zmethod)
        case 'CG'
                y = [u0 - eta_0; u1 - eta_1];
                try
                        [xs, info] = qpwls_pcg1(x, D, W, y, Gdiag(zeros(Nx*Ny, 1)), ...
                                'niter', arg.inner_iter, 'stop_grad_tol', 1e-13, 'precon', P);
			x = xs;
                catch
                        display('qpwls failed');
                        keyboard
                end
        case 'FFT'
                rhs = reshape(mu0 * R' * (u0 - eta_0) + mu1 * (u1 - eta_1), Nx, Ny);
                invMat = mu0 * eigvalsrr + mu1;
                x = Q' * col((Q * rhs(:)) ./ col(invMat)) / (Nx*Ny);
        otherwise
                display(sprintf('unknown option for z-update: %s', arg.zmethod));
end
% keyboard
end
