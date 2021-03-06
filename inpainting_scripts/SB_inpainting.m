function [x, xsave, err, costOrig, time] = SB_inpainting(y, D, R, ...
        xinit, niters, lambda, xtrue, varargin)
% 
% implements Split Bregman with option for CG on x-update for noncirc R
% min 1/2 ||y - Dx||^2 + beta ||u||_1, u = Rx 
%
% min 1/2 ||y - Dx||^2 + beta ||u||_1 + mu/2||u - Rx||^2 
%
% 
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
%       mu [3 1] AL tuning parameters, default all 1
%       mask (logical) [Nx Ny] mask for NRMSE calculation
%       zmethod (string) 'CG' for non circulant R or 'FFT' for circulant R
%	pot (potential_fun), default l1
%       inner iter (int) inner CG steps for z-update
% 	debug (boolean)
%		plots all aux vars every iteration
%	precon (boolean) uses preconditioner for inner CG
%	save_x (boolean) if false, saves memory
%	parFFT (boolean) uses parfor for u update
%		
% outputs:
%       x [Nx Ny] reconstructed image
%       xsave [Nx Ny niters] estimated image at each iteration
%       err [niters 1] NRMSE at each iteratiop

%       costOrig [niters 1] objective value of original cost func at each iter
%       time [niters 1] wall time per each iteration
%

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
arg.save_progress = [];
arg = vararg_pair(arg, varargin);

tic;
if arg.betaw ~= 0
	justR = R.arg.blocks{1};
	justW = R.arg.blocks{2};
	RR = justR'*justR + justW'*justW;
else       
	RR = R'*R;
end
Q = Gdft('mask', true(Nx, Ny));
e0 = zeros(Nx, Ny);
e0(1, 1) = 1; %end/2,end/2) = 1;
eigvalsrr = real(reshape(Q * RR * e0(:), Nx, Ny)); % approximate, for use in precon
if arg.betaw ~= 0
%	eigvalsrr = eigvalsrr + sqrt(Nx*Ny)* (arg.betaw/lambda).^2;
end

% eigvals for DD, get mus
DD = D'*D;
eigvalsdd = reshape(DD * ones(D.idim), Nx, Ny); 

if isempty(arg.mu)
	mu = 1; % kapp?
else
	mu = arg.mu{1};
end

y = single(y(:));
x = single(xinit(:));
u = R*x;
eta = zeros(size(u));

% shrinkage (e.g. soft thresholding) function, t is input, a is threshold
if isempty(arg.pot)
	arg.pot = potential_fun('l1', lambda/mu);
end
shrink = @(a, t) arg.pot.shrink(a, t);

if strcmpi(arg.zmethod, 'FFT')
	A_CG = [];
	W_CG = [];
	P_CG = [];
else
% 	Qn = (1/sqrt(Nx * Ny)) * Gdft('mask', true(Nx, Ny));
	A_CG = [D; R];
	if isfield(R.arg, 'dim')
		Rodim = R.dim(1);
	elseif isfield(R, 'odim')
		Rodim = prod(R.odim);
	else
		display('can not figure out odim of R');
		keyboard;
	end
        W_CG = Gdiag([ones(1, D.odim) mu * ones(1, Rodim)]); 
	if arg.precon
		%P_CG = qpwls_precon('circ0', {A_CG, W_CG}, Gdiag(zeros(Nx*Ny,1)), true(Nx, Ny)); 
		T = Q * Gdiag((D.odim / numel(x)) + eigvalsrr) * Q';
		P_CG = qpwls_precon('circ0', {T}, Gdiag(zeros(Nx*Ny,1)), true(Nx, Ny)); 
	else
		P_CG = 1;
	end
end

calc_errcost = 1;
if (calc_errcost)
        calc_orig_cost = @(y, D, R, x, lambda) calc_cost_tridiag_inpaint(y, D, R, 0, x, lambda);
	%norm(col(y - D*x),2)^2/2 + ...
         %       lambda* sum(col(arg.pot.potk(col(R*x))));
        err = calc_NRMSE(x, xtrue, true(size(arg.mask)));
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
       	
	x = x_update(Q, A_CG, W_CG, P_CG, y, u, x, mu, R, D, eta, eigvalsrr, eigvalsdd, Nx, Ny, arg);
	Rx = R*x;
	u = shrink(Rx + eta, lambda/mu);

        eta = eta - (u - Rx);
        time = [time toc(iter_start)];
        if arg.debug && mod(ii, 10) == 0
                subplot(2,2,1); im(reshape(x, Nx, Ny));
                subplot(2,2,2); im(reshape(u, numel(u)/Ny, Ny));
		subplot(2,2,4); plot(costOrig)
                title(sprintf('cost: %2.2f', calc_orig_cost(y, D, R, x, lambda)));
                drawnow;
                keyboard;
        end
        if (calc_errcost)
                err(ii + 1) = calc_NRMSE(x, xtrue, true(size(arg.mask)));
                costOrig(ii + 1) = calc_orig_cost(y, D, R, x, lambda);
        end
        if mod(ii,100) == 0
                printf('%d/%d iterations', ii, niters)
        	if mod(ii,1000) == 0 &  ~isempty(arg.save_progress)
			save(arg.save_progress)
		end
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

% x-update
% problematic for non-circulant R
% x = argmin (mu R'R + D'D)^(-1)(mu R'(0 - eta) + D'y)
% x = argmin mu/2 ||u - Rx + eta||^2 + 1/2 ||y - Dx||^2
function x = x_update(Q, A_CG, W_CG, P_CG, y, u, x, mu, R, D, eta, eigvalsrr, eigvalsdd, Nx, Ny, arg)
switch upper(arg.zmethod)
        case 'CG'
                Y = [y; u - eta];
                try
                        [xs, info] = qpwls_pcg1(double(x), A_CG, W_CG, double(Y), Gdiag(zeros(Nx*Ny, 1)), ...
                                'niter', arg.inner_iter, 'stop_grad_tol', 1e-13, 'precon', P_CG);
			x = xs;
                catch
                        display('qpwls failed');
                        keyboard
                end
        case 'FFT'
                rhs = reshape(mu * R' * (u - eta) + D' * y, Nx, Ny);
                invMat = mu * eigvalsrr + eigvalsdd;
                x = Q' * col((Q * rhs(:)) ./ col(invMat)) / (Nx*Ny);
        otherwise
                display(sprintf('unknown option for z-update: %s', arg.zmethod));
end
% keyboard
end
