function [x, xsave, err, costOrig, time] = AL_P2_gen(y, A, S, R, ...
        xinit, niters, lambda, xtrue, varargin)
% function [x, xsave, err, costOrig, time] = AL_P2_gen(y, A, S, R, ...
%         xinit, niters, lambda, xtrue, varargin)
% 
% implements AL-P2 with option for CG on x-update for noncirc R
% inputs:
%       y [Ns*Nc 1] undersampled data vector
%       A (fatrix2) undersampled Fourier encoding matrix
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
%       maxv (real scalar) debugging, stops program if diverges (i.e. any
%               values in x > maxv
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
% 01/26/2016 Mai Le
% University of Michigan
Nx = size(xinit,1);
Ny = size(xinit,2);

arg.zmethod = 'CG';
arg.maxv = Inf;
arg.mu = [];
arg.mask = true(Nx, Ny);
arg.inner_iter = 1;
arg.debug = false;
arg.precon = true;
arg.save_x = false;
arg.parFFT = false;
arg.pot = [];
arg = vararg_pair(arg, varargin);

Nc = S.arg.Nc;

if arg.parFFT
        tmp = gcp('nocreate');
        if prod(size(tmp)) == 0
                parpool(Nc)
        end
end

tic;
% eigvals for SS, get mus
SS = S'*S;
eigvalsss = SS * ones(Nx * Ny, 1);
RR = R'*R;
Q = Gdft('mask', true(Nx, Ny));
e0 = zeros(Nx, Ny);
e0(1, 1) = 1; %end/2,end/2) = 1;
eigvalsrr = reshape(Q * RR * e0(:), Nx, Ny); % approximate, for use in precon
if isempty(arg.mu)
	arg.mu = get_mu(eigvalsss(arg.mask(:)), eigvalsrr, Nx*Ny, lambda, 'split', 'AL-P2');
end

if length(arg.mu) ~= 3
	display('wrong size for mu convergence parameters');
	keyboard;
end

y = single(y(:));
x = single(xinit(:));
u = S*x;
v = R*x;
z = x;
u_u = arg.mu{1};
u_v = arg.mu{2};
u_z = arg.mu{3};
eta_u = zeros(size(u));
eta_v = zeros(size(v));
eta_z = zeros(size(z));

% shrinkage (e.g. soft thresholding) function, t is input, a is threshold
if isempty(arg.pot)
	arg.pot = potential_fun('l1', lambda/u_v);
end
shrink = @(a, t) arg.pot.shrink(a, t);

if strcmp(upper(arg.zmethod), 'FFT')
	A_CG = [];
	W_CG = [];
	P_CG = [];
else
	Qn = (1/sqrt(Nx * Ny)) * Gdft('mask', true(Nx, Ny));
	A_CG = [R; Gdiag(ones(Nx*Ny,1),'mask', true(Nx, Ny))];
	if isfield(R.arg, 'dim')
		Rodim = R.dim(1);
	elseif isfield(R, 'odim')
		Rodim = prod(R.odim);
	else
		display('can not figure out odim of R');
		keyboard;
	end
        W_CG = Gdiag([u_v * ones(1, Rodim) u_z * ones(1, Nx*Ny)]);
	if arg.precon
		P_CG = qpwls_precon('circ0', {A_CG, W_CG}, Gdiag(zeros(Nx*Ny,1)), true(Nx, Ny)); 
	else
		P_CG = 1;
	end
end

[eigvalsaa, Qbig] = get_eigs(A, Nc); 

calc_errcost = 1;
if (calc_errcost)
        calc_orig_cost = @(y, A, S, R, x, lambda) norm(col(y - A*(S*x)),2)^2/2 + ...
                lambda* sum(col(arg.pot.potk(col(R*x))));
        err = calc_NRMSE_over_mask(x, xtrue, arg.mask);
        costOrig = calc_orig_cost(y, A, S, R, x, lambda);
end

if arg.save_x
	xsave = zeros(Nx, Ny, niters);
else
	xsave = 0;
end
time(1) = toc;
for ii = 1:niters
        iter_start = tic;
        x = x_update(S, u, z, eta_u, eta_z, u_u, u_z, eigvalsss);
        if any(isnan(x)) || any(isinf(x)) || any(abs(x) > arg.maxv)
                keyboard
        end
        u = u_update(Qbig, A, S, y, x, eta_u, u_u, Nx*Ny, eigvalsaa, arg.parFFT);
        v = shrink(R*z + eta_v, lambda/u_v); 
        z = z_update(Q, A_CG, W_CG, P_CG, v, x, z, u_v, u_z, R, eta_v, eta_z, eigvalsrr, Nx, Ny, arg);
        if any(isnan(z)) || any(isinf(z)) || any(abs(z) > arg.maxv)
                keyboard
        end
	if arg.debug
		subplot(2,2,1); im(reshape(abs(x), Nx, Ny));
		subplot(2,2,2); im(reshape(abs(u), Nx, Ny, Nc));
		subplot(2,2,3); im(reshape(abs(v), Nx, Ny, 2));
		subplot(2,2,4); im(reshape(abs(z), Nx, Ny));
		drawnow;
		pause(1);
	end

        eta_u = eta_u - (u - S*x);
        eta_v = eta_v - (v - R*z);
        eta_z = eta_z - (z - x);
        time = [time toc(iter_start)];
        if (calc_errcost)
                err(ii + 1) = calc_NRMSE_over_mask(x, xtrue, arg.mask);
                costOrig(ii + 1) = calc_orig_cost(y, A, S, R, x, lambda);
        end
        if mod(ii,10) == 0
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

function x = x_update(S, u, z, eta_u, eta_z, u_u, u_z, eigvalsss)
Hx = u_u * eigvalsss + u_z;
rhs = u_u * S' * (u - eta_u) + u_z * (z - eta_z);
x = rhs ./ Hx;
end

function u = u_update(Qbig, A, S, y, x, eta_u, u_u, N, eigvalsaa, par)
Hu = eigvalsaa + u_u;
rhs = A' * y + u_u * (S * x + eta_u);
if par
        rhs = reshape(rhs, A.arg.Nx, A.arg.Ny, A.arg.Nc);
        parfor ii = 1:A.arg.Nc
                U(:,:,ii) = fft2(rhs(:,:,ii));
        end
        U = U./reshape(Hu, A.arg.Nx, A.arg.Ny, A.arg.Nc);
        parfor ii = 1:A.arg.Nc
                u(:,:,ii) = ifft2(U(:,:,ii));
        end
        u = col(u);
else
        u = Qbig' * ((Qbig * rhs) ./ Hu) / N;
end
end

% z-update
% problematic for non-circulant R
% z = argmin (u_v R'R + u_z I)^(-1)(u_v R'(v - etav) + u_z (x + etaz))
% z = argmin u_v/2 ||v - Rz - etav||^2 + u_z/2 ||z - x - etaz||^2
function z = z_update(Q, A, W, P, v, x, z, u_v, u_z, R, eta_v, eta_z, eigvalsrr, Nx, Ny, arg)
switch upper(arg.zmethod)
        case 'CG'
                y = [v - eta_v; x + eta_z];
                try
                        z = qpwls_pcg1(z, A, W, y, Gdiag(zeros(Nx*Ny, 1)), ...
                                'niter', arg.inner_iter, 'stop_grad_tol', 1e-13, 'precon', P);
                catch
                        display('qpwls failed');
                        keyboard
                end
        case 'FFT'
                rhs = reshape(u_v * R' * (v - eta_v) + u_z * (x + eta_z), Nx, Ny);
                invMat = u_v * eigvalsrr + u_z;
                z = Q' * col((Q * rhs(:)) ./ col(invMat)) / (Nx*Ny);
        otherwise
                display(sprintf('unknown option for z-update: %s', arg.zmethod));
end
end
