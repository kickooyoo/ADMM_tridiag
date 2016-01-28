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
%       maxv (real scalar) debugging, stops program if diverges (i.e. any
%               values in x > maxv
% outputs:
%       x [Nx Ny] reconstructed image
%       xsave [Nx Ny niters] estimated image at each iteration
%       err [niters 1] NRMSE at each iteration
%       costOrig [niters 1] objective value of original cost func at each iter
%       time [niters 1] wall time per each iteration
%
% 01/26/2016 Mai Le
% University of Michigan
Nx = size(xinit,1);
Ny = size(xinit,2);

arg.zmethod = 'CG';
arg.maxv = Inf;
arg.mu = ones(1,3);
arg.mask = true(Nx, Ny);
arg = vararg_pair(arg, varargin);

Nc = length(S.arg.Nc);

y = single(y(:));
x = single(xinit(:));
u = S*x;
v = R*x;
z = x;
u_u = arg.mu(1);
u_v = arg.mu(2);
u_z = arg.mu(3);
eta_u = zeros(size(u));
eta_v = zeros(size(v));
eta_z = zeros(size(z));

% soft thresholding function, t is input, a is threshold
soft = @(t,a) (t - a * sign(t)) .* (abs(t) > a);
eigvals1 = @(A,Q) Q * A(:,1);

% fatrix objects

SS = S'*S;
eigvalsss = SS * ones(Nx * Ny, 1);
if strcmp(upper(arg.zmethod), 'FFT')
        RR = R'*R;
        Q = Gdft('mask', true(Nx,Ny));
        eigvalsrr = Q*RR(:,1); 
	A_CG = [];
	W_CG = [];
else
        eigvalsrr = [];
        Q = [];
	A_CG = [R; Gdiag(ones(Nx*Ny,1),'mask', true(Nx, Ny))];
        W_CG = Gdiag([u_v * ones(1, prod(R.odim)) u_z * ones(1, Nx*Ny)]);
end
[eigvalsaa, Qbig] = get_eigs(A, Nc); 

calc_errcost = 1;
if (calc_errcost)
        calc_orig_cost = @(y, A, S, R, x, lambda) norm(y - A*(S*x),2)^2/2 + ...
                lambda*norm(R*x,1);
        err = calc_NRMSE_over_mask(x, xtrue, arg.mask);
        costOrig = calc_orig_cost(y, A, S, R, x, lambda);
end

xsave = zeros(Nx, Ny, niters);
time = zeros(niters,1);

for ii=1:niters
        iter_start = tic;
        x = x_update(S, u, z, eta_u, eta_z, u_u, u_z, eigvalsss);
        if any(isnan(x)) || any(isinf(x)) || any(abs(x) > arg.maxv)
                keyboard
        end
        u = u_update(Qbig, A, S, y, x, eta_u, u_u, Nx*Ny, eigvalsaa);
        v = soft(R*z + eta_v, lambda/u_v); 
        z = z_update(Q, A_CG, W_CG, v, x, z, u_v, u_z, R, eta_v, eta_z, eigvalsrr, Nx, Ny, arg);
        if any(isnan(z)) || any(isinf(z)) || any(abs(z) > arg.maxv)
                keyboard
        end
        eta_u = eta_u - (u - S*x);
        eta_v = eta_v - (v - R*z);
        eta_z = eta_z - (z - x);
        time(ii) = toc(iter_start);
        if (calc_errcost)
                err = [err calc_NRMSE_over_mask(x, xtrue, arg.mask)];
                costOrig = [costOrig calc_orig_cost(y, A, S, R, x, lambda)];
        end
        if mod(ii,10) == 0
                printf('%d/%d iterations',ii,niters)
        end
        xsave(:,:,ii) = reshape(x, Nx, Ny);
end
if (~calc_errcost)
        err = zeros(niters+1, 1);
        costOrig = zeros(niters+1, 1);
end
end

function x = x_update(S, u, z, eta_u, eta_z, u_u, u_z, eigvalsss)
Hx = u_u * eigvalsss + u_z;
rhs = u_u * S' * (u - eta_u) + u_z * (z - eta_z);
x = rhs ./ Hx;
end

function u = u_update(Qbig, A, S, y, x, eta_u, u_u, N, eigvalsaa)
Hu = eigvalsaa + u_u;
rhs = A' * y + u_u * (S * x + eta_u);
u = Qbig' * ((Qbig * rhs) ./ Hu) / N;
end

% z-update
% problematic for non-circulant R
% z = argmin (u_v R'R + u_z I)^(-1)(u_v R'(v - etav) + u_z (x + etaz))
% z = argmin u_v/2 ||v - Rz - etav||^2+u_z/2 ||z - x - etaz||^2
function z = z_update(Q, A, W, v, x, z, u_v, u_z, R, eta_v, eta_z, eigvalsrr, Nx, Ny, arg)
switch upper(arg.zmethod)
        case 'CG'
                y = [v - eta_v; x + eta_z];
                try
                        z = qpwls_pcg1(z, A, W, y, Gdiag(zeros(Nx*Ny, 1)), ...
                                'niter', 10, 'stop_grad_tol', 1e-11, 'precon', A'*A);
                catch
                        display('qpwls failed');
                        keyboard
                end
        case 'FFT'
                rhs = reshape(u_v * R' * (v - eta_v) + u_z * (x + eta_z), Nx, Ny);
                invMat = u_v * eigvalsrr + u_z;
                z = Q' * col((Q * rhs(:)) ./ invMat) / (Nx*Ny);
        otherwise
                display(sprintf('unknown option for z-update: %s', arg.zmethod));
end
end
