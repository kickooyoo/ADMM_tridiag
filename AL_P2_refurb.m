function [x, xsave, err, costOrig, time] = AL_P2_refurb(ylong, A, S, R, ...
        initx, niters, lambda, mu, nx, ny, xtrue, varargin)

arg.method = 'fft';
arg = vararg_pair(arg, varargin);

npix = nx*ny;
nc = length(S.arg.Nc);

y = ylong(:);
x = initx(:);

% soft thresholding function, t is input, a is threshold
soft = @(t,a) (t - a * sign(t)) .* (abs(t) > a);
eigvals1 = @(A,Q) Q*A(:,1);

% fatrix objects
Q = Gdft('mask', true(nx,ny));

% R: regularizer, currently
%R = [C1; C2];
SS = S'*S;
eigvalsss = SS*ones(nx*ny,1);
RR = R'*R;
eigvalsrr = Q*RR(:,1); % NOT CIRCULANT ANYMO'
[eigvalsaa,Qbig] = get_eigs(A,nc); % match, they are okay!
%eigvalsaa = repmat(Q*AA(:,1),nc,1);

% hand picked mask for computing NRMSE over mask
mask = true(nx,ny);%generate_mask('slice67',1,nx,ny);

u = S*x;
v = R*x;
z = x;

u_u = mu(1);
u_v = mu(2);
u_z = mu(3);
eta_u = zeros(size(u));
eta_v = zeros(size(v));
eta_z = zeros(size(z));

calc_errcost = 1;
if (calc_errcost)
        calc_orig_cost = @(y, A, S, R, x, lambda) norm(y - A*(S*x),2)^2/2 + ...
                lambda*norm(R*x,1);
        err = calc_NRMSE_over_mask(x, xtrue(:), mask);
        costAL = calc_AL_cost(A, S, R, y, x, u, v, z, u_u, u_v, u_z, ...
                eta_u, eta_v, eta_z, lambda);
        costOrig = calc_orig_cost(y, A, S, R, x, lambda);
end

xsave = zeros(nx, ny, niters);
% [k_u,k_v,k_x] = diag_cond_numbers(eigvalsmssm, eigvalsrr, npix, Q, R, M, S, A, u_u, u_v, u_z);
time = zeros(niters,1);
%tic
for ii=1:niters
        iter_start = tic;
        x = x_update(S, u, z, eta_u, eta_z, u_u, u_z, eigvalsss);
        if any(isnan(x))
                keyboard
        end
        u = u_update(Qbig, A, S, y, x, eta_u, u_u, npix, eigvalsaa);
        %v = soft(R*x + eta_v, lambda/u_v);
        v = soft(R*z + eta_v, lambda/u_v); % TEST THIS
        z = z_update(Q, v, x, z, u_v, u_z, R, eta_v, eta_z, eigvalsrr, nx, ny, arg);
        eta_u = eta_u - (u - S*x);
        eta_v = eta_v - (v - R*z);
        eta_z = eta_z - (z - x);
        time(ii) = toc(iter_start);
        if (calc_errcost)
                err = [err calc_NRMSE_over_mask(x, xtrue(:), mask)];
                costAL = [costAL calc_AL_cost(A, S, R, y, x, u, v, z, u_u, ...
                        u_v, u_z, eta_u, eta_v, eta_z, lambda)];
                costOrig = [costOrig calc_orig_cost(y, A, S, R, x, lambda)];
        end
        if mod(ii,10) == 0
                printf('%d/%d iterations',ii,niters)
        end
        xsave(:,:,ii) = reshape(x, nx, ny);
end
%toc
if (~calc_errcost)
        err = zeros(niters+1, 1);
        costOrig = zeros(niters+1, 1);
end
% keyboard
end

function x = x_update(S, u, z, eta_u, eta_z, u_u, u_z, eigvalsss)
inv_vals = (u_u*eigvalsss + u_z);
rhs = u_u*S'*(u-eta_u) + u_z*(z-eta_z);
x = rhs./inv_vals;
%xtest = inv(u_u*diag(eigvalsss) + u_z*eye(numel(x)))*rhs;
%display(sprintf('x diff: %d', norm(x - xtest)/numel(x)));
end

function u = u_update(Qbig, A, S, y, x, eta_u, u_u, N, eigvalsaa)
invMat = eigvalsaa + u_u;
latter = (A'*y + u_u*(S*x + eta_u));
u = Qbig'*((Qbig*latter)./invMat)/N;
% utest = inv(u_u*eye(numel(u)) + full(A'*A))*latter;
% display(sprintf('u diff: %d', norm(u - utest)/numel(u)));
end

% Conjugate Gradient Descent!
% z = argmin (u_v R'R + u_z I)^(-1)(u_v R'(v - etav) + u_z (x + etaz))
% z = argmin u_?/2 ||v - Rz - etav||^2+u_?/2 ||z - x - etaz||^2
function z = z_update(Q, v, x, z, u_v, u_z, R, eta_v, eta_z, eigvalsrr, nx, ny, arg)

% arg.method = 'CG'; % CG
% switch arg.method
%         case 'CG'
                n1 = length(v);
                n2 = length(x);
                W = Gdiag([u_v*ones(1,n1) u_z*ones(1,n2)]);
                A = [Gmatrix(R); Gmatrix(Gdiag(ones(1,n2)))]; % nightmare, pass in
                y = [v-eta_v; x+eta_z];
                try
                        [z_pcg, info] = qpwls_pcg1(z, A, W, y, Gdiag(zeros(n2,1)), ...
                                'niter', 20, 'stop_grad_tol', 1e-11, 'precon', A'*A);
                catch
                        display('qpwls failed');
                        keyboard
                end
%         case 'fft'
                rhs = reshape(u_v*R'*(v - eta_v) + u_z*(x + eta_z), nx, ny);
                %z_fft = ifft2((fft2(rhs))./(reshape(u_v*eigvalsrr+u_z*ones(nx*ny,1),nx,ny)));
                %z_fft = z_fft(:);
                invMat = u_v*eigvalsrr + u_z;
                z_fft = Q'*((Q*rhs(:))./invMat)/(nx*ny);
                
%         otherwise
%                 display(sprintf('unknown option for z-update: %s', method));
% end
display(sprintf('z diff: %d', norm(z_pcg - z_fft)/numel(z)));
z = z_fft;
end

function cost = calc_AL_cost(A, S, R, y, x, u, v, z, u_u, u_v, u_z, eta_u, ...
        eta_v, eta_z, lambda)
cost = norm(y - A*u)^2/2 + lambda*norm(v,1) + u_u*norm(u - S*x - eta_u)^2/2 + ...
        u_v*norm(v - R*z - eta_v)^2/2 + u_z*norm(z - x - eta_z)^2/2;
end
