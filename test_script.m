% tests out AL_P2_gen and tridiag_ADMM on simulated brainweb data

nx = 240;
ny = 200;
dims = [nx ny];
nc = 4;

img = make_sim_image(nx,ny);
sense_maps = mri_sensemap_sim('nx', nx, 'ny', ny, 'ncoil', nc, 'rcoil', 700 );
% samp = gen_poisson_sampling_pattern('2D',[nx ny], [floor(nx/6) floor(ny/6)],6); % CAUTION: concatenating real poisson disk samplings
params.Nx = nx;
params.Ny = ny;
params.Nf = 1;
params.h = 1;
params.R = round(nx*ny/3);
samp = logical(gen_new_sampling_pattern(params, 'all_kspace', 0));

imgs = repmat(img,[1 1 nc]).*sense_maps;
% S = construct_sensefat(sense_maps);
S = GsplineS(sense_maps, 1);

% subsampling Fourier encoding matrix F
% F = construct_fourierfat(samp,nc);
F = GsplineF(nx, ny, 1, nc, 'samp', samp);

y = F*(S*img(:));
snr = 40; % specify desired SNR (of sampled values!) in dB
snr = 80;
%sig = 10^(-snr/20) * norm(ytrue(samp)) / sqrt(sum(samp(:)));
sig = 10^(-snr/20) * norm(y) / sqrt(length(y));
noisey = y + sig*randn(size(y)) + i*sig*randn(size(y));
% finite-differencing matrices C1 and C2
[C1, C2] = construct_finite_diff(dims);

% choose parameters

% beta, spatial regularization parameter
niters = 50;

beta = 1.85;

% alpha, tradeoff between S
alph = 0.5;
% mu, convergence parameters
mu(1) = 1;
mu(2) = 1;
mu(3) = 4;%mu2
mu(4) = 4;
mu(5) = 4;
mu(6) = 1;
mu(7) = 1;
mu(8) = 1;
mu(9) = 1;

% [k_u,k_v,k_x] = diag_cond_numbers(eigvalsmssm, eigvalsrr, Nx*Ny, Q, R, M, S, A, u_u, u_v, u_z);

load([home 'Dropbox/fessler/experimental_data/tridiag_vega/xinf_via_ALP2_xinf_only'])

xzfill = S'*(F'*noisey);
xzfill = ir_wls_init_scale((F*S), noisey, xzfill);
xinit = xzfill;%zeros(nx,ny);

%%

% ALP2
[xhat_P2, xsaved_P2, err_P2, costOrig_P2] = AL_P2_gen(noisey, F, S, [C1; C2],...
        xinit, niters, beta, img);

for iter = 1:niters
        dist_to_xinf_P2(iter) = sqrt(sum(abs(col(xsaved_P2(:,:,iter)) - xinf(:)).^2)/(nx*ny));
end

% tridiag solver
[xhat_tri, xsaved_tri, cost_tri] = tridiag_ADMM(noisey, F, S, C1, C2, ...
        alph, beta, xinit, img, niters); % xtrue should be xinf

for iter = 1:niters
        dist_to_xinf_tri(iter) = sqrt(sum(abs(col(xsaved_tri(:,:,iter)) - xinf(:)).^2)/(nx*ny));
end

% % "parallel"
% [xhat_FP, xsaved_FP, cost_tri] = SENSE_FP(noisey, F, S, C1, C2, alph, beta, ...
%         mu, nx, ny, xinit, niters);
% 
% 
% for iter = 1:niters
%         dist_to_xinf_FP(iter) = sqrt(sum(abs(col(xsaved_FP(:,:,iter)) - xinf(:)).^2)/(nx*ny));
% end


return
save('races_results.mat','beta','alph','err','xhat_P2','xsaved_P2','dist_to_xinf_P2','xhat_tri','xsaved_tri','dist_to_xinf_tri','xhat_FP','xsaved_tri','dist_to_xinf_FP','img','noisey','F','S','C1','C2','mu','snr');

