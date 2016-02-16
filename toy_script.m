% tests out AL_P2_gen and tridiag_ADMM on simulated brainweb data

nx = 24;
ny = 20;
dims = [nx ny];
nc = 4;

img = make_sim_image(nx,ny);
sense_maps = mri_sensemap_sim('nx', nx, 'ny', ny, 'ncoil', nc, 'rcoil', 3*ny);
% samp = gen_poisson_sampling_pattern('2D',[nx ny], [floor(nx/6) floor(ny/6)],6); % CAUTION: concatenating real poisson disk samplings
params.Nx = nx;
params.Ny = ny;
params.Nf = 1;
params.h = 1;
params.R = round(nx*ny/2);
samp = logical(gen_new_sampling_pattern(params, 'all_kspace', 0));

imgs = repmat(img,[1 1 nc]).*sense_maps;
% S = construct_sensefat(sense_maps);
S = GsplineS(sense_maps, 1);

% subsampling Fourier encoding matrix F
% F = construct_fourierfat(samp,nc);
F = GsplineF(nx, ny, 1, nc, 'samp', samp);

rng(0);
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
niters = 300;

beta = 22.85;

alph = 0.5;

xzfill = S'*(F'*noisey);
xzfill = ir_wls_init_scale((F*S), noisey, xzfill);
xinit = reshape(xzfill, nx, ny);
pot = potential_fun('quad');

R = [C1; C2];
% A_bs = [sqrt(1/2)*full(F*S); sqrt(beta)*full(R)];
% y_bs = [sqrt(1/2)*col(y); zeros(size(R*xinit(:)))];
% xhat_bs = reshape(A_bs\y_bs, nx, ny);

% ALP2
% [xhat_P2, ~, err_P2, costOrig_P2] = AL_P2_gen_genpot(noisey, F, S, R,...
%         xinit, niters, beta, img, 'pot', pot, 'inner_iter', 3, 'mu', plain_mus(1:3));

% tridiag solver
% [xhat_tri, xsaved_tri, cost_tri] = tridiag_ADMM_genpot(noisey, F, S, C1, C2, ...
%         alph, beta, xinit, img, niters, 'pot', pot, 'fancy_mu', 1); 

[xhat_tri, xsaved_tri, cost_tri] = ADMM_tridiag(noisey, F, S, C1, C2, ...
        alph, beta, xinit, img, niters, 'fancy_mu', 1);

figure; subplot(2,2,1); im(xinit); subplot(2,2,2); im(xhat_bs); subplot(2,2,3); im(xhat_P2); subplot(2,2,4); im(xhat_tri)
calc_NRMSE_over_mask(xhat_bs, xhat_tri, true(nx,ny))

return
save('races_results.mat','beta','alph','err','xhat_P2','xsaved_P2','dist_to_xinf_P2','xhat_tri','xsaved_tri','dist_to_xinf_tri','xhat_FP','xsaved_tri','dist_to_xinf_FP','img','noisey','F','S','C1','C2','mu','snr');

