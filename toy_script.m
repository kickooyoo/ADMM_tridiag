% tests out AL_P2_gen and tridiag_ADMM on simulated brainweb data

Nx = 24;
Ny = 20;
dims = [Nx Ny];
Nc = 4;

img = make_sim_image(Nx, Ny);
sense_maps = mri_sensemap_sim('nx', Nx, 'ny', Ny, 'ncoil', Nc, 'rcoil', 3*Ny);
% samp = gen_poisson_sampling_pattern('2D',[Nx Ny], [floor(Nx/6) floor(Ny/6)],6); 
params.Nx = Nx;
params.Ny = Ny;
params.Nf = 1;
params.h = 1;
params.R = round(Nx*Ny/4);
samp = logical(gen_new_sampling_pattern(params, 'all_kspace', 0));

imgs = repmat(img,[1 1 Nc]).*sense_maps;
% S = construct_sensefat(sense_maps);
S = staticS(sense_maps);

% subsampling Fourier eNcoding matrix F
% F = construct_fourierfat(samp,Nc);
F = staticF(Nx, Ny, Nc, 'samp', samp);

SFFS = full(S'*F'*F*S);
[SFeigs, ~] = get_eigs(S'*F'*F*S, 1, 'Nx', Nx, 'Ny', Ny);
[Feigs, Q] = get_eigs(F'*F, Nc, 'Nx', Nx, 'Ny', Ny);
[V, D] = eig(SFFS);
%%
Rcirc = Cdiffs([Nx Ny],'offsets', [1 Nx], 'type_diff','circshift');
RR = full(Rcirc'*Rcirc);
[Reigs, Q] = get_eigs(Rcirc'*Rcirc, 1, 'Nx', Nx, 'Ny', Ny);
VRRV = V*RR*V';
VLRV = V*diag(Reigs)*V';
figure; im(VLRV)
%% 
params.N = [Nx Ny];
params.A = F*S;
params.W = ones(prod(F.odim), 1);
params.eigtol = eps; % Matlab epsilon
params.eigpowermaxitr = 10000;
params.eigdispitr = 10;
mEAWA = get_MaxEigenValue_1by1(params, 'AWA');

return
%% 

rng('default');
y = F*(S*img(:));
snr = 40; % specify desired SNR (of sampled values!) in dB
snr = 80;
%sig = 10^(-snr/20) * norm(ytrue(samp)) / sqrt(sum(samp(:)));
sig = 10^(-snr/20) * norm(y) / sqrt(length(y));
noisey = y + sig*randn(size(y)) + i*sig*randn(size(y));
% finite-differeNcing matrices C1 and C2
[C1, C2] = construct_finite_diff(dims);

% choose parameters

% beta, spatial regularization parameter
niters = 300;

beta = 22.85;

alph = 0.5;

xzfill = S'*(F'*noisey);
xzfill = ir_wls_init_scale((F*S), noisey, xzfill);
xinit = reshape(xzfill, Nx, Ny);
pot = potential_fun('quad');

R = [C1; C2];
curr_folder = '.';
[xMFIS, CMFIS, TFIS, l2DFIS, RMSEFIS] = MFISTA_wrapper(Nx, Ny, R, noisey, xinit, F, S, beta, niters, curr_folder);

return
%%
% A_bs = [sqrt(1/2)*full(F*S); sqrt(beta)*full(R)];
% y_bs = [sqrt(1/2)*col(y); zeros(size(R*xinit(:)))];
% xhat_bs = reshape(A_bs\y_bs, Nx, Ny);

% ALP2
% [xhat_P2, ~, err_P2, costOrig_P2] = AL_P2_gen_genpot(noisey, F, S, R,...
%         xinit, niters, beta, img, 'pot', pot, 'inner_iter', 3, 'mu', plain_mus(1:3));

% tridiag solver
% [xhat_tri, xsaved_tri, cost_tri] = tridiag_ADMM_genpot(noisey, F, S, C1, C2, ...
%         alph, beta, xinit, img, niters, 'pot', pot, 'faNcy_mu', 1); 

[xhat_tri, xsaved_tri, cost_tri] = ADMM_tridiag(noisey, F, S, C1, C2, ...
        beta, xinit, img, niters, 'mu_args', {'faNcy_mu34', true, 'test', 'edgeRRapproxcheck'});

figure; subplot(2,2,1); im(xinit); subplot(2,2,2); im(xhat_bs); subplot(2,2,3); im(xhat_P2); subplot(2,2,4); im(xhat_tri)
calc_NRMSE_over_mask(xhat_bs, xhat_tri, true(Nx,Ny))

return
save('races_results.mat','beta','alph','err','xhat_P2','xsaved_P2','dist_to_xinf_P2','xhat_tri','xsaved_tri','dist_to_xinf_tri','xhat_FP','xsaved_tri','dist_to_xinf_FP','img','noisey','F','S','C1','C2','mu','snr');

