% experimental data for ADMM_tridiag

exp_setup;
niters = 100;

if issim
	slice_str = 'sim';
else
	slice_str = sprintf('slice%d', slice);
end

load(sprintf('./reviv/curr/x_MFISTA_inf_%s_beta%.*d.mat', slice_str, 3, beta), 'x*MFIS*');
if isvar('xMFIS')
        x_MFISTA = xMFIS;
end

if ~isvar('save_suffix')
	save_suffix = '';
end
fudges = 10.^(-2);%:0.5:4);
ktris = 12;%  10:2:24; % 12 still best

if 0
for ii = 1:length(ktris)
	ktri = ktris(ii);
	[xhat_trif(:,:,ii), ~, nrmsd_trif(ii,:), costOrig_trif(ii,:), time_trif(ii,:)] = ADMM_tridiag(...
		y_noise, F, S, CH, CV, beta, xinit, x_MFISTA, niters, 'mask', mask, ...
		'mu_args', {'ktri', ktri});
end
end
[~, ~, nrmsd_tri_edge, costOrig_tri_edge, time_tri_edge] = ADMM_tridiag(...
	y_noise, F, S, CH, CV, beta, xinit, x_MFISTA, niters, 'mask', mask, ...
	'mu_args', {'test', 'edge'});%, 'fancy_mask', true);


[~, ~, nrmsd_tri_edgeRR, costOrig_tri_edgeRR, time_tri_edgeRR] = ADMM_tridiag(...
	y_noise, F, S, CH, CV, beta, xinit, x_MFISTA, niters, 'mask', mask, ...
	'mu_args', {'test', 'edgeRRapprox'});%, 'fancy_mask', true);

[~, ~, nrmsd_tri_RR, costOrig_tri_RR, time_tri_RR] = ADMM_tridiag(...
	y_noise, F, S, CH, CV, beta, xinit, x_MFISTA, niters, 'mask', mask, ...
	'mu_args', {'test', 'RRapprox'});

[~, ~, nrmsd_tri_01RR, costOrig_tri_01RR, time_tri_01RR] = ADMM_tridiag(...
	y_noise, F, S, CH, CV, beta, xinit, x_MFISTA, niters, 'mask', mask, ...
	'mu_args', {'fancy_mu01', true, 'test', 'RRapproxedge'});

[~, ~, nrmsd_alp2t, costOrig_alp2t, time_alp2t] = AL_P2_gen(y_noise, ...
	F, S, R, xinit, niters, beta, x_MFISTA,'inner_iter', 1, 'mask', mask);
save(sprintf('./reviv/curr/fancy_mu_test_MFISTAinf_%dx%d_%diter_%s%s.mat', ...
	Nx, Ny, niters, slice_str, save_suffix));
%send_mai_text('done with fancy mu exp timing');

plot_timing;

display('DONE');



