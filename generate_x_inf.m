% generate xinf for tridiag experiment
exp_setup;

niters_inf = 50000;

do_tri = 0;
do_circ = 1;
do_MFISTA = 0;

if do_tri
	[xhat_tri, ~, ~, costOrig_tri, time_tri] = ADMM_tridiag(y_noise, F, S, CH, CV, beta, SoS, zeros(size(SoS)), niters_inf, 'mu', plain_mu, 'save_progress', 'tridiag_xinf');
	x_tri_inf = xhat_tri;
	save(sprintf('%s/x_tri_inf_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta), 'x_tri_inf', 'niters_inf');
end
if do_circ
	display('doing MFISTA for periodic boundary conditions!')
	[xMFISc, CMFISc, TFISc, ~, ~] = MFISTA_wrapper(Nx, Ny, R, y_noise, xinit, F, S, beta, niters_inf, curr_folder);
	save(sprintf('%s/x_MFISTA_circ_inf_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta), 'xMFISc', 'niters_inf', 'CMFISc', 'TFISc');
end

if 0%1
	[xhat_alp2, ~, ~, costOrig_alp2, time_alp2] = AL_P2_gen(y_noise, F, S, Rcirc, SoS, niters_inf, beta, zeros(size(SoS)), 'inner_iter', 3);
	x_alp2_inf = xhat_alp2;
	save(sprintf('%s/x_alp2_inf_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta),'x_alp2_inf', 'niters_inf');
end
if wavelets

	[xhat_triw, ~, ~, costOrig_triw, time_triw] = ADMM_tridiag_W(y_noise, F, S, CH, CV, beta, SoS, zeros(size(SoS)), niters_inf, 'mu', mu);
	x_triw_inf = xhat_triw;
	save(sprintf('%s/x_triw_inf_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta), 'x_triw_inf', 'niters_inf');


	[xhat_alp2w, ~, ~, costOrig_alp2w, time_alp2w] = AL_P2_gen(y_noise, F, S, RcircW, SoS, niters_inf, beta, zeros(size(SoS)), 'mu', mu, 'zmethod', 'FFT');
	x_alp2cw_inf = xhat_alp2w;
	save(sprintf('%s/x_alp2cw_inf_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta),'x_alp2cw_inf', 'niters_inf');
end
if do_circ && do_MFISTA
	[xMFIS_circ, CMFIS_circ, TFIS_circ, ~, ~] = MFISTA_wrapper(Nx, Ny, Rcirc, y_noise, xinit, F, S, beta, niters_inf, curr_folder);
	save(sprintf('%s/x_MFISTA_circ_inf_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta), 'xMFIS_circ', 'niters_inf');
end

if do_MFISTA
	[xMFIS, CMFIS, TFIS, ~, ~] = MFISTA_wrapper(Nx, Ny, R, y_noise, xinit, F, S, beta, niters_inf, curr_folder);
	save(sprintf('%s/x_MFISTA_inf_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta), 'xMFIS', 'niters_inf', 'CMFIS', 'TFIS');
end

send_mai_text('done with xinf, now run timing tests')
return;
niters_inf = 10000;
	[xMFIS_sc, CMFIS_sc, TFIS_sc, ~, ~] = MFISTA_wrapper(Nx, Ny, R, y_noise, xinit, F, S_sc, beta_sc, niters_inf, curr_folder);

	save(sprintf('%s/x_MFISTA_inf_%s_beta%.*d_scaled.mat', curr_folder, slice_str, 3, beta), 'xMFIS', 'niters_inf');


