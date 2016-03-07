% generate xinf for tridiag experiment
exp_setup;

niters_inf = 50000;

do_tri = 1;
do_circ = 0;
do_MFISTA = 1;

tri_inf_fname = sprintf('%s/x_tri_inf_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta);
if do_tri && ~exist(tri_inf_fname, 'file')
	[x_tri_inf, ~, nrmsd_tri_inf, costOrig_tri_inf, time_tri_inf] = ADMM_tridiag(y_noise, F, S, CH, CV, beta, xinit, zeros(size(SoS)), niters_inf, 'save_progress', 'tridiag_xinf');
	save(tri_inf_fname, 'x_tri_inf', 'nrmsd_tri_inf', 'costOrig_tri_inf', 'time_tri_inf', 'niters_inf');
end
if do_circ
	display('doing MFISTA for periodic boundary conditions!')
	[xMFISc, CMFISc, TFISc, ~, ~] = MFISTA_wrapper(Nx, Ny, R, y_noise, xinit, F, S, beta, niters_inf, curr_folder);
	save(sprintf('%s/x_MFISTA_circ_inf_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta), 'xMFISc', 'niters_inf', 'CMFISc', 'TFISc');
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

MFISTA_inf_fname = sprintf('%s/x_MFISTA_inf_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta);
if do_MFISTA && ~exist(MFISTA_inf_fname, 'file')
	[xMFIS, CMFIS, TFIS, ~, ~] = MFISTA_wrapper(Nx, Ny, R, y_noise, xinit, F, S, beta, niters_inf, curr_folder);
	save(MFISTA_inf_fname, 'xMFIS', 'niters_inf', 'CMFIS', 'TFIS');
end

send_mai_text(sprintf('done with xinf on %s, now run timing tests', machine(1:3)))


