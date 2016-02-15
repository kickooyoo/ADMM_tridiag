% generate xinf for tridiag experiment
tridiag_exp_setup;

niters = 50000;

plain_mu = num2cell(ones(1,5));
if do_tri
	[xhat_tri, ~, ~, costOrig_tri, time_tri] = tridiag_ADMM(y_noise, F, S, CH, CV, alph, beta, SoS, zeros(size(SoS)), niters, 'mu', plain_mu, 'save_progress', 'tridiag_xinf');
	x_tri_inf = xhat_tri;
	save(sprintf('./reviv/x_tri_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x_tri_inf');
end
do_circ = 0;
if do_circ
	[xhat_alp2, ~, ~, costOrig_alp2, time_alp2] = AL_P2_gen(y_noise, F, S, Rcirc, SoS, niters, beta, zeros(size(SoS)), 'zmethod', 'FFT');
	x_alp2c_inf = xhat_alp2;
	save(sprintf('./reviv/x_alp2c_inf_slice%d_beta%.*d.mat', slice, 3, beta),'x_alp2c_inf');
end

if 0
	[xhat_alp2, ~, ~, costOrig_alp2, time_alp2] = AL_P2_gen(y_noise, F, S, Rcirc, SoS, niters, beta, zeros(size(SoS)), 'inner_iter', 3');
	x_alp2_inf = xhat_alp2;
	save(sprintf('./reviv/x_alp2_inf_slice%d_beta%.*d.mat', slice, 3, beta),'x_alp2_inf');
end
if wavelets

	[xhat_triw, ~, ~, costOrig_triw, time_triw] = tridiag_ADMM_W(y_noise, F, S, CH, CV, alph, beta, SoS, zeros(size(SoS)), niters, 'mu', mu);
	x_triw_inf = xhat_triw;
	send_mai_text('check xhat');
	keyboard;
	save(sprintf('./reviv/x_triw_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x_triw_inf');


	[xhat_alp2w, ~, ~, costOrig_alp2w, time_alp2w] = AL_P2_gen(y_noise, F, S, RcircW, SoS, niters, beta, zeros(size(SoS)), 'mu', mu, 'zmethod', 'FFT');
	x_alp2cw_inf = xhat_alp2w;
	send_mai_text('check xhat');
	keyboard;
	save(sprintf('./reviv/x_alp2cw_inf_slice%d_beta%.*d.mat', slice, 3, beta),'x_alp2cw_inf');
end
if do_circ && do_MFISTA
	[xMFIS_circ, CMFIS_circ, TFIS_circ, ~, ~] = MFISTA_wrapper(Nx, Ny, Rcirc, y_noise, xinit, F, S, beta, niters);
	save(sprintf('./reviv/x_MFISTA_circ_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'xMFIS_circ');
end

if do_MFISTA
	[xMFIS, CMFIS, TFIS, ~, ~] = MFISTA_wrapper(Nx, Ny, R, y_noise, xinit, F, S, beta, niters);
	save(sprintf('./reviv/x_MFISTA_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'xMFIS');
	send_mai_text('done with xinf, now run timing tests')
end

