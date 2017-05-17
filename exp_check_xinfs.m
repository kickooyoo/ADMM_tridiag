% checking if MFISTA, ADMM-tridiag, and AL-P2 converge to the same solution

exp_setup;
niters = 25000;

load(sprintf('./reviv/curr/x_tri_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x_tri_inf');

% cost function for MFISTA lacks 1/2
load(sprintf('./reviv/curr/x_MFISTA_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x*MFIS*');
	if isvar('xMFIS')
		x_MFISTA = xMFIS;
	end
        xMFIS = x_MFISTA;

load(sprintf('./reviv/curr/x_alp2c_inf_slice%d_beta%.*d.mat', slice, 3, beta),'x_alp2c_inf');

load(sprintf('./reviv/curr/x_alp2_inf_slice%d_beta%.*d.mat', slice, 3, beta),'x_alp2_inf');

[xhat_alp2, ~, nrmsd_alp2, costOrig_alp2, time_alp2] = AL_P2_gen(y_noise, F, S, R, x_tri_inf, niters, beta, x_alp2_inf,'inner_iter',3, 'mask', mask);
[xhat_alp2_selfinit, ~, nrmsd_alp2_selfinit, costOrig_alp2_selfinit, time_alp2_selfinit] = AL_P2_gen(y_noise, F, S, R, x_alp2_inf, niters, beta, x_alp2_inf,'inner_iter',3, 'mask', mask);
[xhat_tri, ~, nrmsd_tri, costOrig_tri, time_tri] = ADMM_tridiag(y_noise, F, S, CH, CV, beta, xMFIS, x_tri_inf, niters, 'mask', mask, 'mu', plain_mu);
[xhat_tri_selfinit, ~, nrmsd_tri_self_init, costOrig_tri_selfinit, time_tri_selfinit] = ADMM_tridiag(y_noise, F, S, CH, CV, beta, x_tri_inf, x_tri_inf, niters, 'mask', mask, 'mu', plain_mu);
[xhat_tri_alp2init, ~, nrmsd_tri_alp2init, costOrig_tri_alp2init, time_tri_alp2init] = ADMM_tridiag(y_noise, F, S, CH, CV, beta, x_alp2_inf, x_tri_inf, niters, 'mask', mask, 'mu', plain_mu);
[xhat_MFIS, CMFIS, TFIS, l2DFIS, RMSEFIS] = MFISTA_wrapper(Nx, Ny, R, y_noise, x_tri_inf, F, S, beta, niters, curr_folder);

save(sprintf('./reviv/curr/check_xinfs_%dx%d_%diter_%s_revis.mat', Nx, Ny, niters, slice_str));
send_mai_text('done with mpel8 timing');

display('DONE');
