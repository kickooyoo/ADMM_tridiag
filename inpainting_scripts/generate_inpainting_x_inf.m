% generate xinf for tridiag experiment
niters_inf = 50000;

do_MFISTA = 1;
if do_MFISTA
	[xMFIS, CMFIS, TFIS, ~, ~] = MFISTA_inpainting_wrapper(Nx, Ny, RW, y, xinit, D, beta, niters_inf, curr_folder, slice_str);
	save(MFISTA_inf_fname, 'xMFIS', 'niters_inf', 'CMFIS', 'TFIS');
end

do_ADMM = 1;
if do_ADMM
	[x_ADMM_inf, ~, ~, cost_ADMM_inf, time_ADMM_inf] = ADMM_tridiag_inpaint(y, D, CHW, CVW, beta, xinit, zeros(size(xinit)), niters_inf, 'betaw', betaw, 'alphw', alphw, 'alph', alph, 'kapp', 7, 'pos', 0.01, 'save_progress', ADMM_inf_fname);
	save(ADMM_inf_fname, 'x_ADMM_inf', 'niters_inf', 'cost_ADMM_inf', 'time_ADMM_inf');
end

do_ALP2 = 1;
if do_ALP2
	[x_ALP2_inf, ~, ~, cost_ALP2_inf, time_ALP2_inf] = AL_P2_inpainting(y, D, RW, xinit, niters_inf, beta, zeros(size(xinit)), 'betaw', betaw, 'alphw', alphw, 'kapp', 7, 'inner_iter', 3, 'save_progress', ALP2NC_inf_fname);
	save(ALP2NC_inf_fname, 'x_ALP2_inf', 'niters_inf', 'cost_ALP2_inf', 'time_ALP2_inf');
end

do_SB = 1;
if do_SB
	[x_SB_inf, ~, ~, cost_SB_inf, time_SB_inf] = SB_inpainting(y, D, RW, xinit, niters_inf, beta, zeros(size(xinit)), 'betaw', betaw, 'alphw', alphw, 'kapp', 7, 'inner_iter', 3);%, 'save_progress', SBNC_inf_fname);
	save(SBNC_inf_fname, 'x_SB_inf', 'niters_inf', 'cost_SB_inf', 'time_SB_inf');
end

send_mai_text(sprintf('done with xinf on %s, now run timing tests', machine(1:3)))


