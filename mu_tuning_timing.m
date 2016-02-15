% experimental data for tridiag_ADMM

tridiag_exp_setup;
niters = 1000;

load(sprintf('./reviv/x_MFISTA_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x*MFIS*');
if isvar('xMFIS')
        x_MFISTA = xMFIS;
end

load(sprintf('./reviv/x_tri_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x_tri_inf');

for ii = 1:2
	if ii == 1
		x_tri_inf = x_MFISTA;
	end
	nthread = int32(160);
	[xhat_trit, ~, nrmsd_trit, costOrig_trit, time_trit] = tridiag_ADMM(y_noise, ...
		F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask, 'nthread', nthread);
	[xhat_trit2, ~, nrmsd_trit2, costOrig_trit2, time_trit2] = tridiag_ADMM(y_noise, ...
		F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask, 'fancy_mu34', 1, 'nthread', nthread);
	[xhat_trit2_80, ~, nrmsd_trit2_80, costOrig_trit2_80, time_trit2_80] = tridiag_ADMM(y_noise, ...
		F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask, 'fancy_mu34', 1, 'nthread', int32(80));

	[xhat_tri, ~, nrmsd_tri, costOrig_tri, time_tri] = tridiag_ADMM(y_noise, ...
		F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask, 'mu', plain_mu, 'nthread', nthread);

	[xhat_alp2t, ~, nrmsd_alp2t, costOrig_alp2t, time_alp2t] = AL_P2_gen(y_noise, ...
		F, S, R, xinit, niters, beta, x_tri_inf,'inner_iter', 1, 'mask', mask);

	[xhat_trit2p, ~, nrmsd_trit2p, costOrig_trit2p, time_trit2p] = tridiag_ADMM(y_noise, ...
		F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask, 'parFFT', true, 'fancy_mu34', 1, 'nthread', nthread);
	[xhat_alp2tp, ~, nrmsd_alp2tp, costOrig_alp2tp, time_alp2tp] = AL_P2_gen(y_noise,  ...
		F, S, R, xinit, niters, beta, x_tri_inf,'inner_iter', 1, 'mask', mask, 'parFFT', true);
	if ii == 1 
		save(sprintf('./reviv/mu_test_MFISTAinf_%dx%d_%diter_%dslice.mat', Nx, Ny, niters, slice));
	else
		save(sprintf('./reviv/mu_test_triinf_%dx%d_%diter_%dslice.mat', Nx, Ny, niters, slice));
	end
end
send_mai_text('done with mu exp timing');

plot_timing;

display('DONE');
