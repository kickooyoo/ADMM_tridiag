% experimental data for ADMM_tridiag

exp_setup;
niters = 700;

if issim
	slice_str = 'sim';
else
	slice_str = sprintf('slice%d', slice);
end

load(sprintf('./reviv/curr/x_MFISTA_inf_%s_beta%.*d.mat', slice_str, 3, beta), 'x*MFIS*');
if isvar('xMFIS')
        x_MFISTA = xMFIS;
end

load(sprintf('./reviv/curr/x_tri_inf_%s_beta%.*d.mat', slice_str, 4, beta), 'x_tri_inf');

if ~isvar('save_suffix')
	save_suffix = '';
end
fancy = 100; % fudge factor in get_mu
for ii = 1:2
	if ii == 1
		x_tri_inf = x_MFISTA;
	end
	nthread = int32(160);
	[xhat_trit2m, ~, nrmsd_trit2m, costOrig_trit2m, time_trit2m] = ADMM_tridiag(y_noise, ...
		F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask, 'fancy_mu34', 1, 'fancy_mask', fancy, 'nthread', nthread);
%	[xhat_trit, ~, nrmsd_trit, costOrig_trit, time_trit] = ADMM_tridiag(y_noise, ...
%		F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask, 'nthread', nthread);
	[xhat_trit2, ~, nrmsd_trit2, costOrig_trit2, time_trit2] = ADMM_tridiag(y_noise, ...
		F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask, 'fancy_mu34', fancy, 'nthread', nthread);
	[xhat_trit2_80, ~, nrmsd_trit2_80, costOrig_trit2_80, time_trit2_80] = ADMM_tridiag(y_noise, ...
		F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask, 'fancy_mu34', fancy, 'nthread', int32(80));
	[xhat_trit2_40, ~, nrmsd_trit2_40, costOrig_trit2_40, time_trit2_40] = ADMM_tridiag(y_noise, ...
		F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask, 'fancy_mu34', fancy, 'nthread', int32(40));

	[xhat_alp2t, ~, nrmsd_alp2t, costOrig_alp2t, time_alp2t] = AL_P2_gen(y_noise, ...
		F, S, R, xinit, niters, beta, x_tri_inf,'inner_iter', 1, 'mask', mask);

	if ii == 1 
		save(sprintf('./reviv/curr/mu_test_MFISTAinf_%dx%d_%diter_%s%s.mat', Nx, Ny, niters, slice_str, save_suffix));
	else
		save(sprintf('./reviv/curr/mu_test_triinf_%dx%d_%diter_%s%s.mat', Nx, Ny, niters, slice_str, save_suffix));
	end
end
send_mai_text('done with mu exp timing');

plot_timing;

display('DONE');
