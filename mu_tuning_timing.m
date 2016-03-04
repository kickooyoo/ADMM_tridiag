% experimental data for ADMM_tridiag


exp_setup;
niters = 700;

load(sprintf('%s/x_MFISTA_inf_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta), 'x*MFIS*');
if isvar('xMFIS')
	x_MFISTA = xMFIS;
end
if ~isempty(strfind(machine, 'iv1'))
	save_suffix = '_iv1';
end

if ~isvar('save_suffix')
	save_suffix = '';
end
fancy = 100; % fudge factor in get_mu
nthread = int32(160);
[xhat_trit2m, ~, nrmsd_trit2m, costOrig_trit2m, time_trit2m] = ADMM_tridiag(y_noise, ...
	F, S, CH, CV, beta, xinit, x_true, niters, 'mask', mask, 'fancy_mu34', 1, 'fancy_mask', fancy, 'nthread', nthread);
[xhat_trit2, ~, nrmsd_trit2, costOrig_trit2, time_trit2] = ADMM_tridiag(y_noise, ...
	F, S, CH, CV, beta, xinit, x_true, niters, 'mask', mask, 'fancy_mu34', fancy, 'nthread', nthread);
[xhat_trit2_jf, ~, nrmsd_trit2_jf, costOrig_trit2_jf, time_trit2_jf] = ADMM_tridiag(y_noise, ...
	F, S, CH, CV, beta, xinit, x_true, niters, 'mask', mask, 'fancy_mu34', fancy);
[xhat_trit2_80, ~, nrmsd_trit2_80, costOrig_trit2_80, time_trit2_80] = ADMM_tridiag(y_noise, ...
	F, S, CH, CV, beta, xinit, x_true, niters, 'mask', mask, 'fancy_mu34', fancy, 'nthread', int32(80));
[xhat_trit2_40, ~, nrmsd_trit2_40, costOrig_trit2_40, time_trit2_40] = ADMM_tridiag(y_noise, ...
	F, S, CH, CV, beta, xinit, x_true, niters, 'mask', mask, 'fancy_mu34', fancy, 'nthread', int32(40));

[xhat_alp2t, ~, nrmsd_alp2t, costOrig_alp2t, time_alp2t] = AL_P2_gen(y_noise, ...
	F, S, R, xinit, niters, beta, x_true,'inner_iter', 1, 'mask', mask);

save(sprintf('%s/mu_test_%sinf_%dx%d_%diter_%s%s.mat', curr_folder, inf_str, Nx, Ny, niters, slice_str, save_suffix));
send_mai_text('done with mu exp timing');

plot_timing;

display('DONE');
