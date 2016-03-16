% experimental data for tridiag_ADMM

exp_setup;
niters = 1000;


nthread_vals = int32([1 2 4 8 20 40 80 160]);
max_thread = jf('ncore');
nthread_vals = nthread_vals(nthread_vals <= 2*max_thread);
if ~isvar('xhat_tri')
	for ii = 1:length(nthread_vals)
		nthread = nthread_vals(ii);
		[xhat_tri(:,:,ii), ~, nrmsd_tri(ii,:), costOrig_tri(ii,:), time_tri(ii,:)] = ADMM_tridiag(y_noise, F, S, CH, CV, beta, xinit, xinf, niters, 'mu_args', mu_args, 'nthread', nthread, 'timing', 'tridiag');
	end
end
if ~isvar('xhat_AL')
	for ii = 1:length(nthread_vals)
		nthread = nthread_vals(ii);
		[xhat_AL(:,:,ii), ~, nrmsd_AL(ii,:), costOrig_AL(ii,:), time_AL(ii,:)] = AL_tridiag(y_noise, F, S, CH, CV, beta, xinit, xinf, niters, 'mu_args', mu_args, 'nthread', nthread, 'timing', 'tridiag');
	end
end

save(sprintf('%s/%s_timing_%dx%d_%diter_%s_%strue_tridiag_only.mat', curr_folder, machine(1:3), Nx, Ny, niters, slice_str, true_opt));
send_mai_text(sprintf('done with %s timing', machine(1:3)));



display('DONE');
