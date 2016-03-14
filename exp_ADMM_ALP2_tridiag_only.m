% experimental data for tridiag_ADMM

exp_setup;
niters = 1000;

%for ii = 1:1%3%5
%        [xhat_alp2t(:,:,ii), ~, nrmsd_alp2t(ii,:), costOrig_alp2t(ii,:), time_alp2t(ii,:)] = AL_P2_gen(y_noise, F, S, R, xinit, niters, beta, xinf,'inner_iter', ii, 'mask', mask);
%	if (length(time_alp2t) ~= niters + 1), keyboard; end
%end
if ~isvar('xhat_alp2')
	[xhat_alp2, ~, nrmsd_alp2, costOrig_alp2, time_alp2] = AL_P2_gen(y_noise, F, S, R, xinit, niters, beta, xinf,'inner_iter', 1);
end
save(sprintf('%s/%s_timing_%dx%d_%diter_%s_%strue.mat', curr_folder, machine(1:3), Nx, Ny, niters, slice_str, true_opt));

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
save(sprintf('%s/%s_timing_%dx%d_%diter_%s_%strue.mat', curr_folder, machine(1:3), Nx, Ny, niters, slice_str, true_opt));

%if ~isvar('xhat_ALsm')
%	[xhat_ALsm, ~, nrmsd_ALsm, costOrig_ALsm, time_ALsm] = AL_tridiag(y_noise, F, S, CH, CV, beta, xinit, xinf, niters, 'mu_args', [mu_args {'fancy_mu34', false}]);
%end
save(sprintf('%s/%s_timing_%dx%d_%diter_%s_%strue_tridiag_only.mat', curr_folder, machine(1:3), Nx, Ny, niters, slice_str, true_opt));
send_mai_text(sprintf('done with %s timing', machine(1:3)));



display('DONE');