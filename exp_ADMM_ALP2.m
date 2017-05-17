% experimental data for tridiag_ADMM

exp_setup;
if ~isvar('niters')
	niters = 5000;
end

%for ii = 1:1%3%5
%        [xhat_alp2t(:,:,ii), ~, nrmsd_alp2t(ii,:), costOrig_alp2t(ii,:), time_alp2t(ii,:)] = AL_P2_gen(y_noise, F, S, R, xinit, niters, beta, xinf,'inner_iter', ii, 'mask', mask);
%	if (length(time_alp2t) ~= niters + 1), keyboard; end
%end
if ~isvar('xhat_alp2')
	[xhat_alp2, ~, nrmsd_alp2, costOrig_alp2, time_alp2] = AL_P2_gen(y_noise, F, S, R, xinit, niters, beta, xinf,'inner_iter', 1);
end
save(sprintf('%s/%s_timing_%dx%d_%diter_%s_%strue_%dbeta.mat', curr_folder, machine(1:3), Nx, Ny, niters, slice_str, true_opt, beta));

nthread_vals = int32([1 2 4 8 20 40 80 160]);
max_thread = jf('ncore');
nthread_vals = nthread_vals(nthread_vals <= 2*max_thread);
nthread_vals = 2;
if ~isvar('x_circ')
	[x_circ, ~, err_circ, costOrig_circ, time_circ] = AL_P2_gen(y_noise, F, S, Rcirc, xinit, niters, beta, xinf, 'inner_iter', 1);
end
save(sprintf('%s/%s_timing_%dx%d_%diter_%s_%strue_%dbeta.mat', curr_folder, machine(1:3), Nx, Ny, niters, slice_str, true_opt, beta));
if ~isvar('xhat_tri')
	for ii = 1:length(nthread_vals)
		nthread = nthread_vals(ii);
		[xhat_tri(:,:,ii), ~, nrmsd_tri(ii,:), costOrig_tri(ii,:), time_tri(ii,:)] = ADMM_tridiag(y_noise, F, S, CH, CV, beta, xinit, xinf, niters, 'mu_args', mu_args, 'nthread', nthread);
	end
end
save(sprintf('%s/%s_timing_%dx%d_%diter_%s_%strue_%dbeta.mat', curr_folder, machine(1:3), Nx, Ny, niters, slice_str, true_opt, beta));
if ~isvar('xhat_AL')
	for ii = 1:length(nthread_vals)
		nthread = nthread_vals(ii);
		[xhat_AL(:,:,ii), ~, nrmsd_AL(ii,:), costOrig_AL(ii,:), time_AL(ii,:)] = AL_tridiag(y_noise, F, S, CH, CV, beta, xinit, xinf, niters, 'mu_args', mu_args, 'nthread', nthread');
	end
end
save(sprintf('%s/%s_timing_%dx%d_%diter_%s_%strue_%dbeta.mat', curr_folder, machine(1:3), Nx, Ny, niters, slice_str, true_opt, beta));

if ~isvar('xhat_tri_FP')
	[xhat_tri_FP, ~, nrmsd_tri_FP, costOrig_tri_FP, time_tri_FP] = ADMM_FP_tridiag(y_noise, F, S, CH, CV, beta, xinit, xinf, niters, 'mu_args', mu_args);
end
if ~isvar('nrmsd_MFIS')
	xinf_norm = norm(col(xinf), 2);
	[xMFIS, costOrig_MFIS, time_MFIS, nrmsd_MFIS, ~] = MFISTA_wrapper(Nx, Ny, R, y_noise, xinit, F, S, beta, round(niters*1.5), curr_folder, 'xinf', xinf, 'xinfnorm', xinf_norm);
end

if ~isvar('x_circ')
	[x_circ, ~, err_circ, costOrig_circ, time_circ] = AL_P2_inpainting(y, D, RcircW, xinit, niters, beta, xtrue, 'alphw', alphw, 'betaw', betaw, 'zmethod', 'FFT');
end

save(sprintf('%s/%s_timing_%dx%d_%diter_%s_%strue_%dbeta.mat', curr_folder, machine(1:3), Nx, Ny, niters, slice_str, true_opt, beta));
send_mai_text(sprintf('done with %s timing', machine(1:3)));



display('DONE');
