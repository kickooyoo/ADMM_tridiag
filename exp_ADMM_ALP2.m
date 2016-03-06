% experimental data for tridiag_ADMM

exp_setup;
niters = 20000;

%for ii = 1:1%3%5
%        [xhat_alp2t(:,:,ii), ~, nrmsd_alp2t(ii,:), costOrig_alp2t(ii,:), time_alp2t(ii,:)] = AL_P2_gen(y_noise, F, S, R, xinit, niters, beta, xinf,'inner_iter', ii, 'mask', mask);
%	if (length(time_alp2t) ~= niters + 1), keyboard; end
%end
if ~isvar('xhat_alp2')
	[xhat_alp2, ~, nrmsd_alp2, costOrig_alp2, time_alp2] = AL_P2_gen(y_noise, F, S, R, xinit, niters, beta, xinf,'inner_iter', 1);
end
if ~isvar('xhat_tri')	
	[xhat_tri, ~, nrmsd_tri, costOrig_tri, time_tri] = ADMM_tridiag(y_noise, F, S, CH, CV, beta, xinit, xinf, niters, 'mu_args', mu_args);
end
if ~isvar('xhat_AL')
	[xhat_AL, ~, nrmsd_AL, costOrig_AL, time_AL] = AL_tridiag(y_noise, F, S, CH, CV, beta, xinit, xinf, niters, 'mu_args', mu_args);
end
%if ~isvar('xhat_ALsm')
%	[xhat_ALsm, ~, nrmsd_ALsm, costOrig_ALsm, time_ALsm] = AL_tridiag(y_noise, F, S, CH, CV, beta, xinit, xinf, niters, 'mu_args', [mu_args {'fancy_mu34', false}]);
%end
if ~isvar('xhat_tri_FP')
	[xhat_tri_FP, ~, nrmsd_tri_FP, costOrig_tri_FP, time_tri_FP] = ADMM_FP_tridiag(y_noise, F, S, CH, CV, beta, xinit, xinf, niters, 'mu_args', mu_args);
end
if ~isvar('nrmsd_MFIS')
	xinf_norm = norm(col(xinf), 2);
	[xMFIS, costOrig_MFIS, time_MFIS, nrmsd_MFIS, ~] = MFISTA_wrapper(Nx, Ny, R, y_noise, xinit, F, S, beta, round(niters*1.5), curr_folder, 'xinf', xinf, 'xinfnorm', xinf_norm);
end
save(sprintf('%s/%s_timing_%dx%d_%diter_%s_%strue.mat', curr_folder, machine(1:3), Nx, Ny, niters, slice_str, true_opt));
send_mai_text(sprintf('done with %s timing', machine(1:3)));



display('DONE');
