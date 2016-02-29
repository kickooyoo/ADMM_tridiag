% experimental data for tridiag_ADMM

exp_setup;
niters = 500;

xtrue = load_x_inf(slice, beta);

%for ii = 1:1%3%5
%        [xhat_alp2t(:,:,ii), ~, nrmsd_alp2t(ii,:), costOrig_alp2t(ii,:), time_alp2t(ii,:)] = AL_P2_gen(y_noise, F, S, R, xinit, niters, beta, xtrue,'inner_iter', ii, 'mask', mask);
%	if (length(time_alp2t) ~= niters + 1), keyboard; end
%end
[xhat_alp2, ~, nrmsd_alp2, costOrig_alp2, time_alp2] = AL_P2_gen(y_noise, F, S, R, xinit, niters, beta, xtrue,'inner_iter', 1, 'mask', mask);
if issim
	mu_args = {'noise', 2};
else
	mu_args = {};
end
	
[xhat_tri, ~, nrmsd_tri, costOrig_tri, time_tri] = ADMM_tridiag(y_noise, F, S, CH, CV, beta, xinit, xtrue, niters, 'mask', mask, 'mu_args', mu_args);
[xhat_AL, ~, nrmsd_AL, costOrig_AL, time_AL] = AL_tridiag(y_noise, F, S, CH, CV, beta, xinit, xtrue, niters, 'mask', mask, 'mu_args', mu_args);
[xhat_ALsm, ~, nrmsd_ALsm, costOrig_ALsm, time_ALsm] = AL_tridiag(y_noise, F, S, CH, CV, beta, xinit, xtrue, niters, 'mask', mask, 'mu_args', [mu_args {'fancy_mu34', false}]);

%[xhat_tri_max, ~, nrmsd_tri_max, costOrig_tri_max, time_tri_max] = ADMM_tridiag(y_noise, F, S, CH, CV, beta, xinit, xtrue, niters, 'mask', mask, 'nthread', int32(maxNumCompThreads('automatic')),'mu', ones(1,5));
[xhat_tri_FP, ~, nrmsd_tri_FP, costOrig_tri_FP, time_tri_FP] = ADMM_FP_tridiag(y_noise, F, S, CH, CV, beta, xinit, xtrue, niters, 'mask', mask, 'mu_args', mu_args);

save(sprintf('./reviv/curr/%s_timing_%dx%d_%diter_%s_MFISTAtrue.mat', machine(1:3), Nx, Ny, niters, slice_str));
%send_mai_text('done with mpel8 timing');

display('DONE');
