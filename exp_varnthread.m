% timing for varied number of threads in ADMM tridiag
exp_setup;
niters = 300;

% load xinfs
load(sprintf('%s/x_tri_inf_slice%d_beta%.*d.mat', curr_folder, slice, 3, beta), 'x_tri_inf');
if truncate
        x_tri_inf = x_tri_inf(3:end-2, 3:end-2);
end

nthread_vals = int32([1 2 4 10 20 40]);% 80 160]);
for ii = 1:length(nthread_vals);
	nthread = nthread_vals(ii);
	[xhat_tri(:,:,ii), ~, nrmsd_tri(:,ii), costOrig_tri(:,ii), time_tri(:,ii)] = ADMM_tridiag(y_noise, F, S, CH, CV, alph, beta, xinit, x_tri_inf, niters, 'mask', mask, 'nthread', nthread, 'timing', 'tridiag');
end

save(sprintf('%s/timing_tridiag_only_%dx%d_%diter_varthread.mat', curr_folder, Nx, Ny, niters));

send_mai_text(sprintf('done with %s timing', machine(1:3));

display('DONE');
