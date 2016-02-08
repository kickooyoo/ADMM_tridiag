% experimental data for tridiag_ADMM

tridiag_exp_setup;
niters = 300;

% load xinfs
%load([home 'Documents/mai_code/ADMM_tridiag/reviv/tri_chcv_5000iter.mat'],'xhat_tri')
load(sprintf('./reviv/x_tri_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x_tri_inf');
if truncate
        x_tri_inf = reshape(x_tri_inf, 256, 144);
        x_tri_inf = x_tri_inf(3:end-2, 3:end-2);
else
        x_tri_inf = x_tri_inf;
end

nthread_vals = int32([1 2 4 10 20 40]);% 80 160]);
for ii = 1:length(nthread_vals);
	nthread = nthread_vals(ii);
	[xhat_tri(:,:,ii), ~, nrmsd_tri(:,ii), costOrig_tri(:,ii), time_tri(:,ii)] = tridiag_ADMM(y_noise, F, S, CH, CV, alph, beta, xinit, x_tri_inf, niters, 'mask', mask, 'mu', ones(1,5), 'nthread', nthread, 'timing', 'tridiag');
%	[xhat_tri_mu(:,:,ii), ~, nrmsd_tri_mu(:,ii), costOrig_tri_mu(:,ii), time_tri_mu(:,ii)] = tridiag_ADMM(y_noise, F, S, CH, CV, alph, beta, xinit, x_tri_inf, niters, 'mask', mask, 'nthread', nthread, ,'timing', 'tridiag');
end

save(sprintf('./reviv/mpel8_timing_%dx%d_%diter_varthread_ir63.mat', Nx, Ny, niters));
send_mai_text('done with mpel8 timing');

display('DONE');