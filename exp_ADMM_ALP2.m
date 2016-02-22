% experimental data for tridiag_ADMM

exp_setup;
niters = 350;

x_tri_inf = load_x_inf(slice, beta);

if wavelets 
[xhat_alp2c, ~, nrmsd_alp2c, costOrig_alp2c, time_alp2c] = AL_P2_gen(...
        y_noise, F, S, RcircW, SoS, niters, beta, x_alp2c_inf, 'zmethod', 'fft');

[xhat_alp2t, ~, nrmsd_alp2t, costOrig_alp2t,time_alp2t] = AL_P2_gen(...
        y_noise, F, S, RW, SoS, niters, beta, x_tri_inf, 'inner_iter', 3);
end

for ii = 1:1%3%5
        [xhat_alp2t(:,:,ii), ~, nrmsd_alp2t(:,ii), costOrig_alp2t(:,ii), time_alp2t(:,ii)] = AL_P2_gen(y_noise, F, S, R, xinit, niters, beta, x_tri_inf,'inner_iter', ii, 'mask', mask);
	if (length(time_alp2t) ~= niters + 1), keyboard; end
end

%[xhat_alp2c, ~, nrmsd_alp2c, costOrig_alp2c, time_alp2c] = AL_P2_gen(y_noise, F, S, Rcirc, xinit, niters, beta, x_alp2c_inf, 'zmethod','fft', 'mask', mask);
[xhat_tri, ~, nrmsd_tri, costOrig_tri, time_tri] = ADMM_tridiag(y_noise, F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask, 'mu', ones(1,5));
%[xhat_tri_max, ~, nrmsd_tri_max, costOrig_tri_max, time_tri_max] = ADMM_tridiag(y_noise, F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask, 'nthread', int32(maxNumCompThreads('automatic')),'mu', ones(1,5));
%[xhat_tri_2max, ~, nrmsd_tri_2max, costOrig_tri_2max, time_tri_2max] =  ADMM_tridiag(y_noise, F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask, 'nthread', int32(160),'mu', ones(1,5));
%[xhat_tri_mu, ~, nrmsd_tri_mu, costOrig_tri_mu, time_tri_mu] = ADMM_tridiag(y_noise, F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask);
%[xhat_tri_max_mu, ~, nrmsd_tri_max_mu, costOrig_tri_max_mu, time_tri_max_mu] = ADMM_tridiag(y_noise, F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask, 'nthread', int32(maxNumCompThreads('automatic')));

save(sprintf('./reviv/curr/mpel8_timing_%dx%d_%diter_%s_ADMMtrue.mat', Nx, Ny, niters, slice_str));
send_mai_text('done with mpel8 timing');

display('DONE');
