% experimental data for tridiag_ADMM

tridiag_exp_setup;
niters = 25000;

load(sprintf('./reviv/x_tri_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x_tri_inf');
if truncate
        x_tri_inf = reshape(x_tri_inf, 256, 144);
        x_tri_inf = x_tri_inf(3:end-2, 3:end-2);
else
        x_tri_inf = x_tri_inf;
end

% cost function for MFISTA lacks 1/2
load(sprintf('./reviv/x_MFISTA_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x*MFIS*');
	if isvar('xMFIS')
		x_MFISTA = xMFIS;
	end
if truncate
        xMFIS = reshape(x_MFISTA, 256, 144);
        xMFIS = xMFIS(3:end-2, 3:end-2);
else
        xMFIS = x_MFISTA;
end

if 0
load(sprintf('./reviv/x_alp2c_inf_slice%d_beta%.*d.mat', slice, 3, beta),'x_alp2c_inf');
if truncate
        x_alp2c_inf = reshape(x_alp2c_inf, 256, 144);
        x_alp2c_inf = x_alp2c_inf(3:end-2, 3:end-2);
else
        x_alp2c_inf = x_alp2c_inf;
end
end
load(sprintf('./reviv/x_alp2_inf_slice%d_beta%.*d.mat', slice, 3, beta),'x_alp2_inf');
if truncate
        x_alp2_inf = reshape(x_alp2_inf, 256, 144);
        x_alp2_inf = x_alp2_inf(3:end-2, 3:end-2);
else
        x_alp2_inf = x_alp2_inf;
end

[xhat_alp2, ~, nrmsd_alp2, costOrig_alp2, time_alp2] = AL_P2_gen(y_noise, F, S, R, x_tri_inf, niters, beta, x_alp2_inf,'inner_iter',3, 'mask', mask);
[xhat_alp2_selfinit, ~, nrmsd_alp2_selfinit, costOrig_alp2_selfinit, time_alp2_selfinit] = AL_P2_gen(y_noise, F, S, R, x_alp2_inf, niters, beta, x_alp2_inf,'inner_iter',3, 'mask', mask);
[xhat_tri, ~, nrmsd_tri, costOrig_tri, time_tri] = tridiag_ADMM(y_noise, F, S, CH, CV, alph, beta, xMFIS, x_tri_inf, niters, 'mask', mask, 'mu', ones(1,5));
[xhat_tri_selfinit, ~, nrmsd_tri_self_init, costOrig_tri_selfinit, time_tri_selfinit] = tridiag_ADMM(y_noise, F, S, CH, CV, alph, beta, x_tri_inf, x_tri_inf, niters, 'mask', mask, 'mu', ones(1,5));
[xhat_tri_alp2init, ~, nrmsd_tri_alp2init, costOrig_tri_alp2init, time_tri_alp2init] = tridiag_ADMM(y_noise, F, S, CH, CV, alph, beta, x_alp2_inf, x_tri_inf, niters, 'mask', mask, 'mu', ones(1,5));
[xhat_MFIS, CMFIS, TFIS, l2DFIS, RMSEFIS] = MFISTA_wrapper(Nx, Ny, R, y_noise, x_tri_inf, F, S, beta, niters);

save(sprintf('./reviv/check_xinfs_%dx%d_%diter_%dslice.mat', Nx, Ny, niters, slice));
send_mai_text('done with mpel8 timing');

display('DONE');
