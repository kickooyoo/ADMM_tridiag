% generate xinf for tridiag experiment

% for slice 67, beta = 2^20.1
tridiag_exp_setup;

niters = 5000;


[xhat_tri, ~, ~, costOrig_tri, time_tri] = tridiag_ADMM(y, F, S, CH, CV, alph, beta, SoS, 0, niters, 'mu', mu);
x_tri_inf = xhat_tri;
save(sprintf('x_tri_inf_slice%d,beta%.*d.mat', slice, 3, beta), 'x_tri_inf');

[xhat_triw, ~, ~, costOrig_triw, time_triw] = tridiag_ADMM_W(y, F, S, CH, CV, alph, beta, SoS, 0, niters, 'mu', mu);
x_triw_inf = xhat_triw;
save(sprintf('x_triw_inf_slice%d,beta%.*d.mat', slice, 3, beta), 'x_triw_inf');

[xhat_alp2, ~, ~, costOrig_alp2, time_alp2] = AL_P2_gen(y, F, S, Rcirc, SoS, niters, beta, 0, 'mu', mu, 'zmethod', 'FFT');
x_alp2c_inf = xhat_alp2;
save(sprintf('x_alp2c_inf_slice%d,beta%.*d.mat', slice, 3, beta),'x_alp2c_inf');

[xhat_alp2w, ~, ~, costOrig_alp2w, time_alp2w] = AL_P2_gen(y, F, S, RcircW, SoS, niters, beta, 0, 'mu', mu, 'zmethod', 'FFT');
x_alp2cw_inf = xhat_alp2w;
save(sprintf('x_alp2cw_inf_slice%d,beta%.*d.mat', slice, 3, beta),'x_alp2cw_inf');
