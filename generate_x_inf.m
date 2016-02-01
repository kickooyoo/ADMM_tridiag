% generate xinf for tridiag experiment

tridiag_exp_setup;

niters = 5000;


[xhat_tri, ~, ~, costOrig_tri, time_tri] = tridiag_ADMM(y, F, S, CH, CV, alph, beta, SoS, zeros(size(SoS)), niters, 'mu', mu);
x_tri_inf = xhat_tri;
save(sprintf('./reviv/x_tri_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x_tri_inf');

[xhat_alp2, ~, ~, costOrig_alp2, time_alp2] = AL_P2_gen(y, F, S, Rcirc, SoS, niters, beta, zeros(size(SoS)), 'mu', mu, 'zmethod', 'FFT');
x_alp2c_inf = xhat_alp2;
save(sprintf('./reviv/x_alp2c_inf_slice%d_beta%.*d.mat', slice, 3, beta),'x_alp2c_inf');
send_mai_text('done with xinf, now time test');

return;

[xhat_triw, ~, ~, costOrig_triw, time_triw] = tridiag_ADMM_W(y, F, S, CH, CV, alph, beta, SoS, zeros(size(SoS)), niters, 'mu', mu);
x_triw_inf = xhat_triw;
send_mai_text('check xhat');
keyboard;
save(sprintf('./reviv/x_triw_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x_triw_inf');


[xhat_alp2w, ~, ~, costOrig_alp2w, time_alp2w] = AL_P2_gen(y, F, S, RcircW, SoS, niters, beta, zeros(size(SoS)), 'mu', mu, 'zmethod', 'FFT');
x_alp2cw_inf = xhat_alp2w;
send_mai_text('check xhat');
keyboard;
save(sprintf('./reviv/x_alp2cw_inf_slice%d_beta%.*d.mat', slice, 3, beta),'x_alp2cw_inf');
