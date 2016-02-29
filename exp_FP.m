issim = 0;
exp_setup;
xtrue = load_x_inf(slice, beta);
niters = 100;

[x_smpd, ~, nrmsd_smpd, ~, time_smpd] = ADMM_FP_tridiag(y, F, S, CH, CV, beta, xinit, ...
        xtrue, niters, 'attempt_par', true, 'pmethod', 'spmd');
[x_pfor, ~, nrmsd_pfor, ~, time_pfor] = ADMM_FP_tridiag(y, F, S, CH, CV, beta, xinit, ...
        xtrue, niters, 'attempt_par', true, 'pmethod', 'pfor');
[x_fevl, ~, nrmsd_fevl, ~, time_fevl] = ADMM_FP_tridiag(y, F, S, CH, CV, beta, xinit, ...
        xtrue, niters, 'attempt_par', true, 'pmethod', 'feval');
[x, ~, nrmsd, ~, time] = ADMM_FP_tridiag(y, F, S, CH, CV, beta, xinit, ...
        xtrue, niters);
[x_warm, ~, nrmsd_warm, ~, time_warm] = ADMM_FP_tridiag(y, F, S, CH, CV, beta, xinit, ...
        xtrue, niters, 'warmup', 10);

save(sprintf('./reviv/curr/%s_PF_timing_%dx%d_%diter_%s.mat', machine(1:3), Nx, Ny, niters, slice_str));
