% experimental data for ADMM_tridiag
issim = 0;
exp_setup;

niters = 100;

xtrue = load_x_inf(slice, beta);

if ~isvar('save_suffix')
	save_suffix = '';
end
fudges = 10.^(-2:0.5:4);
clear nrmsd*
clear time*
for ii = 1:length(ktris)
        fudge = fudges(ii);
	[x, xsaved, err, cost, time] = ADMM_FP_tridiag(y, F, S, CH, CV, beta, xinit, xtrue, niters, 'mu_args', {'mu0_fudge', fudge});
end

save(sprintf('./reviv/curr/mu_sweep_FP_MFISTAinf_%dx%d_%diter_%s%s.mat', ...
	Nx, Ny, niters, slice_str, save_suffix));
send_mai_text('done with fancy mu exp timing');

plot_timing;

display('DONE');



