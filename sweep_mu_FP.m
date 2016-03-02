% experimental data for ADMM_tridiag
issim = 0;
exp_setup;

niters = 100;

xtrue = load_x_inf(slice, beta, curr_folder);

if ~isvar('save_suffix')
	save_suffix = '';
end
fudges = 10.^(-2:0.5:4);
ktris = 1:10;
clear nrmsd*
clear time*
for ii = 1:length(fudges)
	for jj = 1:length(ktris)
		fudge = fudges(ii);
		ktri = ktris(jj);
		[x_FP(:,:,ii,jj), ~, nrmsd_FP(ii,jj,:), costOrig_FP(ii,jj,:), time_FP(ii,jj,:)] = ADMM_FP_tridiag(y, F, S, CH, CV, beta, xinit, xtrue, niters, 'mu_args', {'mu0_fudge', fudge, 'ktri', ktri});
	end
end

save(sprintf('%s/%s_mu_sweep_FP_MFISTAinf_%dx%d_%diter_%s%s.mat', curr_folder,...
	machine(1:3), Nx, Ny, niters, slice_str, save_suffix));
send_mai_text('done with fancy mu exp timing');


display('DONE');



