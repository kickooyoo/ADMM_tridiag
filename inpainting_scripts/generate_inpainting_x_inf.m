% generate xinf for tridiag experiment
niters_inf = 50000;

if ~alphw == 1
	MFISTA_inf_fname = sprintf('%s/x_MFISTA_inf_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta);
else
	MFISTA_inf_fname = sprintf('%s/x_MFISTA_inf_%s_beta%.*d_%1.1dalphw.mat', curr_folder, slice_str, 3, beta, alphw);
end
[xMFIS, CMFIS, TFIS, ~, ~] = MFISTA_inpainting_wrapper(Nx, Ny, RW, y, xinit, D, beta, niters_inf, curr_folder, slice_str);
save(MFISTA_inf_fname, 'xMFIS', 'niters_inf', 'CMFIS', 'TFIS');

send_mai_text(sprintf('done with xinf on %s, now run timing tests', machine(1:3)))


