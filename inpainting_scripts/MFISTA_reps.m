% test tridiag inpaint
niters = 5000;


MFISTA_inf_fname = sprintf('%s/x_MFISTA_inf_reps_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta);
for aa = 1:25
	if aa > 1
		x_prev = x_MFIS(:,:,aa-1);
		[x_MFIS(:,:,aa), C_MFIS(:,aa), T_MFIS(:,aa), err_MFIS(:,aa), ~] = MFISTA_inpainting_wrapper(Nx, Ny, RW, y, x_prev, D, beta, niters, curr_folder, slice_str, 'xinf', x_prev, 'xinfnorm', norm(col(x_prev),2));
	else
		[x_MFIS(:,:,aa), C_MFIS(:,aa), T_MFIS(:,aa), ~, ~] = MFISTA_inpainting_wrapper(Nx, Ny, RW, y, xinit, D, beta, niters, curr_folder, slice_str);
	end
	save(MFISTA_inf_fname);
end

send_mai_text(sprintf('done with reps of MFISTA on %s', machine(1:3)))




