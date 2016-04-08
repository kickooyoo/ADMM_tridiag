% generate xinf for tridiag experiment
exp_setup;

niters = 30000;
xinf_norm = norm(col(xinf), 2);
[xMFIS, costOrig_MFIS, time_MFIS, nrmsd_MFIS, ~] = MFISTA_wrapper(Nx, Ny, R, y_noise, xinit, F, S, beta, niters, curr_folder, 'xinf', xinf, 'xinfnorm', xinf_norm);
save(sprintf('%s/%s_x_MFISTA_timing_%s_beta%.*d.mat', curr_folder, machine(1:3), slice_str, 3, beta), 'xMFIS', 'niters', 'costOrig_MFIS', 'time_MFIS', 'nrmsd_MFIS');

send_mai_text('done with MFISTA timing')
% why is there such a big discrepancy between my NRMSD and MFISTA/Sathish's??
