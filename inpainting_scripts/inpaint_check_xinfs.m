% checking if MFISTA, ADMM-tridiag, and AL-P2 converge to the same solution

inpainting_setup;
D = Ginpaint(samp);
load(ALP2NC_inf_fname);

niters = 5000;

if isvar('xMFIS')
	x_MFISTA = xMFIS;
end

%load(sprintf('%s/x_alp2c_inf_slice%d_beta%.*d.mat', curr_folder, slice, 3, beta),'x_alp2c_inf');
%load(sprintf('%s/x_alp2_inf_slice%d_beta%.*d.mat', curr_folder, slice, 3, beta),'x_alp2_inf');

[xhat_alp2, ~, nrmsd_alp2, costOrig_alp2, time_alp2] = AL_P2_inpainting(y, D, RW, x_ADMM_inf, niters, beta, x_ALP2_inf,  'betaw', betaw, 'alphw', alphw, 'kapp', 7, 'inner_iter',3);
[xhat_alp2_MFISinit, ~, nrmsd_alp2_MFISinit, costOrig_alp2_MFISinit, time_alp2_MFISinit] = AL_P2_inpainting(y, D, RW, xMFIS, niters, beta, x_ALP2_inf,  'betaw', betaw, 'alphw', alphw, 'kapp', 7, 'inner_iter',3);
[xhat_alp2_selfinit, ~, nrmsd_alp2_selfinit, costOrig_alp2_selfinit, ti0ime_alp2_selfinit] = AL_P2_inpainting(y, D, RW, x_ALP2_inf, niters, beta, x_ALP2_inf,  'betaw', betaw, 'alphw', alphw, 'kapp', 7, 'inner_iter',3);
%[xhat_tri, ~, nrmsd_tri, costOrig_tri, time_tri] = ADMM_tridiag(y, F, S, CH, CV, beta, xMFIS, x_ADMM_inf, niters);
%[xhat_ADMM_selfinit, ~, nrmsd_ADMM_self_init, costOrig_ADMM_selfinit, time_ADMM_selfinit] = ADMM_tridiag_inpaint(y, D, CHW, CVW, beta, x_ADMM_inf, x_ADMM_inf, niters, 'betaw', betaw, 'alphw', alphw, 'alph', alph, 'kapp', 7);
[xhat_ADMM_alp2init, ~, nrmsd_ADMM_alp2init, costOrig_ADMM_alp2init, time_ADMM_alp2init] = ADMM_tridiag_inpaint(y, D, CHW, CVW, beta, x_ALP2_inf, x_ADMM_inf, niters, 'betaw', betaw, 'alphw', alphw, 'alph', alph, 'kapp', 7);
[xhat_MFIS_alp2init, cost_MFIS_alp2init, time_MFIS_alp2init, ~, rmsd_MFIS_alp2init] = MFISTA_inpainting_wrapper(Nx, Ny, RW, y, x_ALP2_inf, D, beta, niters, curr_folder, slice_str);



[xhat_sb, ~, nrmsd_sb, costOrig_sb, time_sb] = SB_inpainting(y, D, RW, x_ADMM_inf, niters, beta, x_ALP2_inf,  'betaw', betaw, 'alphw', alphw, 'kapp', 7, 'inner_iter',3);

send_mai_text('done with mpel8 timing');

display('DONE');
