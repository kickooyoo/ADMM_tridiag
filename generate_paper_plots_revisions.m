% script for making plots and images on vega
db_path = '~/Dropbox/fessler/writing/tridiag/';
FontSize = 12;
FP_incl = 0;
do_all = [0 0 0 0 1 0 0];
preload = true;
%% sagittal slice 69 --------------------------------------------------------
clearvars -except db_path FontSize FP_incl
eval(load_except([db_path 'sagittal_slice69_r6/ir7_timing_256x128_5000iter_sagittal69_avgtrue_8388608beta.mat'], 'db_path'))
clear time*inf
orn_ndx = 800;
plot_axes = {[0 200 -88 0];[0 3000 -88 0]};
if FP_incl
	short_slice_str = 'ax_38';
	order = [5 2 3 4 1];
else
	clear *FP
	short_slice_str = 'ax_38_noFP';
	colors = 'mgkbry'; % for w/o FP
	markers = 'o+*sd.^v><ph'; % for w/o FP
	order = [2 3 4 1];
end
eval(load_except([db_path 'axial_slice38_truncate/ir7_timing_256x128_5000iter_slice38_truetrue.mat'], 'db_path'))
clear time*inf
orn_ndx = 50;
plot_axes = {[0 8 -20 0];[0 150 -20 0]};
if FP_incl
	short_slice_str = 'ax_38_err';
	order = [5 2 3 4 1];
else
	clear *FP
	short_slice_str = 'ax_38_err_noFP';
	colors = 'mgkbry'; % for w/o FP
	markers = 'o+*sd.^v><ph'; % for w/o FP
	order = [2 3 4 1];
end
ystr = 'NRMSE to x_{SENSE} (dB)';
lstring = {'AL-tridiag';'MFISTA-5';'AL-P2-NC';'ADMM-tridiag';'ADMM-FP-tridiag'};
plot_timing

%% timing plot for parallelization methods
clearvars -except db_path
eval(load_except([db_path 'curr/ir7_PF_timing_256x144_100iter_slice38.mat'], 'db_path'))
load([db_path 'curr/ir7_x_MFISTA_timing_slice38_beta1048576.mat'])
clear time_comp; clear *_tri*; clear *alp2*; clear *warm; plot_timing; close;
lh = legend('ADMM-tridiag,nopar', 'ADMM-tridiag,feval', 'ADMM-tridiag,pfor', 'ADMM-tridiag,spmd'); set(lh, 'FontSize', FontSize);
title('');
save_im(db_path, 'ax_38_timing_FP_parmethods_ir72')
end
%% sim sagittal  --------------------------------------------------------
clearvars -except db_path FontSize order
eval(load_except([db_path 'sim/sim_r6/ir7_timing_240x200_5000iter_sim_avgtrue.mat'], 'db_path'))
if 1
	orn_ndx = 250;
	plot_axes = {[0 200 -90 0];[0 1500 -90 0]};
	short_slice_str = 'sim';
else
	orn_ndx = 90;
	plot_axes = {[0 9 -55 0];[0 150 -55 0]};
	short_slice_str = 'sim_zoom';
end
eval(load_except([db_path 'sim/ir7_timing_240x200_500iter_sim_truetrue.mat'], 'db_path'))
orn_ndx = 55;
plot_axes = {[0 8 -22 0];[0 125 -22 0]};
short_slice_str = 'sim_err';
ystr = 'NRMSE to x_{true} (dB)';
lstring = {'AL-tridiag'; 'MFISTA-5'%'AL-tridiag,svt';
        'AL-P2-NC';'ADMM-tridiag';'ADMM-FP-tridiag'};
plot_timing

clearvars -except db_path FontSize
eval(load_except([db_path 'sim/ir7_timing_240x200_500iter_sim_truetrue.mat'], 'db_path'))
orn_ndx = 15;
plot_axes = {[0 9 -22 0];[0 150 -22 0]};
lstring = {'AL-tridiag'; 'MFISTA-5'%'AL-tridiag,svt';
        'AL-P2';'ADMM-tridiag';'ADMM-FP-tridiag'};
plot_timing
short_slice_str = 'sim_err';
ylabel('NRMSE to x_{true} (dB)')


%% --------------------------- save timing ---------------------------
title('');
save_im(db_path, sprintf('figs/%s_iter_%s', short_slice_str, machine(1:3)))
close;
title('');
save_im(db_path, sprintf('figs/%s_timing_%s', short_slice_str, machine(1:3)))

%% --------------------------- save images ---------------------------

pad_samp = zeros(Nx + 1, Ny + 1);
pad_samp(1 : end-1, 1 : end-1) = fftshift(samp);
figure; im(1-pad_samp)
axis off; title('');
save_im(db_path, sprintf('figs/%s_samp', short_slice_str))

if ~isvar('xinf') && isvar('xtrue')
        brain_max = max(abs([xhat_AL(:); xtrue(:); xtrue(:); xhat_alp2(:)]));
else
        brain_max = max(abs([xhat_AL(:); xinit(:); xinf(:); xhat_alp2(:)]));
end
figure; im_brain(xhat_AL, 'clim', [0 brain_max], 'orient', orient)
axis off; title(''); colorbar;
save_im(db_path, sprintf('figs/%s_x_AL', short_slice_str));

figure; im_brain(xhat_tri, 'clim', [0 brain_max], 'orient', orient)
axis off; title(''); colorbar;
save_im(db_path, sprintf('figs/%s_x_tri', short_slice_str));

if isvar('body_coil') && isvar('xinf')
        [bc_scaled, bc_scale] = ir_wls_init_scale(1, xinf, body_coil);
        figure; im_brain(bc_scaled, 'clim', [0 brain_max], 'orient', orient)
        axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
        save_im(db_path, sprintf('figs/%s_body_coil', short_slice_str));
end

figure; im_brain(xinit, 'clim', [0 brain_max], 'orient', orient)
axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
save_im(db_path, sprintf('figs/%s_SoS_init', short_slice_str));

if isvar('xinf')
        xinf_tmp = xinf;
else
        xinf_tmp = xMFIS;
end
figure; im_brain(xMFIS, 'clim', [0 brain_max], 'orient', orient)
axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
save_im(db_path, sprintf('figs/%s_x_MFISTA', short_slice_str));

figure; im_brain(xhat_AL - xMFIS, 'orient', orient)
axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
save_im(db_path, sprintf('figs/%s_diff_AL_MFISTA', short_slice_str));

figure; im_brain(xhat_tri - xMFIS, 'orient', orient)
axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
save_im(db_path, sprintf('figs/%s_diff_tri_MFISTA', short_slice_str));

figure; im_brain(body_coil - xMFIS, 'orient', orient)
axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
save_im(db_path, sprintf('figs/%s_diff_BC_MFISTA', short_slice_str));

figure; im_brain(body_coil - xhat_tri, 'orient', orient)
axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
save_im(db_path, sprintf('figs/%s_diff_BC_tri', short_slice_str));

if strcmp(true_opt, 'true')
        
      figure; im_brain(xinf, 'clim', [0 brain_max], 'orient', orient)
      axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
      save_im(db_path, sprintf('figs/%s_x_SENSE', short_slice_str));
      
      figure; im_brain(xinf - xhat_tri, 'orient', orient)
      axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
      save_im(db_path, sprintf('figs/%s_diff_SENSE_tri', short_slice_str));
end
%%
if isvar('xtrue')
        figure; im_brain(xtrue, 'clim', [0 brain_max], 'orient', orient)
        axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
        save_im(db_path, sprintf('figs/%s_x_true', short_slice_str));
        
        figure; im_brain(xhat_AL - xtrue, 'orient', orient)
        axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
        save_im(db_path, sprintf('figs/%s_diff_AL_true', short_slice_str));
end