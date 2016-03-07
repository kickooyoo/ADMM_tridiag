% script for making plots and images on vega
db_path = '~/Dropbox/fessler/experimental_data/tridiag_vega/';
FontSize = 16;

%% axial slice 38 --------------------------------------------------------
close all; clearvars -except db_path
eval(load_except([db_path 'axial_slice38/ir7_timing_256x144_20000iter_slice38_avgtrue.mat'], 'db_path'))
plot_axes = {[0 200 -80 0];[0 2000 -80 0]};
plot_timing
short_slice_str = 'ax_38';
% lstring = {'AL-tridiag';'AL-tridiag,svt';'MFISTA-5';'AL-P2';'ADMM-tridiag';'ADMM-FP-tridiag'};
lstring = {'AL-tridiag';'MFISTA-5';'AL-P2';'ADMM-tridiag';'ADMM-FP-tridiag'};

%% axial slice 38 RAW --------------------------------------------------------
close all; clearvars -except db_path
eval(load_except([db_path 'axial_slice38_raw/ir6_timing_256x144_20000iter_slice38_avgtrue.mat'], 'db_path'))
plot_timing
short_slice_str = 'ax_38';
% lstring = {'AL-tridiag';'AL-tridiag,svt';'MFISTA-5';'AL-P2';'ADMM-tridiag';'ADMM-FP-tridiag'};
lstring = {'AL-tridiag';'MFISTA-5';'AL-P2';'ADMM-tridiag';'ADMM-FP-tridiag'};

%% timing plot for parallelization methods
clearvars -except db_path
eval(load_except([db_path 'curr/ir7_PF_timing_256x144_100iter_slice38.mat'], 'db_path'))
load([db_path 'curr/ir7_x_MFISTA_timing_slice38_beta1048576.mat'])
clear time_comp; clear *_tri*; clear *alp2*; clear *warm; plot_timing; close;
lh = legend('ADMM-tridiag,nopar', 'ADMM-tridiag,feval', 'ADMM-tridiag,pfor', 'ADMM-tridiag,spmd'); set(lh, 'FontSize', FontSize);
title('');
save_im(db_path, 'ax_38_timing_FP_parmethods_ir72')

%% sim --------------------------------------------------------
clearvars -except db_path
load([db_path 'sim/ir7_x_MFISTA_timing_sim_beta8192.mat'])
eval(load_except([db_path 'sim/ir7_timing_240x200_10000iter_sim_MFISTAtrue.mat'], 'db_path'))
plot_timing
short_slice_str = 'sim';
orient = 'sim';
% lstring = {'AL-tridiag';'AL-tridiag,svt';'MFISTA-5';'AL-P2'};

lstring = {'AL-tridiag';'AL-tridiag,svt';'AL-P2';'ADMM-tridiag';'ADMM-FP-tridiag'};
% reminder: want NRMSE!
%% coronal slice 135 --------------------------------------------------------
clearvars -except db_path
eval(load_except([db_path 'coronal_135/ir7_timing_144x128_6000iter_coronal135_MFISTAtrue.mat'], 'db_path'))
if ~isvar('xMFIS')
        load([db_path 'coronal_135/ir7_x_MFISTA_timing_coronal135_beta1048576.mat'], '*MFIS')
end
plot_timing
lstring = {'AL-tridiag';'MFISTA-5';'AL-P2';'ADMM-tridiag';'ADMM-FP-tridiag'};
short_slice_str = 'co_135';

%% sagittal slice 69 ----------------------------------------------------
clearvars -except db_path
eval(load_except([db_path 'sagittal_slice69/iv1_timing_256x128_6000iter_sagittal69_MFISTAtrue.mat'], 'db_path'))
plot_timing
short_slice_str = 'sa_69';
lstring = {'AL-tridiag'; 'AL-tridiag,svt'; 'MFISTA-5';'AL-P2'; 'ADMM-tridiag'; 'ADMM-FP-tridiag'};

%% axial slice 90 --------------------------------------------------------
close all; clearvars -except db_path
eval(load_except([db_path 'axial_slice90/iv1_timing_256x144_20000iter_slice90_avgtrue.mat'], 'db_path'))

plot_timing
short_slice_str = 'ax_90';
lstring = {'AL-tridiag';'MFISTA-5';'AL-P2';'ADMM-tridiag';'ADMM-FP-tridiag'};
%% axial slice 90 RAW--------------------------------------------------------
close all; clearvars -except db_path
eval(load_except([db_path 'axial_slice90_raw/iv1_timing_256x144_20000iter_slice90_avgtrue.mat'], 'db_path'))

plot_timing
short_slice_str = 'ax_90';
lstring = {'AL-tridiag';'MFISTA-5';'AL-P2';'ADMM-tridiag';'ADMM-FP-tridiag'};

%% axial slice 67 --------------------------------------------------------
clearvars -except db_path FontSize
eval(load_except([db_path 'axial_slice67/iv1_timing_256x144_20000iter_slice67_avgtrue.mat'], 'db_path'))

plot_timing
short_slice_str = 'ax_67';
lstring = {'AL-tridiag'; 'MFISTA-5'; 'AL-P2'; 'ADMM-tridiag'; 'ADMM-FP-tridiag'};

%% --------------------------- save timing ---------------------------
lh = legend(lstring); set(lh, 'FontSize', FontSize);
title('');
save_im(db_path, sprintf('%s_iter_%s', short_slice_str, machine(1:3)))

close;

lh = legend(lstring); set(lh, 'FontSize', FontSize);
title('');
save_im(db_path, sprintf('%s_timing_%s', short_slice_str, machine(1:3)))


%% --------------------------- save images ---------------------------

figure; im(fftshift(1-samp))
axis off; title('');
save_im(db_path, sprintf('%s_samp', short_slice_str))

if ~isvar('xinf') && isvar('xtrue')
        brain_max = max(abs([xhat_AL(:); xtrue(:); xtrue(:); xhat_alp2(:)]));
else
        brain_max = max(abs([xhat_AL(:); xinit(:); xinf(:); xhat_alp2(:)]));
end
figure; im_brain(xhat_AL, 'clim', [0 brain_max], 'orient', orient)
axis off; title(''); colorbar;
save_im(db_path, sprintf('%s_x_AL', short_slice_str));

figure; im_brain(xhat_tri, 'clim', [0 brain_max], 'orient', orient)
axis off; title(''); colorbar;
save_im(db_path, sprintf('%s_x_tri', short_slice_str));

if isvar('body_coil') && isvar('xinf')
        [bc_scaled, bc_scale] = ir_wls_init_scale(1, xinf, body_coil);
        figure; im_brain(bc_scaled, 'clim', [0 brain_max], 'orient', orient)
        axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
        save_im(db_path, sprintf('%s_body_coil', short_slice_str));
end

figure; im_brain(xinit, 'clim', [0 brain_max], 'orient', orient)
axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
save_im(db_path, sprintf('%s_SoS_init', short_slice_str));

if isvar('xinf')
        xinf_tmp = xinf;
else
        xinf_tmp = xMFIS;
end
figure; im_brain(xinf_tmp, 'clim', [0 brain_max], 'orient', orient)
axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
save_im(db_path, sprintf('%s_x_MFISTA', short_slice_str));

figure; im_brain(xhat_AL - xinf_tmp, 'orient', orient)
axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
save_im(db_path, sprintf('%s_diff_AL_MFISTA', short_slice_str));

figure; im_brain(xhat_tri - xinf_tmp, 'orient', orient)
axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
save_im(db_path, sprintf('%s_diff_tri_MFISTA', short_slice_str));

figure; im_brain(bc_scaled - xinf_tmp, 'orient', orient)
axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
save_im(db_path, sprintf('%s_diff_BC_MFISTA', short_slice_str));

%%
if isvar('xtrue')
        figure; im_brain(xtrue, 'clim', [0 brain_max], 'orient', orient)
        axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
        save_im(db_path, sprintf('%s_x_true', short_slice_str));
        
        figure; im_brain(xhat_AL - xtrue, 'orient', orient)
        axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
        save_im(db_path, sprintf('%s_diff_AL_true', short_slice_str));
end