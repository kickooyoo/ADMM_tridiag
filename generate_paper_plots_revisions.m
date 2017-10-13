% script for making plots and images on vega
db_path = '~/Dropbox/fessler/writing/tridiag/';
dbr_path = '~/Dropbox/fessler/experimental_data/tridiag_vega/';
FontSize = 12;
FP_incl = 0;
%% sagittal slice 69 --------------------------------------------------------
% clearvars -except db_path dbr_path FontSize FP_incl
eval(load_except([dbr_path 'sagittal_slice69/ir7_timing_256x128_5000iter_sagittal69_avgtrue_8388608beta.mat'], 'db_path'))
clear time*inf
orn_ndx = 800;
plot_axes = {[0 150 -88 0];[0 2500 -88 0]};
clear *circ
if FP_incl
	short_slice_str = 'sag_69';
	order = [1 2 3 4];
else
	clear *FP
	short_slice_str = 'sag_69_noFP';
	colors = 'mgkbry'; % for w/o FP
	markers = 'o+*sd.^v><ph'; % for w/o FP
	order = [1 2 3 4];
end

ystr = 'NRMSE to x_{SENSE} (dB)';
lstring = {'AL-tridiag';'MFISTA-5';'AL-P2-NC';'ADMM-tridiag'; 'ADMM-FP-tridiag'};
plot_timing
%%
eval(load_except([dbr_path 'sagittal_slice69/ir7_timing_256x128_5000iter_sagittal69_avgtrue_8388608beta.mat'], 'db_path'))
clear time*inf
orn_ndx = 50;
plot_axes = {[0 8 -20 0];[0 150 -20 0]};
if FP_incl
	short_slice_str = 'sag_69_err';
	order = [5 2 3 4 1];
else
	clear *FP
	short_slice_str = 'sag_69_err_noFP';
	colors = 'mgkbry'; % for w/o FP
	markers = 'o+*sd.^v><ph'; % for w/o FP
	order = [2 3 4 1];
end
ystr = 'NRMSE to x_{SENSE} (dB)';
lstring = {'AL-tridiag';'MFISTA-5';'AL-P2-NC';'ADMM-tridiag';'ADMM-FP-tridiag'};
plot_timing

%% sim sagittal  --------------------------------------------------------
clearvars -except db_path dbr_path FontSize order
eval(load_except([dbr_path 'sim/ir7_timing_240x200_5000iter_sim_truetrue_8192beta.mat'], 'db_path'))
slice_str = [slice_str '_sag'];
clear time_circ

orn_ndx = 40;
plot_axes = {[0 10 -33 0];[0 200 -33 0]};
short_slice_str = 'sim_sag_err_r2';
order = [5 2 3 4 1];

ystr = 'NRMSE to x_{true} (dB)';
lstring = {'AL-tridiag'; 'MFISTA-5'%'AL-tridiag,svt';
        'AL-P2-NC';'ADMM-tridiag';'ADMM-FP-tridiag'};
plot_timing

%%

clearvars -except db_path FontSize
eval(load_except([db_path 'sim/ir7_timing_240x200_500iter_sim_truetrue.mat'], 'db_path'))
orn_ndx = 15;
plot_axes = {[0 9 -22 0];[0 150 -22 0]};
lstring = {'AL-tridiag'; 'MFISTA-5'%'AL-tridiag,svt';
        'AL-P2';'ADMM-tridiag';'ADMM-FP-tridiag'};
plot_timing
short_slice_str = 'sim_sag_err';
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
if strcmp(orient, 'sagittal')
        brain_max = brain_max*0.5;
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

if isvar('x_circ')
        figure; im_brain(x_circ, 'clim', [0 brain_max], 'orient', orient)
        axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
        save_im(db_path, sprintf('figs/%s_x_circ', short_slice_str));
        figure; im_brain(xhat_AL - x_circ, 'orient', orient)
        axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
        save_im(db_path, sprintf('figs/%s_diff_AL_circ', short_slice_str));
        ylim([128-50+1 128]); xlim([101 256-100])
        save_im(db_path, sprintf('figs/%s_diff_AL_circ_zoom', short_slice_str));
        if strcmp('true_opt', 'true')
                
                if strcmp(orient, 'sim')
                        SENSE_str = 'true';
                else
                        SENSE_str = 'SENSE';
                end
                figure; im_brain(xhat_AL - xinf, 'orient', orient)
                axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
                save_im(db_path, sprintf('figs/%s_diff_AL_%s', short_slice_str, SENSE_str));
               
                figure; im_brain(x_circ - xinf, 'orient', orient)
                axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
                save_im(db_path, sprintf('figs/%s_diff_circ_%s', short_slice_str, SENSE_str));
                
                figure; im_brain(abs(x_circ - xinf) - abs(xhat_tri - xinf) , 'orient', orient, 'abs', 0)
                title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
                axis on
                set(gca,'ytick',[])
                set(gca,'yticklabel',[])
                set(gca,'xtick',[])
                set(gca,'xticklabel',[])
                save_im(db_path, sprintf('figs/%s_diff_of_circ_tri_diffs_%s', short_slice_str, SENSE_str));
                axis([50 120 180 200 ])
                axis on
                set(gca,'ytick',[])
                set(gca,'yticklabel',[])
                save_im(db_path, sprintf('figs/%s_diff_of_circ_tri_diffs_%s_zoom', short_slice_str, SENSE_str));
                
                figure; im_brain(x_circ - xinf, 'orient', orient)
                axis off; title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
                save_im(db_path, sprintf('figs/%s_diff_circ_%s', short_slice_str, SENSE_str));
                
                figure; im_brain(abs((x_circ - xinf) - (xhat_tri - xinf)) , 'orient', orient, 'abs', 0)
                title(''); ch = colorbar; set(ch, 'FontSize', FontSize);
                axis on
                set(gca,'ytick',[])
                set(gca,'yticklabel',[])
                set(gca,'xtick',[])
                set(gca,'xticklabel',[])
                save_im(dbr_path, sprintf('figs/%s_diff_of_complex_circ_tri_diffs_%s', short_slice_str, SENSE_str));
                axis([50 120 180 200 ])
                axis on
                set(gca,'ytick',[])
                set(gca,'yticklabel',[])
                save_im(dbr_path, sprintf('figs/%s_diff_of_complex_circ_tri_diffs_%s_zoom', short_slice_str, SENSE_str));
                
        end
end

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
