curr_path = '~/Dropbox/fessler/writing/tridiag/figs/';

%%
if size(x,3) > 1
        x = x(:,:,ceil(size(x,3)/2));
%         time_tri = time(:,3);
%         err_tri = err(:,3);
        time = time(:,3);
        err = err(:,3);
end
time_MFIS = time_FIS;
clear time_FIS

%% timing
orn_ndx = 400;
plot_axes = {[0 130 -84 0];[0 1400 -84 0]};

short_slice_str = 'ax_90_noFP';
colors = 'mgkbry'; % for w/o FP
markers = 'o+*sd.^v><ph'; % for w/o FP
order = [2 3 4 1];
y_val = 'err';

lstring = {'AL-tridiag';'MFISTA-5';'AL-P2-NC';'ADMM-tridiag'};%;'ADMM-FP-tridiag'};
plot_timing
%%

%%
figure; im(xtrue.', [0 1])
axis off; title(''); %colorbar
save_im(curr_path, 'inpaint_xtrue')
figure; im(xinit.', [0 1])
axis off; title(''); %colorbar
save_im(curr_path, sprintf('%s/inpaint_xinit_SNR%d_r%1.2d', obj, SNR, reduce))
%%
full_err = calc_NRMSE_over_mask(x_circ, xtrue);
circ_diff = abs(x_circ - xtrue);
figure; imshow(circ_diff, [0 0.1]);
colormap(flipud(colormap));
cbarh = colorbar;
set(cbarh, 'FontSize', 20);
delete(cbarh);
axis off; title(''); %title(sprintf('NRMSE: %2.2f', full_err)); %colorbar
save_im(curr_path, sprintf('circ_err_SNR%d_r%1.2d', SNR, reduce))
xlim([470 540]); ylim([280 432]);
save_im(curr_path, sprintf('%s/circ_err_SNR%d_r%1.2d_detail1', obj, SNR, reduce))

%%
figure; im(x_circ.', [0 1])
axis off; title(''); %colorbar
save_im(curr_path, sprintf('%s/inpaint_circ_SNR%d_r%1.2d', obj, SNR, reduce))

%%
full_err = calc_NRMSE_over_mask(x, xtrue);
diff = abs(x - xtrue);
figure; imshow(diff,[0 0.1]);
colormap(flipud(colormap));
cbarh = colorbar;
set(cbarh, 'FontSize', 20);
delete(cbarh)
axis off; title('');%title(sprintf('NRMSE: %2.2f', full_err)); %colorbar
save_im(curr_path, sprintf('%s/err_SNR%d_r%1.2d', obj, SNR, reduce))
xlim([470 540]); ylim([280 432]);
save_im(curr_path, sprintf('%s/err_SNR%d_r%1.2d_detail1', obj, SNR, reduce))

%%
figure; im(x.', [0 1])
axis off; title(''); %colorbar
save_im(curr_path, sprintf('%s/inpaint_SNR%d_r%1.2d', obj, SNR, reduce))

%%
figure; im(abs(x_circ - xtrue).' - abs(x - xtrue).', [-0.1 0.1]);

%%

figure; plot(cumsum(time(:,1),1), err(:,1), 'b--')
hold on; plot(cumsum(time_P2), err_P2, 'k--')
