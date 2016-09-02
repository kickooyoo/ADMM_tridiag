curr_path = '~/Dropbox/fessler/writing/tridiag/figs/';
%%
figure; im(xtrue.', [0 1])
axis off; title(''); %colorbar
save_im(curr_path, 'inpaint_xtrue')
figure; im(xinit.', [0 1])
axis off; title(''); %colorbar
save_im(curr_path, sprintf('inpaint_xinit_SNR%d_r%1.2d', SNR, reduce))
%%
full_err = calc_NRMSE_over_mask(x_circ, xtrue);
figure; im(abs(x_circ - xtrue).',[0 0.1]);
axis off; title(sprintf('NRMSE: %2.2f', full_err)); %colorbar
save_im(curr_path, 'circ_err_SNR10_r2')
xlim([450 540]); ylim([300 432]);
save_im(curr_path, sprintf('err_SNR%d_r%1.2d_detail1', SNR, reduce))

%%
figure; im(x_circ.', [0 1])
axis off; title(''); %colorbar
save_im(curr_path, sprintf('inpaint_circ_SNR%d_r%1.2d', SNR, reduce))

%%
full_err = calc_NRMSE_over_mask(x(:,:,3), xtrue);
figure; im(abs(x(:,:,3) - xtrue).',[0 0.1]);
axis off; title(sprintf('NRMSE: %2.2f', full_err)); %colorbar
save_im(curr_path, sprintf('err_SNR%d_r%1.2d', SNR, reduce))
xlim([450 540]); ylim([300 432]);
save_im(curr_path, sprintf('err_SNR%d_r%1.2d_detail1', SNR, reduce))

%%
figure; im(x(:,:,3).', [0 1])
axis off; title(''); %colorbar
save_im(curr_path, sprintf('inpaint_SNR%d_r%1.2d', SNR, reduce))

%%
figure; im(abs(x_circ - xtrue).' - abs(x(:,:,3) - xtrue).', [-0.1 0.1]);

%%

figure; plot(cumsum(time(:,1),1), err(:,1), 'b--')
hold on; plot(cumsum(time_P2), err_P2, 'k--')
