function [sense_maps, body_coil, Sxtrue] = invivo_exp(home_path, slice);
%[smaps, body_coil, xtrue] = function invivo_exp(home_path, slice);
% sets up slice 38 experiment from retrospectively sampled in vivo data

fn = [home_path 'Documents/data/2010-07-06-fessler-3d/raw/p23040-3dspgr-8ch.7'];

nc = 8;
[im4d, d] = recon3dft(fn,nc);
[body_coil_images, d] = recon3dft([home_path 'Documents/data/2010-07-06-fessler-3d/raw/p24064-3dspgr-body.7'],1);

mapped_im = squeeze(im4d(:,:,slice,:));

smap_fname = sprintf('%sDocuments/data/2010-07-06-fessler-3d/slice38/ramani/Smaps%d.mat', home_path, slice);
if exist(smap_fname)
	load(smap_fname, 'Smap_QPWLS'); 
else
	display(sprintf('cannot load sense maps for slice %d', slice));
	keyboard
end
sense_maps = Smap_QPWLS;
body_coil = body_coil_images(:,:,slice);
Sxtrue = mapped_im;

