function [sense_maps, body_coil, Sxtrue] = invivo_exp(home_path, slice, varargin);
%[smaps, body_coil, xtrue] = function invivo_exp(home_path, slice);
% sets up slice 38 experiment from retrospectively sampled in vivo data
%
% varargin:
%	orient (string) 'axial', 'sagittal', 'coronal'
%

arg.orient = 'axial';
arg = vararg_pair(arg, varargin);

fn = [home_path 'Documents/data/2010-07-06-fessler-3d/raw/p23040-3dspgr-8ch.7'];

nc = 8;
[im4d, d] = recon3dft(fn,nc);
[body_coil_images, d] = recon3dft([home_path 'Documents/data/2010-07-06-fessler-3d/raw/p24064-3dspgr-body.7'],1);

[Nx, Ny, Nz, Nc] = size(im4d);

switch arg.orient
case 'axial'
	assert(slice < Nz, sprintf('slice choice %d not possible, only %d slices axially', slice, Nz));
	mapped_im = squeeze(im4d(:,:,slice,:));
	body_coil = body_coil_images(:,:,slice);

case 'sagittal'
	assert(slice < Ny, sprintf('slice choice %d not possible, only %d slices sagitally', slice, Ny));
	mapped_im = squeeze(im4d(:,slice,:,:));
	body_coil = body_coil_images(:,slice,:);

case 'coronal'
	assert(slice < Nx, sprintf('slice choice %d not possible, only %d slices coronally', slice, Nx));
	mapped_im = squeeze(im4d(slice,:,:,:));
	body_coil = body_coil_images(slice,:,:);

otherwise
	display(sprintf('unknown orientation: %s', arg.orient))
	keyboard;
end

if strcmp(arg.orient, 'axial') && (slice == 38)
	smap_fname = sprintf('%sDocuments/data/2010-07-06-fessler-3d/slice38/ramani/Smaps%d.mat', home_path, slice);
elseif strcmp(arg.orient, 'sagittal')
	smap_fname = sprintf('./reviv/sag_slice%d/sag_slice%d_smap.mat', slice, slice);
end
if exist(smap_fname)
	load(smap_fname, '*map*'); 
else
	display(sprintf('cannot load sense maps for slice %d', slice));
	display('try est_S_reg')
	keyboard
end
if isvar('Smap_QPWLS')
	sense_maps = Smap_QPWLS;
elseif isvar('smap')
	sense_maps = smap;
else
	display('not sure which var is smap');
	keyboard
end
Sxtrue = mapped_im;

