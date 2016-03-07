function [sense_maps, body_coil, Sxtrue, y] = invivo_exp(home_path, slice, varargin);
%function [smaps, body_coil, Sxtrue] = invivo_exp(home_path, slice, varargin);
% sets up slice 38 experiment from retrospectively sampled in vivo data
%
% varargin:
%	orient (string) 'axial', 'sagittal', 'coronal'
% 	base (string) './reviv', path to saved .mat file tree

arg.orient = 'axial';
arg.base = './reviv';
arg.force_smap = false;
arg = vararg_pair(arg, varargin);

fn = [home_path 'Documents/data/2010-07-06-fessler-3d/raw/p23040-3dspgr-8ch.7'];

Nc = 8;
[im4d, d] = recon3dft(fn, Nc);
[body_coil_images, dbc] = recon3dft([home_path 'Documents/data/2010-07-06-fessler-3d/raw/p24064-3dspgr-body.7'],1);
[Nx, Ny, Nz, Nc] = size(im4d);

switch arg.orient
case 'axial'
	assert(slice < Nz, sprintf('slice choice %d not possible, only %d slices axially', slice, Nz));
	mapped_im = squeeze(im4d(:,:,slice,:));
	body_coil = squeeze(body_coil_images(:,:,slice));
	d = ifft(fftshift(d),[],3);
        shift_slice = fftshift(1:Nz);
	y = squeeze(d(:,:, shift_slice(slice),:));
case 'sagittal'
	assert(slice < Ny, sprintf('slice choice %d not possible, only %d slices sagitally', slice, Ny));
	mapped_im = squeeze(im4d(:,slice,:,:));
	body_coil = squeeze(body_coil_images(:,slice,:));
	y = squeeze(d(:,slice,:,:));
case 'coronal'
	assert(slice < Nx, sprintf('slice choice %d not possible, only %d slices coronally', slice, Nx));
	mapped_im = squeeze(im4d(slice,:,:,:));
	body_coil = squeeze(body_coil_images(slice,:,:));
	y = squeeze(d(slice,:,:,:));
otherwise
	display(sprintf('unknown orientation: %s', arg.orient))
	keyboard;
end

if strcmp(arg.orient, 'axial') && (slice == 38)
	smap_fname = sprintf('%sDocuments/data/2010-07-06-fessler-3d/slice38/ramani/Smaps%d.mat', home_path, slice);
elseif strcmp(arg.orient, 'axial')
	smap_fname = sprintf('%s/axial/axial_slice%d_smap.mat', arg.base, slice);
elseif strcmp(arg.orient, 'sagittal')
	smap_fname = sprintf('%s/sagittal/sag_slice%d_smap.mat', arg.base, slice);
elseif strcmp(arg.orient, 'coronal')
	smap_fname = sprintf('%s/coronal/cor_slice%d_smap.mat', arg.base, slice);
end
if exist(smap_fname) && ~arg.force_smap
	load(smap_fname, '*map*'); 
else
        centersamp = logical(coverDC_SamplingMask(zeros(Nx, Ny), 16, 16));
        F_center = staticF(Nx, Ny, 1, 'samp', centersamp);
        mapped_im_lp = F_center'*F_center*mapped_im;
        body_coil_lp = F_center'*F_center*body_coil;
	display(sprintf('cannot load sense maps for slice %d', slice));
	display('try est_S_reg')
	sense_maps = est_S_reg(mapped_im, 'bodycoil', body_coil, 'l2b', 4);
	keyboard
	save(smap_fname, 'sense_maps');
end
if ~isvar('sense_maps') && isvar('Smap_QPWLS')
	sense_maps = Smap_QPWLS;
elseif  ~isvar('sense_maps') && isvar('smap')
	sense_maps = smap;
elseif ~isvar('sense_maps')
	display('not sure which var is smap');
	keyboard
end
Sxtrue = mapped_im;

