function [sense_maps, body_coil, Sxtrue] = sim_setup();
%function [sense_maps, body_coil, Sxtrue] = sim_setup();
% setup file for simulation image
%

Nx = 240;
Ny = 200;
Nc = 6;%12;

img = make_sim_image(Nx, Ny);


%sense_maps = mri_sensemap_sim('nx', Nx, 'ny', Ny, 'ncoil', Nc, 'rcoil', 3*min(Nx, Ny));

% for sagittal orient
sense_maps = ir_mri_sensemap_sim('nx', Nx, 'ny', Nx, 'nz', Ny, 'ncoil', 2*Nc, 'rcoil', 3*min(Nx,Ny), 'nring', round(2*Nc/4));
sense_maps = squeeze(sense_maps(:,round(end/2),:,:));
if Nc == 6
sense_maps = sense_maps(:,:,1:2:end);
end

mask = generate_mask('sim', 0, Nx, Ny);
sense_maps = truncate_sense_maps(sense_maps, mask);

Sxtrue = repmat(img, [1 1 size(sense_maps,3)]).*sense_maps;
body_coil = img;

