% test tridiag inpaint
if ~isvar('wavelets')
	wavelets = 1;
end
% ------ construct true image ------
fname = 'textbook_contrast.jpg';
xtrue = imread(fname);
if size(xtrue, 3) > 1
	xtrue = rgb2gray(xtrue);
end
xtrue = double(xtrue);
xtrue = downsample2(xtrue, 7);
% make even dimensions
if mod(size(xtrue, 1), 2) == 1
	xtrue = xtrue(2:end, :);
end
if mod(size(xtrue, 2), 2) == 1
	xtrue = xtrue(:, 2:end);
end
xtrue = xtrue(:,1:540);
xtrue = xtrue./max(xtrue(:));
[Nx, Ny] = size(xtrue);

% ------ take measurements ------
rng(0);
if ~isvar('reduce')
	reduce = 2;
end
samp = (rand(Nx, Ny) <= 1/reduce);
D = Ginpaint(samp);
[CH, CV] = construct_finite_diff([Nx Ny]);
if ~isvar('SNR')
	SNR = 10;
end
y_noiseless = D * xtrue;
sig = 10^(-SNR/20) * norm(y_noiseless) / sqrt(length(y_noiseless));
y = y_noiseless + sig*randn(size(y_noiseless));
% ------ initialize x with nearest neighbors ------
[xx, yy] = ndgrid(1:Nx, 1:Ny);
xx_D = xx(samp);
yy_D = yy(samp);
xinit = griddata(xx_D, yy_D, y, xx, yy, 'nearest');
% ------ construct regularizers ------- 
R = [CH; CV];
Rcirc = Cdiffs([Nx Ny],'offsets', [1 Nx], 'type_diff','circshift');
if wavelets
	W = Godwt1(true(Nx, Ny));
end
% ------ optimization params ------
niters = 100;
mu0 = 1;
mu1 = 1;
mu2 = 1;
alph = 0.5;
alphw = 0.5;
% ------ regularization parameters for SNR = 20, reduce = 1.5 ------
beta_search_fname = sprintf('wavelet%d_SNR%d_reduce%1.2d.mat', wavelets, SNR, reduce);
d = dir('inpainting_mat');
for ii = 1:length(d)
	if ~isempty(strfind(d(ii).name, beta_search_fname))
		load(d(ii).name, 'betas', 'betaws', 'betas_circ', 'betaws_circ');
	end
end
if ~isvar('betas')
	display(sprintf('no file found matching %s, if doing reg search no prob', beta_search_fname));
else
	beta = betas(ceil(length(betas)/2));
	betaw = betaws(ceil(length(betaws)/2));
	beta_circ = betas_circ(ceil(length(betas_circ)/2));
	betaw_circ = betaws_circ(ceil(length(betaws_circ)/2));
end
	%beta = 0.02743;
	%beta_circ = 0.02743;
% ------
% wavelets, CH, CV, R, Rcirc, RcircW, D, y, xinit, niters, mu0, mu1, mu2

