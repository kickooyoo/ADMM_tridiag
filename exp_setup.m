% setup tridiag exp

truncate = 0;
wavelets = 0;
if ~isvar('home_path')
	tridiag_setup;
end

if ~strcmp(orient, 'sim')
	switch orient
		case 'axial'
			slice = 38;
			slice = 90;
		case 'sagittal'
			slice = 69;
		case 'coronal'
			slice = 135; 
		otherwise
			display(sprintf('need to choose slice for %s orientation', orient))
			keyboard
	end
	if ~isvar('Sxtrue')
		[sense_maps, body_coil, Sxtrue, y_full] = invivo_exp(home_path, slice, 'orient', orient);
	end
else
	slice = 0;
	[sense_maps, body_coil, Sxtrue] = sim_setup();
end

[Nx, Ny, Nc] = size(Sxtrue);

if truncate
	Sxtrue = Sxtrue(3:end-2, 3:end-2, :);
	sense_maps = sense_maps(3:end-2, 3:end-2, :);
	body_coil = body_coil(3:end-2, 3:end-2);
	[Nx, Ny] = size(body_coil);
end

% make sampling pattern
if ~isvar('reduction')
	reduction = 6;%6;
	display(sprintf('sampling factor set to: %d', reduction));
end
[slice_str, curr_folder] = get_exp_labels(orient, slice, reduction);

if ~isvar('samp')
	% generate sampling pattern
	if 1	
		samp_fname = sprintf('./reviv/%s/PD_sampling_%dx%d_R%d_center8.mat', orient, Nx, Ny, reduction);
		if exist(samp_fname)
			load(samp_fname);
		else
			params.Nx = Nx;
			params.Ny = Ny;
			params.Nf = 1;
			params.h = 1;
			params.R = round(Nx*Ny/reduction);
			samp = logical(gen_new_sampling_pattern(params, 'all_kspace', 0, 'ellipse', 0, 'ncenter', 8, 'vardensity', 0, 'fudge', 0.75));
			save(samp_fname,'samp');
		end
	else
		load('al-p2/PDiskR256x256L1R2.252.25B0.8R1D1pctg20.2972.mat', 'SP', 'sampname');
		if size(SP, 1) > Nx
			SP = SP(1:Nx,:);
		elseif size(SP, 1) < Nx
			display('need bigger Poisson disk sampling pattern');
			keyboard;
		end
		if size(SP, 2) > Ny
			SP = SP(:,1:Ny);
		elseif size(SP, 2) < Ny
			display('need bigger Poisson disk sampling pattern');
			keyboard;
		end
		Ncentx = 8; % sample a window of 2*N around DC along x
		Ncenty = 8; % sample a window of 2*N around DC along y
		samp = logical(coverDC_SamplingMask(SP, Ncentx, Ncenty));
	end
end

% construct fatrices
S = staticS(sense_maps);
F = staticF(Nx, Ny, Nc, 'samp', samp);
[CH, CV] = construct_finite_diff([Nx Ny]); 
R = [CH; CV];
Rcirc = Cdiffs([Nx Ny],'offsets', [1 Nx], 'type_diff','circshift');
if wavelets
        W = Godwt1(true(Nx, Ny));
        betaw = beta;
        alphw = 0.5;
        RcircW = [Rcirc; W];
        CHW = [CH; betaw * alphw / beta * W];
        CVW = [CV; betaw * (1-alphw) / beta * W];     
	RW = [CHW; CVW];
end

% generate data
if strcmp(orient, 'sim')
	y = F * Sxtrue(:);
	SNR = 40;
	sig = 10^(-SNR/20) * norm(y) / sqrt(length(y));
	rng(0, 'twister')
	y_noise = y + sig*randn(size(y)) + 1i*sig*randn(size(y));
elseif isvar('y_full')
%	y_noise = masker(y_full, repmat(fftshift(samp), [1 1 Nc]));
	y_noise = masker(y_full, repmat(samp, [1 1 Nc]));
else
	display('should have y_noise from invivo_exp');
	keyboard
end
% initialize with SoS zero-fill solution
zero_fill = reshape(F'*y_noise, Nx, Ny, Nc)/(Nx * Ny);
SoS = sqrt(sum(abs(zero_fill).^2,3));
xinit = SoS;
if 0
        [xinit, scale] = ir_wls_init_scale(A, y_center_noise, SoS);
        [xinit_MFIS, scale_MFIS] = ir_wls_init_scale(A, y_center_noise, xMFIS);
        figure; im(cat(1, cat(2, SoS, xinit), cat(2, xMFIS, xinit_MFIS)))
        title(sprintf('SoS scale: %.2d + %.2di, MFISTA scale: %.2d + %.2di',real(scale),imag(scale),real(scale_MFIS),imag(scale_MFIS)))
        [~, testscale] = ir_wls_init_scale(1, xMFIS, SoS)
        % testscale = 0.1600 - 0.4581i
end

% build mask
switch orient
case 'axial'
	if slice == 67
		mask = generate_mask('slice67',1,Nx,Ny);
	elseif slice == 38
		mask = generate_mask('slice38',1,Nx,Ny);
	else
		mask = true(Nx, Ny);
	end
case 'coronal'
	mask = cat(1, zeros(5, 128), ones(144-13,128), zeros(8,128));
otherwise
	mask = true(Nx, Ny);
end


% parameters
beta = choose_beta(orient, slice, reduction);
plain_mu = num2cell(ones(1,5));

xinf = load_x_inf(slice, beta, curr_folder, slice_str);
if strcmp(orient, 'sim')
	mu_args = {'noise', 0.07*max([col(abs(CH*xinf)); col(abs(CV*xinf))])};
else
	mu_args = {};
end



