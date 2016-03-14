% setup tridiag exp

if ~isvar('use_raw')
	use_raw = 0;
end
wavelets = 0;
if ~isvar('home_path')
	tridiag_setup;
end
force_smap = 0;
if ~isvar('truncate')
	truncate = 0;
end
if ~strcmp(orient, 'sim')
	if ~isvar('slice')
		switch orient
			case 'axial'
				slice = 38;
				%slice = 90;
				slice = 67;
			case 'sagittal'
				slice = 69;
			case 'coronal'
				slice = 135; 
				slice = 155; 
			otherwise
				display(sprintf('need to choose slice for %s orientation', orient))
				keyboard
		end
	end
	if ~isvar('Sxtrue')
		[sense_maps, body_coil, Sxtrue, y_full] = invivo_exp(home_path, slice, 'orient', orient, 'force_smap', force_smap);
	end
else
	slice = 0;
	[sense_maps, body_coil, Sxtrue] = sim_setup();
end


if truncate == 1
	Sxtrue = Sxtrue(3:end-2, 3:end-2, :);
	sense_maps = sense_maps(3:end-2, 3:end-2, :);
	body_coil = body_coil(3:end-2, 3:end-2);
	if use_raw
		display('cannot use raw data if truncate');
		use_raw = 0;
	end
elseif truncate == 2
	Sxtrue = Sxtrue(:, 9:end-8, :);
	sense_maps = sense_maps(:, 9:end-8, :);
	body_coil = body_coil(:, 9:end-8, :);
	if use_raw
		display('cannot use raw data if truncate');
		use_raw = 0;
	end
end
[Nx, Ny, Nc] = size(Sxtrue);

% make sampling pattern
if ~isvar('reduction')
	reduction = 6;%6;
	display(sprintf('sampling factor set to: %d', reduction));
end
[slice_str, curr_folder] = get_exp_labels(orient, slice, reduction, use_raw, truncate);

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
if strcmp(orient, 'sim') || ~use_raw
	F = staticF(Nx, Ny, Nc, 'samp', samp, 'shift_img', false);
else
	F = staticF(Nx, Ny, Nc, 'samp', samp, 'shift_img', true);
end
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
	if use_raw
		y_noise = masker(y_full, repmat(samp, [1 1 Nc]));
	else
		y_noise = F*Sxtrue(:);
	end
else
	display('should have y_noise from invivo_exp');
	keyboard
end

% build mask
mask = generate_mask(orient, slice, Nx, Ny);

% initialize with SoS zero-fill solution
zero_fill = reshape(F'*y_noise, Nx, Ny, Nc)/(Nx * Ny);
SoS = sqrt(sum(abs(zero_fill).^2,3));
[xinit, scale] = ir_wls_init_scale(F*S, y_noise, SoS);
xinit = xinit.*mask;
if 0
        [xinit, scale] = ir_wls_init_scale(A, y_noise, SoS);
        [xinit_MFIS, scale_MFIS] = ir_wls_init_scale(A, y_center_noise, xMFIS);
        figure; im(cat(1, cat(2, SoS, xinit), cat(2, xMFIS, xinit_MFIS)))
        title(sprintf('SoS scale: %.2d + %.2di, MFISTA scale: %.2d + %.2di',real(scale),imag(scale),real(scale_MFIS),imag(scale_MFIS)))
        [~, testscale] = ir_wls_init_scale(1, xMFIS, SoS)
        % testscale = 0.1600 - 0.4581i
end



% parameters
beta = choose_beta(orient, slice, reduction);
plain_mu = num2cell(ones(1,5));

true_opt = 'avg';
xinf = load_x_inf(slice, beta, curr_folder, slice_str, 'method', true_opt);
if strcmp(orient, 'sim') && ~isscalar(xinf)
	mu_args = {'noise', 0.07*max([col(abs(CH*xinf)); col(abs(CV*xinf))])};
else
	mu_args = {};
end



