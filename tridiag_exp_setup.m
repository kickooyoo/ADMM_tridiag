% setup tridiag exp
if ~issim
	slice = 38;
	% load in in vivo data
	if ~isvar('Sxtrue')
		[sense_maps, body_coil, Sxtrue] = invivo_exp(home_path, slice);
	end
else
	slice = 0;
	[sense_maps, body_coil, Sxtrue] = sim_setup();
end

[Nx, Ny, Nc] = size(Sxtrue);
truncate = 0;

wavelets = 0;

if truncate
	Sxtrue = Sxtrue(3:end-2, 3:end-2, :);
	sense_maps = sense_maps(3:end-2, 3:end-2, :);
	body_coil = body_coil(3:end-2, 3:end-2);
	[Nx, Ny] = size(body_coil);
end

% make sampling pattern
if ~isvar('samp')
	% generate sampling pattern
	if 1	
		samp_fname = sprintf('PD_sampling_%dx%d_R6_center8.mat', Nx, Ny);
		if exist(samp_fname)
			load(samp_fname);
		else
			reduction = 6;
			params.Nx = Nx;
			params.Ny = Ny;
			params.Nf = 1;
			params.h = 1;
			params.R = round(Nx*Ny/reduction);
			samp = logical(gen_new_sampling_pattern(params, 'all_kspace', 0, 'ellipse', 0, 'ncenter', 8, 'vardensity', 0));
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
S = GsplineS(sense_maps, 1);
F = GsplineF(Nx, Ny, 1, Nc, 'samp', samp);
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
y = F * Sxtrue(:);
SNR = 40;
sig = 10^(-SNR/20) * norm(y) / sqrt(length(y));
rng(0, 'twister')
y_noise = y + sig*randn(size(y)) + 1i*sig*randn(size(y));


% initialize with SoS zero-fill solution
center_samp = logical(coverDC_SamplingMask(zeros(Nx, Ny), 16, 16));
F_center = GsplineF(Nx, Ny, 1, Nc, 'samp', center_samp);
y_center = F_center * Sxtrue(:);
y_center_noise = y_center + sig*randn(size(y_center)) + 1i*sig*randn(size(y_center));
zero_fill = reshape(F_center'*y_center_noise, Nx, Ny, Nc)/(Nx * Ny);
SoS = sqrt(sum(abs(zero_fill).^2,3));
% [xinit, scale] = ir_wls_init_scale(F*S, y_noise, SoS);
%SoS_compensate = sum(conj(sense_maps).*Sxtrue,3)./(sum(abs(sense_maps).^2,3));
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
if slice == 67
	mask = generate_mask('slice67',1,Nx,Ny);
elseif slice == 38
	mask = generate_mask('slice38',1,Nx,Ny);
else
	mask = true(Nx, Ny);
end

% parameters
if slice == 67
	beta = 2^19;
elseif slice == 38
	beta = 2^24; % for slice 38 and l2b = 12 samp
	beta = 2^28; % for slice 38 and l2b = 16 samp
	beta = 2^25; % for slice 38 and sathish samp and sathish samp
	beta = 2^20; % for slice 38 and new samp R=6, sathish smap
else	
	display(sprintf('unknown beta for slice %d', slice));
	keyboard;
end
% convergence parameters
plain_mu = num2cell(ones(1,5));

