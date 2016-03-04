function mu = get_mu(SS, RR, Nr, lambda, varargin)
%function mu = get_mu(eigvalss)
% set mu for nice condition numbers, based on AL-P2 splits

% vals for x_tri_inf with beta = 2^19 and SNR 40 on slice 67
arg.noise = 24000; % pixel diff for noise (default val for experiment with reduction of 6)
arg.split = 'AL-P2'; % 'ADMM-tridiag'
arg.alph = 0.5; % tridiag design param
arg.mask = true(size(SS));
arg.fancy_mu34 = true;
arg.mu0_fudge = 0.32;% 0.01;% found to be best for slice 38 []; % for ramani, tridiag
arg.fancy_mu01 = false;
arg.test = '';
arg.ktri = 12;
arg.kapu_tri = 24;
arg.CH = [];
arg.CV = [];
arg.set_FP = [];
arg = vararg_pair(arg, varargin);

SSmax = max(col(abs(SS)));
SSmin = min(col(abs(SS)));
SSmaxm = SSmax;%max(col(abs(SS(arg.mask(:)))));
SSminm = SSmin;%min(col(abs(SS(arg.mask(:)))));
[Nx Ny] = size(SS);
if isempty(RR)
        if ~isempty(strfind(arg.test, 'RRapprox'))
                CCx = sparse(diag(-1*ones(1,Nx-1),-1) + diag(cat(2, 1, 2*ones(1,Nx - 2), 1)) + diag(-1*ones(1,Nx-1),1));
                CCy = sparse(diag(-1*ones(1,Ny-1),-1) + diag(cat(2, 1, 2*ones(1,Ny - 2), 1)) + diag(-1*ones(1,Ny-1),1));
                [V,D] = eigs(CCx);
                RRmaxx = max(D(:));
                [V,D] = eigs(CCy);
                RRmaxy = max(D(:));
        else
                RRmaxx = 4;
                RRmaxy = 4;
        end
        RRmax = 4;
        RRmin = 0;
else
        RRmax = max(col(abs(RR)));
        RRmin = min(col(abs(RR)));
end


kapx = 0.9*SSmaxm/SSminm;
kapu = 24; % condition number of (FF + mu*I)
kapz = 12; % condition number of (RR + nu2/n1*I)`
switch arg.split
        case 'AL-P2'
                % mu -> mu_u
                % mu*nu2 -> mu_z
                % mu*nu1 -> mu_v
                nu2 = (SSmaxm - SSminm * kapx) / (kapx - 1);
                mu_P2 = Nr / (kapu - 1);
                nu1 = (kapz - 1) * nu2 / (RRmax - RRmin * kapz);
                
                mu = {mu_P2; mu_P2*nu1; mu_P2*nu2};
        case 'ADMM-tridiag'
                mu_2 = Nr / (arg.kapu_tri - 1);
                if arg.fancy_mu34
                        fudge = arg.mu0_fudge;
                        if arg.fancy_mu01
                                fudge = fudge.*arg.mask;
                        end
                        mu_0 = mu_2*(1-arg.alph).^2*SSmax*(arg.ktri - 1)/RRmaxx + fudge;
                        mu_1 = mu_2*(1-arg.alph).^2*SSmax*(arg.ktri - 1)/RRmaxy + fudge;
                        %                         if arg.fancy_mu01
                        %                                 mu_0(~arg.mask(:)) = mu_0(~arg.mask(:))./10;
                        %                                 mu_1(~arg.mask(:)) = mu_1(~arg.mask(:))./10;
                        %                         end
                        
                        mu_3 = RRmaxy*mu_1/(arg.ktri - 1) - (mu_2*(1-arg.alph).^2*SS);
                        mu_4 = RRmaxx*mu_0/(arg.ktri - 1) - (mu_2*(arg.alph).^2*SS);
                        if arg.fancy_mu01
                                figure; subplot(2,2,1); im(mu_3);
                                subplot(2,2,2); im(mu_4)
                                mu_3 = max(mu_3, 1e-3);
                                mu_4 = max(mu_4, 1e-3);
                                subplot(2,2,3); im(mu_3);
                                subplot(2,2,4); im(mu_4)
                        end
                        if ~isempty(strfind(arg.test, 'edge'))
                                vert_edge = zeros(size(SS));
                                vert_edge(:, [1 end]) = 1;
                                horz_edge = zeros(size(SS));
                                horz_edge([1 end], :) = 1;
                                mu_3 = mu_3 + mu_1.*vert_edge;
                                mu_4 = mu_4 + mu_0.*horz_edge;
                        end
                        %c = mu_2*arg.alph.^2*SS + mu_3;
                        %Hx = sparse(diag(-mu0*ones(1,Nx-1), -1) + diag(-mu0
                else
                        mu_0 = lambda/arg.noise;
                        mu_1 = mu_0;
                        SSaprox = mean(col(SS));
                        mu_3 = kapz*mu_0 - (mu_2*(1-arg.alph).^2*SSaprox);
                        mu_4 = kapz*mu_1 - (mu_2*(arg.alph).^2*SSaprox);
                        if (mu_3 < 0) || (mu_4 < 0)
                                display('invalid mu_3 or mu_4');
                                keyboard;
                        end
                end
                
                mu = {mu_0; mu_1; mu_2; mu_3; mu_4};
                
                if any(col(mu_0) < 0) || any(col(mu_1)< 0) || mu_2 < 0 || any(col(mu_3) < 0) || any(col(mu_4) < 0)
                        %if any(col(mu_0(arg.mask(:))) < 0) || any(col(mu_1(arg.mask(:)))< 0) || mu_2 < 0 || any(col(mu_3(arg.mask(:))) < 0) || any(col(mu_4(arg.mask(:))) < 0)
                        display('invalid negative mu');
                        keyboard;
                end
                
                if any(imag(col(mu_0)) > 0) || any(imag(col(mu_1)) > 0) ||  imag(mu_2) > 0 || any(imag(col(mu_3)) > 0) || any(imag(col(mu_4)) < 0)
                        display('invalid complex mu');
                        keyboard;
                end
                
                if ~isempty(strfind(arg.test, 'check')) && prod(Nx, Ny) < 10000
                        sub3 = col(-1*mu_1*(cat(1, ones(Ny - 1, Nx), zeros(1, Nx))));
                        H3 = diag(sub3(1:end-1), -1) + ...
                                diag(col(mu_1*(cat(1, ones(1, Nx), 2*ones(Ny - 2, Nx), ones(1, Nx))) + mu_2*(1-arg.alph).^2*SS.' + mu_3.'))+ ...
                                diag(sub3(1:end-1), 1);
                        subx = col(-1*mu_0*(cat(1, ones(Nx - 1, Ny), zeros(1, Ny))));
                        Hx = diag(subx(1:end-1), -1) + ...
                                diag(col(mu_0*(cat(1, ones(1, Ny), 2*ones(Nx - 2, Ny), ones(1, Ny))) + mu_2*arg.alph.^2*SS + mu_4))+ ...
                                diag(subx(1:end-1), 1);
                        [~, D3] = eigs(sparse(double(H3)));
                        [~, Dx] = eigs(sparse(double(Hx)));
                        keyboard
                end
        case 'ADMM-FP-tridiag'
                mu_0 = lambda/arg.noise;
                mu_2 = mu_0;                
                mu_4 = Nr / (arg.kapu_tri - 1);
               	if isempty(arg.set_FP) 
			mu_3 = mu_0;
			c3 = mu_3 * RRmaxy / arg.ktri;
			mu_5 = min(mu_4, c3/((1-arg.alph).^2* SSmax));
			mu_7 = c3 - mu_5 * (1-arg.alph).^2 * SS + arg.mu0_fudge; 
			
			mu_1 = mu_0;
			cx = mu_1 * RRmaxx / arg.ktri;
			mu_6 = min(mu_4, cx/(arg.alph.^2* SSmax));
			mu_8 = cx - mu_6 * arg.alph.^2 * SS + arg.mu0_fudge; 
		else
			mu_0 = arg.set_FP(1);
			keyboard;	
		end
                if any(cat(1, mu_0, mu_1, mu_2, mu_3, mu_4, mu_5, mu_6, col(mu_7), col(mu_8)) < 0)
                        display('invalid negative mu');
                        keyboard;
                end
                mu = {mu_0; mu_1; mu_2; mu_3; mu_4; mu_5; mu_6; mu_7; mu_8};
                
        otherwise
                display(sprintf('unknown splitting scheme %s', arg.split));
                keyboard
end

