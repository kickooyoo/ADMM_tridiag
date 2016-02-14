function mu = get_mu(SS, RR, Nr, lambda, varargin)
%function mu = get_mu(eigvalss)
% set mu for nice condition numbers, based on AL-P2 splits

% vals for x_tri_inf with beta = 2^19 and SNR 40 on slice 67
arg.edge = 13600; % pixel diff for an edge
arg.noise = 9300; % pixel diff for noise
%arg.edge = 13600*10; % pixel diff for an edge
%arg.noise = 9300*10; % pixel diff for noise
arg.split = 'AL-P2'; % 'ADMM-tridiag'
arg.author = 'ramani';
arg.alph = 0.5; % tridiag design param
arg = vararg_pair(arg, varargin);

SSmax = max(col(abs(SS)));
SSmin = min(col(abs(SS)));
if isempty(RR)
        display('what to do for noncirc tridiag?');
        RRmax = 4;
        RRmin = 0;
else
        RRmax = max(col(abs(RR)));
        RRmin = min(col(abs(RR)));
end

switch arg.author
        case 'mai'
                mu_v = mean(lambda./[arg.edge arg.noise]);
                
                % kappa_u = (Nr + mu_u)/mu_u
                %ku = @(mu, Nr) (Nr + mu(1)) / mu(1);
                
                % kappa_z = (mu_v * RRmax + mu_z)/(mu_v * RRmin + mu_z) for only horizontal and vertical finite diff
                %kz = @(mu) (mu(2) * RRmax + mu(3)) / (mu(2) * RRmin + mu(3);
                
                % kappa_x = (mu_u * SSmax + mu_z)/(mu_u * SSmin + mu_z)
                %kx = @(mu, SSmax, SSmin) (mu(1) * SSmax + mu(3)) / (mu(1) * SSmin + mu(3));
                
                % x = [mu_u, mu_z]
                kappas = @(x) double([ (Nr + x(1)) / x(1); (mu_v * RRmax + x(2)) / (mu_v * RRmin + x(2)); (x(1) * SSmax + x(2)) / (x(1) * SSmin + x(2))]);
                old_kappas = @(x) [ (Nr + x(1)) / x(1); (4 + x(2)) / (x(2)); (x(1) * SSmax + x(2)) / (x(1) * SSmin + x(2))];
                
                
                x0 = ones(1,2);
                %for ii = 1:10
                x = lsqnonlin(kappas, x0, zeros(1,2), Inf(1,2));
                %end
                %display('there are some vals');
                %keyboard
                %kappas(x)
                
                
                switch arg.split
                        case 'AL-P2'
                                mu = [x(1); mu_v; x(2)];
                        case 'ADMM-tridiag'
                                mu = [mu_v; mu_v; x(1); x(2); x(2)];
                        otherwise
                                display(sprintf('unknown splitting scheme %s', arg.split));
                                keyboard
                end
                
        case 'ramani'
                kapx = 0.9*SSmax/SSmin;
                kapu = 24; % condition number of (FF + mu*I)
                kapz = 12; % condition number of (RR + nu2/n1*I)`
                switch arg.split
                        case 'AL-P2'
                                % mu -> mu_u
                                % mu*nu2 -> mu_z
                                % mu*nu1 -> mu_v
                                nu2 = (SSmax - SSmin * kapx) / (kapx - 1);
                                mu_P2 = Nr / (kapu - 1); 
                                nu1 = (kapz - 1) * nu2 / (RRmax - RRmin * kapz);
                                
                                mu = mu_P2*[1 nu1 nu2];
                        case 'ADMM-tridiag'
                                % brainstorm:
                                % enforce mu_0 = mu_1 for
                                % symmetry
                                
                                nu2 = (SSmax - SSmin * kapx) / (kapx - 1);
                                mu_P2 = Nr / (kapu - 1); 
                                nu1 = (kapz - 1) * nu2 / (RRmax - RRmin * kapz);
                                
                                
                                mu_2 = Nr / (kapu - 1); 
                                
                                SSaprox = mean(col(SS));
                                mu_0 = nu1*mu_P2; % same as AL-P2
                                mu_1 = mu_0;
                                mu_3 = kapz*mu_0 - (mu_2*(1-arg.alph).^2*SSaprox);
                                mu_4 = kapz*mu_1 - (mu_2*(arg.alph).^2*SSaprox);

                                mu = [mu_0 mu_1 mu_2 mu_3 mu_4];
                        otherwise
                                display(sprintf('unknown splitting scheme %s', arg.split));
                                keyboard
                end
        otherwise
                display('unknown author option for get_mu.m')
end