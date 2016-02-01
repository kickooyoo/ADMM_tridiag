% mai_setup.m
% run in /mai_code

[~, machine] = system('hostname');

switch machine(1:end-1) % last char is some sort of new line
    case {'eecs-IV1', 'iv1'}
        addpath('/n/ire/Volumes/s2/fessler/web/irt/irt')
        addpath(genpath('~/Documents/mai_code/ADMM_tridiag'))
        addpath(genpath('~/Documents/mai_code/util'))
        addpath(genpath('~/Documents/mai_code/spline_basis'))
        addpath(genpath('~/Documents/mai_code/pthread_tutor'))
        addpath(genpath('~/Documents/data/2010-07-06-fessler-3d/code'))
	home = '~/';
    case 'mpel8.eecs.umich.edu'
        addpath('/n/ire/Volumes/s2/fessler/web/irt/irt')
        addpath(genpath('~/iv1h/Documents/mai_code/ADMM_tridiag'))
        addpath(genpath('~/iv1h//Documents/mai_code/util'))
        addpath(genpath('~/iv1h//Documents/mai_code/spline_basis'))
        addpath(genpath('~/iv1h//Documents/mai_code/pthread_tutor'))
        addpath(genpath('~/iv1h//Documents/data/2010-07-06-fessler-3d/code'))
	home = '~/iv1h/';
    case 'ir63.eecs.umich.edu'
        addpath('/n/ire/Volumes/s2/fessler/web/irt/irt')
        addpath(genpath('~'))
	home = '~/iv1h';
	case 'mpel8.eecs.umich.edu'
        addpath('/n/ire/Volumes/s2/fessler/web/irt/irt')
        addpath(genpath('~/iv1h/Documents/mai_code/util'))
        addpath(genpath('~/iv1h/Documents/mai_code/spline_basis'))
        addpath(genpath('~/iv1h/Documents/mai_code/pthread_tutor'))
        addpath(genpath('~/iv1h/Documents/data/2010-07-06-fessler-3d/code'))
	home = '~/iv1h/';
    case 'vega'
        addpath('/Users/mai/Documents/irt')
        addpath(genpath('/Users/mai/Documents/mai_code'))
        addpath(genpath('/Users/mai/Documents/mai_code/ADMM_tridiag'))
        addpath(genpath('/Users/mai/Documents/contrib'))
        addpath(genpath('/Users/mai/Documents/mai_code/util'))
        addpath(genpath('/Users/mai/Documents/mai_code/spline_basis'))
	home = '/Users/mai/';
    otherwise
        display('unknown machine');
end
setup;
