% mai_setup.m
% run in /mai_code

[~, machine] = system('hostname');

switch machine(1:end-1) % last char is some sort of new line
    case {'eecs-IV1', 'iv1'}
        addpath('/n/ire/Volumes/s2/fessler/web/irt/irt')
	addpath('/n/ire/Volumes/s2/fessler/web/irt/irt/contrib/ramani/al-p2');
        addpath(genpath('~/Documents/mai_code/ADMM_tridiag'))
        addpath(genpath('~/Documents/mai_code/util'))
        addpath(genpath('~/Documents/mai_code/spline_basis'))
        addpath(genpath('~/Documents/mai_code/pthread_tutor'))
        addpath(genpath('~/Documents/data/2010-07-06-fessler-3d/code'))
	addpath('~/Documents/contrib/ramani_MFISTA');
	home_path = '~/';
	db_path = '~/Dropbox/fessler/experimental_data/tridiag/';
    case 'mpel8.eecs.umich.edu'
        addpath('/n/ire/Volumes/s2/fessler/web/irt/irt')
	addpath('/n/ire/Volumes/s2/fessler/web/irt/irt/contrib/ramani/al-p2');
        addpath(genpath('~/iv1h/Documents/mai_code/ADMM_tridiag'))
        addpath(genpath('~/iv1h//Documents/mai_code/util'))
        addpath(genpath('~/iv1h//Documents/mai_code/spline_basis'))
        addpath(genpath('~/iv1h//Documents/mai_code/pthread_tutor'))
        addpath(genpath('~/iv1h//Documents/data/2010-07-06-fessler-3d/code'))
	addpath('~/iv1h/Documents/contrib/ramani_MFISTA');
	home_path = '~/iv1h/';
	db_path = '~/iv1h/Dropbox/fessler/experimental_data/tridiag/';
    case {'ir63.eecs.umich.edu', 'ir72.eecs.umich.edu'}
        addpath('/n/ire/Volumes/s2/fessler/web/irt/irt')
	addpath('/n/ire/Volumes/s2/fessler/web/irt/irt/contrib/ramani/al-p2');
        addpath(genpath('~/iv1h/Documents/mai_code/util'))
        addpath(genpath('~/iv1h/Documents/mai_code/spline_basis'))
        addpath(genpath('~/iv1h/Documents/mai_code/pthread_tutor'))
        addpath(genpath('~/iv1h/Documents/data/2010-07-06-fessler-3d/code'))
	addpath('~/iv1h/Documents/contrib/ramani_MFISTA');
	home_path = '~/iv1h/';
	db_path = '~/iv1h/Dropbox/fessler/experimental_data/tridiag/';
    case 'vega'
        addpath('/Users/mai/Documents/irt')
	addpath('/Users/mai/Documents/irt/contrib/ramani/al-p2');
        addpath(genpath('/Users/mai/Documents/mai_code'))
        addpath(genpath('/Users/mai/Documents/mai_code/ADMM_tridiag'))
        addpath(genpath('/Users/mai/Documents/contrib'))
        addpath(genpath('/Users/mai/Documents/mai_code/util'))
        addpath(genpath('/Users/mai/Documents/mai_code/spline_basis'))
        addpath(genpath('/Users/mai/Documents/data/2010-07-06-fessler-3d/code'))
	addpath('/Users/mai/Documents/irt/reproduce/ramani-12-asb/2011-12-19');
	addpath('/Users/mai/Documents/contrib/ramani_MFISTA');
	home_path = '/Users/mai/';
	db_path = '~/Dropbox/fessler/experimental_data/tridiag_vega/';
    otherwise
        display('unknown machine');
end
if ~exist('im','file')
	setup;
end
