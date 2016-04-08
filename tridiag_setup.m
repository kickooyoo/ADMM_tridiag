% mai_setup.m
% run in /mai_code
nfs_broken = true;
[~, machine] = system('hostname');

nfs_broken = ~exist('~/iv1h', 'dir') && ...
        isempty([strfind(lower(machine), 'iv1') strfind(lower(machine), 'vega')]);
if nfs_broken
	machine = 'iv1.';
end

switch machine(1:end-1) % last char is some sort of new line
    case {'eecs-IV1', 'iv1', sprintf('\niv1')}
        addpath('/n/ire/Volumes/s2/fessler/web/irt/irt')
	addpath('/n/ire/Volumes/s2/fessler/web/irt/irt/contrib/ramani/al-p2');
        addpath(genpath('~/Documents/mai_code/ADMM_tridiag'))
        addpath(genpath('~/Documents/mai_code/util'))
        addpath(genpath('~/Documents/mai_code/spline_basis'))
        addpath(genpath('~/Documents/mai_code/pthread_tutor'))
	addpath('~/Documents/contrib/ramani_MFISTA');
	addpath('~/Documents/contrib/ramani_fbrain');
	home_path = '~/';
	db_path = '~/Dropbox/fessler/experimental_data/tridiag/';
	if nfs_broken
		[~, machine] = system('hostname');
		if isempty(strfind(lower(machine), 'iv1'))
			db_path = '';
		end
		display('warning NFS broken, assuming iv1 paths');
	end
    case 'mpel8.eecs.umich.edu'
        addpath('/n/ire/Volumes/s2/fessler/web/irt/irt')
	addpath('/n/ire/Volumes/s2/fessler/web/irt/irt/contrib/ramani/al-p2');
        addpath(genpath('~/iv1h/Documents/mai_code/ADMM_tridiag'))
        addpath(genpath('~/iv1h//Documents/mai_code/util'))
        addpath(genpath('~/iv1h//Documents/mai_code/spline_basis'))
        addpath(genpath('~/iv1h//Documents/mai_code/pthread_tutor'))
	addpath('~/iv1h/Documents/contrib/ramani_MFISTA');
	addpath('~/iv1h/Documents/contrib/ramani_fbrain');
	home_path = '~/iv1h/';
	db_path = '~/iv1h/Dropbox/fessler/experimental_data/tridiag/';
    case {'ir63.eecs.umich.edu', 'ir72.eecs.umich.edu'}
        addpath('/n/ire/Volumes/s2/fessler/web/irt/irt')
	addpath('/n/ire/Volumes/s2/fessler/web/irt/irt/contrib/ramani/al-p2');
        addpath(genpath('~/iv1h/Documents/mai_code/util'))
        addpath(genpath('~/iv1h/Documents/mai_code/spline_basis'))
        addpath(genpath('~/iv1h/Documents/mai_code/pthread_tutor'))
	addpath('~/iv1h/Documents/contrib/ramani_MFISTA');
	addpath('~/iv1h/Documents/contrib/ramani_fbrain');
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
