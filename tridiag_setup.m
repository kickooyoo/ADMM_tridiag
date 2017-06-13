function machine = tridiag_setup(varargin)

if (nargin == 2) && strcmp(varargin{1}, 'irt_choice')
	irt_choice = varargin{2};
else
	irt_choice = 'dev';
end

addpath('../')
[id, machine, home_path] = mai_setup(0, 0, irt_choice);
if length(machine) > 4	
	machine = machine(1:4);
end

switch id
    case 'iv1'
	%addpath('/n/ire/Volumes/s2/fessler/web/irt/irt/contrib/ramani/al-p2');
	db_path = '~/Dropbox/fessler/experimental_data/tridiag/';
    case 'mpel8'
	db_path = '~/iv1h/Dropbox/fessler/experimental_data/tridiag/';
    case {'ir72'; 'ir'}
	db_path = '~/iv1h/Dropbox/fessler/experimental_data/tridiag/';
    case 'vega'
        addpath(genpath('/Users/mai/Documents/data/2010-07-06-fessler-3d/code'))
	addpath('/Users/mai/Documents/irt/reproduce/ramani-12-asb/2011-12-19');
	db_path = '~/Dropbox/fessler/experimental_data/tridiag_vega/';
    otherwise
        display('unknown machine');
end


addpath(genpath(sprintf('%s/Documents/mai_code/spline_basis', home_path)))
addpath(genpath(sprintf('%s/Documents/mai_code/pthread_tutor', home_path)))
addpath(genpath(sprintf('%s/Documents/contrib/ramani_MFISTA', home_path)))
addpath(genpath(sprintf('%s/Documents/contrib/ramani_fbrain', home_path)))
addpath(genpath(sprintf('%s/Documents/mai_code/ADMM_tridiag', home_path)))
