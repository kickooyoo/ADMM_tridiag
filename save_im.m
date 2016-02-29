function save_im(path, name, varargin)
% function save_im(path, name, varargin)
arg.FontSize = [];
arg = vararg_pair(arg, varargin);

% if ~isempty(arg.FontSize)
	

print([path name], '-depsc');
savefig([path name])


