function confirm_compile(file, varargin)
% function confirm_compile(file)
% 
% recompiles mex file 

arg.fpath = './pthread_tutor/';
arg = vararg_pair(arg, varargin);


if ~exist(sprintf('./pthread_tutor/%s.c', file))
        display(sprintf('cannot find ./pthread_tutor/%s.c', file))
        keyboard;
end

% to do: make sure mex compiled
curr = cd;
if ~strcmp(curr(end-11:end), 'ADMM_tridiag')
        display(sprintf('must be in ADMM_tridiag dir to compile %s.c', file));
        keyboard;
end
% mex -O CFLAGS="\$CFLAGS -std=c99 -DMmex" -I./pthread_tutor/def/ ./pthread_tutor/tridiag_inv_mex_noni.c
mex('-O', '-DMmex', sprintf('-I%s/def/', fpath), sprintf('%s/%s.c', fpath, file))
