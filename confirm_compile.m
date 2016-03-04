function confirm_compile(file)
% function confirm_compile(varargin)
% recompiles mex file


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
mex('-O', '-DMmex', '-I./pthread_tutor/def/', sprintf('./pthread_tutor/%s.c', file))
