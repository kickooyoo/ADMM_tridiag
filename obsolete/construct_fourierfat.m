function F = construct_fourierfat(samp,ncoils)

% assumer that sampling pattern is 2D, can change to 3D later
assert(size(samp,3)==1,'must have 2D sampling pattern');

nx = size(samp,1);
ny = size(samp,2);
A = Gdft('mask',true(nx,ny),'samp',samp);
Acells = repmat({A},ncoils,1);
F = block_fatrix(Acells);
%Ftest = block_fatrix({Gdft('mask',true(nx,ny),'samp',samp)},'type','kron','Mkron',nc);

% result is BCCB with <ncoils> blocks
