function [C1, C2] = construct_finite_diff(dims);

% assume dims is 1x2 vector
% what if for dynamic? Then need to add option so don't regularize in 3rd
% dim


% output is sparse matrices
%c0 = C2sparse('tight',true(dims(1),dims(2)),4,1);

%C1 = c0(1:dims(1)*dims(2),:);
%C2 = c0(dims(1)*dims(2)+1:end,:);

% FIND BETTER WAY TO MAKE FATRICES
% gives bad inner idim, odim
%C1test = Gmatrix(C1);
%C2test = Gmatrix(C2);

% 
%C1 = Cdiff(dims, 'offsets',[1],'edge_type','tight');
%C2 = Cdiff(dims, 'offsets',dims(1),'edge_type','tight');

R1 = Reg1(true(dims(1),dims(2)),'offsets',1,'type_diff','mex');
C1 = R1.C;


R2 = Reg1(true(dims(1),dims(2)),'offsets',dims(1),'type_diff','mex');
C2 = R2.C;
