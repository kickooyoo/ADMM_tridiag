function [x, xsaved, err, cost, time] = SENSE_ADMM_ALP2_wavelets(y,F,S,C1,C2,alph,beta,mu,nx,ny,xinit,xtrue,niters)

nc = size(S,1)/size(S,2);
% nc = S.odim/S.idim;
assert(mod(nc,1)==0,'nc not integer');

%x = zeros(nx*ny,1);
x = xinit(:);

iter = 0;


alph2 = 0.5;
beta2 = beta;

% construct wavelets
%W = Gwave1('mask',true(nx,ny),'wname','haar'); % how do i make sure it is undecimated?
W = Godwt1(true(nx,ny));
CHW = [C1; beta2*alph2/beta*W];
CVW = [C2; beta2*alph2/beta* W];

u0 = CHW*x;
u3 = x;
u1 = CVW*u3;
u2 = (1-alph)*S*u3+alph*S*x;
v3 = zeros(size(u3)); % size of u3 and eta3, different initializations?
v4 = -v3;
eta0 = zeros(size(u0));
eta1 = zeros(size(u1));
eta2 = zeros(size(u2));
eta3 = zeros(size(u3));
eta4 = zeros(size(x));


[eig_FF,Qbig] = get_eigs(F,nc);

%display('check inner dims of Q');
%keyboard;
SS = S'*S;
eig_SS = SS*ones(S.idim,1); % diagonal, not BCCB
time = zeros(niters,1);
mask = generate_mask('slice67',1,nx,ny);
err(1) = calc_NRMSE_over_mask(x,xtrue(:),mask);
cost(1) = calc_cost(beta,beta2,C1,C2,F,S,W,y,x);
while(iter < niters)
    
	iter_start = tic;
    % u0 update
    u0 = soft(CHW*double(x)-eta0,beta/mu(1)); %mu0

    % u1 update
    u1 = soft(CVW*double(u3)-eta1,beta/mu(2));
    
    % u2 update
    u2 = u2_update(mu(3),alph,eig_FF,Qbig,F,S,y,u3,x,eta2,nx,ny,nc);

    % u3 update
    [u3, u3_extra_time] = u3_update(mu,alph,alph2,beta,beta2,eig_SS,C1,CVW,S,u1,u2,x,v3,eta1,eta2,eta3,nx,ny);

    % x update
    [x, x_extra_time] = x_update(mu,alph,alph2,beta,beta2,eig_SS,CHW,S,u0,u2,u3,v4,eta0,eta2,eta4,nx,ny);
    
    % skip v0, v1, v2 because they are constrained to be zero
    
    % v3 update
    v3 = (mu(4)*(-u3-eta3)+mu(5)*(-x+eta4))./(mu(4)+mu(5));
    
    % v4 update
    v4 = -v3;
    
% 	eta updates % these are the correct ones
	eta0 = eta0 - (-u0 + CHW*x);
	eta1 = eta1 - (-u1 + CVW*u3);
	eta2 = eta2 - (-u2 + (1-alph)*S*u3 + alph*S*x);
	eta3 = eta3 - (-u3 - v3);
	eta4 = eta4 - (x - v4);


    % these are the experimental ones
% 	eta0 = eta0 - (-u0 + C1*x);
% 	eta1 = eta1 - (-u1 + C2*u3);
% 	eta2 = eta2 - (-u2 + (1-alph)*S*u3 + alph*S*x);
% 	eta3 = eta3 - (-u3 - v3);
% 	eta4 = eta4 + (x - v4);

    
    iter = iter+1;
	iter_time = toc(iter_start);
	time(iter) = iter_time - u3_extra_time - x_extra_time;
	err(iter) = calc_NRMSE_over_mask(x,xtrue(:),mask);
    
    if mod(iter,10) == 0
        printf('%d/%d iterations',iter,niters)
    end
%    u0saved(:,:,iter) = reshape(u0,nx,ny);
%    u1saved(:,:,iter) = reshape(u1,nx,ny);
%    u2saved(:,:,:,iter) = reshape(u2,nx,ny,nc);
%    u3saved(:,:,iter) = reshape(u3,nx,ny);
%    v3saved(:,:,iter) = reshape(v3,nx,ny);
    xsaved(:,:,iter) = reshape(x,nx,ny);

%	AL_cost(iter) = calc_AL_cost(mu,alph,beta,C1,C2,F,S,y,u0,u1,u2,u3,x,v3,eta0,eta1,eta2,eta3,eta4);
	cost(iter+1) = calc_cost(beta,beta2,C1,C2,F,S,W,y,x);


end
%figure; %subplot(1,3,1); 
%im(xsaved); colorbar;
%title('xsaved');
%subplot(1,3,2); 
%figure;
%im(u0saved); colorbar; title('u0saved');
%subplot(1,3,3); 
%figure;
%im(u1saved); colorbar; title('u1saved');
%figure; im(reshape(u2,nx,ny,nc));
%figure; im(reshape(S*x,nx,ny,nc));
%display('STOP');
%keyboard;


end

function AL_cost = calc_AL_cost(mu,alph,beta,C1,C2,F,S,y,u0,u1,u2,u3,x,v3,eta0,eta1,eta2,eta3,eta4)
	AL_cost = norm(y-F*u2,2)^2/2 + beta*norm(u0,1)+beta*norm(u1,1)+mu(1)*norm(-u0+C1*x-eta0,2)^2/2+mu(2)*norm(-u1+C2*u3-eta1,2)^2/2+mu(3)*norm(-u2+(1-alph)*S*u3+alph*S*x-eta2,2)^2/2 +mu(4)*norm(-u3-v3-eta3,2)^2/2+mu(5)*norm(x+v3-eta4,2)^2/2;
end

function cost = calc_cost(beta,beta2,C1,C2,F,S,W,y,x)
	cost = norm(y-F*(S*x),2)^2/2 + beta*norm(C1*x,1) + beta*norm(C2*x,1)+beta2*norm(W*x,1);
end

function out = soft(in,thresh)
	out = (in - thresh*sign(in)).*(abs(in) > thresh);
end

function u2 = u2_update(mu2,alph,eig_FF,Q,F,S,y,u3,x,eta2,nx,ny,nc)
	arg_u2 = F'*y+mu2*((1-alph)*S*u3+alph*S*x-eta2);
	new_u2 = Q'*((Q*arg_u2)./(eig_FF+mu2))/(nx*ny);
	

%	A = full(F'*F);
%	A = A+mu2*eye(length(A));
%	test = A\arg_u2;
%
%	display('test u2');
%	figure; plot(abs(new_u2),'k'); 
%	hold on; plot(abs(test),'r');
%	title('u2 update');
%	legend('my soln','backslash soln');
%	keyboard; % passed test! :)
	u2 = new_u2;
end

function [u3, time_subtract] = u3_update(mu,alph,alph2,beta1,beta2,eig_SS,C1,CVW,S,u1,u2,x,v3,eta1,eta2,eta3,nx,ny)
	u3 = mu(2)*CVW'*(u1+eta1)+mu(3)*(1-alph)*S'*(u2-alph*S*x+eta2)+mu(4)*(-v3-eta3);

	% NOTE: C1 CURRENTLY NOT USED!

	% take argument of C2'C2 
	% reshape to rect
	% transpose,  stretch, apply C1'C1
	flipu3 = reshape(u3,nx,ny);
	flipu3 = flipu3.';
	flipu3 = flipu3(:);
	flipSS = reshape(eig_SS,nx,ny);
	flipSS = flipSS.';
	flip_eig_SS = flipSS(:);

	% construct diagonal entries
	subdiag = -mu(2)*ones(nx*ny-1,1);
	subdiag(ny:ny:end) = 0;

	supdiag = subdiag;

	diagvals = 2*ones(nx*ny,1);
	diagvals(1:ny:end) = 1;
	diagvals(ny:ny:end) = 1;
	diagvals = mu(2)*diagvals+mu(3)*(1-alph)^2.*flip_eig_SS+(mu(4)+mu(2)*(beta2/beta1)^2*(1-alph2)^2);
   

	%u3out = apply_tridiag_inv(subdiag,diagvals,supdiag,flipu3);
	
	u3out = zeros(nx*ny,1);
	tridiag_start = tic;
	for block_ndx=1:nx
		u3out((block_ndx-1)*ny+1:block_ndx*ny) = apply_tridiag_inv(subdiag((block_ndx-1)*ny+1:block_ndx*ny-1),diagvals((block_ndx-1)*ny+1:block_ndx*ny),supdiag((block_ndx-1)*ny+1:block_ndx*ny-1),flipu3((block_ndx-1)*ny+1:block_ndx*ny));
	end
	tridiag_time = toc(tridiag_start);
	tridiag_time_per_block = tridiag_time/nx;
	time_subtract = ceil(nx*7/8)*tridiag_time_per_block; %numcores = 4

	% FLIP SOLUTION BACK
	flipu3 = reshape(u3out,ny,nx); % NOT nx ny
	flipu3 = flipu3.';
	u3final = flipu3(:);
 	
	u3 = u3final;
%	u3 = test;
end

function [x, time_subtract] = x_update(mu,alph,alph2,beta1,beta2,eig_SS,CHW,S,u0,u2,u3,v4,eta0,eta2,eta4,nx,ny)
	x = mu(1)*CHW'*(u0+eta0)+mu(3)*alph*S'*(u2-(1-alph)*S*u3+eta2)+mu(5)*(v4+eta4);

	subdiag = -1*ones(nx*ny-1,1);
	subdiag(nx:nx:end) = 0;
	
	supdiag = subdiag;
	
	diagvals = 2*ones(nx*ny,1);
	diagvals(1:nx:end) = 1;
	diagvals(nx:nx:end) = 1;
	diagvals = mu(1)*diagvals+mu(3)*alph^2.*eig_SS+(mu(5)+mu(1)*(beta2/beta1)^2*alph2^2);
	
	%xout = apply_tridiag_inv(subdiag,diagvals,supdiag,x);
	xout = zeros(nx*ny,1);
	tridiag_start = tic;
	for block_ndx = 1:ny
		xout((block_ndx-1)*nx+1:block_ndx*nx) = apply_tridiag_inv(subdiag((block_ndx-1)*nx+1:block_ndx*nx-1),diagvals((block_ndx-1)*nx+1:block_ndx*nx),supdiag((block_ndx-1)*nx+1:block_ndx*nx-1),x((block_ndx-1)*nx+1:block_ndx*nx));
	end
	tridiag_time = toc(tridiag_start);
	tridiag_time_per_block = tridiag_time/nx;
	time_subtract = ceil(ny*7/8)*tridiag_time_per_block; %numcores = 4
	
	x = xout;
	%x = A*x;
end





