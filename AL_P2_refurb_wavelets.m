function [x,xsave,err,costOrig,time] = AL_P2_refurb_wavelets(ylong,A,S,R,initx,niters,lambda,mu,nx,ny,xtrue)

npix = nx*ny;
nc = length(A.arg.blocks); 

y = ylong(:);
x = initx(:);

% soft thresholding function, t is input, a is threshold
soft = @(t,a) (t - a * sign(t)) .* (abs(t) > a);
eigvals1 = @(A,Q) Q*A(:,1);  

% fatrix objects
Q = Gdft('mask', true(nx,ny));

% R: regularizer, currently 
%R = [C1; C2];
SS = S'*S;
eigvalsss = SS*ones(nx*ny,1);
RR = R'*R;
eigvalsrr = Q*RR(:,1); % NOT CIRCULANT ANYMO'
[eigvalsaa,Qbig] = get_eigs(A,nc); % match, they are okay!
%eigvalsaa = repmat(Q*AA(:,1),nc,1);

W = Godwt1(true(nx,ny));

% hand picked mask for computing NRMSE over mask
mask = generate_mask('slice67',1,nx,ny);


u = S*x;
v = [R*x; W*x]; %v = [R; W]z
z = x;

u_u = mu(1);
u_v = mu(2);
u_z = mu(3);
eta_u = zeros(size(u));
eta_v = zeros(size(v));
eta_z = zeros(size(z));

calc_errcost = 1;
if (calc_errcost)
    calc_orig_cost = @(y,A,S,R,W,x,lambda) norm(y-A*(S*x),2)^2/2 + ...
        lambda*norm(R*x,1) + lambda*norm(W*x,1);
    err = calc_NRMSE_over_mask(x,xtrue(:),mask);
    costOrig = calc_orig_cost(y,A,S,R,W,x,lambda);
end

xsave = zeros(nx,ny,niters);
% [k_u,k_v,k_x] = diag_cond_numbers(eigvalsmssm,eigvalsrr,npix,Q,R,M,S,A,u_u,u_v,u_z);
time = zeros(niters,1);
%tic
for ii=1:niters
	iter_start = tic;
    x = x_update(S,u,z,eta_u,eta_z,u_u,u_z,eigvalsss);
    u = u_update(Qbig,A,S,y,x,eta_u,u_u,npix,eigvalsaa);
    %v = soft(R*x+eta_v,lambda/u_v);
    v = soft([R*z; W*z] + eta_v,lambda/u_v); % TEST THIS
    z = z_update(Q,v,x,z,u_v,u_z,R,W,eta_v,eta_z,eigvalsrr,nx,ny);
    eta_u = eta_u - (u - S*x);
    eta_v = eta_v - (v - [R*z; W*z]);
    eta_z = eta_z - (z - x);
    time(ii) = toc(iter_start);
    if (calc_errcost)
    	err = [err calc_NRMSE_over_mask(x,xtrue(:),mask)];
        costOrig = [costOrig calc_orig_cost(y,A,S,R,W,x,lambda)];
    end
    if mod(ii,10) == 0
        printf('%d/%d iterations',ii,niters)
    end
    xsave(:,:,ii) = reshape(x,nx,ny);
end
%toc
if (~calc_errcost)
   err = zeros(niters+1,1);
   costOrig = zeros(niters+1,1);nor
end
end

function x = x_update(S,u,z,eta_u,eta_z,u_u,u_z,eigvalsss)
inv_vals = (u_u*eigvalsss+u_z);
x = (u_u*S'*(u-eta_u)+u_z*(z-eta_z))./inv_vals;
end

function u = u_update(Qbig,A,S,y,x,eta_u,u_u,N,eigvalsaa)
invMat = eigvalsaa + u_u;
latter = (A'*y+u_u*(S*x+eta_u));
u = Qbig'*((Qbig*latter)./invMat)/N;
end

% Conjugate Gradient Descent!
% z = argmin (u_v R'R+u_z I)^(-1)(u_v R'(v-etav)+u_z (x+etaz))
% z = argmin u_?/2 ||v-Rz-etav||^2+u_?/2 ||z-x-etaz||^2
function z = z_update(Q,v,x,z,u_v,u_z,R,W,eta_v,eta_z,eigvalsrr,nx,ny)
circ = 1;

if (~circ) % need to update for wavelet version!
	n1 = length(v); 
	n2 = length(x);
	W = Gdiag([(u_v/2)*ones(1,n1) (u_z/2)*ones(1,n2)]);
	A = [Gmatrix(R); Gmatrix(Gdiag(-1*ones(1,n2)))]; % nightmare, pass in 
	y = [v-eta_v;-(x+eta_z)];
	[z_pcg, info] = qpwls_pcg1(z,A,W,y,Gdiag(zeros(n2,1)),'niter',120,'stop_grad_tol',1e-11,'precon',A'*A);
	z = z_pcg;
else 
	rhs = reshape(u_v*[R; W]'*(v-eta_v)+u_z*(x+eta_z),nx,ny);
	%z_fft = ifft2((fft2(rhs))./(reshape(u_v*eigvalsrr+u_z*ones(nx*ny,1),nx,ny)));
	%z_fft = z_fft(:);
	invMat = u_v*eigvalsrr+u_v+u_z;
	z_fft2 = Q'*((Q*rhs(:))./invMat)/(nx*ny);
	%keyboard;
	z = z_fft2;
end

end

