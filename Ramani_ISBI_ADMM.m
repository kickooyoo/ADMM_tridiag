function [x, xsaved, cost] = Ramani_ISBI_ADMM(y,F,S,W,beta,mu,nu1,nu2,nx,ny,xinit,niters)

nc = size(S,1)/size(S,2);
% nc = S.odim/S.idim;
assert(mod(nc,1)==0,'nc not integer');

%x = zeros(nx*ny,1);
x = xinit(:);

iter = 0;

u0 = S*x;
u2 = x;
u1 = W*u2;
v0 = zeros(size(u0));
v1 = zeros(size(S*x));
v2 = zeros(size(u1));
v3 = zeros(size(W*u2));
v4 = zeros(size(u2));
v5 = zeros(size(x));
eta0 = zeros(size(v0));
eta1 = zeros(size(v1));
eta2 = zeros(size(v2));
eta3 = zeros(size(v3));
eta4 = zeros(size(v4));
eta5 = zeros(size(v5));


[eig_FF,Qbig] = get_eigs(F,nc);

%display('check inner dims of Q');
%keyboard;
SS = S'*S;
eig_SS = SS*ones(S.idim,1); % diagonal, not BCCB

while(iter < niters)
    
    % u0 update
    u0arg = F'*y + mu*(v0-eta0);
    u0 = Qbig'*((Qbig*u0arg)./(eig_FF + mu));

    % u1 update
    u1 = soft((v2-eta2)/sqrt(nu1), beta/(mu*nu1));
    
    % u2 update
    u2arg = sqrt(nu2)*(v4+eta4)-sqrt(nu1*W'*(v3+eta3));
    % if using undecimated wavelets (W'W circulant)
%    u2 =  Qbig'*((Qbig*u2arg)./(nu1*eig_WW + nu2));
    % if using orthogonal wavelets (W'W = I)
    u2 = u2arg./(nu1 + nu2); 

    % x update
    xarg = -sqrt(nu2)*(v5+eta5)-S'*(v1+eta1);
    x = xarg./(eig_SS+nu2);


	if(isnan(sum(x(:))))
		keyboard;
	end
    
    % v0 update
    %v0 = (mu(1)*(u0-eta0)+mu(2)*(S*x+eta1))/(mu(1)+mu(2));
    v0 = (mu*(u0-eta0)+mu*(S*x+eta1))/(2*mu);
	if(isnan(sum(v0(:))))
		keyboard;
	end

    % v1 update
    v1 = -v0;

    % v2 update
    %v2 = (mu(3)*(sqrt(nu1)*u1-eta2)+mu(4)*(sqrt(nu1)*W*u2+eta3))/(mu(3)+mu(4));
    v2 = (mu*(sqrt(nu1)*u1-eta2)+(sqrt(nu1)*W*u2+eta3))/(2*mu);

    % v3 update
    v3 = -v2;

    % v4 update
    %v4 = (mu(5)*(sqrt(nu2)*u2-eta4)+mu(6)*(sqrt(nu2)*x-eta5))/(mu(5)+mu(6));
    v4 = (mu*(sqrt(nu2)*u2-eta4)+mu*(sqrt(nu2)*x-eta5))/(2*mu);

    % v5 update
    v5 = -v4;

% 	eta updates % these are the correct ones
	eta0 = eta0 - (u0 - v0);
	eta1 = eta1 - (-S*x - v1);
	eta2 = eta2 - (sqrt(nu1)*u1 -v2);
	eta3 = eta3 - (-sqrt(nu1)*W*u2 - v3);
	eta4 = eta4 - (sqrt(nu2)*u2 - v4);
	eta5 = eta5 - (-sqrt(nu2)*x - v5);


    
    iter = iter+1;
    
    %u0saved(:,:,iter) = reshape(u0,nx,ny);
    %u1saved(:,:,iter) = reshape(u1,nx,ny);
    %u2saved(:,:,:,iter) = reshape(u2,nx,ny,nc);
    %u3saved(:,:,iter) = reshape(u3,nx,ny);
    %v3saved(:,:,iter) = reshape(v3,nx,ny);
    %xsaved(:,:,iter) = reshape(x,nx,ny);

%	AL_cost(iter) = calc_AL_cost(mu,alph,beta,C1,C2,F,S,y,u0,u1,u2,u3,x,v3,eta0,eta1,eta2,eta3,eta4);
%	cost(iter) = calc_cost(beta,W,F,S,y,x);



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

function cost = calc_cost(beta,W,F,S,y,x)
	cost = norm(y-F*(S*x),2)^2/2 + beta*norm(W*x,1);
end

function out = soft(in,thresh)
	out = (in - thresh*sign(in)).*(abs(in) > thresh);
end



