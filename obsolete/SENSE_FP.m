function [x, xsaved, cost] = SENSE_FP(y,F,S,C1,C2,alph,beta,mu,nx,ny,xinit,niters,attempt_par)

	nc = size(S,1)/size(S,2);
	% nc = S.odim/S.idim;
	assert(mod(nc,1)==0,'nc not integer');

	x = zeros(nx*ny,1);
	%x = xinit(:);


	u0 = C1*x;
	u3 = x;
	u1 = C2*u3;
	u2 = (1-alph)*S*u3+alph*S*x;
	v0 = zeros(size(u0));
	v1 = v0;
	v2 = zeros(size(u1));
	v3 = v2;
	v4 = zeros(size(u2));
	v5 = v4;
	v6 = v4;
	v7 = zeros(size(u3));
	v8 = v7;
	eta0 = zeros(size(v0));
	eta1 = zeros(size(v1));
	eta2 = zeros(size(v2));
	eta3 = zeros(size(v3));
	eta4 = zeros(size(v4));
	eta5 = zeros(size(v5));
	eta6 = zeros(size(v6));
	eta7 = zeros(size(v7));
	eta8 = zeros(size(v8));


	[eig_FF,Qbig] = get_eigs(F,nc);

	%display('check inner dims of Q');
	%keyboard;
	SS = S'*S;
	eig_SS = SS*ones(S.idim,1); % diagonal, not BCCB

	%u_update_list = {@u0_update, @u1_update, @u2_update, @u3_update, @x_update};
	%u_in_list = 
	%iter = 0;
	%while(iter < niters)

	attempt_tridiag_par = 0;

	if(attempt_par || attempt_tridiag_par)
		%if matlabpool('size') == 0
			%matlabpool('open',ncore);
	%		pool = parpool(ncore); % only in 2013b
%		en
		ncore = 4;
		pool = parpool(ncore); % only in 2013b
	end

	pmethod = 'feval'; %'pfor'; % or 'spmd' or 'feval'
	profile on;	
	tic
	%matlabpool open 
	for iter=1:niters 
		if (attempt_par)
			if strcmp(pmethod,'pfor')
				parfor ui = 0:4
					switch ui
					case 0 
					%u0 = soft(-v0-eta0,beta/mu(1)); %mu0
						u0 = u0_update(v0,eta0,beta,mu(1));
					case 1
						%u1 = soft(-v2-eta2,beta/mu(3));
						u1 = u1_update(v2,eta2,beta,mu(3));
					case 2
						u2 = u2_update(mu(5),eig_FF,Qbig,F,y,v4,eta4,nx,ny,nc);
					case 3
						u3 = u3_update(mu,alph,eig_SS,C1,C2,S,v3,v5,v7,eta3,eta5,eta7,nx,ny);
					case 4
						x = x_update(mu,alph,eig_SS,C1,S,v1,v6,v8,eta1,eta6,eta8,nx,ny);
					otherwise
						display('no such case for ui loop');
					end
				end
			elseif strcmp(pmethod,'spmd'); 
				spmd;
					if labindex == 1	
						u0 = u0_update(v0,eta0,beta,mu(1));
					elseif labindex == 2
						u1 = u1_update(v2,eta2,beta,mu(3));
					elseif labindex == 3
						u2 = u2_update(mu(5),eig_FF,Qbig,F,y,v4,eta4,nx,ny,nc);
					elseif labindex == 4
						u3 = u3_update(mu,alph,eig_SS,C1,C2,S,v3,v5,v7,eta3,eta5,eta7,nx,ny);
					elseif labindex == 5
						x = x_update(mu,alph,eig_SS,C1,S,v1,v6,v8,eta1,eta6,eta8,nx,ny);
				end
			end
			% extract vectors from composites
			u0 = u0{1};	
			u1 = u1{2};
			u2 = u2{3};
			u3 = u3{4};
			% x
			elseif strcmp(pmethod,'feval')
				u0 = parfeval(pool,@u0_update,1,v0,eta0,beta,mu(1));
				u1 = parfeval(pool,@u1_update,1,v2,eta2,beta,mu(3)); 
				u2 = parfeval(pool,@u2_update,1,mu(5),eig_FF,Qbig,F,y,v4,eta4,nx,ny,nc);
				u3 = parfeval(pool,@u3_update,1,mu,alph,eig_SS,C1,C2,S,v3,v5,v7,eta3,eta5,eta7,nx,ny);
				x = parfeval(pool,@x_update,1,mu,alph,eig_SS,C1,S,v1,v6,v8,eta1,eta6,eta8,nx,ny);
				u0 = fetchOutputs(u0);
				u1 = fetchOutputs(u1);
				u2 = fetchOutputs(u2);
				u3 = fetchOutputs(u3);
				x = fetchOutputs(x);
			end
		else %not parallel
			u0 = u0_update(v0,eta0,beta,mu(1));
			u1 = u1_update(v2,eta2,beta,mu(3));
			u2 = u2_update(mu(5),eig_FF,Qbig,F,y,v4,eta4,nx,ny,nc);
			if (attempt_tridiag_par)
				u3 = u3_update(mu,alph,eig_SS,C1,C2,S,v3,v5,v7,eta3,eta5,eta7,nx,ny,pool);
				x = x_update(mu,alph,eig_SS,C1,S,v1,v6,v8,eta1,eta6,eta8,nx,ny,pool);
			else
				u3 = u3_update(mu,alph,eig_SS,C1,C2,S,v3,v5,v7,eta3,eta5,eta7,nx,ny,0);
				x = x_update(mu,alph,eig_SS,C1,S,v1,v6,v8,eta1,eta6,eta8,nx,ny,0);
				
			end
		end
		if (attempt_par)
			if strcmp(pmethod,'pfor')
				parfor vi = 0:4
					switch vi
					case 0
						v0 = (mu(1)*(-u0-eta0)+mu(2)*(-C1*x+eta1))/(mu(1)+mu(2));
					case 1
						v2 = (mu(3)*(-u1-eta2)+mu(4)*(-C2*u3+eta3))/(mu(3)+mu(4));
					case 2
						v4 = v4_update(u2,u3,x,alph,S,eta4,eta5,eta6,mu); 
					case 3
						v5 = v5_update(u2,u3,x,alph,S,eta4,eta5,eta6,mu); 
					case 4
						v7 = (mu(8)*(-u3-eta7)+mu(9)*(-x+eta8))/(mu(8)+mu(9));
					otherwise
					end
				end
			elseif strcmp(pmethod,'spmd') 
				spmd;
					if labindex == 1	
						v0 = (mu(1)*(-u0-eta0)+mu(2)*(-C1*x+eta1))/(mu(1)+mu(2));
					elseif labindex == 2	
						v2 = (mu(3)*(-u1-eta2)+mu(4)*(-C2*u3+eta3))/(mu(3)+mu(4));
					elseif labindex == 3	
						v4 = v4_update(u2,u3,x,alph,S,eta4,eta5,eta6,mu); 
					elseif labindex == 4	
						v5 = v5_update(u2,u3,x,alph,S,eta4,eta5,eta6,mu); 
					elseif labindex == 5	
						v7 = (mu(8)*(-u3-eta7)+mu(9)*(-x+eta8))/(mu(8)+mu(9));
					end
				end
		
				% extract vectors from composite
				v0 = v0{1};
				v2 = v2{2};
				v4 = v4{3};
				v5 = v5{4};
				% v7
			elseif strcmp(pmethod,'feval')
				v0 = parfeval(pool,@v0_update,1,mu(1),mu(2),u0,eta0,C1,x,eta1);
				v2 = parfeval(pool,@v2_update,1,mu(3),mu(4),u1,eta2,C2,u3,eta3);
				v4 = parfeval(pool,@v4_update,1,u2,u3,x,alph,S,eta4,eta5,eta6,mu); 
				v5 = parfeval(pool,@v5_update,1,u2,u3,x,alph,S,eta4,eta5,eta6,mu); 
				v7 = parfeval(pool,@v7_update,1,mu(8),mu(9),u3,eta7,x,eta8);
				v0 = fetchOutputs(v0);
				v2 = fetchOutputs(v2);
				v4 = fetchOutputs(v4);
				v5 = fetchOutputs(v5);
				v7 = fetchOutputs(v7);

			end
		else % not parallel
			v0 = (mu(1)*(-u0-eta0)+mu(2)*(-C1*x+eta1))/(mu(1)+mu(2));
			v2 = (mu(3)*(-u1-eta2)+mu(4)*(-C2*u3+eta3))/(mu(3)+mu(4));
			v4 = v4_update(u2,u3,x,alph,S,eta4,eta5,eta6,mu,v5); 
			v5 = v5_update(u2,u3,x,alph,S,eta4,eta5,eta6,mu,v4); 
			v7 = (mu(8)*(-u3-eta7)+mu(9)*(-x+eta8))/(mu(8)+mu(9));
		end
		v1 = -v0; 
		v3 = -v2;
		v6 = -v4 - v5;
		v8 = -v7;

	% 	eta updates 
	%	eta  = eta - (Au-v)
		eta0 = eta0 - (-u0 - v0);
		eta1 = eta1 - (C1*x - v1);
		eta2 = eta2 - (-u1 - v2);
		eta3 = eta3 - (C2*u3 - v3);
		eta4 = eta4 - (-u2 - v4);
		eta5 = eta5 - ((1-alph)*S*u3 - v5); 
		eta6 = eta6 - (alph*S*x - v6);
		eta7 = eta7 - (-u3 - v7);
		eta8 = eta8 - (x - v8);
	    
		%iter = iter+1; % for while loop
	    
		u0saved(:,:,iter) = reshape(u0,nx,ny);
		u1saved(:,:,iter) = reshape(u1,nx,ny);
		u2saved(:,:,:,iter) = reshape(u2,nx,ny,nc);
		u3saved(:,:,iter) = reshape(u3,nx,ny);
		v3saved(:,:,iter) = reshape(v3,nx,ny);
		xsaved(:,:,iter) = reshape(x,nx,ny);

		%AL_cost(iter) = calc_AL_cost(mu,alph,beta,C1,C2,F,S,y,u0,u1,u2,u3,x,v3,eta0,eta1,eta2,eta3,eta4);
		cost(iter) = calc_cost(beta,C1,C2,F,S,y,x);

	end
	toc
	profile viewer
	keyboard;

	if(attempt_par || attempt_tridiag_par)
		matlabpool close
	end

end

function AL_cost = calc_AL_cost(mu,alph,beta,C1,C2,F,S,y,u0,u1,u2,u3,x,v3,eta0,eta1,eta2,eta3,eta4)
	AL_cost = norm(y-F*u2,2)^2/2 + beta*norm(u0,1)+beta*norm(u1,1)+mu(1)*norm(-u0+C1*x-eta0,2)^2/2+mu(2)*norm(-u1+C2*u3-eta1,2)^2/2+mu(3)*norm(-u2+(1-alph)*S*u3+alph*S*x-eta2,2)^2/2 +mu(4)*norm(-u3-v3-eta3,2)^2/2+mu(5)*norm(x+v3-eta4,2)^2/2;
end

function cost = calc_cost(beta,C1,C2,F,S,y,x)
	cost = norm(y-F*(S*x),2)^2/2 + beta*norm(C1*x,1) + beta*norm(C2*x,1);
end

function out = soft(in,thresh)
	out = (in - thresh*sign(in)).*(abs(in) > thresh);
end

function u0 = u0_update(v0,eta0,beta,mu);
	u0 = soft(-v0-eta0,beta/mu);
end

function u1 = u1_update(v2,eta2,beta,mu)
	u1 = soft(-v2-eta2,beta/mu);
end

function u2 = u2_update(mu4,eig_FF,Q,F,y,v4,eta4,nx,ny,nc);
	arg_u2 = F'*y+mu4*(-v4-eta4);
	new_u2 = Q'*((Q*arg_u2)./(eig_FF+mu4))/(nx*ny);
	u2 = new_u2;
end

function u3 = u3_update(mu,alph,eig_SS,C1,C2,S,v3,v5,v7,eta3,eta5,eta7,nx,ny,pool);
	arg = mu(4)*C2'*(v3+eta3)+mu(6)*(1-alph)*S'*(v5+eta5)+mu(8)*(-v7-eta7);

	% take argument of C2'C2, reshape to rect, transpose,  stretch, apply C1'C1
	flipu3 = reshape(arg,nx,ny);
	flipu3 = flipu3.';
	flipu3 = flipu3(:);
	flipSS = reshape(eig_SS,nx,ny);
	flipSS = flipSS.';
	flip_eig_SS = flipSS(:);

	% construct diagonal entries
	subdiag = -mu(4)*ones(nx*ny-1,1);
	subdiag(ny:ny:end) = 0;

	supdiag = subdiag;

	diagvals = 2*ones(nx*ny,1);
	diagvals(1:ny:end) = 1;
	diagvals(ny:ny:end) = 1;
	diagvals = mu(4)*diagvals+mu(6)*(1-alph)^2.*flip_eig_SS+mu(8);
		
	if pool == 0
		u3out = apply_tridiag_inv(subdiag,diagvals,supdiag,flipu3);

	else
	       	%par_u3 = zeros(size(flipu3));
		par_u3 = cell(nx,1);
		parfor ii = 1:nx
			par_u3{ii} = apply_tridiag_inv(subdiag((ii-1)*ny+1:ii*ny-1),diagvals((ii-1)*ny+1:ii*ny),supdiag((ii-1)*ny+1:ii*ny-1),flipu3((ii-1)*ny+1:ii*ny));
		end
		u3out = cat(1,par_u3{:});
	end

	% FLIP SOLUTION BACK
	flipu3 = reshape(u3out,ny,nx); % NOT nx ny
	flipu3 = flipu3.';
	u3final = flipu3(:);
	
 	u3 = u3final;
end

function x = x_update(mu,alph,eig_SS,C1,S,v1,v6,v8,eta1,eta6,eta8,nx,ny,pool);
	x = mu(2)*C1'*(v1+eta1)+mu(7)*alph*S'*(v6+eta6)+mu(9)*(v8+eta8);

	subdiag = -1*ones(nx*ny-1,1);
	subdiag(nx:nx:end) = 0;
	
	supdiag = subdiag;
	
	diagvals = 2*ones(nx*ny,1);
	diagvals(1:nx:end) = 1;
	diagvals(nx:nx:end) = 1;
	diagvals = mu(2)*diagvals+mu(7)*alph^2.*eig_SS+mu(9);

	if (pool == 0)	
		xout = apply_tridiag_inv(subdiag,diagvals,supdiag,x);
	else
		%par_x = zeros(size(x));
		par_x = cell(ny,1);
		parfor ii = 1:ny
			par_x{ii} = apply_tridiag_inv(subdiag((ii-1)*nx+1:ii*nx-1),diagvals((ii-1)*nx+1:ii*nx),supdiag((ii-1)*nx+1:ii*nx-1),x((ii-1)*nx+1:ii*nx));
		end
		xout = cat(1,par_x{:});
	end
	x = xout;
end

function v0 = v0_update(mu0,mu1,u0,eta0,C1,x,eta1);
	v0 = (mu0*(-u0-eta0)+mu1*(-C1*x+eta1))/(mu0+mu1);
end

function v2 = v2_update(mu2,mu3,u1,eta2,C2,u3,eta3);
	v2 = (mu2*(-u1-eta2)+mu3*(-C2*u3+eta3))/(mu2+mu3);
end

function v4 = v4_update(u2,u3,x,alph,S,eta4,eta5,eta6,mu,v5); 
	mu4 = mu(5);
	mu5 = mu(6);
	mu6 = mu(7);
	v4 = mu4*(mu5+mu6)*(-u2-eta4)-mu5*mu6*((1-alph)*S*u3-eta5)+mu5*mu6*(-alph*S*x+eta6);
	v4 = v4/(2*(mu4*mu5+mu5*mu6+mu4*mu6));

	% not parallel
	%v4 = (mu4*(-u2-eta4)+mu6*(-alph*S*x-v5+eta6))./(mu4+mu6);
end

function v5 = v5_update(u2,u3,x,alph,S,eta4,eta5,eta6,mu,v4); 
	mu4 = mu(5);
	mu5 = mu(6);
	mu6 = mu(7);
	v5 = -mu4*mu6*(-u2-eta4)+mu5*(mu4+mu6)*((1-alph)*S*u3-eta5)+mu4*mu6*(-alph*S*x+eta6);
	v5 = v5/(2*(mu4*mu5+mu5*mu6+mu4*mu6));

	% not parallel
	%v5 = (mu5*((1-alph)*S*u3-eta5)+mu6*(-alph*S*x-v4+eta6))./(mu5+mu6);
end

function v7 = v7_update(mu7,mu8,u3,eta7,x,eta8);
	v7 = (mu7*(-u3-eta7)+mu8*(-x+eta8))/(mu7+mu8);
end
