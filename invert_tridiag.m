function A_inv = invert_tridiag(a,b,c,o)
% assume tridiagonal matrix is of form diagonal + tri-diagonal circulant
% A = diag(a)+diag(b,1)+diag(c,-1);
% a: center diagonal
% b: above, superdiag
% c: under, subdiag
% offset of diagonals

% WARNING: theta values scale exponentially with values of a, so for larger
% mu, can cause overflow (Inf) error!

% figure; plot(a); title('center diag');
assert((length(b) == length(c)) && (length(a) == length(b)+o),'incompatible vector lengths');

n = length(a);

if (o == 1)
	theta = NaN(n+1,1);
	phi = NaN(n+1,1);
	
	theta(1) = 1;
	theta(2) = a(1);
	phi(n+1) = 1;
	phi(n) = a(n);
	for ii = 2:n
	    theta(ii+1) = double(a(ii)*theta(ii-1+1)-b(ii-1)*c(ii-1)*theta(ii-2+1));
	    if isnan(theta(ii+1)) || isinf(theta(ii+1))
	        display('overflow! theta too big');
	        keyboard;
	    end
	end
	for ii = n-1:-1:1
	    phi(ii) = a(ii)*phi(ii+1)-b(ii)*c(ii)*phi(ii+2);  
	    if isnan(phi(ii)) || isinf(phi(ii))
	        display('overflow! phi too big');
	        keyboard;
	    end
	end
	
	for ii = 1:n
	    for jj = 1:n
	        if (ii <= jj)
	            A_inv(ii,jj) = (-1)^(ii+jj)*prod(b(ii:jj-1))*theta(ii)*phi(jj+1)/theta(n+1);
	        else
	            A_inv(ii,jj) = (-1)^(ii+jj)*prod(c(jj:ii-1))*theta(jj)*phi(ii+1)/theta(n+1);
	        end
	    end
	end
	
	A = diag(a)+diag(b,1)+diag(c,-1);
	true_inv = inv(A);
	
else % stripes are not on sub/superdiagonal
	
	A = diag(a)+diag(b,o)+diag(c,-o);

	% size of blocks is mxm, which is same as nx san e as o?
	m = o;
	n = length(a)/m; % number of blocks

	F = ones(n*m,m);
	S = eye(m);
	G = ones(m,m*n);
	D = ones(m,m*n);
	D = construct_D(A,m,n);
	R = eye(m);
	%R = construct_R(D,A,m,n);
	E = ones(m*n,m);

	A_inv = zeros(m*n);

	for ii = 1:n
	    for jj = 1:n
	        if (ii <= jj)
	            %A_inv((ii-1)*m+1:ii*m,(jj-1)*m+1:jj*m) = -F((ii-1)*m+1:ii*m,:)*inv(S)*G(:,(jj-1)*m+1:jj*m);
	            A_inv = set_block( -F((ii-1)*m+1:ii*m,:)*inv(S)*G(:,(jj-1)*m+1:jj*m),ii,jj,m);;
	        else
	            %A_inv((ii-1)*m+1:ii*m,(jj-1)*m+1:jj*m) = -E((ii-1)*m+1:ii*m,:)*inv(R)*D(:,(jj-1)*m+1:jj*m);
	            A_inv = set_block( -E((ii-1)*m+1:ii*m,:)*inv(R)*D(:,(jj-1)*m+1:jj*m),ii,jj,m);

	        end
	    end
	end
	display('got here');
	keyboard;	
	
	
end


% norm(A_inv-true_inv) % 3 e-5

 display('pause! start work here');
 keyboard;

end

function block = get_block(A,ii,jj,m)
	block = A((ii-1)*m+1:ii*m,(jj-1)*m+1:jj*m);
end

function A = set_block(block,ii,jj,m)
	A((ii-1)*m+1:ii*m,(jj-1)*m+1:jj*m) = block;
end

function D = construct_D(A,m,n)
	D = zeros(m,m*n);
	D(:,1:m) = eye(m);
	A11 = get_block(A,1,1,m);
	A21 = get_block(A,2,1,m);
	% is inv(A21) always going to be -I?
	D(:,m+1:2*m) = -A11*inv(A21);
	for ii = 3:n
		A21 = get_block(A,ii-2,ii-1,m);
		A11 = get_block(A,ii-1,ii-1,m);
		A01 = get_block(A,ii,ii-1,m);
		D2 = D(:,(ii-2-1)*m+1:(ii-2)*m);
		D1 = D(:,(ii-1-1)*m+1:(ii-1)*m);
	% is inv(A01) always going to be -I?
		D(:,(ii-1)*m+1:ii*m) = -(D2*A21+D1*A11)*inv(A01);
	end
%	display('done with D');
%	keyboard;
end










