x = ones(5);
x(:,end) = x(:,end)+1;
x(end,:) = x(end,:)+1;
x(end-1,:) = x(end-1,:)+1;
x(2,2) = x(2,2)+2;

y = [zeros(5,2) ones(5,2) 2*ones(5,1)];

[c1, c2] = construct_finite_diff([5 5]);

y1 = reshape(c1*x(:),5,5);
y2 = reshape(c2*x(:),5,5);


%c11 = diag(ones(5,1));
%c11 = c11(2:end,:);
%c11 = c11*c1;


% this is not a fatrix but a sparse matrix!
%c0 = C2sparse('tight',true(5),4,1);
