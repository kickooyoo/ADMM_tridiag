function [cost, costDT, costRT] = calc_cost_tridiag_inpaint(y, D, CH, CV, x, beta, varargin)
%function [cost, costDT, costRT] = calc_cost_tridiag_inpaint(y, D, CH, CV, x, beta, varargin)

arg.potx = [];
arg.poty = [];
arg.mu0 = 1;
arg.mu1 = 1;
arg = vararg_pair(arg, varargin);

if isempty(arg.potx) || isempty(arg.poty)
	arg.potx = potential_fun('l1', beta/arg.mu0);
	arg.poty = potential_fun('l1', beta/arg.mu1);
end

costDT = norm(double(col(y) - col(D * x)),2)^2/2;

costRT = sum(col(beta(:) .* arg.potx.potk(col(CH * x))), 'double') + sum(col(beta(:) .* arg.poty.potk(col(CV * x))), 'double');


cost = costDT + costRT;
