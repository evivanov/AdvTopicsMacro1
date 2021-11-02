function value = euler_equation(k_pr, k, alpha, beta, delta)
% The Euler equation from which we have to find the roots (in terms of
% k_pr):
value = (k^alpha - (1-delta)*k - k_pr)^(-1) - beta*(1 + alpha*k_pr^(alpha-1) - delta) * (k_pr^alpha - (1-delta)*k_pr - k_pr)^(-1);
end