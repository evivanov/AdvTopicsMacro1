var  c k z;

varexo eps;

parameters alpha beta rho sigma A B C D E;

beta = 0.99;
alpha = 0.33;
rho = 0.95;
sigma = 0.01;

A = alpha*beta;
B = 2*alpha*beta -alpha^2*beta - 1;
C = -alpha;
D = -alpha*beta;
E = -1;


model(linear);

%1/(exp(z)*k(-1)^alpha - k) = beta*( exp(z(+1)) * alpha * k^(alpha-1) / (exp(z(+1))*k^alpha - k(+1)));
%0 = A*k(+1) + B*k + C*k(-1) + D*z(+1) + E*z;

0 = z(+1) + (alpha-1)*k - c(+1) - c;
c = 1/(1-alpha*beta) * (z + alpha*k(-1)) - (alpha*beta)/(1-alpha*beta) *k;


z = rho*z(-1) + eps;

end;



%% Shocks
shocks;
var eps = sigma^2;
end;

model_diagnostics;

%% Results
stoch_simul(order = 1);
