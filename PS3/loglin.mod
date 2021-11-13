close all;

var c k z;

varexo eps;

parameters alpha beta rho sigma;

beta = 0.99;
alpha = 0.33;
rho = 0.95;
sigma = 0.01;

model;

1/c = beta*( exp(z(+1)) * alpha * k^(alpha-1) / c(+1));
c = exp(z)*k(-1)^alpha - k;

%1/(exp(z)*k(-1)^alpha - k) = beta*( exp(z(+1)) * alpha * k^(alpha-1) /(exp(z(+1))*k^alpha - k(+1) ));
log(z) = rho*log(z(-1)) + sigma*eps;

end;

initval;

k = (alpha*beta)^(1/(1-alpha));
c = k^alpha - k;
z = 1;

end;

%% Shocks
shocks;
var eps = 1;
end;

model_diagnostics;

%% Results
stoch_simul(loglinear, order = 1);
