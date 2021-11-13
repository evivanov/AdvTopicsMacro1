beta = 0.99;
alpha = 0.33;
rho = 0.95;
sigma = 0.01;

k_ss = (alpha*beta)^(1/(1-alpha));

A = alpha*beta;
B = -alpha^2*beta - 1;
C = alpha;
D = -alpha*beta;
E = 1;


H1 = 1/(2*A) * (-B + sqrt(B^2 - 4*A*C));
H2 = 1/(2*A) * (-B - sqrt(B^2 - 4*A*C));

%use H2
H = H2

G = (-D*rho-E)/(A*H+A*rho+B)



rng(1);
epsilon = normrnd(0,1, 200, 1);

z = zeros(200,1);
z(1) = sigma * epsilon(1);
for t = 2:200
    z(t) = rho*z(t-1) + sigma* epsilon(t);
end

k_analytical = zeros(200,1);
k_analytical(1) = k_ss;

for t = 2:200
    k_analytical(t) = alpha*beta*exp(z(t))*k_analytical(t-1)^alpha;
end

figure('Name', 'Analytical path')
plot(k_analytical)

k_loglin = zeros(200,1);

for t = 2:200
    k_loglin(t) = H*k_loglin(t-1) + G*z(t);
end

k_loglin = exp(k_loglin + log(k_ss));

figure('Name', 'Loglin')
plot(k_loglin)

