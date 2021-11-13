beta = 0.99;
alpha = 0.33;
rho = 0.95;
sigma = 0.01;


A = alpha*beta;
B = 2*alpha*beta -alpha^2*beta - 1;
C = -alpha;
D = -alpha*beta;
E = -1;


H1 = 1/(2*A) * (-B + sqrt(B^2 - 4*A*C));
H2 = 1/(2*A) * (-B - sqrt(B^2 - 4*A*C));

%use H2

G = (-D*rho-E)/(A*H2+A*rho+B);