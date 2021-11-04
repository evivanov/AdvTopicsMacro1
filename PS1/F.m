function F = F(k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
z0 = 1;
z1 = 2;
alpha = 0.5;
beta = 0.8;
delta = 1;

k_star = ( 1/(alpha*z0) * ( 1/beta - 1 + delta ) )^(1/(alpha-1));
k_star_new = ( 1/(alpha*z1) * ( 1/beta - 1 + delta ) )^(1/(alpha-1));


F = zeros(201,1);
F(1) = k_star - k(1);
for i = 2:199
    F(i) = 1/(z1*k(i)^alpha + (1-delta)*k(i) - k(i+1)) - beta*(z1*alpha*k(i+1)^(alpha-1)/(z1*(k(i+1))^alpha + (1-delta)*k(i+1) - k(i+2)));
end
F(200) = k_star_new - k(200);
F(201) = k_star_new - k(201);
end