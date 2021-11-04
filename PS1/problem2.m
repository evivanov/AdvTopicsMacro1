clear all

% Parameters
alpha = 0.3;
beta = 0.97;
delta = 0.01;

z0 = 1;
z1 = 2*z0;

k_star = ( 1/(alpha*z0) * ( 1/beta - 1 + delta ) )^(1/(alpha-1));
k_star_new = ( 1/(alpha*z1) * ( 1/beta - 1 + delta ) )^(1/(alpha-1));

T = 200;

%%

k_sys = fsolve(@(x) F(x), k_star_new*ones(201,1));
plot(1:T, k_sys(1:T))

%%
k = zeros(T,1);
k(1) = k_star;
c = zeros(T,1);

lower_bound_guess_k2 = k(1);
upper_bound_guess_k2 = z1*k(1)^alpha + (1-delta) * k(1);
k(2) = (lower_bound_guess_k2+upper_bound_guess_k2)/2;

c(1) = z1*k(1)^alpha + (1-delta) * k(1) - k(2);

ki = k(2);
ki1 = 0;

while abs(ki1 - ki) > 1e-45
    ki = k(2);
    k(3:T) = zeros(T-3+1,1);
    c(1) = z1*k(1)^alpha + (1-delta) * k(1) - k(2);
    c(2) = beta*c(1)*( z1*alpha*k(2)^(alpha-1) + 1 - delta);
    for t = 3:T
       k(t) = z1*k(t-1)^alpha + (1-delta)*k(t-1) - c(t-1);
       if k(t) < 0
           lower_bound_guess_k2 = ki;
           break;
       end
       if k(t) > k_star_new
           upper_bound_guess_k2 = ki;
           break;
       end
       c(t) = beta*c(t-1)*(z1*alpha*k(t)^(alpha-1) + 1 - delta);
    end
    if (k(T) < k_star_new && k(T) > 0)
           lower_bound_guess_k2 = ki;
    end
    if (k(T) > k_star_new)
           upper_bound_guess_k2 = ki;
    end
    ki1 = (lower_bound_guess_k2+upper_bound_guess_k2)/2;
    k(2) = ki1;
end



plot(1:T, k(1:T))


%%

k = zeros(T,1);
k(T) = k_star_new;
c = zeros(T,1);

lower_bound_guess_k_t_1 = k_star;
upper_bound_guess_k_t_1 = k(T);
%z1*k(1)^alpha + (1-delta) * k(1);
k(T-1) = k_star_new;%(lower_bound_guess_k_t_1+upper_bound_guess_k_t_1)/2;

c(T) = z1*k(T)^alpha + (1-delta) * k(T) - k(T);
c(T-1) = (1/beta)*c(T)/( z1*alpha*k(T)^(alpha-1) + 1 - delta);


ki = k(T-1);
ki1 = 0;

while abs(lower_bound_guess_k_t_1 - upper_bound_guess_k_t_1) > 1e-3
    ki = k(T-1);
    k(1:T-2) = zeros(T-3+1,1);
    c(1:T-2) = zeros(T-3+1,1);
    for t = T-2:-1:1
        c(t) = z1*k(t)^alpha + (1-delta) * k(t) - k(t+1);
        if c(t) < 0
            break;
        end
        func = @(x) 1/(z1*x^alpha + (1-delta)*x - k(t+1)) - beta*(z1*alpha*k(t+1)^(alpha-1)/(z1*(k(t+1))^alpha + (1-delta)*k(t+1) - k(t+2)));
        k(t) = lsqnonlin(func, k_star_new, k_star, k_star_new);
        if k(t) < k_star
           lower_bound_guess_k_t_1 = ki;
           break;
       end
       if k(t) > k_star_new
           upper_bound_guess_k_t_1 = ki;
           break;
       end
       if (k(t) > k(t+1))
           lower_bound_guess_k_t_1 = ki;
           break;
       end
    end
%     if (k(1) < k_star && k(1) > 0)
%            lower_bound_guess_k_t_1 = ki;
%     end
%     if (k(1) > k_star_new)
%            upper_bound_guess_k_t_1 = ki;
%     end
    ki1 = (lower_bound_guess_k_t_1+upper_bound_guess_k_t_1)/2;
    k(T-1) = ki1;
end



plot(1:T, k(1:T))


%%

alpha = 0.3;
beta = 0.9;
delta = 0.01;

z0 = 1;
z1 = 2*z0;

k_star = ( 1/(alpha*z0) * ( 1/beta - 1 + delta ) )^(1/(alpha-1));
k_star_new = ( 1/(alpha*z1) * ( 1/beta - 1 + delta ) )^(1/(alpha-1));


N = 500;
k_grid = linspace(k_star, k_star_new, N);

k = zeros(T+1,1);
c = zeros(T,1);

k(T+1) = k_star_new;
k(T) = k_star_new;

for i = N:-1:1
    k(T-1) = k_grid(i);
    c(T-1) = z1*k(T-1)^alpha + (1-delta) * k(T-1) - k(T);
    if c(T-1) < 0
        continue;
    end
    for t = T-2:-1:1
        c(t) = z1*k(t)^alpha + (1-delta) * k(t) - k(t+1);
        if c(t) < 0
            break;
        end
        func = @(x) 1/(z1*x^alpha + (1-delta)*x - k(t+1)) - beta*(z1*alpha*k(t+1)^(alpha-1)/(z1*(k(t+1))^alpha + (1-delta)*k(t+1) - k(t+2)));
        k(t) = lsqnonlin(func, k_star_new, k_star, k_star_new);
    end
    if abs(k(1) - k_star) < 1e-5
        break
    end
end

plot(1:T, k(1:T))