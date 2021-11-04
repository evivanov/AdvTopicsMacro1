clear all

% Parameters
alpha = 0.3;
beta = 0.96;
delta = 1;
M = 100;

tol = 1e-7;
dif = 1;

kmin_grid = 0.001;
kmax_grid = 1;
k_grid = linspace(kmin_grid, kmax_grid, M)';

%% Value function iteration: brute force
v0 = ones(M,1);
v1 = zeros(M,1);
k11 = zeros(M,1);
w = zeros(M,1);

while dif > tol
    for i = 1:M
        ki = k_grid(i);
        for a = 1:M
            c = ki^alpha + (1-delta)*ki - k_grid(a);
            if c <= 0
                w(a) = -Inf;
            else
                w(a) = log(c) + beta*v0(a);
            end
        end
        [q, j] = max(w);
        k1 = j;
        v1(i) = w(j);
        k11(i) = k1;
    end
    dif = max(abs(v1-v0));
    v0 = v1;
end

%Analytical results
v_true = alpha/(1-alpha*beta)*log(k_grid) + 1/(1-beta)*log(1-alpha*beta)+1/(1-beta)*beta*alpha/(1-alpha*beta)*log(alpha*beta);
g_true = alpha*beta*k_grid.^alpha;

%5
figure('Name', 'Value Function wo')
plot(k_grid, v1);
hold on;
plot(k_grid, v_true, 'g');
title('Value Function without Interpolation');
% legend('without interpolation','true value function', 'Location','southeast');
xlabel('current capital');
ylabel('value function');
% saveas(gcf,'1_Value_function_without','epsc')


g = k_grid(k11);
figure('Name','Policy function w/o')
plot(k_grid, g);
hold on;
plot(k_grid, g_true, 'g');
title('Policy function without Interpolation');
% legend('without interpolation','true policy function','Location','southeast');
xlabel('current capital');
ylabel('next period');
% saveas(gcf,'1_Policy_Function_without','epsc')


%% Value function iteration: exploiting monotonicity
v0 = ones(M,1);
v1 = zeros(M,1);
k11 = zeros(M,1);
w = zeros(M,1);
dif = 1;
monoton_iter = 0;

tic
while dif > tol
    for i = 1:M
        ki = k_grid(i);
        if i == 1
            amin = 1;
        else
            amin = k11(i-1);
        end
        for a = amin:M
            c = ki^alpha + (1-delta)*ki - k_grid(a);
            if c <= 0
                w(a) = -Inf;
            else
                w(a) = log(c) + beta*v0(a);
            end
        end
        [q, j] = max(w);
        k1 = j;
        v1(i) = w(j);
        k11(i) = k1;
    end
    dif = max(abs(v1-v0));
    v0 = v1;
    monoton_iter = monoton_iter + 1;
end
toc
monoton_iter

%% Value function iteration: exploiting concavity
v0 = ones(M,1);
v1 = zeros(M,1);
k11 = zeros(M,1);
w = zeros(M,1);
dif = 1;
conc_iter = 0;

tic
while dif > tol
    for i = 1:M
        ki = k_grid(i);
        for a = 1:M
            c = ki^alpha + (1-delta)*ki - k_grid(a);
            if c <= 0
                w(a) = -Inf;
            else
                w(a) = log(c) + beta*v0(a);
                if a > 1
                    if w(a) < w(a-1)
                        w(a:M) = -Inf*ones(M-a+1,1);
                        break;
                    end
                end
            end
        end
        [q, j] = max(w);
        k1 = j;
        v1(i) = w(j);
        k11(i) = k1;
    end
    dif = max(abs(v1-v0));
    v0 = v1;
    conc_iter = conc_iter + 1;
end
toc
conc_iter

%% Howard's policy function iteration
v0 = ones(M,1);
v1 = zeros(M,1);
k11 = zeros(M,1);
w = zeros(M,1);
dif = 1;
iter_pol = 10;
howard_iter = 0;

tic
while dif > tol
    for i = 1:M
        ki = k_grid(i);
        for a = 1:M
            c = ki^alpha + (1-delta)*ki - k_grid(a);
            if c <= 0
                w(a) = -Inf;
            else
                w(a) = log(c) + beta*v0(a);
            end
        end
        [q, j] = max(w);
        v1(i) = w(j);
        k11(i) = j;
    end
    for i = 1:iter_pol-1
        v1(i+1) = w(k11(i)) + beta*v1(i);
    end
    dif = max(abs(v1-v0));
    v0 = v1;
    howard_iter = howard_iter + 1;
end
toc
howard_iter







