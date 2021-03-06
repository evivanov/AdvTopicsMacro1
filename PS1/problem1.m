clear all

% Parameters
alpha = 0.3;
beta = 0.96;
delta = 1;
alpha = 0.5;
beta = 0.8;
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

brute_iter=0;
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
        k1 = j;
        v1(i) = w(j);
        k11(i) = k1;
    end
    dif = max(abs(v1-v0));
    v0 = v1;
    brute_iter=brute_iter+1;
end
toc
brute_iter

%Analytical results
v_true = alpha/(1-alpha*beta)*log(k_grid) + 1/(1-beta)*log(1-alpha*beta)+1/(1-beta)*beta*alpha/(1-alpha*beta)*log(alpha*beta);
g_true = alpha*beta*k_grid.^alpha;

%5
figure('Name', 'Value Function wo')
plot(k_grid, v1);
hold on;
plot(k_grid, v_true, 'g');
title('Value Function');
legend('computed value function','true value function', 'Location','southeast');
xlabel('current capital');
ylabel('value function');
saveas(gcf,'1_Value_function_without','epsc')


g = k_grid(k11);
figure('Name','Policy function w/o')
plot(k_grid, g);
hold on;
plot(k_grid, g_true, 'g');
title('Policy function');
legend('computed policy function','true policy function','Location','southeast');
xlabel('current capital');
ylabel('next period');
saveas(gcf,'1_Policy_Function_without','epsc')


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


%% Howard's policy function (and using Umat and Vmat...)
tic
Umat = zeros(M,M); % rows are K and cols are K'
V0 = zeros(M,1);
how_iter = 30;
V_iter = zeros(M,how_iter);
dif = 1;

% Create Umat:
for i=1:M
    for j=1:M
        c = k_grid(i)^alpha + (1-delta)*k_grid(i) - k_grid(j); %i=K, j=K'
        if (c >0)
            Umat(i,j) = log(c);
        else
            Umat(i,j) = -Inf;
        end
    end
end

iterations = 0;
% Perform operation:
while dif > tol
    Vmat = repmat(V0',M,1); %transpose is very important for some reason...
    W = Umat + beta*Vmat;
    [V,pol] = max(W,[],2);
    
    % Policy iteration:
    V_iter(:,1) = V;
    for i=2:how_iter % number of iters to do for howard's policy iteration
        for j = 1:M % iter over all rows
            V_iter(j, i) = Umat(j,pol(j)) + beta*V_iter(pol(j),i-1);
        end
    end
    V = V_iter(:,how_iter);

    dif = max(abs(V-V0));
    V0 = V;
    iterations = iterations + 1;
end
iterations
toc




