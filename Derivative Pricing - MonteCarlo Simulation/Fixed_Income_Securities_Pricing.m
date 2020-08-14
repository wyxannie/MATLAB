
clear; close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Question 1 (Vasicek Model) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r0 = 0.05; sigma = 0.1; k = 0.82; r_bar = 0.05;
%%%%%%%%%%%
%%%% Part a
%%%%%%%%%%%
rng('default');
F = 1000;
T = 0.5; dt = 1/252;
Ndt = round(T/dt); dt = T/Ndt;
N = 10000;

Z = normrnd(0, 1,[Ndt+1, N]);
r = zeros([Ndt+1, N]); r(1, :) = r0;
for i = 2:(Ndt+1)
    r(i, :) = r(i-1, :) + k*(r_bar - r(i-1, :))*dt + sigma*sqrt(dt)*Z(i, :);
end
Price_Q1_a = F* mean(exp(-sum(r)*dt))

%%%%%%%%%%%
%%%% Part b
%%%%%%%%%%%
rng('default');
F = 1000;
T = 4; dt = 1/252;
Ndt = round(T/dt); dt = T/Ndt;
N = 10000;

Z = normrnd(0, 1,[Ndt+1, N]);
r = zeros([Ndt+1, N]); r(1, :) = r0;
for i = 2:(Ndt+1)
    r(i, :) = r(i-1, :) + k*(r_bar - r(i-1, :))*dt + sigma*sqrt(dt)*Z(i, :);
end
CF = repmat(30, [1, T*2]); CF(end) = CF(end) + F;
T_cf = [0.5:0.5:4] ./ dt + 1;
DiscountFactor = exp(-cumsum(r)*dt);
Price_Q1_b = mean(CF*DiscountFactor(T_cf, :))

%%%%%%%%%%%
%%%% Part c
%%%%%%%%%%%
rng('default');
F = 1000;
T_bond= 0.5; dt = 1/252;

K = 980; T_option = 3/12;
Ndt = round(T_option/dt); dt = T_option/Ndt;
N = 10000;

Z = normrnd(0, 1,[Ndt+1, N]);
r = zeros([Ndt+1, N]); r(1, :) = r0;
for i = 2:(Ndt+1)
    r(i, :) = r(i-1, :) + k*(r_bar - r(i-1, :))*dt + sigma*sqrt(dt)*Z(i, :);
end
r_To = r(end, :); % interest rate at the time of option maturity

%%%% Vasicek Explicit Formula for bond price P(T,S)
B = (1/k) * (1-exp(-k*(T_bond - T_option)));
A = exp((r_bar - ((sigma^2)/(2*k^2)))*(B - (T_bond - T_option)) - ((sigma^2)/(4*k)*(B^2)));
Price_To_Tb = A*exp(-B*r_To);

Price_Q1_c = mean(exp(-sum(r)*dt) .* max(Price_To_Tb*F - K, 0))


%%%%%%%%%%%
%%%% Part d
%%%%%%%%%%%

%%%% Use Jamshidian method + MC simulation

rng('default');

F = 1000;
T_bond= 4; dt = 1/252;

K = 980; T_option = 3/12;
Ndt = round(T_option/dt); dt = T_option/Ndt;
N = 100; % N = 10000;

%%%% Solve r*, then Kseq(Ki)
CF = repmat(30, [1, T_bond*2]); CF(end) = CF(end) + F;
T_cf = [T_option:0.5:(T_bond - T_option)]';
r_star = fsolve(@(x) CF*exp(-x*T_cf) - K, 0);
Kseq = exp(-r_star*T_cf)';

%%%% call option on ZCB - MC simulation
Z = normrnd(0, 1,[Ndt+1, N]);
r = zeros([Ndt+1, N]); r(1, :) = r0;
for i = 2:(Ndt+1)
    r(i, :) = r(i-1, :) + k*(r_bar - r(i-1, :))*dt + sigma*sqrt(dt)*Z(i, :);
end
r_To = r(end, :); % interest rate at the time of option maturity

%%%% restart and simulate the ZCB
Ndt_new = round((T_bond - T_option)/dt); dt_new = (T_bond - T_option)/Ndt_new;
N_new = 1000;
Price_ZCB = zeros(length(Kseq), length(r_To));
for r0_new = r_To
    
    Z_new = normrnd(0, 1, [Ndt_new+1, N_new]);
    r_new = zeros(size(Z_new)); r_new(1, :) = r0_new;
    for i = 2:(Ndt_new+1)
        r_new(i, :) = r_new(i-1, :) + k*(r_bar - r_new(i-1, :))*dt_new + sigma*sqrt(dt_new)*Z_new(i, :);
    end
    
    T_cf_new = (T_option:0.5:(T_bond - T_option)) ./ dt_new + 1;
    DiscountFactor_new = exp(-cumsum(r_new)*dt_new);
    Price_ZCB(:, r_To == r0_new) = mean(DiscountFactor_new(T_cf_new, :), 2);
    
end
OptionPrice_ZCB = repmat(exp(-sum(r)*dt), [8, 1]) .* max(Price_ZCB - repmat(Kseq', [1, N]), 0);
Price_Q1_d = mean(CF * OptionPrice_ZCB)

%%%%%%%%%%%
%%%% Part e
%%%%%%%%%%%
%%%% Use Jamshidian method + explicit formula

rng('default');

F = 1000;
T_bond= 4; dt = 1/252;

K = 980; T_option = 3/12;
Ndt = round(T_option/dt); dt = T_option/Ndt;
N = 10000; 

%%%% Solve r*, then Kseq(Ki)
CF = repmat(30, [1, T_bond*2]); CF(end) = CF(end) + F;
T_cf = [T_option:0.5:(T_bond - T_option)]';
r_star = fsolve(@(x) CF*exp(-x*T_cf) - K, 0);
Kseq = exp(-r_star*T_cf)';

%%%% call option on ZCB - explicit formula
Z = normrnd(0, 1,[Ndt+1, N]);
r = zeros([Ndt+1, N]); r(1, :) = r0;
for i = 2:(Ndt+1)
    r(i, :) = r(i-1, :) + k*(r_bar - r(i-1, :))*dt + sigma*sqrt(dt)*Z(i, :);
end
r_To = r(end, :); % interest rate at the time of option maturity

%%%% Vasicek Explicit Formula for bond price P(T,S)
B = (1/k) * (1-exp(-k*((0.5:0.5:T_bond) - T_option)));
A = exp((r_bar - ((sigma^2)/(2*k^2)))*(B - ((0.5:0.5:T_bond) - T_option)) - ((sigma^2)/(4*k)*(B .^ 2)));
Price_ZCB = repmat(A', [1, length(r_To)]) .* exp(-repmat(B', [1, length(r_To)]) .* repmat(r_To, [length(B), 1]));

OptionPrice_ZCB = repmat(exp(-sum(r)*dt), [8, 1]) .* max(Price_ZCB - repmat(Kseq', [1, N]), 0);
Price_Q1_e = mean(CF * OptionPrice_ZCB)


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Question 2 (CIR Model) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r0 = 0.05; sigma = 0.12; k = 0.92; r_bar = 0.055;
%%%%%%%%%%%
%%%% Part a
%%%%%%%%%%%

%%%% MC simulation + CIR model

F = 1000;
T_bond = 1; dt = 1/252;

T_option = 0.5; K = 980;
Ndt = round(T_option/dt); dt = T_option/Ndt;
N = 100;

Z = normrnd(0, 1,[Ndt+1, N]);
r = zeros([Ndt+1, N]); r(1, :) = r0;
for i = 2:(Ndt+1)
    r(i, :) = r(i-1, :) + k*(r_bar - r(i-1, :))*dt + sigma*sqrt(r(i-1, :)).*sqrt(dt).*Z(i, :);
end
r_To = r(end, :); % interest rate at the time of option maturity

%%%% restart and simulate the ZCB
Ndt_new = round((T_bond - T_option)/dt); dt_new = (T_bond - T_option)/Ndt_new;
N_new = 100; %N_new = 100;
Price_ZCB = zeros(size(r_To));
for r0_new = r_To

    Z_new = normrnd(0, 1, [Ndt_new+1, N_new]);
    r_new = zeros(size(Z_new)); r_new(1, :) = r0_new;
    for i = 2:(Ndt_new+1)
        r_new(i, :) = r_new(i-1, :) + k*(r_bar - r_new(i-1, :))*dt_new + sigma*sqrt(r_new(i-1, :)).*sqrt(dt_new).*Z_new(i, :);
    end
    Price_ZCB(r_To == r0_new) = mean(F*exp(-sum(r_new)*dt_new));
    
end
Price_Q2_a = mean(exp(-sum(r)*dt).* max(Price_ZCB-K, 0))

%%%%%%%%%%%
%%%% Part b
%%%%%%%%%%%

%%%% Implicit Finite-Difference Method + CIR model

F = 1000;
T_bond = 1; dt = 1/252;

T_option = 0.5; K = 980; 
Ndt = round(T_option/dt); dt = T_option/Ndt;

%%%% Implicit Finite-Difference Method
dr = 0.01; rMax = 0.10; rseq = rMax:-dr:0;
OptionGrid = zeros(length(rseq), Ndt+1);
[nrow, ncol] = size(OptionGrid);

%%%% Payoff P(T,S) - explicit formula
h1 = sqrt(k^2 + 2*(sigma^2)); h2 = (k+h1)/2; h3 = (2*k*r_bar)/(sigma^2);
A_cir = ((h1*exp(h2*(T_bond-T_option))) / (h2*(exp(h1*(T_bond-T_option))-1)+h1))^h3;
B_cir = (exp(h1*(T_bond-T_option))-1) / (h2*(exp(h1*(T_bond-T_option))-1)+h1);
Price_ZCB = A_cir*exp(-B_cir*rseq);

OptionGrid(:, end) = max(Price_ZCB*F - K, 0);
OptionGrid(end, :) = A_cir*F - K;
OptionGrid(1, :) = 0;
 
j_param = (length(rseq)-2):-1:1;
p_u = (-1/2/dr) * dt* ((sigma^2)*j_param + k*(r_bar - dr*j_param));
p_m = 1 + (dt/dr) * (sigma^2) * j_param + dr*j_param*dt;
p_d = (-1/2/dr) * dt* ((sigma^2)*j_param - k*(r_bar - dr*j_param));
A = [diag(p_u,0),zeros(length(j_param), 2)] +...
    [zeros(length(j_param), 1), diag(p_m,0),zeros(length(j_param), 1)] +...
    [zeros(length(j_param), 2), diag(p_d,0)];
B = zeros(nrow, ncol-1);
B(:, end) = OptionGrid(:, ncol);
B(1, :) = 0;
B(end, :) = OptionGrid(end, 1:(end-1));
for i = ncol-1:-1:1
    OptionGrid(:, i) = A \ B(2:(end-1), i);
    OptionGrid(end, i) = A_cir*F - K;
    OptionGrid(1, i) = 0;
    if i > 1
        B(2:(end-1), i-1) = OptionGrid(2:(end-1), i);
    end
end
Price_Q2_b = interp1(rseq, OptionGrid(:, 1), 0.05)

%%%%%%%%%%%
%%%% Part c
%%%%%%%%%%%

%%%% Explicit formula + CIR model
F = 1000;
T_bond = 1;
T_option = 0.5; K = 980;
t = 0; T = T_option; S = T_bond;

h1 = sqrt(k^2 + 2*(sigma^2)); h2 = (k+h1)/2; h3 = (2*k*r_bar)/(sigma^2);
A_tS = ((h1*exp(h2*(S-t))) / (h2*(exp(h1*(S-t))-1)+h1))^h3;
B_tS = (exp(h1*(S-t))-1) / (h2*(exp(h1*(S-t))-1)+h1);
A_tT = ((h1*exp(h2*(T-t))) / (h2*(exp(h1*(T-t))-1)+h1))^h3;
B_tT = (exp(h1*(T-t))-1) / (h2*(exp(h1*(T-t))-1)+h1);
A_TS = ((h1*exp(h2*(S-T))) / (h2*(exp(h1*(S-T))-1)+h1))^h3;
B_TS = (exp(h1*(S-T))-1) / (h2*(exp(h1*(S-T))-1)+h1);

P_tS = A_tS*exp(-B_tS*r0);
P_tT = A_tT*exp(-B_tT*r0);

theta = sqrt(k^2 + 2*sigma^2);
phi = (2*theta) / ((sigma^2)*(exp(theta*(T-t))-1));
psy = (k + theta) / (sigma^2);
r_star = log(A_TS/K*F) / B_TS;

Price_Q2_c = P_tS*F * ncx2cdf(2*r_star*(phi+psy+B_TS), (4*k*r_bar)/(sigma^2), (2*(phi^2)*r0*exp(theta*(T-t)))/(phi+psy+B_TS))...
    - K*P_tT * ncx2cdf(2*r_star*(phi+psy), (4*k*r_bar)/(sigma^2), (2*(phi^2)*r0*exp(theta*(T-t)))/(phi+psy))


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Question 3 (G2++ Model) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = 0; y0 = 0; phi0 = 0.03; r0 = 0.03; phit = 0.03;
rho = 0.7; a = 0.1; b = 0.3; sigma = 0.03; eta = 0.08;

F = 1000; T_bond = 1;
K = 950; T_option = 0.5; dt = 1/252;

%%%%%%%%%%%%%%%%%%
%%%% MC simulation
%%%%%%%%%%%%%%%%%%
Ndt = round(T_option/dt); dt = T_option/Ndt;
N = 1000;

%%%% Dependent Brownian Motion
Z_1 = normrnd(0, 1,[Ndt+1, N]); Z_2 = normrnd(0, 1,[Ndt+1, N]);
W_1 = 0+1*Z_1; W_2 = 0 + 1*rho*Z_1 + 1*sqrt(1-rho^2)*Z_2;

r = zeros([Ndt+1, N]); r(1, :) = r0;
x = zeros(size(r)); x(1, :) = x0;
y = zeros(size(r)); y(1, :) = y0;
for i = 2:(Ndt+1)
    x(i, :) = x(i-1, :) + (-a*x(i-1, :))*dt + sigma*sqrt(dt)*W_1(i, :);
    y(i, :) = y(i-1, :) + (-b*y(i-1, :))*dt + eta*sqrt(dt)*W_2(i, :);
    r(i, :) = x(i, :) + y(i, :) + phit;
end
r_To = r(end, :); % interest rate at the time of option maturity

%%%% restart and simulate the ZCB
Ndt_new = round((T_bond - T_option)/dt); dt_new = (T_bond - T_option)/Ndt_new;
N_new = 1000;
Price_ZCB = zeros(size(r_To));
for r0_new = r_To

    Z_1 = normrnd(0, 1,[Ndt_new+1, N_new]); Z_2 = normrnd(0, 1,[Ndt_new+1, N_new]);
    W_1_new = 0+1*Z_1; W_2_new = 0 + 1*rho*Z_1 + 1*sqrt(1-rho^2)*Z_2;
    r_new = zeros([Ndt_new+1, N_new]); r_new(1, :) = r0;
    x_new = zeros(size(r_new)); x_new(1, :) = x0;
    y_new = zeros(size(r_new)); y_new(1, :) = y0;
    for i = 2:(Ndt_new+1)
        x_new(i, :) = x_new(i-1, :) + (-a*x_new(i-1, :))*dt_new + sigma*sqrt(dt_new)*W_1_new(i, :);
        y_new(i, :) = y_new(i-1, :) + (-b*y_new(i-1, :))*dt_new + eta*sqrt(dt_new)*W_2_new(i, :);
        r_new(i, :) = x_new(i, :) + y_new(i, :) + phit;
    end
    Price_ZCB(r_To == r0_new) = mean(F*exp(-sum(r_new)*dt_new));
    
end
Price_Q3_1 = mean(exp(-sum(r)*dt).* max(K - Price_ZCB, 0))

%%%%%%%%%%%%%%%%%%%%%
%%%% Explicit formula
%%%%%%%%%%%%%%%%%%%%%
t = 0; T = T_option; S = T_bond;

V_tT = (sigma^2)/(a^2) * (T-t + (2/a)*exp(-a*(T-t)) - (1/(2*a))*exp(-2*a*(T-t)) - 3/(2*a))...
    + (eta^2)/(b^2) * (T-t + (2/b)*exp(-b*(T-t)) - (1/(2*b))*exp(-2*b*(T-t)) - 3/(2*b))...
    + 2*rho*(sigma*eta)/(a*b) * (T-t + (exp(-a*(T-t))-1)/a + (exp(-b*(T-t))-1)/b - (exp(-(a+b)*(T-t))-1)/(a+b));
V_tS = (sigma^2)/(a^2) * (S-t + (2/a)*exp(-a*(S-t)) - (1/(2*a))*exp(-2*a*(S-t)) - 3/(2*a))...
    + (eta^2)/(b^2) * (S-t + (2/b)*exp(-b*(S-t)) - (1/(2*b))*exp(-2*b*(S-t)) - 3/(2*b))...
    + 2*rho*(sigma*eta)/(a*b) * (S-t + (exp(-a*(S-t))-1)/a + (exp(-b*(S-t))-1)/b - (exp(-(a+b)*(S-t))-1)/(a+b));

P_tT = exp(-T*phit - ((1-exp(-a*(T-t)))/a)*x0 - ((1-exp(-b*(T-t)))/b)*y0 + (1/2)*V_tT);
P_tS = exp(-S*phit - ((1-exp(-a*(S-t)))/a)*x0 - ((1-exp(-b*(S-t)))/b)*y0 + (1/2)*V_tS);

SIGMA2 = (sigma^2)/(2*a^3)*((1-exp(-a*(S-T)))^2) * (1-exp(-2*a*(T-t)))...
    + (eta^2)/(2*b^3) * ((1-exp(-b*(S-T)))^2) * (1-exp(-2*b*(T-t)))...
    + 2*rho* (sigma*eta)/(a*b*(a+b)) * (1-exp(-a*(S-T))) * (1-exp(-b*(S-T))) * (1-exp(-(a+b)*(T-t)));
SIGMA = sqrt(SIGMA2);

Price_Q3_2 = -P_tS*F * normcdf((log(K*P_tT/(P_tS*F)) / SIGMA) - (1/2)*SIGMA)...
    + P_tT*K * normcdf((log(K*P_tT/(P_tS*F)) / SIGMA) + (1/2)*SIGMA)