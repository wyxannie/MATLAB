
clear; clc;

r = 0.06; sigma = 0.2; 
S0 = 40; X = 40; 

N  =100000;

Tvalues = [0.5, 1, 2];

%%%% 1.a Laguerre Polynomials
for T = Tvalues
    disp(['Laguerre Polynomials Method:  ', 'When T = ', num2str(T)])
    result = Project_5(S0, X, T, r, sigma, N, "Laguerre")
end


%%%% 1.b Hermite Polynomials
for T = Tvalues
    disp(['Hermite Polynomials Method:  ', 'When T = ', num2str(T)])
    result = Project_5(S0, X, T, r, sigma, N, "Hermite")
end


%%%% 1.c Simple Monomials
for T = Tvalues
    disp(['Simple Monomials Method:  ', 'When T = ', num2str(T)])
    result = Project_5(S0, X, T, r, sigma, N, "Monomials")
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = Project_5(S0, X, T, r, sigma, N, method)

% To ensure the stability of the algo, use stock price divided by strike
% price to do the simulation
const = X;
S0 = S0/const; X = X/const;

kvalues = 2:4; nobj = length(kvalues);
result = table(zeros(nobj,1), zeros(nobj,1), zeros(nobj,1));
result.Properties.VariableNames = {'T', 'k', 'Price'};

count = 0;
for k = kvalues
    count = count +1;
    dt = T / round(T*sqrt(N)); % recalculate the dt to ensure equally divided time intervals

    StockPrice = SimStockPrice(S0, N, T, dt, r, sigma); % simulate stock price
    OptionPrice = LSMC(StockPrice, dt, X, r, method, k) * const; % use LSMC to price American put option

    result{count, 'T'} = T;
    result{count, 'k'} = k;
    result{count, 'Price'} = OptionPrice;

end

end


%%%%%%%%%%%%%%%%%%%%%%%%% Simulate Stock Price %%%%%%%%%%%%%%%%%%%%%%%%%%%

function StockPrice = SimStockPrice(S0, N, T, dt, r, sigma)
% Simulate stock price paths

StockPrice = zeros(N, round(T/dt) + 1); StockPrice(:, 1) = S0;
Z = randn(N/2, round(T/dt));
Z = [Z; -Z]; % anthetic

for i = 2: round(T/dt) + 1
    StockPrice(:, i) = StockPrice(:, i-1) +...
        r*StockPrice(:, i-1)*dt + sigma*StockPrice(:, i-1) .* sqrt(dt) .* Z(:, i-1);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LSMC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function OptionPrice = LSMC(StockPrice, dt, K, r, method, k)

% create early exercise flag matrix
Flag = zeros(size(StockPrice));

Payoff = zeros(size(StockPrice));
ECV = zeros(size(StockPrice)); EV = zeros(size(StockPrice));

% Payoff at time to maturity
Payoff(:, end) = max(K - StockPrice(:, end), 0);
Flag(:, end) = Payoff(:, end) > 0;

for t = size(StockPrice, 2):-1:3

ITM_index = find(K > StockPrice(:, t-1)); % select ITM path at time t-1
Y = Payoff(ITM_index, t) * exp(-r*dt); % discount Y back to time t-1; create Y variable
X = StockPrice(ITM_index, t-1); % create X variable

ECV(ITM_index, t-1) = LSFun(X, Y, method, k); % compute ECV by using least square
EV(ITM_index, t-1) = K - StockPrice(ITM_index, t-1); % compare EV and ECV

% update flag matrix
Flag(ITM_index, t-1) = EV(ITM_index, t-1) > ECV(ITM_index, t-1);
Flag(ITM_index, t:end) = 0;

% update payoff matrix
Payoff(Flag(:, t-1) == 1, t-1) = EV(Flag(:, t-1) == 1, t-1);
Payoff(Payoff(:, t)>0 & Flag(:, t-1)==0, t-1) = Payoff(Payoff(:, t)>0 & Flag(:, t-1)==0, t) * exp(-r*dt);

end

OptionPrice = mean(Payoff(:, 2)*exp(-r*dt));

end


%%%%%%%%%%%% Basis Function of Least-Square Estimation %%%%%%%%%%%%%%%%%%%%

function ECV = LSFun(x, y, method, kterm)
% kterm: max=4

if strcmp(method, 'Laguerre')
    Lterms = exp(-x/2) .* [ones(size(x)), 1-x, 1-2*x+(x.^2)/2, 1-3*x+3*(x.^2)/2-(x.^3)/6];
    X = Lterms(:, 1:kterm);
    beta = (X' * X) \ (X' * y);
    ECV = X * beta;
end

if strcmp(method, 'Hermite')
    Lterms = [ones(size(x)), 2*x, 4*(x.^2)-2, 8*(x.^3)-12*x];
    X = Lterms(:, 1:kterm);
    beta = (X' * X) \ (X' * y);
    ECV = X * beta;
end

if strcmp(method, 'Monomials')
    Lterms = [ones(size(x)), x, x.^2, x.^3];
    X = Lterms(:, 1:kterm);
    beta = (X' * X) \ (X' * y);
    ECV = X * beta;
end

end