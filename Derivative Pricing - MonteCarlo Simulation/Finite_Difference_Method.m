
clear; close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Question 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% input parameters
T = 0.5; K = 10;
sigma = 0.2;
dt = 0.002;
dXValues = [sigma*sqrt(dt), sigma*sqrt(3*dt), sigma*sqrt(4*dt)];

S0Values = 4:16;
Price_BS = zeros(size(S0Values));
for s0 = S0Values
   Price_BS(S0Values == s0) = BSOptionPrice('Put', s0, K, 0.04, sigma, T, 0); 
end

for i = 1:3
dX = dXValues(i);
Pa = Proj7_Q1(K, sigma, T, dt, dX, 'EFD');
Pb = Proj7_Q1(K, sigma, T, dt, dX, 'IFD');
Pc = Proj7_Q1(K, sigma, T, dt, dX, 'C-NFD');
Result = array2table([Pa; Pb; Pc],...
    'RowNames', {'Pa', 'Pb', 'Pc'},...
    'VariableNames', {'4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'});
eval(['ResultQ1_', num2str(i), '=','Result', ';']);

ErrorTable = array2table([Pa-Price_BS; Pb-Price_BS; Pc-Price_BS],...
    'RowNames', {'EFD', 'IFD', 'C-NFD'},...
    'VariableNames', {'4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'});
eval(['ErrorQ1_', num2str(i), '=','ErrorTable', ';']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Question 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% input parameters
T = 0.5; K = 10;
sigma = 0.2;
dt = 0.002;

Pa = Proj7_Q2(K, sigma, T, dt, 'Put', 'EFD');
Pb = Proj7_Q2(K, sigma, T, dt, 'Put', 'IFD');
Pc = Proj7_Q2(K, sigma, T, dt, 'Put', 'C-NFD');
Ca = Proj7_Q2(K, sigma, T, dt, 'Call', 'EFD');
Cb = Proj7_Q2(K, sigma, T, dt, 'Call', 'IFD');
Cc = Proj7_Q2(K, sigma, T, dt, 'Call', 'C-NFD');
Result_Q2 = array2table([Ca; Cb; Cc; Pa; Pb; Pc],...
    'RowNames', {'Ca', 'Cb', 'Cc', 'Pa', 'Pb', 'Pc'},...
    'VariableNames', {'4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'});

%%%% Call Option Figure
figure;
plot(4:16, Ca, 'yo-', 'LineWidth', 4); hold on;
plot(4:16, Cb, 'c*--', 'LineWidth', 2); hold on;
plot(4:16, Cc, 'r^:');
xlabel('Stock Price'); ylabel('Call Option Price'); legend('EFD', 'IFD', 'C-NFD', 'Location', 'northwest');
%%%% Put Option Figure
figure;
plot(4:16, Pa, 'yo-', 'LineWidth', 4); hold on;
plot(4:16, Pb, 'c*--', 'LineWidth', 2); hold on;
plot(4:16, Pc, 'r^:');
xlabel('Stock Price'); ylabel('Put Option Price'); legend('EFD', 'IFD', 'C-NFD', 'Location', 'northeast');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Price = Proj7_Q1(K, sigma, T, dt, dX, method)

%%%% default setting parameters
r = 0.04; 

StockMax = 16*(1+r+3*sigma); StockMin = 4*(1+r-3*sigma);
LogStockSeq = [log(StockMax):-dX:log(StockMin), log(StockMin)];
StockSeq = exp(LogStockSeq); 

OptionGrid = zeros(length(StockSeq), T/dt+1);
[nrow, ncol] = size(OptionGrid);
OptionGrid(:, end) = max(K - StockSeq, 0);
OptionGrid(end, :) = K * exp(-r*dt*((ncol-1) - (0:(ncol-1))))-StockMin;
OptionGrid(1, :) = 0;


if strcmp(method, 'EFD')
    %%%% Explicit Finite-Difference Method
    p_u = dt* ( (sigma^2)/(2*(dX^2)) + (r - ((sigma^2)/2))/(2*dX) );
    p_m = 1 - dt * ((sigma^2)/(dX^2)) - r*dt;
    p_d = dt* ( (sigma^2)/(2*(dX^2)) - (r - ((sigma^2)/2))/(2*dX) );

    A = zeros(nrow-2, nrow);
    A_0 = [p_u, p_m, p_d, zeros(1, nrow-3)];
    for i = 1:(nrow-2)
        A(i, :) = circshift(A_0, i-1);
    end
    
    for i = (ncol-1):-1:1
        OptionGrid(2:(end-1), i) = A * OptionGrid(:, i+1);
    end

end


if strcmp(method, 'IFD')
    %%%% Implicit Finite-Difference Method
    p_u = (-1/2) * dt* ( (sigma^2)/(dX^2) + (r - ((sigma^2)/2))/dX );
    p_m = 1 + dt * ((sigma^2)/(dX^2)) + r*dt;
    p_d = (-1/2) * dt* ( (sigma^2)/(dX^2) - (r - ((sigma^2)/2))/dX );

    A = zeros(nrow, nrow);
    A_0 = [p_u, p_m, p_d, zeros(1, nrow-3)];
    A(1, :) = [1, zeros(1, nrow-1)]; A(end, :) = [zeros(1, nrow-1), 1];
    for i = 2:(nrow-1)
        A(i, :) = circshift(A_0, i-2);
    end

    B = zeros(nrow, ncol-1);
    B(:, end) = OptionGrid(:, ncol);
    B(1, :) = 0;
    B(end, :) = OptionGrid(end, 1:(end-1));
    for i = ncol-1:-1:1
        OptionGrid(:, i) = A \ B(:, i);
        if i > 1
            B(2:(end-1), i-1) = OptionGrid(2:(end-1), i);
        end
    end

end


if strcmp(method, 'C-NFD')
    %%%% Crank-Nicolson Finite-Difference Method
    p_u = -1/4 * dt* ( (sigma^2)/(dX^2) + (r - ((sigma^2)/2))/dX );
    p_m = 1 + dt * ((sigma^2)/(2*(dX^2))) + r*dt/2;
    p_d = -1/4 * dt* ( (sigma^2)/(dX^2) - (r - ((sigma^2)/2))/dX );

    A = zeros(nrow, nrow);
    A_0 = [p_u, p_m, p_d, zeros(1, nrow-3)];
    A(1, :) = [1, zeros(1, nrow-1)]; A(end, :) = [zeros(1, nrow-1), 1];
    for i = 2:(nrow-1)
        A(i, :) = circshift(A_0, i-2);
    end

    Z = zeros(nrow, ncol-1);
    Z_A = zeros(nrow, nrow);
    Z_A_0 = [p_u, p_m-2, p_d, zeros(1, nrow-3)];
    for i = 2:(nrow-1)
        Z_A(i, :) = -circshift(Z_A_0, i-2);
    end
    Z_B = [zeros(nrow-1, ncol-1); OptionGrid(end, 1:(end-1))];
    Z(:, end) = Z_A * OptionGrid(:, end) + Z_B(:, end);

    for i = ncol-1:-1:1
        OptionGrid(:, i) = A \ Z(:, i);
        if i > 1
            Z(:, i-1) = Z_A * OptionGrid(:, i) + Z_B(:,i-1);
        end
    end
    
end

Price = interp1(StockSeq, OptionGrid(:, 1), 4:16);

end


function Price = BSOptionPrice(type, S0, K, r, sigma, T, delta)

d1 = (log(S0 ./ K) + (r - delta + (1/2)*(sigma .^ 2))*T ) ./ (sigma .* sqrt(T));
d2 = d1 - sigma * sqrt(T);

if strcmp(type, 'Call')
    Price = S0*exp(-delta*T) .* normcdf(d1, 0, 1) - K*exp(-r*T) .* normcdf(d2, 0, 1);
end

if strcmp(type, 'Put')
    Price = -(S0*exp(-delta*T) .* normcdf(-d1, 0, 1) - K*exp(-r*T) .* normcdf(-d2, 0, 1));
end

end



function Price = Proj7_Q2(K, sigma, T, dt, type, method)

%%%% default setting parameters
r = 0.04; 
dS = 1; 

StockMax = ceil(16*(1+r+3*sigma)); StockMin = 0;
StockSeq = StockMax:-dS:StockMin;

OptionGrid = zeros(length(StockSeq), T/dt+1);
[nrow, ncol] = size(OptionGrid);

if strcmp(type, 'Call')
    OptionGrid(:, end) = max(StockSeq - K, 0);
    OptionGrid(1, :) = StockMax - K * exp(-r*dt*((ncol-1) - (0:(ncol-1))));
    OptionGrid(end, :) = 0;
end

if strcmp(type, 'Put')
    OptionGrid(:, end) = max(K - StockSeq, 0);
    OptionGrid(1, :) = 0;
    OptionGrid(end, :) = K * exp(-r*dt*((ncol-1) - (0:(ncol-1))))-StockMin;
end

if strcmp(method, 'EFD')
    %%%% Explicit Finite-Difference Method
    j_param = (length(StockSeq)-2):-1:1; % when StockMin == 0
    p_u = dt* ( (r .* j_param) ./ 2 + (sigma^2) .* (j_param.^2) ./ 2 );
    p_m = 1 - dt * ( (sigma^2) .* (j_param.^2) + r );
    p_d = dt* ( -(r * j_param) ./ 2 + (sigma^2) * (j_param.^2) ./ 2 );

    A = [diag(p_u,0),zeros(length(j_param), 2)] +...
        [zeros(length(j_param), 1), diag(p_m,0),zeros(length(j_param), 1)] +...
        [zeros(length(j_param), 2), diag(p_d,0)];

    for i = (ncol-1):-1:1
        OptionGrid(2:(end-1), i) = A * OptionGrid(:, i+1);
        %----------American Put-----------
        if strcmp(type, 'Put')
            OptionGrid(:, i) = max([OptionGrid(:, i)'; K - StockSeq]);
        end
        %----------American Call----------
        if strcmp(type, 'Call')
            OptionGrid(:, i) = max([OptionGrid(:, i)'; StockSeq - K]);
        end
    end

end


if strcmp(method, 'IFD')
    %%%% Implicit Finite-Difference Method
    j_param = (length(StockSeq)-2):-1:1; % when StockMin == 0
    p_u = (-1/2) * dt* ( (sigma^2)*(j_param.^2) + r*j_param );
    p_m = 1 + dt * (sigma^2) * (j_param.^2) + r*dt;
    p_d = (-1/2) * dt* ( (sigma^2)*(j_param.^2) - r*j_param );

    A = [diag(p_u,0),zeros(length(j_param), 2)] +...
        [zeros(length(j_param), 1), diag(p_m,0),zeros(length(j_param), 1)] +...
        [zeros(length(j_param), 2), diag(p_d,0)];
    B = zeros(nrow, ncol-1);
    B(:, end) = OptionGrid(:, ncol);
    B(1, :) = 0;
    B(end, :) = OptionGrid(end, 1:(end-1));
    for i = ncol-1:-1:1
        OptionGrid(:, i) = A \ B(2:(end-1), i);

        if strcmp(type, 'Put')
            %------- Put option---------
            OptionGrid(1, i) = 0;
            OptionGrid(end, i) = K * exp(-r*dt*(ncol - i))-StockMin;
            %----------American Put-----------
            OptionGrid(:, i) = max([OptionGrid(:, i)'; K - StockSeq]);
        end

        if strcmp(type, 'Call')
            %------- Call option---------
            OptionGrid(1, i) = StockMax - K * exp(-r*dt*(ncol - i));
            OptionGrid(end, i) = 0;
            %----------American Call-----------
            OptionGrid(:, i) = max([OptionGrid(:, i)'; StockSeq - K]);
        end

        if i > 1
            B(2:(end-1), i-1) = OptionGrid(2:(end-1), i);
        end

    end
    
end


if strcmp(method, 'C-NFD')
    %%%% Crank-Nicolson Finite-Difference Method
    j_param = (length(StockSeq)-2):-1:1; % when StockMin == 0
    p_u = -1/4 * dt* ( (sigma^2)*(j_param.^2) + r*j_param );
    p_m = 1 + dt * (sigma^2) * (j_param.^2) ./ 2 + r*dt/2;
    p_d = -1/4 * dt* ( (sigma^2)*(j_param.^2) - r*j_param );

    A = [diag(p_u,0),zeros(length(j_param), 2)] +...
        [zeros(length(j_param), 1), diag(p_m,0),zeros(length(j_param), 1)] +...
        [zeros(length(j_param), 2), diag(p_d,0)];
    Z = zeros(nrow-2, ncol-1);
    Z_A = -([diag(p_u,0),zeros(length(j_param), 2)] +...
        [zeros(length(j_param), 1), diag(p_m-2,0),zeros(length(j_param), 1)] +...
        [zeros(length(j_param), 2), diag(p_d,0)]);
    Z(:, end) = Z_A * OptionGrid(:, end);

    for i = ncol-1:-1:1

        OptionGrid(:, i) = A \ Z(:, i);

        if strcmp(type, 'Put')
            %------- Put option---------
            OptionGrid(1, i) = 0;
            OptionGrid(end, i) = K * exp(-r*dt*(ncol - i))-StockMin;
            %----------American Put-----------
            OptionGrid(:, i) = max([OptionGrid(:, i)'; K - StockSeq]);
        end

        if strcmp(type, 'Call')
            %------- Call option---------
            OptionGrid(1, i) = StockMax - K * exp(-r*dt*(ncol - i));
            OptionGrid(end, i) = 0;
            %----------American Call-----------
            OptionGrid(:, i) = max([OptionGrid(:, i)'; StockSeq - K]);
        end

        if i > 1
            Z(:, i-1) = Z_A * OptionGrid(:, i);
        end

    end
end

Price = interp1(StockSeq, OptionGrid(:, 1), 4:16);

end