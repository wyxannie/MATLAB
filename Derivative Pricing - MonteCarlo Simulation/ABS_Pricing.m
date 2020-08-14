
clear; close all;

T = 30; WAC = 0.08; L = 100000;
r0 = 0.078; k = 0.6; rbar = 0.08; sigma = 0.12;
N = 10000; dt = 1/12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Question 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Part a - Compute the price of the MBS

[Price_Q1_a, CFMatrix, rMatrix, ~] = MBSPrice_byCIR(T, WAC, L, N, k, rbar, sigma, r0);
Price_Q1_a

%%%% Part b - Compute and plot
kValues = 0.3:0.1:0.9; 
Price_Q1_b = zeros(size(kValues));
for k = kValues
    [Price_Q1_b(kValues==k), ~] = MBSPrice_byCIR(T, WAC, L, N, k, rbar, sigma, r0);
end

figure; plot(kValues, abs(Price_Q1_b), 'o-');
xlabel("k"); ylabel("Price of MBS");

%%%% Part c - Compute and plot
rbarValues = 0.03:0.01:0.09; k = 0.6;
Price_Q1_c = zeros(size(rbarValues));
IOPrice = zeros(size(rbarValues)); POPrice = zeros(size(rbarValues));
for rbar = rbarValues
    [Price_Q1_c(rbarValues==rbar), ~, ~, IOPrice(rbarValues==rbar), POPrice(rbarValues==rbar)]...
        = MBSPrice_byCIR(T, WAC, L, N, k, rbar, sigma, r0);
end

figure; plot(rbarValues, abs(Price_Q1_c), 'o-');
xlabel("r_{bar}"); ylabel("Price of MBS");


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Question 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Compute OAS
rbar = 0.08; k = 0.6; P0 = 110000;
OAS = real(fsolve(@(x) mean(sum(CFMatrix .* exp(-cumsum(rMatrix(2:end, :)+x)*dt))) - P0, 0))


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Question 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Compute OAS-Adjusted duration and convexity
[Duration, Convexity] = Get_DurationAndConvexity_OAS(OAS, P0, CFMatrix, rMatrix)


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Question 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Price IO and PO
figure; plot(rbarValues, IOPrice, 'o-');
hold on; plot(rbarValues, POPrice, '*-');
xlabel("r_{bar}"); ylabel("Price");
legend('IO', 'PO');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   FUNCTION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Simulate CIR model
function [DiscountFactor, r] = CIR_Model(T, N, k, r_bar, sigma, r0)

% N = 10000;

dt = 1/12; % Simulate monthly rate
Ndt = round(T/dt); dt = T/Ndt;

Z = normrnd(0, 1,[Ndt+1, N]);
r = zeros([Ndt+1, N]); r(1, :) = r0;
for i = 2:(Ndt+1)
    r(i, :) = r(i-1, :) + k*(r_bar - r(i-1, :))*dt + sigma*sqrt(r(i-1, :)).*sqrt(dt).*Z(i, :);
end

DiscountFactor = exp(-cumsum(r(2:end, :))*dt);
end

%%%% Compute 10 tr rate
function r10 = Get_r10(k, sigma, r_bar, r0)

t = 0;
T = 10;

h1 = sqrt(k^2 + 2*(sigma^2));
h2 = (k+h1)/2; 
h3 = (2*k*r_bar)/(sigma^2);
A = ((h1*exp(h2*(T - t))) / (h2*(exp(h1*(T - t))-1)+h1))^h3;
B = (exp(h1*(T - t))-1) / (h2*(exp(h1*(T - t))-1) + h1);
Price_ZCB = A*exp(-B*r0);

r10 = (-1/T)*log(Price_ZCB);

end

%%%% CPR model from Numerix 4-factor
function CPR = CPR_Numerix(R, r10, t, PV_pre, PV0)
% t: Month
% PV_pre: PV_t-1

RI = 0.28 + 0.14*atan(-8.57 + 430*(R - r10));
BU = 0.3 + 0.7*PV_pre/PV0;
SG = min(1, t/30);
SY = [0.98, 0.94, 0.76, 0.74, 0.95, 0.98, 0.92, 0.98, 1.10, 1.18, 1.22, 1.23]; % start with December

CPR = RI.*BU * SG * SY(mod(t, 12)+1);

end

%%%% Compute payment-relevant matrix
function [CashFlowt, PVt, IPt, TPPt] = Get_CtAndPVt(PV_pre, r, N, t, CPRt)
% PV_pre: PV_t-1
% CPRt: CPR_t
% N: #of payment

CashFlowt = (PV_pre*r)/(1-(1+r)^(-N+(t-1))) +...
    (PV_pre - PV_pre*r*(1/(1-(1+r)^(-N+(t-1))) - 1)) .* (1 - (1-CPRt).^(1/12));

IPt = PV_pre*r;
TPPt = CashFlowt - IPt;
PVt = PV_pre - TPPt;

end

%%%% Compute MBS price, IO, PO
function [MBSPrice, CFMatrix, rMatrix, IOPrice, POPrice] = MBSPrice_byCIR(T, WAC, L, N, k, rbar, sigma, r0)

% L: Notional amount of the loan

[DiscountFactor, rMatrix] = CIR_Model(T, N, k, rbar, sigma, r0);
r10Matrix = Get_r10(k, sigma, rbar, rMatrix);

r = WAC/12; R = WAC; PV0 = L;
PVMatrix = zeros(size(rMatrix)); PVMatrix(1, :) = PV0;
CFMatrix = zeros(size(rMatrix)-[1, 0]); % CahFlow and CPR begin at time=1
CPRMatrix = zeros(size(CFMatrix));

IOMatrix = zeros(size(CFMatrix)); POMatrix = zeros(size(CFMatrix));

for t = 1:(12*T) % Monthly payment
    
    % CF and CPR begin at time=1
    % r10 and PV start at time=0
    CPRMatrix(t, :) = CPR_Numerix(R, r10Matrix(t, :), t, PVMatrix(t, :), PV0);
    [CFMatrix(t, :), PVMatrix(t+1, :), IOMatrix(t, :), POMatrix(t, :)] = Get_CtAndPVt(PVMatrix(t, :), r, 30*12, t, CPRMatrix(t, :));
    
end

MBSPrice = abs(mean(sum(DiscountFactor .* CFMatrix)));

IOPrice = abs(mean(sum(DiscountFactor .* IOMatrix)));
POPrice = abs(mean(sum(DiscountFactor .* POMatrix)));

end

%%%% Compute Duration and Convexity
function [Duration, Convexity] = Get_DurationAndConvexity_OAS(OAS, P0, CFMatrix, rMatrix)

y = 5/10000;

    function P = Get_Price(CFMatrix, rMatrix, x)
        dt = 1/12;
        P = mean(sum(CFMatrix .* exp(-cumsum(rMatrix(2:end, :)+x)*dt)));
    end

Duration = real((Get_Price(CFMatrix, rMatrix, OAS-y) - Get_Price(CFMatrix, rMatrix, OAS+y))/(2*y*P0));
Convexity = real((Get_Price(CFMatrix, rMatrix, OAS+y) + Get_Price(CFMatrix, rMatrix, OAS-y) - 2*P0)/(2*(y^2)*P0));

end