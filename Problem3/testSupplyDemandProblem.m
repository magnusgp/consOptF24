clear
clc
close all

% Load the problem data from the file LP_Test.mat
load('LP_Test.mat'); % This will load C, U, Pd_max, Pg_max

% The objective function coefficients need to be constructed based on your U and C vectors
% Since we are maximizing, we negate the objective function for linprog
g = [-U; C]; 

% The equality constraint (power balancing) seems to be the sum of p_d - sum of p_g = 0
% Assuming the number of p_d variables is equal to the length of U and the same for p_g and C
A = [ones(1, length(U)), -ones(1, length(C))];
b = 0;

% Lower bounds are zeros for both p_d and p_g, upper bounds are Pd_max and Pg_max respectively
lb = [zeros(length(U), 1); zeros(length(C), 1)]; 
ub = [Pd_max; Pg_max];

% Solve the linear program using the dual-simplex algorithm
options = optimoptions('linprog','Algorithm','dual-simplex');
tic;
[x, fval, exitflag, output, lambda] = linprog(g, [], [], A, b, lb, ub, options);
t = toc;
% Extract the solution for p_d and p_g from the result x
p_d = x(1:length(U));
p_g = x(length(U)+1:end);

% get number of iterations 
num_iter = output.iterations;

% Display the results
disp('The optimal demand prices, p_d, are:');
disp(p_d);
disp('The optimal generation prices, p_g, are:');
disp(p_g);

% display cpu time and number of iterations
disp('CPU time: ');
disp(t);
disp('Number of iterations: ');
disp(num_iter);

% Post-processing steps for supply-demand curve and market clearing price can be added here
% sum p_d and p_g to check if market is cleared
assert(sum(p_d) == sum(p_g), 'Market not cleared');
% display 
disp('sum of p_d: ');
disp(sum(p_d));
disp('sum of p_g: ');
disp(sum(p_g));

% determine the market clearing price which is where p_d and p_g intersect 
% sort p_d in descending order
[sorted_p_d, ~] = sort(p_d, 'descend');
cumulative_demand = cumsum(ones(size(sorted_p_d)));

% sort p_g in ascending order
[sorted_p_g, ~] = sort(p_g, 'ascend');
cumulative_supply = cumsum(ones(size(sorted_p_g)));

% Plotting the stairs plot for supply
figure;
stairs(cumulative_supply, sorted_p_g, 'LineWidth', 2);
hold on;

% Plotting the stairs plot for demand
stairs(cumulative_demand, sorted_p_d, 'LineWidth', 2);

% Adding labels and title
xlabel('Quantity');
ylabel('Price');
title('Supply-Demand Curve linprog');
legend('Supply', 'Demand');
grid on

% find the intersection of the two curves
% find the index of the first element in sorted_p_d that is greater than sorted_p_g just look in the first 10 elements
% this is the market clearing price
for i = 1:15
    if sorted_p_d(i) < sorted_p_g(i)
        market_clearing_price = sorted_p_d(i);
        break;
    end
end

disp('Market clearing price: ');
disp(market_clearing_price);


% now solve with function [x,info,mu,lambda,iter] = LPippd(g,A,b,x) 
% where g = f, A = [A_ineq; A_eq], b = [b_ineq; b_eq], x = [p_d; p_g]

[m, n] = size(A);

Abar = [A -A zeros(m,n) zeros(m,n);
        eye(n) -eye(n) -eye(n) zeros(n);
        eye(n) -eye(n) zeros(n) eye(n)];

bbar = [b;lb;ub];

gbar = [g;-g;zeros(n,1);zeros(n,1)];

xtilde = Abar'*inv(Abar*Abar')*bbar;
lambdatilde = inv(Abar*Abar')*Abar*gbar;
stilde = gbar - Abar'*lambdatilde;

deltax = max([-(3/2)*min(xtilde),0]);
deltas = max([-(3/2)*min(stilde),0]);

e = ones(length(xtilde),1);

xhat = xtilde + deltax*e;
shat = stilde + deltas*e;

deltaxhat = 0.5*(xhat'*shat)/(e'*shat);
deltashat = 0.5*(xhat'*shat)/(e'*xhat);

x0 = xhat + deltaxhat*e;
lambda0 = lambdatilde;
s0 = shat + deltashat*e;

tol = 1.0e-9;
[xout,lambda,s,iter,info,rcres,rbres,mures] = lpsolverInteriorPoint(gbar,Abar,bbar,x0,lambda0,s0,tol);

xhat = xout(1:n)-xout((n+1):(2*n));

disp(max(abs(xhat-x)))

disp("IP")
disp(iter)
disp(g'*xhat)
disp(A*xhat)
disp(prod(xhat>=0))

disp("linprog")
disp(g'*x)
disp(A*x)
disp(prod(x>=0))

figure
subplot(1,3,1)
plot(vecnorm(rcres,2),'-o',LineWidth=1.5)
xlabel("Iterations")
ylabel("$||r_{c}||_2$")
grid on
title(sprintf("$r_{c}$"))

subplot(1,3,2)
plot(vecnorm(rbres,2,1),'-o',LineWidth=1.5)
xlabel("Iterations")
ylabel("$||r_{b}||_2$")
grid on
title("$r_b$")

subplot(1,3,3)
plot(mures,'-o',LineWidth=1.5)
xlabel("Iterations")
ylabel("$\mu$")
grid on
title("$\mu$")
sgtitle("Residuals w.r.t. iterations")

% Extract the solution for p_d and p_g from the result x
p_d = xhat(1:length(U));
p_g = xhat(length(U)+1:end);

[sorted_p_d, ~] = sort(p_d, 'descend');
cumulative_demand = cumsum(ones(size(sorted_p_d)));

% sort p_g in ascending order
[sorted_p_g, ~] = sort(p_g, 'ascend');
cumulative_supply = cumsum(ones(size(sorted_p_g)));

% Plotting the stairs plot for supply
figure;
stairs(cumulative_supply, sorted_p_g, 'LineWidth', 2);
hold on;

% Plotting the stairs plot for demand
stairs(cumulative_demand, sorted_p_d, 'LineWidth', 2);

% Adding labels and title
xlabel('Quantity');
ylabel('Price');
title('Supply-Demand Curve Interior point');
legend('Supply', 'Demand');
grid on
