% Load the problem data from the file LP_Test.mat
load('LP_Test.mat'); % This will load C, U, Pd_max, Pg_max

% The objective function coefficients need to be constructed based on your U and C vectors
% Since we are maximizing, we negate the objective function for linprog
f = [-U; C]; 

% There are no inequality constraints apart from bounds in this problem
A_ineq = [];  
b_ineq = [];

% The equality constraint (power balancing) seems to be the sum of p_d - sum of p_g = 0
% Assuming the number of p_d variables is equal to the length of U and the same for p_g and C
A_eq = [ones(1, length(U)), -ones(1, length(C))];
b_eq = 0;

% Lower bounds are zeros for both p_d and p_g, upper bounds are Pd_max and Pg_max respectively
lb = [zeros(length(U), 1); zeros(length(C), 1)]; 
ub = [Pd_max; Pg_max];

% Solve the linear program using the dual-simplex algorithm
options = optimoptions('linprog','Algorithm','dual-simplex');
tic;
[x, fval, exitflag, output, lambda] = linprog(f, A_ineq, b_ineq, A_eq, b_eq, lb, ub, options);
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

% Assuming each p_d and p_g corresponds to a quantity of 1 unit
% Sorting demand prices in descending order for the demand curve
[sorted_p_d, ~] = sort(p_d, 'descend');
cumulative_demand = cumsum(ones(size(sorted_p_d)));

% Sorting generation prices in ascending order for the supply curve
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
title('Supply-Demand Curve');
legend('Supply', 'Demand');

% find the market clearing price i.e. where supply is equal to demand
