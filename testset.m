% Parameters
n = 100;           % Size of x, choose according to your requirement
beta = 0.15;       % Adjust beta as per the problem statement
alpha = 1e-5;      % Small constant to ensure positive definiteness
density = 0.15;    % Density of the nonzero elements

% Calculate m
m = round(beta * n);

% Generate x
x = randn(n, 1);

% Generate sparse A with 15% non-zero elements
A = sprandn(n, m, density);

% Generate sparse M with 15% non-zero elements
M = sprandn(n, n, density);

% Ensure M is full for the next operation
M_full = full(M);

% Construct H
H = M_full * M_full' + alpha * eye(n);

% Generate b (using lambda in the problem description)
b = randn(m, 1);

% You now have H, A, x, and b as per the problem statement and can proceed with the solver

