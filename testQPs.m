% import solver functions from EqualityQPSolver.m



function [x, lambda] = testQPs(n, beta, alpha, solver)
    % Calculate m
    m = round(beta * n);
    
    % Generate sparse random matrices A and M
    A = sprandn(n, m, 0.15);
    M = sprandn(n, n, 0.15);
    
    % Generate H
    H = M * M' + alpha * eye(n);
    
    % Generate x and lambda
    x_init = randn(n, 1);
    lambda_init = randn(m, 1);
    
    % Generate g and b
    g = H * x_init + A * lambda_init;
    b = A' * x_init;
    
    % Call the solver
    [x, lambda] = EqualityQPSolver(H, g, A, b, solver);
    
    % Display the solution
    % disp('Solution x:');
    % disp(x_sol);
    % disp('Lagrange multipliers lambda:');
    % disp(lambda_sol);
end

% Test the solver
testQPs(100, 0.5, 0.1, 'LUdense')