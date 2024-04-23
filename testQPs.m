function [x, lambda] = testQPs(n, beta, alpha, solver)
    % Calculate m
    m = round(beta * n);
    
    % Generate sparse random matrices A and M
    A = sprand(n, m, 0.15, 0.9);
    M = sprand(n, n, 0.15, 0.9);
    
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