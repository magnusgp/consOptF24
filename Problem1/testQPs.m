function [x, lambda] = testQPs(n, beta, alpha, solver)
    % beta: amount of constraints
    % alpha: scales diagonal of H

    % KKT matrix condition tolerance
    tol = 1e+05;

    % Calculate m
    m = round(beta * n);
    
    % Generate sparse random matrices A and M
    A = sprandn(n, m, 0.15, 0.5);
    
    M = sprandn(n, n, 0.15, 0.5);
    
    % Generate H
    H = M * M' + alpha * eye(n);
    
    KKT_matrix = [H A;A' sparse(size(A', 1), size(A, 2))];

    while cond(full(KKT_matrix)) > tol

        % Generate sparse random matrices A and M
        A = sprandn(n, m, 0.15);
        
        M = sprandn(n, n, 0.15);
        
        % Generate H
        H = M * M' + alpha * eye(n);
        
        KKT_matrix = [H A;A' sparse(size(A', 1), size(A, 2))];

    end

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