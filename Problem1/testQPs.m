function [x, lambda, t] = testQPs(n, beta, alpha, solver)
    % beta: amount of constraints
    % alpha: scales diagonal of H

    % KKT matrix condition tolerance
    tol = 1e+03;

    % Calculate m
    m = round(beta * n);

    % Sparsity
    s = 0.15;

    % Reciprocal cond
    r = 1;
    
    % Generate sparse random matrices A and M
    A = sprandn(n, m, s, r);
    
    M = sprandn(n, n, s, r);
    
    % Generate H
    H = M * M' + alpha * eye(n);
    
    KKT_matrix = [H A;A' sparse(size(A', 1), size(A, 2))];

    while cond(full(KKT_matrix)) > tol

        % Generate sparse random matrices A and M
        A = sprandn(n, m, s, r);
        
        M = sprandn(n, n, s, r);
        
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

    switch solver
        case 'LUdense'
            H = full(H);
            A = full(A);
        case 'LDLdense'
            H = full(H);
            A = full(A);
        case 'range-space'
            H = full(H);
            A = full(A);
    end
    
    % Call the solver
    tic;
    [x, lambda] = EqualityQPSolver(H, g, A, b, solver);
    t = toc;
    % Display the solution
    % disp('Solution x:');
    % disp(x_sol);
    % disp('Lagrange multipliers lambda:');
    % disp(lambda_sol);
end

