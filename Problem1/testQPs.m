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
function [x_opt, fval, exitflag, output] = SQP_solver(x0, options)
    % SQP_solver implements the SQP algorithm with line search.
    % Inputs:
    % x0 - initial guess
    % options - structure with fields:
    %   'objective' - objective function handle
    %   'constraints' - constraints function handle
    %   'hessian' - Hessian function handle or 'BFGS'
    %   'tol' - tolerance for convergence
    %   'maxIter' - maximum number of iterations
    % Outputs:
    % x_opt - optimal solution
    % fval - value of the objective function at the optimal solution
    % exitflag - convergence flag
    % output - additional output information

    % Initialize variables
    xk = x0;
    lambda_k = options.lambda0;
    Bk = eye(length(x0));  % Initial Hessian approximation (identity matrix)
    tol = options.tol;
    maxIter = options.maxIter;
    objective = options.objective;
    constraints = options.constraints;
    hessian = options.hessian;

    for k = 1:maxIter
        % Evaluate objective, constraints, and their gradients
        [fk, grad_f] = objective(xk);
        [c_eq, c_ineq, grad_ceq, grad_cineq] = constraints(xk);

        % Form the Hessian matrix
        if strcmp(hessian, 'BFGS')
            if k > 1
                sk = xk - xk_prev;
                yk = grad_f - grad_f_prev;
                if sk' * yk > 1e-10 * norm(sk) * norm(yk)
                    Bk = Bk - (Bk * (sk * sk') * Bk) / (sk' * Bk * sk) + (yk * yk') / (yk' * sk);
                end
            end
        else
            Bk = hessian(xk, lambda_k);
        end

        % Solve the QP subproblem
        [pk, lambda_k] = solve_qp_subproblem(grad_f, Bk, c_eq, c_ineq, grad_ceq, grad_cineq);

        % Perform line search
        alpha_k = line_search(xk, pk, objective, constraints, grad_f, c_eq, c_ineq, lambda_k);

        % Update variables
        xk_prev = xk;
        grad_f_prev = grad_f;
        xk = xk + alpha_k * pk;

        % Check for convergence
        if norm(pk) < tol
            exitflag = 1; % Converged
            break;
        end
    end

    % Prepare output
    x_opt = xk;
    fval = objective(xk);
    output.iterations = k;
    output.lambda = lambda_k;
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

