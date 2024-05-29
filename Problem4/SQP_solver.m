% Implements the SQP procedure for solving the Himmelblau Test Problem NLP

function [x_opt, fval, x_steps, exitflag, output] = SQP_solver(x0, options)
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
    % x_steps - history of iterates used for plotting
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
    exitflag = 0;
    x_steps = zeros(length(x0), maxIter);

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
        [pk, fval, exitflag, output, lambda_vec] = quadprog(Bk, grad_f, grad_cineq, -c_ineq, grad_ceq, -c_eq, [], [], [], optimoptions('quadprog', 'Display', 'off'));

        lambda_k = lambda_vec.eqlin;
        mu_k = lambda_vec.ineqlin;

        % Perform line search
        alpha_k = line_search(xk, fk, grad_f, c_eq, c_ineq, grad_ceq, grad_cineq, lambda_k, mu_k, grad_f, objective, constraints);
        % disp(alpha_k);

        % Update variables
        xk_prev = xk;
        grad_f_prev = grad_f;
        xk = xk + alpha_k * pk;

        % Save the iterates for plotting
        x_steps(:, k) = xk;

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

function alpha = backtrack_line_search(x, p, lagrangian, lambda, mu, g)
    % Backtracking line-search
    alpha = 1;
    c1 = 0.1;
    rho = 0.5;
    while lagrangian(x + alpha * p, lambda + alpha * constraint(x), mu + alpha * max(0, constraint(x))) > ...
            lagrangian(x, lambda, mu) + c1 * alpha * g' * p
        alpha = rho * alpha;
    end
end

function grad = gradient(func, x, lambda, mu)
    % Numerical approximation of gradient
    epsilon = 1e-6;
    grad = zeros(length(x), 1);
    for i = 1:length(x)
        x_plus = x;
        x_plus(i) = x_plus(i) + epsilon;
        grad(i) = epsilon./(func(x_plus, lambda, mu) - func(x, lambda, mu));
    end
end