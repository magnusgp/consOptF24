function [x_opt, f_opt] = sqp_line_search(x0, lambda0, mu0)
    % Objective function
    objective = @(x) (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2;

    % Inequality constraint function
    constraint = @(x) [(x(1) + 2).^2 - x(2); -4*x(1) + 10*x(2)];

    % Lagrangian function
    lagrangian = @(x, lambda, mu) objective(x) + lambda'*constraint(x) + mu'*max(0, constraint(x));

    % Optimization options
    options = optimoptions('fminunc', 'Display', 'off');

    % Maximum number of iterations
    max_iter = 100;

    % Tolerance for convergence
    tol = 1e-6;

    % Initial point
    x = x0;
    lambda = lambda0;
    mu = mu0;

    for iter = 1:max_iter
        % Compute gradients
        grad = @(x) gradient(lagrangian, x, lambda, mu);
        g = grad(x);

        % Compute Hessian
        hess = @(x) hessian(lagrangian, x, lambda, mu);

        % Solve quadratic subproblem
        [p, ~, ~] = fminunc(@(p) g'*p + 0.5*p'*hess(x)*p, zeros(2,1), options);

        % Line-search
        alpha = backtrack_line_search(x, p, lagrangian, lambda, mu, g);

        % Update x
        x = x + alpha * p;

        % Update Lagrange multipliers
        lambda = lambda + alpha * constraint(x);
        mu = max(0, mu + alpha * constraint(x));

        % Check convergence
        if norm(g) < tol
            break;
        end
    end

    x_opt = x;
    f_opt = objective(x);
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

function hess = hessian(func, x, lambda, mu)
    % Numerical approximation of Hessian
    epsilon = 1e-6;
    n = length(x);
    hess = zeros(n, n);
    for i = 1:n
        for j = 1:n
            x_plus_i = x;
            x_plus_j = x;
            x_plus_i(i) = x_plus_i(i) + epsilon;
            x_plus_j(j) = x_plus_j(j) + epsilon;
            hess(i, j) = (func(x_plus_i, lambda, mu) - func(x, lambda, mu) - ...
                func(x_plus_j, lambda, mu) + func(x, lambda, mu)) / epsilon^2;
        end
    end
end
