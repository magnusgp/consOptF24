function alpha_k = line_search(xk, pk, objective, constraints, grad_f, c_eq, c_ineq, lambda_k)
    % Perform a line search to find step size alpha_k
    alpha_k = 1;
    beta = 0.5;
    sigma = 1e-4;
    max_line_search_iters = 20;

    for i = 1:max_line_search_iters
        x_new = xk + alpha_k * pk;
        [f_new, ~] = objective(x_new);
        [c_eq_new, c_ineq_new] = constraints(x_new);
        
        % Check Wolfe conditions
        if f_new <= objective(xk) + sigma * alpha_k * grad_f' * pk && ...
                all(c_eq_new >= 0) && all(c_ineq_new >= 0)
            break;
        end

        alpha_k = beta * alpha_k;
    end
end
