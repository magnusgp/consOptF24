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
        
        % Check Armijo condition
        if f_new <= grad_f' * pk * sigma * alpha_k
            % Check feasibility
            if all(c_eq_new <= 1e-6) && all(c_ineq_new <= 1e-6)
                break;
            end
        end

        alpha_k = beta * alpha_k;
    end
end
