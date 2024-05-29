function alpha_k = line_search(xk, fxk, grad_fxk, hxk, gxk, grad_ceq, grad_cineq, lambda, mu, gradxk, objective, constraints)
    % Perform a line search to find step size alpha_k
    alpha_k = 20;
    max_line_search_iters = 20;
    c = fxk + lambda*abs(fxk) + mu*abs(min(0, gxk));
    b = grad_fxk.'*gradxk - lambda.'*abs(hxk) - mu.'*abs(min(0, gxk));

    for i = 1:max_line_search_iters
        fprintf('Iteration no: %d\n', i);
        x_new = xk + alpha_k * gradxk;
        [f_new, ~] = objective(x_new);
        [c_eq_new, c_ineq_new] = constraints(x_new);

        phialpha = f_new + lambda.'*abs(c_eq_new) + mu.'*abs(min(0, c_ineq_new));
        % disp(phialpha);
        % disp(c + 0.1*alpha_k*b);
        % disp(phialpha <= c + 0.1*alpha_k*b);
        % Check Wolfe conditions
        if phialpha <= c + 0.1*alpha_k*b
            break;
        end
        % fprintf('alpha_k: %f\n', alpha_k);
        % fprintf('c: %f\n', c);
        % fprintf('b: %f\n', b);
        % fprintf('phialpha: %f\n', phialpha);
        a = (phialpha - (c + b*alpha_k))/(alpha_k^2);
        % fprintf('a: %f\n', a);
        alpha_k_min = -b/(2*a);
        fprintf('alpha_k_min: %f\n', alpha_k_min);
        max_var = max(max(0.1*alpha_k, alpha_k_min));
        % disp(max_var);
        % disp(min(0.9*alpha_k, max_var));
        alpha_k = min(0.9*alpha_k, max_var);
    end
end