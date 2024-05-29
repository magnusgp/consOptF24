function alpha_k = line_search2(xk, fxk, grad_fxk, hxk, gxk, lambda, mu, gradxk)
    % Perform a line search to find step size alpha_k
    alpha_k = 1;
    max_line_search_iters = 20;
    c = fxk + lambda*abs(fxk) + mu*abs(min(0, gxk));
    b = fxk.'*gradxk - lambda.'*abs(fxk).*abs(gradxk) - mu.'*abs(min(0, gxk));

    for i = 1:max_line_search_iters
        x_new = xk + alpha_k * gradxk;
        [f_new, ~] = objective(x_new);
        [c_eq_new, c_ineq_new] = constraints(x_new);

        phialpha = f_new + lambda*abs(f_new) + mu*abs(min(0, c_eq_new));
        
        % Check Wolfe conditions
        if phialpha <= c + alpha_k*b
            break;
        end
        
        if f_new <= fxk + sigma * alpha_k * grad_fxk' * gradxk && ...
                all(c_eq_new >= 0) && all(c_ineq_new >= 0)
            break;
        end

        a = (phialpha - (c + b*alpha_k))/(alpha_k^2);
        alpha_k_min = -b/(2*a);
        alpha_k = min(0.9*alpha_k, max(0.1*alpha_k, alpha_k_min));
    end
end