function [pk, lambda] = solve_qp_subproblem(grad_f, Bk, c_eq, c_ineq, grad_ceq, grad_cineq)
    % Solve the QP subproblem
    % Inputs are the gradient of the objective, Hessian approximation,
    % constraints, and their gradients.
    H = Bk;
    f = grad_f;
    Aeq = grad_ceq;
    beq = -c_eq;
    Aineq = [grad_cineq; -grad_cineq];
    bineq = [-c_ineq; c_ineq];

    % Solve the QP problem using quadprog
    options = optimoptions('quadprog', 'Display', 'off');
    [pk, ~, ~, ~, lambda] = quadprog(H, f, Aineq, bineq, Aeq, beq, [], [], [], options);
end
