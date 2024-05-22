% Initial guess
x0 = [0, 0.5];
options.objective = @objective_function;
options.constraints = @constraints_function;
% options.hessian = 'BFGS'; % or @hessian_function for analytical Hessian
options.hessian = @hessian_function;
options.tol = 1e-6;
options.maxIter = 100;
options.lambda0 = [1,1];

% Call the SQP solver
[x_opt, fval, exitflag, output] = SQP_solver(x0, options);

% Display results
disp('Optimal Solution:');
disp(x_opt);
disp('Function Value at Optimum:');
disp(fval);
disp('Exit Flag:');
disp(exitflag);
disp('Output:');
disp(output);

function [f, grad_f] = objective_function(x)
    % Define your objective function and its gradient here
    f = (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2; % Himmelblau's function
    grad_f = [4*x(1)*(x(1)^2 + x(2) - 11) + 2*(x(1) + x(2)^2 - 7);
              2*(x(1)^2 + x(2) - 11) + 4*x(2)*(x(1) + x(2)^2 - 7)];
end

function [c_eq, c_ineq, grad_ceq, grad_cineq] = constraints_function(x)
    % Define your constraints and their gradients here
    c_eq = [];
    c_ineq = []; % Inequality constraints
    grad_ceq = [];
    grad_cineq = []; % Gradients of inequality constraints
end

function H = hessian_function(x, lambda)
    % Define the Hessian of the Lagrangian here for the Himmmelblau's test problem
    H = [120*x(1)^2 - 40*x(2) + 2, -40*x(1);
         -40*x(1), 120*x(2)^2 - 40*x(1) + 2];
end
