% Initial guess
x0 = [-4; 1];
options.objective = @objective_function;
options.constraints = @constraints_function;
% options.hessian = 'BFGS'; % or @hessian_function for analytical Hessian
options.hessian = @hessian_function;
options.tol = 1e-6;
options.maxIter = 100;
options.lambda0 = [1, 1];

% Call the SQP solver
[x_opt, fval, x_steps, exitflag, output] = SQP_solver(x0, options);

% Display results
disp('Optimal Solution:');
disp(x_opt);
disp('Function Value at Optimum:');
disp(fval);
disp('Exit Flag:');
disp(exitflag);
disp('Output:');
disp(output);

% Plot the objective function
x1 = linspace(-5, 5, 100);
x2 = linspace(-5, 5, 100);
[X1, X2] = meshgrid(x1, x2);
F = (X1.^2 + X2 - 11).^2 + (X1 + X2.^2 - 7).^2;

% Constraints
C1 = (X1 + 2).^2 - X2;   % C1 >= 0 (feasible region is above this curve)
C2 = -4*X1 + 10*X2;    % C2 >= 0 (feasible region is above this curve)

% Plot the objective function contour
figure;
contour(X1, X2, F, 100);
hold on;

% Plot the optimal point
% Ensure you define x_opt somewhere in your script with the optimal values
plot(x_opt(1), x_opt(2), 'ro');

% Plot the constraint boundaries
contour(X1, X2, C1, [0, 0], 'r', 'LineWidth', 2);
contour(X1, X2, C2, [0, 0], 'b', 'LineWidth', 2);

% Plot the iterates
plot(x_steps(1, :), x_steps(2, :), 'k.-', 'MarkerSize', 10, 'LineWidth', 1);

% Shade the infeasible regions
% For C1 < 0 (infeasible region is below the red curve)
infeasible_C1 = C1 == 0;
% For C2 < 0 (infeasible region is below the blue curve)
infeasible_C2 = C2 < 0;
% Combine both infeasible regions
infeasible = infeasible_C1 | infeasible_C2;

% Convert logical matrix to numeric
infeasible_numeric = double(infeasible);

% % Plot the infeasible region
% h = pcolor(X1, X2, infeasible_numeric);
% set(h, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); % Adjust transparency

% % Set colormap for shading (gray color for infeasible regions)
% colormap([1 1 1; 0.8 0.8 0.8]); % White for feasible, gray for infeasible

% Axis labels and title
xlabel('x1');
ylabel('x2');
title('Himmelblau''s Function with Constraints');
hold off;


function [f, grad_f] = objective_function(x)
    % Define your objective function and its gradient here
    f = (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2; % Himmelblau's function
    grad_f = [4*x(1)*(x(1)^2 + x(2) - 11) + 2*(x(1) + x(2)^2 - 7);
              2*(x(1)^2 + x(2) - 11) + 4*x(2)*(x(1) + x(2)^2 - 7)];
end

function [c_eq, c_ineq, grad_ceq, grad_cineq] = constraints_function(x)
    % Constraints and gradients here
    c_eq = [(x(1) + 2).^2 - x(2)]; % Equality constraints
    c_ineq = [-4*x(1) + 10*x(2)]; % Inequality constraints
    grad_ceq = [2*x(1)+4, -1]; % Gradients of equality constraints
    grad_cineq = [-4, 10]; % Gradients of inequality constraints
end

function H = hessian_function(x, lambda)
    % Hessian of the Lagrangian for the Himmmelblau's test problem
    H = [120*x(1)^2 - 40*x(2) + 2, -40*x(1);
         -40*x(1), 120*x(2)^2 - 40*x(1) + 2];
end
