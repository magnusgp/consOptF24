clear
clc
close all

rng('default') % For reproducibility

% Define the matrices and vectors
H = [5.0000 1.8600 1.2400 1.4800 -0.4600;
     1.8600 3.0000 0.4400 1.1200 0.5200;
     1.2400 0.4400 3.8000 1.5600 -0.5400;
     1.4800 1.1200 1.5600 7.2000 -1.1200;
     -0.4600 0.5200 -0.5400 -1.1200 7.8000];

g = [-16.1000; 
    -8.5000; 
    -15.7000; 
    -10.0200; 
    -18.6800];

A = [16.1000 1.0000; 
    8.5000 1.0000; 
    15.7000 1.0000; 
    10.0200 1.0000; 
    18.6800 1.0000];

b = [15; 
    1];

[x,lambda] = EqualityQPSolver(H,g,A,b,'LUdense');
disp(x)
[x,lambda] = EqualityQPSolver(H,g,A,b,'LUsparse');
disp(x)
[x,lambda] = EqualityQPSolver(H,g,A,b,'LDLdense');
disp(x)
[x,lambda] = EqualityQPSolver(H,g,A,b,'LDLsparse');
disp(x)
[x,lambda] = EqualityQPSolver(H,g,A,b,'range-space');
disp(x)
[x,lambda] = EqualityQPSolver(H,g,A,b,'null-space');
disp(x)

x = quadprog(H,g,[],[],A',b);
disp(x)

%%

% Defining the solvers
solvers = {'LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'range-space', 'null-space'};

% % TASK 1.4/1.5

% Define the problem sizes and beta values (start:step:stop)
n_values = 10:100:800;
% n_values = 50;
beta_values = linspace(0.1,1,5);
% beta_values = 0.5;

% Initialize the solution and runtime matrices
% Define the number of iterations
num_iterations = 8;
current_iteration = 1;

% Initialize matrices to store total runtimes
total_runtimes = zeros(length(n_values), length(beta_values), length(solvers));

% Loop over the number of iterations
for iter = 1:num_iterations
    disp("Current iteration")
    disp(current_iteration)
    % Loop over the problem sizes
    for i = 1:length(n_values)
        disp("n:")
        disp(n_values(i))
        % Loop over the beta values
        for j = 1:length(beta_values)
            % Loop over the solvers
            for k = 1:length(solvers)
                
                % Call the testProblem function
                [x, ~, t] = testQPs(n_values(i), beta_values(j), 10, solvers{k});

                % Add to total runtime
                total_runtimes(i, j, k) = total_runtimes(i, j, k) + t;

                % Benchmark with quadprog
            end
        end
    end
    current_iteration = current_iteration + 1;
end

%%

% Calculate average runtimes
average_runtimes = total_runtimes / num_iterations;

figure;
subplot(1,2,1)
hold on;
for k = 1:length(solvers)
    plot(n_values, squeeze(mean(average_runtimes(:,:,k), 2)), 'o-');
end
hold off;
xlabel('Problem Size (n)');
ylabel('Average Runtime (s)');
title('Average Runtime vs. Problem Size');
legend(solvers, 'Location', 'best');
grid on;

subplot(1,2,2)
hold on;
for k = 1:length(solvers)
    plot(beta_values, squeeze(mean(average_runtimes(:,:,k), 1)), 'o-');
end
hold off;
xlabel('Beta Values');
ylabel('Average Runtime (s)');
title('Average Runtime vs. Beta Values');
grid on;

%%

[n, m] = size(A);

dFdz = [H -A;A' zeros(m,m)];

dFdb = [zeros(length(g),length(b));-eye(length(b),length(b))];

dzdb = -dFdz\dFdb;

b1_range = 8.5:0.1:18.68;

bs = [b1_range;b(2)*ones(1,length(b1_range))];

dz = dzdb*bs;

% TASK 1.6
% Define the range for b(1)
b1_range = 8.5:0.1:18.68;
solver = 'range-space';

% Initialize the solution matrix
solutions = zeros(length(b1_range), length(g)+2);
phi_solutions = zeros(length(b1_range), 1);

btemp = b;
[xOG,lambdaOG] = EqualityQPSolver(H,g,A,b,solver);

disp(dzdb*[-1;0])

% Loop over the range of b(1)
for i = 1:length(b1_range)
    % Update b(1)
    btemp(1) = b1_range(i);
    
    % Solve the problem using the current solver
    [x,lambda] = EqualityQPSolver(H,g,A,btemp,solver);

    % Compute phi and store it
    phi = 0.5 * x' * H * x + g' * x;
    
    % Store the solution
    phi_solutions(i) = phi;

    solutions(i,:) = [x;lambda];
end

% Plot the solutions
figure;
subplot(1,3,1)
plot(b1_range, phi_solutions, '-',LineWidth=1.5);
xlim([min(b1_range),max(b1_range)])
xlabel('$b(1)$', 'Interpreter', 'latex');
ylabel('$\phi(x)$', 'Interpreter', 'latex');
legend('$\phi$')
title("Numerical objective function")

subplot(1,3,2)
plot(b1_range, solutions, '-',LineWidth=1.5);
xlim([min(b1_range),max(b1_range)])
xlabel('$b(1)$', 'Interpreter', 'latex');
ylabel('$x$', 'Interpreter', 'latex');
legend('$x_1$','$x_2$','$x_3$','$x_4$','$x_5$','$\lambda_1$','$\lambda_2$')
title("Numerical solutions")

subplot(1,3,3)
bar({'$x_1$','$x_2$','$x_3$','$x_4$','$x_5$','$\lambda_1$','$\lambda_2$'},dzdb*[1;0])
title("Analytical derivatives")

% subplot(1,3,3)
% plot(b1_range, dz(1:5,:)', '-',LineWidth=1.5);
% xlim([min(b1_range),max(b1_range)])
% xlabel('$b(1)$', 'Interpreter', 'latex');
% ylabel('$x$', 'Interpreter', 'latex');
% legend('$x_1$','$x_2$','$x_3$','$x_4$','$x_5$')
% title("Analytical results")

sgtitle('Objective function $\phi(x)$ and solution for $b(1) \in [8.5, 18.68]$', 'Interpreter', 'latex');