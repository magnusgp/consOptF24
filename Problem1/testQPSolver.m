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

% Defining the solvers
% WORKING SOLVERS: LUsparse, LUdense, LDLdense, LDLsparse, range-space
% MISSING SOLVERS: null-space
% solvers = {'LUsparse', 'LUsparse'};
solvers = {'LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'range-space'};

% % TASK 1.4/1.5

% Define the problem sizes and beta values (start:step:stop)
n_values = 10:100:500;
beta_values = linspace(0.1,1,5);

% Initialize the solution and runtime matrices
% Define the number of iterations
num_iterations = 10;
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
        disp(i)
        % Loop over the beta values
        for j = 1:length(beta_values)
            % Loop over the solvers
            for k = 1:length(solvers)
                % Start the timer
                tic;
                
                if strcmp(solvers{k}, 'benchmark')
                    
                    % Call the standard solver function
                    [x, ~] = EqualityQPSolverLUdense(H, g, A, b);
                    
                    % Stop the timer and add to total runtime
                    total_runtimes(i, j, k) = total_runtimes(i, j, k) + toc;
                else
                    
                    % Call the testProblem function
                    [x, ~] = testQPs(n_values(i), beta_values(j), 100, solvers{k});

                    % Stop the timer and add to total runtime
                    total_runtimes(i, j, k) = total_runtimes(i, j, k) + toc;
                end
            end
        end
    end
    current_iteration = current_iteration + 1;
end

%%

% Calculate average runtimes
average_runtimes = total_runtimes / num_iterations;

figure;
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

figure;
hold on;
for k = 1:length(solvers)
    plot(beta_values, squeeze(mean(average_runtimes(:,:,k), 1)), 'o-');
end
hold off;
xlabel('Beta Values');
ylabel('Average Runtime (s)');
title('Average Runtime vs. Beta Values');
legend(solvers, 'Location', 'best');
grid on;

% TASK 1.6
% % Define the range for b(1)
% b1_range = 8.5:0.1:18.68;

% % Initialize the solution matrix
% x_solutions = zeros(length(b1_range), length(solvers));

% % Loop over the range of b(1)
% for i = 1:length(b1_range)
%     % Update b(1)
%     b(1) = b1_range(i);
    
%     % Loop over the solvers
%     for j = 1:length(solvers)
%         % Solve the problem using the current solver
%         [x,lambda] = EqualityQPSolver(H,g,A,b,solvers{j});

%         % Compute phi and store it
%         phi = 0.5 * x' * H * x + g' * x;
        
%         % Store the solution
%         x_solutions(i,j) = phi;
%     end
% end

% % Plot the solutions
% figure;
% hold on;
% for j = 1:length(solvers)
%     plot(b1_range, x_solutions(:,j), '-o', 'DisplayName', solvers{j});
% end
% hold off;
% legend('Location', 'best', 'Interpreter', 'latex');
% xlabel('$b(1)$', 'Interpreter', 'latex');
% ylabel('$\phi(x)$', 'Interpreter', 'latex');
% title('Solution $\phi(x)$ for $b(1) \in [8.5, 18.68]$ using LUdense solver', 'Interpreter', 'latex');