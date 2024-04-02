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
% WORKING SOLVERS: LUsparse, LDLdense, range-space
% MISSING SOLVERS: LUdense, LDLsparse, null-space
solvers = {'LUsparse', 'LDLdense'};
% solvers = {'LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'range-space', 'null-space'};

% TASK 1.4/1.5

% Define the problem sizes and beta values
n_values = 10:10:100;
beta_values = 0.1:0.1:1.0;

% Initialize the solution and runtime matrices
x_solutions = zeros(length(n_values), length(beta_values), length(solvers));
runtimes = zeros(length(n_values), length(beta_values), length(solvers));

% Loop over the problem sizes
for i = 1:length(n_values)
    % Loop over the beta values
    for j = 1:length(beta_values)
        % Loop over the solvers
        for k = 1:length(solvers)
            % Start the timer
            tic;
            
            % Call the testProblem function
            [x, lambda] = testQPs(n_values(i), beta_values(j), 0.1, solvers{k});
            
            % Stop the timer and store the runtime
            runtimes(i, j, k) = toc;
            
            % Store the solution
            x_solutions(i, j, k) = x(1);
        end
    end
end

% Plot the solutions and runtimes
for k = 1:length(solvers)
    figure;
    surf(n_values, beta_values, squeeze(x_solutions(:,:,k)));
    title(['Solutions for ', solvers{k}]);
    xlabel('Problem size n');
    ylabel('Number of constraints m = beta * n');
    zlabel('Solution x(1)');
    
    figure;
    surf(n_values, beta_values, squeeze(runtimes(:,:,k)));
    title(['Runtimes for ', solvers{k}]);
    xlabel('Problem size n');
    ylabel('Number of constraints m = beta * n');
    zlabel('Runtime (seconds)');
end


clear
clc

% TASK 1.6
% Define the range for b(1)
b1_range = 8.5:0.1:18.68;

% Initialize the solution matrix
x_solutions = zeros(length(b1_range), length(solvers));

% Loop over the range of b(1)
for i = 1:length(b1_range)
    % Update b(1)
    b(1) = b1_range(i);
    
    % Loop over the solvers
    for j = 1:length(solvers)
        % Solve the problem using the current solver
        [x,lambda] = EqualityQPSolver(H,g,A,b,solvers{j});
        
        % Store the solution
        x_solutions(i,j) = x(1);
    end
end

% Plot the solutions
figure;
hold on;
for j = 1:length(solvers)
    plot(b1_range, x_solutions(:,j), '-o', 'DisplayName', solvers{j});
end
hold off;
legend('Location', 'best');
xlabel('b(1)');
ylabel('x(1)');
title('Solution as a function of b(1)');