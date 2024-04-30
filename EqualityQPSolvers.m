function [x,lambda] = EqualityQPSolver(H,g,A,b,solver)
    switch solver
        case 'LUdense'
            [x,lambda] = EqualityQPSolverLUdense(H,g,A,b);
        case 'LUsparse'
            [x,lambda] = EqualityQPSolverLUsparse(H,g,A,b);
        case 'LDLdense'
            [x,lambda] = EqualityQPSolverLDLdense(H,g,A,b);
        case 'LDLsparse'
            [x,lambda] = EqualityQPSolverLDLsparse(H,g,A,b);
        case 'range-space'
            [x,lambda] = EqualityQPSolverRangeSpace(H,g,A,b);
        case 'null-space'
            [x,lambda] = EqualityQPSolverNullSpace(H,g,A,b);
        otherwise
            error('Invalid solver');
    end
end

% Maybe remove checking eigenvalues to save time
% Pretty expensive for large system probably
% Maybe just assume H is symmetric and pos def? 
% Or cheaper way to check

% Define each solver function here
function [x, lambda] = EqualityQPSolverLUdense(H, g, A, b)
    % Solves the Equality Constrained Convex QP problem:
    %   min_x phi = 1/2 * x' * H * x + g' * x
    %   subject to: A' * x = b
    %
    % Inputs:
    %   H: Positive definite matrix (n x n)
    %   g: Column vector of size n
    %   A: Matrix defining the constraints (m x n)
    %   b: Column vector of size m
    %
    % Outputs:
    %   x: Solution vector (n x 1)
    %   lambda: Lagrange multiplier (m x 1)
    
    % Check if H is positive definite
    % if ~isequal(H, H') || any(eig(H) <= 0)
    %     error('H must be a positive definite matrix');
    % end
    
    % Check dimensions
    [n, m] = size(A);
    if size(H, 1) ~= n || size(H, 2) ~= n || numel(g) ~= n || numel(b) ~= size(A', 1)
        error('Dimensions of inputs are inconsistent');
    end
    
    % Formulate the KKT system
    KKT_matrix = [H, A; A', zeros(m,m)];
    rhs = [-g; b];

    % disp("LUdense")
    % disp(issparse(KKT_matrix))

    % Solve the system using dense LU factorization
    [L, U, P] = lu(KKT_matrix);
    y = L \ (P * rhs);
    sol = U \ y;
    
    % Extract solution and Lagrange multiplier
    x = sol(1:n);
    lambda = sol(n+1:end);
    
    end    

function [x, lambda] = EqualityQPSolverLUsparse(H, g, A, b)
    % Solves the Equality Constrained Convex QP problem:
    %   min_x phi = 1/2 * x' * H * x + g' * x
    %   subject to: A' * x = b
    %
    % Inputs:
    %   H: Positive definite matrix (n x n)
    %   g: Column vector of size n
    %   A: Matrix defining the constraints (m x n)
    %   b: Column vector of size m
    %
    % Outputs:
    %   x: Solution vector (n x 1)
    %   lambda: Lagrange multiplier (m x 1)
    
    % Check if H is positive definite
    % if ~isequal(H, H') || any(eig(H) <= 0)
    %     error('H must be a positive definite matrix');
    % end
    
    % Check dimensions
    [n, m] = size(A);
    if size(H, 1) ~= n || size(H, 2) ~= n || numel(g) ~= n || numel(b) ~= size(A', 1)
        error('Dimensions of inputs are inconsistent');
    end
    
    % Formulate the KKT system
    KKT_matrix = [H, A; A', sparse(m, m)];
    rhs = [-g; b];

    % disp("LUsparse")
    % disp(issparse(KKT_matrix))

    % Solve the system using sparse LU factorization
    [L, U, P, Q] = lu(KKT_matrix);
    y = Q * (U \ (L \ (P * rhs)));
    sol = y;
    
    % Extract solution and Lagrange multiplier
    x = sol(1:n);
    lambda = sol(n+1:end);
    
    end
        
function [x, lambda] = EqualityQPSolverLDLdense(H, g, A, b)
    % Solves the Equality Constrained Convex QP problem:
    %   min_x phi = 1/2 * x' * H * x + g' * x
    %   subject to: A' * x = b
    %
    % Inputs:
    %   H: Positive definite matrix (n x n)
    %   g: Column vector of size n
    %   A: Matrix defining the constraints (m x n)
    %   b: Column vector of size m
    %
    % Outputs:
    %   x: Solution vector (n x 1)
    %   lambda: Lagrange multiplier (m x 1)
    
    % Check if H is positive definite
    % if ~isequal(H, H') || any(eig(H) <= 0)
    %     error('H must be a positive definite matrix');
    % end
    
    % Check dimensions
    [n, m] = size(A);
    if size(H, 1) ~= n || size(H, 2) ~= n || numel(g) ~= n || numel(b) ~= size(A', 1)
        error('Dimensions of inputs are inconsistent');
    end
    
    % Formulate the KKT system
    KKT_matrix = [H, A; A', zeros(m, m)];
    rhs = [-g; b];

    % disp("LDLdense")
    % disp(issparse(KKT_matrix))

    % Solve the system using dense LDL factorization
    [L, D] = ldl(KKT_matrix);
    y = (L' \ (D \ (L \ rhs)));
    sol = y;
    
    % Extract solution and Lagrange multiplier
    x = sol(1:n);
    lambda = sol(n+1:end);
    
    end            

function [x, lambda] = EqualityQPSolverLDLsparse(H, g, A, b)
    % Solves the Equality Constrained Convex QP problem:
    %   min_x phi = 1/2 * x' * H * x + g' * x
    %   subject to: A' * x = b
    %
    % Inputs:
    %   H: Positive definite matrix (n x n)
    %   g: Column vector of size n
    %   A: Matrix defining the constraints (m x n)
    %   b: Column vector of size m
    %
    % Outputs:
    %   x: Solution vector (n x 1)
    %   lambda: Lagrange multiplier (m x 1)
    
    % Check if H is positive definite
    % if ~isequal(H, H') || any(eig(H) <= 0)
    %     error('H must be a positive definite matrix');
    % end
    
    % Check dimensions
    [n, m] = size(A);
    if size(H, 1) ~= n || size(H, 2) ~= n || numel(g) ~= n || numel(b) ~= size(A', 1)
        error('Dimensions of inputs are inconsistent');
    end
    
    % Formulate the KKT system
    KKT_matrix = [H, A; A', sparse(m, m)];
    rhs = [-g; b];
    
    % disp("LDLsparse")
    % disp(issparse(KKT_matrix))

    % LDL factorization
    [L, D, P] = ldl(KKT_matrix);

    y = P * (L' \ (D \ (L \ (P' * rhs))));

    sol = y;
    
    % Extract the solution vectors
    x = sol(1:n);
    lambda = sol(n+1:end);
end
    
function [x, lambda] = EqualityQPSolverRangeSpace(H, g, A, b)
    % Solves the Equality Constrained Convex QP problem:
    %   min_x phi = 1/2 * x' * H * x + g' * x
    %   subject to: A' * x = b
    %
    % Inputs:
    %   H: Positive definite matrix (n x n)
    %   g: Column vector of size n
    %   A: Matrix defining the constraints (m x n)
    %   b: Column vector of size m
    %
    % Outputs:
    %   x: Solution vector (n x 1)
    %   lambda: Lagrange multiplier (m x 1)
    
    % Check if H is positive definite
    % if ~isequal(H, H') || any(eig(H) <= 0)
    %     error('H must be a positive definite matrix');
    % end
    
    % Check dimensions
    [n, m] = size(A);
    if size(H, 1) ~= n || size(H, 2) ~= n || numel(g) ~= n || numel(b) ~= size(A', 1)
        error('Dimensions of inputs are inconsistent');
    end

    L = chol(H);

    K = (L')\A;
    w = (L')\g;

    
    
    % % Formulate the KKT system
    % KKT_matrix = [H, A; A', zeros(m, m)];
    % rhs = [-g; b];
    % 
    % % Compute the range-space decomposition of KKT_matrix
    % [U, S, V] = svd(KKT_matrix);
    % 
    % % Determine the rank of the KKT_matrix
    % rank_KKT = sum(diag(S) > eps(S(1)) * max(size(S)));
    % 
    % % Extract relevant matrices from the decomposition
    % U1 = U(:, 1:rank_KKT);
    % V1 = V(:, 1:rank_KKT);
    % S1_inv = diag(1./diag(S(1:rank_KKT, 1:rank_KKT)));
    % 
    % % Solve the system using the range-space factorization
    % sol = V1 * S1_inv * U1' * rhs;
    % 
    % % Extract solution and Lagrange multiplier
    % x = sol(1:n);
    % lambda = sol(n+1:end);
    
    end
        
function [x, lambda] = EqualityQPSolverNullSpace(H, g, A, b)
    % Solves the Equality Constrained Convex QP problem:
    %   min_x phi = 1/2 * x' * H * x + g' * x
    %   subject to: A' * x = b
    %
    % Inputs:
    %   H: Positive definite matrix (n x n)
    %   g: Column vector of size n
    %   A: Matrix defining the constraints (m x n)
    %   b: Column vector of size m
    %
    % Outputs:
    %   x: Solution vector (n x 1)
    %   lambda: Lagrange multiplier (m x 1)
    
    % Check if H is positive definite
    % if ~isequal(H, H') || any(eig(H) <= 0)
    %     error('H must be a positive definite matrix');
    % end
    
    % Check dimensions
    [n, m] = size(A);
    if size(H, 1) ~= n || size(H, 2) ~= n || numel(g) ~= n || numel(b) ~= size(A', 1)
        error('Dimensions of inputs are inconsistent');
    end
    
    % Formulate the KKT system
    % KKT_matrix = [H, A; A', zeros(size(A', 1))];
    % rhs = [-g; b];

    % Step 1: Find the range space and null space of A
    [Q, ~] = qr(A);

    [m,~] = size(A');

    Y = Q(:,1:m); % Range space assuming full row-rank of A'
    Z = Q(:,(m+1):end); % Null space

    % Step 1: From second row in KKT system
    x_y = (A'*Y)\b;

    % Step 2: From first row in KKT system
    R = chol(Z'*H*Z);

    x_z = R\(R'\(-Z'*H*Y*x_y+Z'*g));
    
    % Step 3: Resolve x
    x = Y*x_y + Z*x_z;

    % Step 4: Find lagrange multipliers
    lambda = (Y'*A)\(Y'*H*x+Y'*g);

    % % Investigate dims of R, Q and g
    % 
    % Qinvg = Q' * g;
    % 
    % % Step 2: Solve for x_R
    % x_R = (R' \ (Qinvg));
    % 
    % % Step 3: Solve for x_N
    % x_N = null_space_basis * (null_space_basis' * (b - H * x_R));
    % 
    % % Step 4: Combine x_R and x_N to get the solution x
    % x = x_R + x_N;
    % 
    % % Step 5: Solve for lambda
    % lambda = (H * x + g - A * x);
    
    end

% function to test solvers
function [x, lambda] = testQPs(n, beta, alpha, solver)
    % Calculate m
    m = round(beta * n);
    
    % Generate sparse random matrices A and M
    A = sprandn(n, m, 0.15);
    M = sprandn(n, n, 0.15);
    
    % Generate H
    H = M * M' + alpha * eye(n);
    
    % Generate x and lambda
    x_init = randn(n, 1);
    lambda_init = randn(m, 1);
    
    % Generate g and b
    g = H * x_init + A * lambda_init;
    b = A' * x_init;
    
    % Call the solver
    [x, lambda] = EqualityQPSolver(H, g, A, b, solver);
    
    % Display the solution
    % disp('Solution x:');
    % disp(x_sol);
    % disp('Lagrange multipliers lambda:');
    % disp(lambda_sol);
end

% Test the solvers
n = 100;
beta = 0.5;
alpha = 0.1;

% Test the LUdense solver
disp('Testing LUdense solver');
[x, lambda] = testQPs(n, beta, alpha, 'LUdense');
% print
disp(x);
disp(lambda);

% Test the LUsparse solver
disp('Testing LUsparse solver');
[x, lambda] = testQPs(n, beta, alpha, 'LUsparse');
% print
disp(x);
disp(lambda);