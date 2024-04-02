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
    if ~isequal(H, H') || any(eig(H) <= 0)
        error('H must be a positive definite matrix');
    end
    
    % Check dimensions
    [n, ~] = size(A);
    if size(H, 1) ~= n || size(H, 2) ~= n || numel(g) ~= n || numel(b) ~= size(A', 1)
        error('Dimensions of inputs are inconsistent');
    end
    
    % Formulate the KKT system
    KKT_matrix = [H, A; A', zeros(size(A', 1))];
    rhs = [-g; b];
    
    % Solve the system using dense LU factorization
    [L, U, P] = lu(KKT_matrix);
    y = L \ (P .* rhs);
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
    if ~isequal(H, H') || any(eig(H) <= 0)
        error('H must be a positive definite matrix');
    end
    
    % Check dimensions
    [n, ~] = size(A);
    if size(H, 1) ~= n || size(H, 2) ~= n || numel(g) ~= n || numel(b) ~= size(A', 1)
        error('Dimensions of inputs are inconsistent');
    end
    
    % Formulate the KKT system
    KKT_matrix = [H, A; A', zeros(size(A', 1))];
    rhs = [-g; b];
    
    % Convert KKT matrix to sparse format
    KKT_sparse = sparse(KKT_matrix);
    
    % Solve the system using sparse LU factorization
    [L, U, P, Q] = lu(KKT_sparse);
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
    if ~isequal(H, H') || any(eig(H) <= 0)
        error('H must be a positive definite matrix');
    end
    
    % Check dimensions
    [n, ~] = size(A);
    if size(H, 1) ~= n || size(H, 2) ~= n || numel(g) ~= n || numel(b) ~= size(A', 1)
        error('Dimensions of inputs are inconsistent');
    end
    
    % Formulate the KKT system
    KKT_matrix = [H, A; A', zeros(size(A', 1))];
    rhs = [-g; b];
    
    % Solve the system using dense LDL factorization
    [L, D, P] = ldl(KKT_matrix);
    y = P * (L' \ (D \ (L \ (P' * rhs))));
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
    if ~isequal(H, H') || any(eig(H) <= 0)
        error('H must be a positive definite matrix');
    end
    
    % Check dimensions
    [n, ~] = size(A);
    if size(H, 1) ~= n || size(H, 2) ~= n || numel(g) ~= n || numel(b) ~= size(A', 1)
        error('Dimensions of inputs are inconsistent');
    end
    
    % Formulate the KKT system
    KKT_matrix = [H, A; A', zeros(size(A', 1))];
    rhs = [-g; b];
    
    % Compute sparse Cholesky factorization
    [L, p, s] = chol(KKT_matrix, 'lower', 'vector');
    
    % Check for factorization failure
    if s ~= 0
        error('Cholesky factorization failed');
    end
    
    % Solve the system using sparse forward and backward substitution
    y = L'\(L\rhs(p));
    sol = y;
    
    % Extract solution and Lagrange multiplier
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
    if ~isequal(H, H') || any(eig(H) <= 0)
        error('H must be a positive definite matrix');
    end
    
    % Check dimensions
    [n, ~] = size(A);
    if size(H, 1) ~= n || size(H, 2) ~= n || numel(g) ~= n || numel(b) ~= size(A', 1)
        error('Dimensions of inputs are inconsistent');
    end
    
    % Formulate the KKT system
    KKT_matrix = [H, A; A', zeros(size(A', 1))];
    rhs = [-g; b];
    
    % Compute the range-space decomposition of KKT_matrix
    [U, S, V] = svd(KKT_matrix);
    
    % Determine the rank of the KKT_matrix
    rank_KKT = sum(diag(S) > eps(S(1)) * max(size(S)));
    
    % Extract relevant matrices from the decomposition
    U1 = U(:, 1:rank_KKT);
    V1 = V(:, 1:rank_KKT);
    S1_inv = diag(1./diag(S(1:rank_KKT, 1:rank_KKT)));
    
    % Solve the system using the range-space factorization
    sol = V1 * S1_inv * U1' * rhs;
    
    % Extract solution and Lagrange multiplier
    x = sol(1:n);
    lambda = sol(n+1:end);
    
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
    if ~isequal(H, H') || any(eig(H) <= 0)
        error('H must be a positive definite matrix');
    end
    
    % Check dimensions
    [n, ~] = size(A);
    if size(H, 1) ~= n || size(H, 2) ~= n || numel(g) ~= n || numel(b) ~= size(A', 1)
        error('Dimensions of inputs are inconsistent');
    end
    
    % Formulate the KKT system
    KKT_matrix = [H, A; A', zeros(size(A', 1))];
    rhs = [-g; b];
    
    % Compute the null space of KKT_matrix
    [~, ~, V] = svd(KKT_matrix');
    null_space_basis = V(:, end - rank(KKT_matrix) + 1:end);
    
    % Compute the null-space factorization
    null_space_factor = null_space_basis * null_space_basis';
    
    % Solve the system using the null-space factorization
    sol = eye(size(KKT_matrix)) - null_space_factor;
    sol = sol \ rhs;
    
    % Extract solution and Lagrange multiplier
    x = sol(1:n);
    lambda = sol(n+1:end);
    
    end