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
    KKT_matrix = [H, -A; A', zeros(m,m)];
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
    KKT_matrix = [H, -A; A', sparse(m, m)];
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
    KKT_matrix = [H, -A; A', zeros(m, m)];
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
    KKT_matrix = [H, -A; A', sparse(m, m)];
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

    L = chol(H)';

    K = L\A;
    w = L\g;

    M = chol(K'*K);

    lambda = M\(M'\(K'*w+b));

    x = (L')\(K*lambda-w);
    
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
    
    end