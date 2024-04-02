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
function [x,lambda] = EqualityQPSolverLUdense(H,g,A,b)
    % Perform LU factorization
    [L,U,P] = lu(A);
    
    % Ensure b is a column vector
    b = b(:);

    disp('Size of P:');
    disp(size(P));
    disp('Size of L:');
    disp(size(L));
    disp('Size of U:');
    disp(size(U));
    disp('Size of b:');
    disp(size(b));
    
    % Solve Ly = Pb for y
    y = L \ (P * b);
    
    % Solve Ux = y for x
    x = U \ y;
    
    % Compute lambda
    lambda = H \ (g - A' * x);
end

function [x,lambda] = EqualityQPSolverLUsparse(H,g,A,b)
    % Implement LU sparse factorization here
    [m, n] = size(A);
    [L, U, P] = lu(A);
    y = P * b;
    z = forwardSubstitution(L, y);
    x = backwardSubstitution(U, z);
    lambda = H * x + g;
end

function [x,lambda] = EqualityQPSolverLDLdense(H,g,A,b)
    % Implement LDL dense factorization here
    [m, n] = size(A);
    [L, D] = ldl(A);
    y = forwardSubstitution(L, b);
    z = D \ y;
    x = backwardSubstitution(L', z);
    lambda = H * x + g;
end

function [x,lambda] = EqualityQPSolverLDLsparse(H,g,A,b)
    % Implement LDL sparse factorization here
    [m, n] = size(A);
    [L, D] = ldl(A);
    y = forwardSubstitution(L, b);
    z = D \ y;
    x = backwardSubstitution(L', z);
    lambda = H * x + g;
end

function [x,lambda] = EqualityQPSolverRangeSpace(H,g,A,b)
    % Implement range-space factorization here
    [m, n] = size(A);
    [Q, R] = qr(A');
    Q1 = Q(:, 1:m);
    Q2 = Q(:, m+1:end);
    R1 = R(1:m, :);
    R2 = R(m+1:end, :);
    y = R1 \ b;
    z = R2' \ (H * Q2 * y + g);
    x = Q1 * y + Q2 * z;
    lambda = H * x + g;
end

function [x,lambda] = EqualityQPSolverNullSpace(H,g,A,b)
    % Implement null-space factorization here
    [m, n] = size(A);
    [Q, R] = qr(A);
    Q1 = Q(:, 1:n);
    Q2 = Q(:, n+1:end);
    R1 = R(1:n, :);
    R2 = R(n+1:end, :);
    y = R1 \ b;
    z = R2' \ (H * Q2 * y + g);
    x = Q1 * y + Q2 * z;
    lambda = H * x + g;
end

function z = forwardSubstitution(L, b)
    [m, n] = size(L);
    z = zeros(n, 1);
    for i = 1:n
        z(i) = b(i) / L(i, i);
        for j = i+1:n
            b(j) = b(j) - L(j, i) * z(i);
        end
    end
end

function x = backwardSubstitution(U, z)
    [m, n] = size(U);
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = z(i) / U(i, i);
        for j = i-1:-1:1
            z(j) = z(j) - U(j, i) * x(i);
        end
    end
end