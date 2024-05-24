function [x,Xmat,iter] = lpsolverActiveSet(g,A,b,x0,Bset,Nset,maxiter)
    %          Primal simplex LP solver
    %
    %          min  g'*x
    %           x
    %          s.t. A x  = b      (Lagrange multiplier: mu)
    %                 x >= 0      (Lagrange multiplier: s)
    %
    
    x = x0;

    Xmat = x;

    [n,m] = size(A);

    tol = 1.0e-16;

    % For handling degenerate states
    eps = 0.001;
    perturbed = false;
    bold = b;

    B = A(:,Bset);
    N = A(:,Nset);

    %% Main loop

    % x(Bset) = B\b;
    
    converged = false;
    iter = 0;
    while iter < maxiter && ~converged
        iter = iter + 1;

        mu = (B')\g(Bset);

        % Compute reduced cost
        lambdaN = g(Nset) - N'*mu;

        % Converged if reduced costs all nonnegative
        if lambdaN >= 0
            converged = true;
        else

            s = find(lambdaN < 0,1,'last');

            h = B\A(:,Nset(s));

            if h <= 0
                disp("Unbounded")
                return
            else
                % Find step size
                idx = find(h > 0);
                [alpha,j] = min(x(Bset(idx))./h(idx));
                j = idx(j);

                % Check if xB_i == 0 and h_i < 0 for any i
                degenerate = nnz(intersect(find(x(Bset) < tol), find(h < 0))) > 0;

                % Handle degeneracy by perturbing rhs
                if degenerate && ~perturbed
                    perturbed = true;
                    b = b + B*eps.^(1:n)';
                    x(Bset) = B\b;
                end

                % Compute step
                x(Bset) = x(Bset) - alpha*h;
                x(Bset(j)) = 0;
                x(Nset(s)) = alpha;

                % Change B and N sets
                idxtemp = Bset(j);

                Bset(j) = Nset(s);

                Nset(s) = idxtemp;
            
                B = A(:,Bset);
                N = A(:,Nset);

                Xmat = [Xmat x];

            end
        end
    end

    if perturbed
        x(Bset) = B\bold;
    end

end
