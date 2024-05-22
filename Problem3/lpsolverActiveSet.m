function [x,Xmat,iter] = lpsolverActiveSet(g,A,b,x0,Bset,Nset)
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
    eps = 0.01;
    perturbed = false;
    bold = b;

    B = A(:,Bset);
    N = A(:,Nset);

    %% Main loop
    maxit = 200;

    x(Bset) = B\b;
    
    converged = false;
    iter = 0;
    while iter < maxit && ~converged
        iter = iter + 1;

        mu = (B')\g(Bset);

        % Compute reduced cost
        lambdaN = g(Nset) - N'*mu;

        % Converged if reduced costs all nonnegative
        if lambdaN >= 0
            converged = true;
        else

            s = find(lambdaN < 0,1,'first');
            
            is = Nset(s);

            h = B\A(:,is);

            if h <= 0
                disp("Unbounded")
                return
            else
                idx = find(h > tol);
            
                [alpha,j] = min(x(Bset(idx))./h(idx));
                j = idx(j);

                % Is xB_i == 0 and d_i < 0 for any i?
                degenerate = nnz(intersect(find(x(Bset) < tol), find(h < 0))) > 0;

                % Handle degeneracy
                if degenerate && ~perturbed
                    perturbed = true;
                    b = b + B*eps.^(1:n)';
                    x(Bset) = B\b;
                end

                x(Bset) = x(Bset) - alpha*h;
                x(Bset(j)) = 0;
                x(Nset(s)) = alpha;

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
