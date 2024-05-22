function [x,Xmat,iter] = lpsolverActiveSet(g,A,b,x0,Bset,Nset)
    %          Primal simplex LP solver
    %
    %          min  g'*x
    %           x
    %          s.t. A x  = b      (Lagrange multiplier: lambda)
    %                 x >= 0      (Lagrange multiplier: s)
    %
    
    x = x0;

    Xmat = x;

    [n,m] = size(A);

    tol = 1.0e-16;

    % For handling degenerate states
    eps = 0.1;
    perturbed = false;
    bold = b;

    % Bset = find(abs(x) > tol)';
    % 
    % if length(Bset) < n
    % 
    %     additionalidx = find(abs(x) < tol)';
    % 
    %     neededidx = n - length(Bset);
    % 
    %     Bset = [Bset additionalidx(1:neededidx)];
    % 
    % end
    % 
    % Nset = setdiff(1:m,Bset);

    B = A(:,Bset);
    N = A(:,Nset);
    
    %% Main loop
    maxit = 100;
    
    converged = false;
    iter = 0;
    while iter < maxit && ~converged
        iter = iter + 1;

        gB = g(Bset);
        gN = g(Nset);

        x(Bset) = B\b;

        lambda = B'\gB;

        sN = gN - N'*lambda;

        if sN >= 0
            converged = true;
        else

            [~,Nidx] = min(sN);
            qidx = Nset(Nidx);

            d = B\A(:,qidx);

            if d <= 0
                disp("Unbounded")
                return
            else
                idx = find(d > tol);
            
                [xqplus,pidx] = min(x(Bset(idx))./d(idx));

                % disp(max(d))
                % disp(min(x(Bset)))

                % Handle degeneracy
                if nnz(intersect(find(x(Bset) < tol), find(d < 0))) > 0 % && ~perturbed
                    disp("Perturb!")
                    perturbed = true;
                    b = b + B*eps.^(1:n)';
                end

                x(Bset) = x(Bset) - xqplus*d;
                x(Bset(pidx)) = 0;
                x(Nset(Nidx)) = xqplus;

                idxtemp = Bset(pidx);

                Bset(pidx) = Nset(Nidx);

                Nset(Nidx) = idxtemp;
            
                B = A(:,Bset);
                N = A(:,Nset);

                % disp(Bset)

                Xmat = [Xmat x];
        
            end
        end
    end

    if perturbed
        x(Bset) = B\bold;
    end

end
