function [xopt,lambdaopt,X,Wset,it] = lpsolverActiveSet(H,g,A,b,x0,tol)
    %          Primal simplex LP solver
    %
    %          min  g'*x
    %           x
    %          s.t. A x  = b      (Lagrange multiplier: lambda)
    %                 x >= 0      (Lagrange multiplier: s)
    %
    
    % Tolerances
    tolLx = tol;
    tolc = tol;
    tollambda = tol;
    tolp = tol;
    
    x = x0;
    s = s0;
    
    % Start with empty working set
    [n,m] = size(A);
    
    Bset = zeros(0,1);
    Nset = [1:m]';
    
    B = A(:,Bset);
    N = A(:,Nset);
        
    xB = B\b;
    xN = zeros(m,1);
    
    %% Main loop
    maxit = 100;
    
    it = 0;
    while ~converged && (it < maxit)
        it = it + 1;
        
        gB = g(Bset);
        gN = g(Nset);
    
        lambda = B'\gB;
    
        sN = gN - N'*lambda;
    
        if sN >= 0
            disp("Optimum found")
            return
        end
    
        for i = 1:length(Nset)
            if s(i) < 0
                qidx = i;
            end
        end
    
        d = B\A(:,qidx);
    
        if d <= 0
            disp("Unbounded")
            return
        end
    
        idx = find(d > 0);
    
        [xqplus,pidx] = min(xB(idx)./d(idx));
    
        xB = xB - d*xqplus;

        xN = zeros(length(Nset),1);
    
        Bset = union(Bset,qidx);
    
        Bset = setdiff(Bset,pidx);
    
        Nset = setdiff(Nset,Bset);
    
    end
    



end
