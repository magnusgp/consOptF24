function [x,lambda,s,iter,info,rcres,rbres,mures] = lpsolverInteriorPoint(g,A,b,x0,lambda0,s0,tol)
    % LPIPPD   Primal-Dual Interior-Point LP Solver
    %
    %          min  g'*x
    %           x
    %          s.t. A x  = b      (Lagrange multiplier: lambda)
    %                 x >= 0      (Lagrange multiplier: s)
    %
    %         info = true   : Converged
    %              = false  : Not Converged
    
    %%
    [m,n]=size(A);
    
    maxit = 100;
    tolc = tol;
    tolb = tol;
    tolmu = tol;
    
    eta = 0.9;
    
    x = x0;
    lambda = lambda0;
    s = s0;
    
    % Compute residuals
    rc = A'*lambda + s - g;         % Lagrangian gradient
    rb = A*x - b;                   % Equality Constraint
    rxs = x.*s;                     % Complementarity
    mu = sum(rxs)/length(x);        % Duality gap
    
    rcres = rc;
    rbres = rb;
    mures = mu;

    % Converged
    Converged = (norm(rc,'inf') <= tolc) && ...
                (norm(rb,'inf') <= tolb) && ...
                (abs(mu) <= tolmu);
    
    iter = 0;
    while ~Converged && (iter<maxit)
        iter = iter+1;
        
        %% Form and Factorize Hessian Matrix
        
        X = diag(x);
        S = diag(s);
        
        % Sinv = diag(1./s);
        % Xinv = diag(1./x);
        % 
        % D = diag(x./s);
        % H = A*D*A';
        % L = chol(H,'lower');

        %% Affine step

        % Solve

        % rhs = -rb - A*X*Sinv*rc + A*Sinv*rxs;
        % 
        % dlambda = L'\(L\rhs);
        % 
        % ds = -rc - A'*dlambda;
        % 
        % dx = -Sinv*rxs - X*Sinv*ds;

        sysmat = [zeros(n) A' eye(n);
                  A zeros(m) zeros(m,n);
                  S zeros(n,m) X];

        rhs = [-rc;-rb;-rxs];

        sol = sysmat\rhs;

        dx = sol(1:length(x));
        dlambda = sol((length(x)+1):(length(x)+length(lambda)));
        ds = sol((length(x)+length(lambda)+1):end);

        % Step length
        idx = find(dx < 0.0);
        alphaaffpri = min([1.0; -x(idx)./dx(idx)]);
        
        idx = find(ds < 0.0);
        alphaaffdual = min([1.0; -s(idx)./ds(idx)]);
        
        %% Center step

        % Center Parameter
           
        muaff = (x + alphaaffpri*dx)'*(s + alphaaffdual*ds)/length(x);

        sigma = (muaff/mu)^3;

        % Center-Corrector Step
        
        % % rhs = -rb - (A*x*rc + A*rxs)./s;
        % temp = -x.*s - dx.*ds + sigma*mu;
        % 
        % rhs = -rb - A*X*Sinv*rc + A*Sinv*temp;
        % 
        % dlambda = L'\(L\rhs);
        % 
        % ds = -rc - A'*dlambda;
        % 
        % dx = -Sinv*temp - X*Sinv*ds;
        
        rhs = [-rc;-rb;-rxs - dx.*ds + sigma*mu];

        sol = sysmat\rhs;

        dx = sol(1:length(x));
        dlambda = sol((length(x)+1):(length(x)+length(lambda)));
        ds = sol((length(x)+length(lambda)+1):end);

        % Step length
        idx = find(dx < 0.0);
        alphaprimax = min(-x(idx)./dx(idx));

        idx = find(ds < 0.0);
        alphadualmax = min(-s(idx)./ds(idx));

        alphapri = min([1,eta*alphaprimax]);
        alphadual = min([1,eta*alphadualmax]);

        % Take step 
        
        x = x + alphapri*dx;
        lambda = lambda + alphadual*dlambda;
        s = s + alphadual*ds;

        % Residuals and Convergence

        % Compute residuals
        rc = A'*lambda + s - g;         % Lagrangian gradient
        rb = A*x - b;                   % Equality Constraint
        rxs = x.*s;                     % Complementarity
        mu = sum(rxs)/length(x);        % Duality gap
        
        rcres = [rcres rc];
        rbres = [rbres rb];
        mures = [mures mu];

        % Converged
        Converged = (norm(rc,'inf') <= tolc) && ...
                    (norm(rb,'inf') <= tolb) && ...
                    (abs(mu) <= tolmu);

    end
    
    %%
    % Return solution
    info = Converged;
    
end
