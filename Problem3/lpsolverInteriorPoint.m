function [x,lambda,s,iter,info,rcres,rbres,mures,Xmat] = lpsolverInteriorPoint(g,A,b,tol)
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
    
    %% Find initial point

    xtildetemp = (A*A')\b;
    xtilde = A'*xtildetemp;
    lambdatilde = (A*A')\A*g;
    stilde = g - A'*lambdatilde;
    
    deltax = max([-(3/2)*min(xtilde),0]);
    deltas = max([-(3/2)*min(stilde),0]);
    
    e = ones(length(xtilde),1);
    
    xIP = xtilde + deltax*e;
    shat = stilde + deltas*e;
    
    deltaxhat = 0.5*(xIP'*shat)/(e'*shat);
    deltashat = 0.5*(xIP'*shat)/(e'*xIP);
    
    x0 = xIP + deltaxhat*e;
    lambda0 = lambdatilde;
    s0 = shat + deltashat*e;

    %% Initialize algorithm

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
    Xmat = x;

    % Converged
    Converged = (norm(rc,'inf') <= tolc) && ...
                (norm(rb,'inf') <= tolb) && ...
                (abs(mu) <= tolmu);

    % When close to optimum, use more stable solver
    tol1 = 1.0e+05;
    % solver1 = (norm(rc,'inf') <= tol1) && ...
    %           (norm(rb,'inf') <= tol1) && ...
    %           (abs(mu) <= tol1);
    solver1 = (min(x) < tol1) && (min(s) < tol1);

    iter = 0;
    while ~Converged && (iter<maxit)
        iter = iter+1;
        
        %% Form and Factorize Hessian Matrix
        
        X = diag(x);
        S = diag(s);

        Sinv = diag(1./s);
        Xinv = diag(1./x);
        
        %% Affine step

        if solver1
            sysmat = [zeros(n) A' eye(n);
                      A zeros(m) zeros(m,n);
                      S zeros(n,m) X];
    
            rhs = [-rc;-rb;-rxs];
    
            sol = sysmat\rhs;
    
            dx = sol(1:length(x));
            % dlambda = sol((length(x)+1):(length(x)+length(lambda)));
            ds = sol((length(x)+length(lambda)+1):end);
        else 
            disp("here")
            D = diag(x./s);

            sysmat = [D A';
                      A zeros(m,m)];

            rhs = [-rc+Xinv*rxs;-rb];

            sol = sysmat\rhs;

            dx = sol(1:length(x));
            % dlambda = sol((length(x)+1):(length(x)+length(lambda)));

            ds = -Xinv*rxs - Xinv*S*dx;
        end

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
        
        if solver1
            rhs = [-rc;-rb;-rxs - dx.*ds + sigma*mu];

            sol = sysmat\rhs;
    
            dx = sol(1:length(x));
            dlambda = sol((length(x)+1):(length(x)+length(lambda)));
            ds = sol((length(x)+length(lambda)+1):end);
        else 
            temp = -x.*s - dx.*ds + sigma*mu;

            rhs = [-rc+Xinv*temp;-rb];

            sol = sysmat\rhs;

            dx = sol(1:length(x));
            dlambda = sol((length(x)+1):(length(x)+length(lambda)));

            ds = -Xinv*temp - Xinv*S*dx;
        end
        
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
        Xmat = [Xmat x];

        % Converged
        Converged = (norm(rc,'inf') <= tolc) && ...
                    (norm(rb,'inf') <= tolb) && ...
                    (abs(mu) <= tolmu);

        % solver1 = (norm(rc,'inf') <= tol1) && ...
        %           (norm(rb,'inf') <= tol1) && ...
        %           (abs(mu) <= tol1);
        solver1 = (min(x) < tol1) && (min(s) < tol1);

    end
    
    %%
    % Return solution
    info = Converged;
    
end
