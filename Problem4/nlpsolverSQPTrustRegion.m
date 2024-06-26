function [x,X,iter] = nlpsolverSQPTrustRegion(mu,Delta0,x0,y0,z0,cE,cI,AE,AI,f,df,ddL,QPsolver,BFGS)
    
    iter = 0;

    x = x0;
    y = y0;
    z = z0;
    Delta = Delta0;

    if BFGS
        Bk = eye(length(x0));
    else
        Bk = ddL(x,y,z);
    end

    nE = length(cE(x));
    nI = length(cI(x));
    nx = length(x);

    % Merit function
    mk = @(xk,pk) sum(abs(cE(xk)+AE(xk)*pk)) + sum(min(0,-(cI(xk)+AI(xk)*pk)));
    
    % Penalty function
    phi = @(xk) f(xk) + sum(abs(cE(xk))) + mu*sum(min(0,-cI(xk)));

    % Penalized objective function
    qmu = @(xk,pk,Bk) f(xk) + df(xk)'*pk + 0.5*pk'*Bk*pk + mu*mk(xk,pk);

    converged = false;
    maxiter = 200;
    tol = 1.0e-6;
    eta = 0.1;
    gam = 0.9;
    iterreset = 5; % Reset Hessian every 5th iteration

    X = x;

    while ~converged && iter < maxiter

        iter = iter + 1;

        quadterm = [Bk zeros(nx,2*nE+nI);
                    zeros(2*nE+nI,nx) zeros(2*nE+nI,2*nE+nI)];

        linterm = [df(x);
                   mu*ones(nE,1);
                   mu*ones(nE,1);
                   mu*ones(nI,1)];

        Aeq = [AE(x) -eye(nE) eye(nE) zeros(nE,nI)]';

        beq = -cE(x);

        Aineq = [AI(x) zeros(nI,2*nE) eye(nI);
                 eye(nx) zeros(nx,nE) zeros(nx,nE) zeros(nx,nI);
                 -eye(nx) zeros(nx,nE) zeros(nx,nE) zeros(nx,nI);
                 zeros(nE,nx) eye(nE) zeros(nE,nE) zeros(nE,nI);
                 zeros(nE,nx) zeros(nE,nE) eye(nE) zeros(nE,nI);
                 zeros(nI,nx) zeros(nI,nE) zeros(nI,nE) eye(nI)]';

        bineq = [-cI(x);
                 -Delta*ones(nx,1);
                 -Delta*ones(nx,1);
                 zeros(nE,1);
                 zeros(nE,1);
                 zeros(nI,1)];

        init = [x;ones(nE,1);ones(nE,1);ones(nI,1)];

        % [testsol] = quadprog(quadterm,linterm,-Aineq',-bineq,Aeq',beq);

        % Solve subproblem QP
        [sol,lambda,~,~] = QPsolver(init,quadterm,linterm,Aeq,beq,Aineq,bineq,200,1.0e-9,true);

        px = sol(1:nx);
        y = lambda(1:nE);
        z = lambda((nE+1):(nE+nI));

        rho = (phi(x) - phi(x+px))/(qmu(x,zeros(nx,1),Bk)-qmu(x,px,Bk));

        if rho > eta 
            
            xtemp = x;
            x = x + px;

            % Update Hessian if using BFGS
            if BFGS
                
                dLx = df(xtemp)-AE(xtemp)'*y-AI(xtemp)'*z;
    
                dLxplus = df(x)-AE(x)'*y-AI(x)'*z;
    
                p = x - xtemp;
    
                q = dLxplus - dLx;
    
                if p'*q >= 0.2*p'*Bk*p
                    theta = 1;
                else
                    theta = (0.8*p'*Bk*p)/(p'*Bk*p-p'*q);
                end
    
                r = theta*q + (1-theta)*(Bk*p);
    
                Bk = Bk - ((Bk*p)*(Bk*p)')/(p'*(Bk*p)) + (r*r')/(p'*r);
    
                % Reset if NaNs produced
                if nnz(isnan(Bk)) > 0 || mod(iter,iterreset) == 0
                    Bk = eye(length(x0));
                end
            
            else % Otherwise use analytical Hessian
    
                Bk = ddL(x,y,z);
    
            end

            % Reset Delta
            Delta = Delta0;
        else 
            Delta = gam*max(abs(px));
        end

        X = [X x];

        % Check for convergence
        dLx = df(x)-AE(x)'*y-AI(x)'*z;

        converged = norm(dLx,'inf') < tol;

    end
end