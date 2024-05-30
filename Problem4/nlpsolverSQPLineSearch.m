function [x,X,iter] = nlpsolverSQPLineSearch(mu,x0,y0,z0,cE,cI,AE,AI,f,df,ddL,QPsolver,BFGS)
    
    iter = 0;

    x = x0;
    y = y0;
    z = z0;

    if BFGS
        Bk = eye(length(x0));
    else
        Bk = ddL(x,y,z);
    end

    nE = length(cE(x));
    nI = length(cI(x));
    nx = length(x);

    converged = false;
    maxiter = 100;
    tol = 1.0e-6;

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
                 zeros(nE,nx) eye(nE) zeros(nE,nE) zeros(nE,nI);
                 zeros(nE,nx) zeros(nE,nE) eye(nE) zeros(nE,nI);
                 zeros(nI,nx) zeros(nI,nE) zeros(nI,nE) eye(nI)]';

        bineq = [-cI(x);
                 zeros(nE,1);
                 zeros(nE,1);
                 zeros(nI,1)];

        init = [x;ones(nE,1);ones(nE,1);ones(nI,1)];

        % testsol = quadprog(quadterm,linterm,-Aineq',-bineq,Aeq',beq);

        % Solve subproblem QP
        [sol,lambda,~,~] = QPsolver(init,quadterm,linterm,Aeq,beq,Aineq,bineq,200,1.0e-9,true);

        px = sol(1:nx);
        yhat = lambda(1:nE);
        zhat = lambda((nE+1):(nE+nI));

        % % Solve subproblem QP
        % [px,lambda,~,~] = QPsolver(x,Bk,df(x),AE(x)',-cE(x),AI(x)',-cI(x),200,1.0e-9,true);
        % 
        % yhat = lambda(1:nE);
        % zhat = lambda((nE+1):end);

        % Compute step directions
        py = yhat - y;
        pz = zhat - z;
        
        % Compute optimal step lengths
        alpha = linesearch(x,px,f,df,cE,cI,y,z);
        
        xtemp = x;

        % Take steps
        x = x + alpha*px;
        y = y + alpha*py;
        z = z + alpha*pz;

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
            if nnz(isnan(Bk)) > 0
                Bk = eye(length(x0));
            end
        
        else % Otherwise use analytical Hessian

            Bk = ddL(x,y,z);

        end

        X = [X x];

        % Check for convergence
        dLx = df(x)-AE(x)'*y-AI(x)'*z;

        converged = norm(dLx,'inf') < tol;

    end
end