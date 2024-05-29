function [x] = nlpsolverLineSearch(x0,s0,cE,cI,AE,AI,f,df,dL,ddL)
    
    % len(s) = ni
    % len(y) = ne
    % len(z) = ni

    % Getting dimensions
    [ne,me] = size(AE(0));
    [ni,mi] = size(AI(0));

    x = x0;
    s = s0;

    % Compute initial y and z (19.37)

    Ahat = [AE(x) zeros(ne,mi);
            AI(x) -diag(s)];

    sol = (Ahat*Ahat')\Ahat*[df(x);-0.9*ones(ni,1)];
    
    y = sol(1:ne);
    z = sol((ne+1):end);

    % Initial penalty parameters
    lambda = max(abs(y));
    mu = max(abs(z));

    % Set parameters
    tau = 0.5;
    k = 0;

    E = @(xk,sk,yk,zk,tauk) max([norm(df(xk)-AE(xk)'*yk-AI(xk)'*zk,1); ...
                           norm(sk.*zk-mu*ones(ni,1),1); ...
                           norm(cE(xk),1); ...
                           norm(cI(xk)-sk,1)]);

    while E(x,s,y,z,0) > epsilon
        while E(x,s,y,z,tau) > tau

            Bk = ddL(x);

            % Compute primal-dual direction
            
            sysmat = [Bk zeros(ni,mi) AE(x) AI(x);
                     zeros(ni,mi) diag(z./s) zeros(ni,me) zeros(ni,mi);
                     AE(x)' zeros(ne,mi) zeros(ne,me) zeros(ne,mi);
                     AI(x)' -eye(ni) zeros(ni,me) zeros(ni,mi)];

            rhs = -[dL(x,y,z);
                    z-tau./s;
                    cE(x);
                    cI(x)-s];

            sol = sysmat\rhs;

            px = sol(1:me);
            ps = sol((me+1):(me+ni));
            py = -sol((me+ni+1):(me+ni+ne));
            pz = -sol((me+ni+ne+1):end);

            % Compute step lengths

            alphamaxs = 1;
    
            while s + alphamaxs*ps > (1-tau)*s
                alphamaxs = alphamaxs*0.99;
            end
            
            alphamaxz = 1;
    
            while z + alphamaxz*pz > (1-tau)*z
                alphamaxz = alphamaxz*0.99;
            end

            % Powells update of penalty parameters
            lambda = max(abs(y),0.5*(lambda+abs(y)));
            mu = max(abs(z),0.5*(mu+abs(z)));

            alphas = linesearch(x,p,f,df,cE,cI,lambda,mu,alphamaxs);
            alphaz = linesearch(x,p,f,df,cE,cI,lambda,mu,alphamaxz);

            % Take steps
            x = x + alphas*px;
            y = y + alphaz*py;
            s = s + alphas*ps;
            z = z + alphaz*pz;

            k = k + 1;
        end
    end


end