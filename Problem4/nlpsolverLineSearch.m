function [x] = nlpsolverLineSearch(x0,s0,cE,cI,AE,AI,f,df,dL,ddL)
    
    % len(s) = ni
    % len(y) = ne
    % len(z) = ni

    x = x0;
    s = s0;

    % Getting dimensions
    [ne,me] = size(AE(x0));
    [ni,mi] = size(AI(x0));

    nx = length(x);
    ns = length(s);

    % Compute initial y and z (19.37)

    % disp([AE(x) zeros(ne,ne)])
    % disp([AI(x) -diag(s)])

    % Ahat = [AE(x) zeros(ne,ne);
    %         AI(x) -diag(s)];
    % 
    % disp([df(x);-0.9*ones(ni,1)])
    % disp(Ahat)
    % 
    % sol = (Ahat*Ahat')\(Ahat*[df(x);-0.9*ones(ni,1)]);
    % 
    % y = sol(1:ne);
    % z = sol((ne+1):end);

    y = zeros(ne,1);
    z = zeros(ni,1);

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
    
    epstol = 1.0e-5;
    epsmu = 1.0e-5;
    sigma = 0.5;

    while E(x,s,y,z,0) > epstol
        while E(x,s,y,z,tau) > epsmu

            Bk = ddL(x,y,z,s);

            % Compute primal-dual direction

            disp([Bk zeros(nx,ns) AE(x)' AI(x)'])
            disp([zeros(ns,nx) diag(z./s) zeros(ns,ne) -eye(ns)])
            disp(ne)
            
            sysmat = [Bk zeros(nx,ns) AE(x)' AI(x)';
                     zeros(ns,nx) diag(z./s) zeros(ns,ne) -eye(ns);
                     AE(x) zeros(ne,ns) zeros(ne,ne) zeros(ne,ni);
                     AI(x) eye(ni) zeros(ni,ne) zeros(ni,ni)];

            rhs = -[dL(x,y,z,s);
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
        mu = sigma*mu;
        epsmu = mu;
    end


end