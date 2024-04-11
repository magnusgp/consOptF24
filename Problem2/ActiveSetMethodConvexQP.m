function [X] = ActiveSetMethodConvexQP(z0,W0,FdF,G,c,A,b,Eps,I)

    % Assume form
    % min f(x) = 1/2 x' G x + c' x
    % s.t.
    % ai' x = bi,  i in Epsilon
    % ai' x >= bi, i in I

    pktol = 10^(-15);

    N = 5;

    xidx = 1:size(G,1);

    % Working set
    Wk = W0;

    xk = z0(xidx);
    Lambdak = zeros(size(A,1),1);
    
    X = xk;
    
    for k = 1:N

        % sprintf("Iter: %.0f",k)
       
        % Constraints from working set
        Ak = A(Wk,:);
        bk = b(Wk);

        % RHS of constraints
        b0 = bk*0;

        zk = [xk;zeros(length(Wk),1)];
        
        gk = G*xk+c;

        zk = NewtonMethod(1, @(z) FdF(z,G,gk,Ak,b0), zk);

        pk = zk(xidx);
        Lambdak(Wk) = zk((xidx(end)+1):end);
        
        if norm(pk,'inf') < pktol
            
            % Wk intersection with I
            WkIntersectI = intersect(Wk,I);

            lambdak = Lambdak(WkIntersectI);

            if min(lambdak) >= 0
                X = [X xk];
                return;
            else
                [~,j] = min(lambdak);

                % Remove index j from working set
                Wk(Wk == j) = [];
            end

        else

            alphas = [];
            blockingConstraints = [];

            % Find i not in Wk and for which ai'pk <0
            for i = 1:(length(Eps)+length(I))
                if sum(i == Wk) ~= 1 % Check that i is not in Wk
                    if A(i,:)*pk < 0 % Check if ai'*pk < 0
                        blockingConstraints = [blockingConstraints i];
                        alphas = [alphas (b(i)-A(i,:)*xk)/(A(i,:)*pk)];
                    end
                end
            end

            alphak = min([1 alphas]);
            
            xk = xk + alphak*pk;

            % Check if any blocking constraints (p. 469)
            if alphak < 1
                % Adding one of the blocking constraints to working set
                Wk = [Wk min(blockingConstraints)];
            end
        end
        
        X = [X xk];

    end
end