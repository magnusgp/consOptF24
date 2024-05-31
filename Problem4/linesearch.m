function alpha = linesearch(xk,pk,f,df,cE,cI,lambda,mu)
    % h = cE
    % g = cI
    % Perform a line search to find step size alpha
    alpha = 1;
    stop = false;

    % c:
    phi0 = f(xk) + lambda'*sum(abs(cE(xk))) ... 
                + mu'*abs(min(0,cI(xk)));
    % b:
    dphi0 = df(xk)'*pk - lambda'*abs(cE(xk)) ... 
                - mu'*abs(min(0,cI(xk)));

    while ~stop
        x = xk + alpha * pk;
        
        phia = f(x) + lambda'*sum(abs(cE(x))) ... 
                    + mu'*abs(min(0,cI(x)));

        % Check Armijo condition
        if phia <= phi0 + 1.0e-4*dphi0*alpha
            stop = true;
        else
            a = (phia - (phi0 + dphi0*alpha))/alpha^2;
            alphamin = -dphi0/(2*a);

            alpha = min(0.9*alpha,max(alphamin,0.1*alpha));
        end
    end
end
