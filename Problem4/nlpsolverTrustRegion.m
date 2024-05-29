function [x] = nlpsolverTrustRegion(x0,s0,cE,cI,AE,AI,f,df,dL,ddL)
    
    % Compute initial y and z (19.37)

    Ahat = [AE(x) zeros(ne,mi);
            AI(x) -diag(s)];

    sol = (Ahat*Ahat')\Ahat*[df(x);-mu*ones(ni,1)];
    
    y = sol(1:ne);
    z = sol((ne+1):end);

end