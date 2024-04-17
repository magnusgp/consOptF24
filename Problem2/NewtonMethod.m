function [z,stat] = NewtonMethod(alpha,FdF,z0,varargin)
    
    % max iterations
    maxit = 20;
    tol = 1.0e-5;

    stat.Z = z0;

    % Initialization
    z = z0;
    it = 0;

    % Computing first step
    [F, dF] = FdF(z); 
    H = -dF\F;

    converged = (norm(F,'inf') <= tol);
 
    while ~converged && (it < maxit)
        it = it+1;

        z = z + alpha*H;
        stat.Z = [stat.Z z];
        
        [F, dF] = FdF(z);
        H = -dF\F;

        converged = (norm(F,'inf') <= tol);
    end

end