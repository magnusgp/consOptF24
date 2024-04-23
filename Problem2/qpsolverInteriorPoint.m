function [xopt,lambdaopt,X,it] = qpsolverInteriorPoint(x0,y0,z0,s0,H,g,A,b,C,d,maxIter,tol,predictorCorrector)

    % Assume form
    % min f(x) = 1/2 x' H x + g' x
    % s.t.
    % A' x = b
    % C' x >= d
    %
    % We assume C is non-empty, whereas A may be empty, i.e. A = [];
    % 
    % x: variable we want to determine
    % y: equality constraint lagrange multiplier
    % z: inequality constraint lagrange multiplier
    % s: slack variable
    
    % Initialize
    x = x0;
    y = y0;
    z = z0;
    s = s0;

    mc = 0.5*length(s);

    e = ones(length(s),1);

    %% Compute starting values

    % Compute residuals

    if isempty(A) % If no equality constraints
        rL = H*x+g-C*z;
        rA = [];
    else
        rL = H*x+g-A*y-C*z;
        rA = b-A'*x;
    end
    
    rC = s+d-C'*x;
    rsz = s.*z;

    D = diag(z./s);

    % One step in affine direction

    Hbar = H + C*D*C';

    sysmat = [Hbar -A;-A' zeros(size(A,2),size(A,2))];

    [Lsys,Dsys] = ldl(sysmat);

    rLbar = rL - C*D*(rC - rsz./z);

    sysrhs = -[rLbar;rA];

    [deltaxyaff] = Lsys'\(Dsys\(Lsys\sysrhs));

    deltaxaff = deltaxyaff(1:length(x));
    % deltayaff = deltaxyaff((length(x)+1):end);

    deltazaff = -D*C'*deltaxaff + D*(rC-rsz./z);
    deltasaff = -rsz./z - (s./z).*deltazaff;

    for i = 1:length(z)
        z(i) = max(1,abs(z(i)+deltazaff(i)));
        s(i) = max(1,abs(s(i)+deltasaff(i)));
    end

    %% Start iterations
    
    % Compute residuals

    if isempty(A) % If no equality constraints
        rL = H*x+g-C*z;
        rA = [];
    else
        rL = H*x+g-A*y-C*z;
        rA = b-A'*x;
    end
    
    rC = s+d-C'*x;
    rsz = s.*z;

    D = diag(z./s);

    mu = (z'*s)/mc;

    % Start iterations

    X = x;

    k = 0;

    converged = norm(blkdiag(rL,rA,rC,mu),'inf') < tol;

    while ~converged && k < maxIter

        k = k + 1;
        
        Hbar = H + C*D*C';
        
        sysmat = [Hbar -A;-A' zeros(size(A,2),size(A,2))];

        [Lsys,Dsys] = ldl(sysmat);

        %% Affine direction
        
        rLbar = rL - C*D*(rC - rsz./z);
        
        sysrhs = -[rLbar;rA];

        [deltaxyaff] = Lsys'\(Dsys\(Lsys\sysrhs));

        deltaxaff = deltaxyaff(1:length(x));
        deltayaff = deltaxyaff((length(x)+1):end);

        deltazaff = -D*C'*deltaxaff + D*(rC-rsz./z);
        deltasaff = -rsz./z - (s./z).*deltazaff;

        alphaaff = 1;

        while min(z + alphaaff*deltazaff) < 0 && min(s + alphaaff*deltasaff) < 0
            alphaaff = alphaaff*0.9;
        end

        %% Duality gap and centering parameter
        
        if predictorCorrector
            muaff = ((z + alphaaff*deltazaff)'*(s + alphaaff*deltasaff))/mc;
    
            sigma = (muaff/mu)^3;
    
            %% Affine-centering-correction direction
    
            rszbar = rsz + diag(deltasaff)*diag(deltazaff)*e - sigma*mu*e;
    
            rLbar = rL - C*D*(rC - rszbar./z);
    
            sysrhs = -[rLbar;rA];
    
            [deltaxy] = Lsys'\(Dsys\(Lsys\sysrhs));
    
            deltax = deltaxy(1:length(x));
            deltay = deltaxy((length(x)+1):end);
            
            deltaz = -D*C'*deltax + D*(rC-rszbar./z);
            deltas = -rszbar./z - diag(s./z)*deltaz;

            alpha = 1;

            while min(z + alpha*deltaz) < 0 && min(s + alpha*deltas) < 0
                alpha = alpha*0.9;
            end
            
        end
        %%% Update iteration %%%
        
        if ~predictorCorrector
            alpha = alphaaff;
            deltax = deltaxaff;
            deltay = deltayaff;
            deltaz = deltazaff;
            deltas = deltasaff;
        end

        eta = 0.995;
        alphabar = eta*alpha;
        
        x = x + alphabar*deltax;
        y = y + alphabar*deltay;
        z = z + alphabar*deltaz;
        s = s + alphabar*deltas;
        
        % Compute residuals

        if isempty(A) % If no equality constraints
            rL = H*x+g-C*z;
            rA = [];
        else
            rL = H*x+g-A*y-C*z;
            rA = b-A'*x;
        end

        rC = s+d-C'*x;
        rsz = s.*z;

        D = diag(z./s);

        mu = (z'*s)/mc;

        X = [X x];

        converged = norm(blkdiag(rL,rA,rC,mu),'inf') < tol;

    end

    it = k;
    xopt = x;
    lambdaopt = [y;z];
end