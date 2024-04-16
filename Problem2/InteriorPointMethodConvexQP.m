function [X] = InteriorPointMethodConvexQP(x0,y0,z0,s0,H,g,A,b,C,d,maxIter,tol)

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

    if isempty(A) % If no equality constraints
        rL = H*x+g-C*z;
        rA = [];
    else
        rL = H*x+g-A*y-C*z;
        rA = b-A'*x;
    end
    
    rC = s+d-C'*x;
    rsz = s.*z;

    mc = length(s);

    %% Compute starting values

    mu = (z'*s)/mc;

    % Compute S^-1 Z
    D = diag(z./s);

    Hbar = H + C*D*C';

    sysmat = [Hbar -A;-A' zeros(size(A,2),size(A,2))];

    % [Lsys,Dsys] = ldl(sysmat);

    %%% Affine direction %%%

    rLbar = rL - C*D*(rC - rsz./z);

    sysrhs = -[rLbar;rA];

    % [deltaxyaff] = Lsys'\(Dsys\(Lsys\sysrhs));
    [deltaxyaff] = sysmat\sysrhs;

    deltaxaff = deltaxyaff(1:length(x));
    deltayaff = deltaxyaff((length(x)+1):end);

    deltazaff = -D*C'*deltaxaff + D*(rC-rsz./z);
    deltasaff = -rsz./z - (deltazaff.*s)./z;

    for i = 1:length(y0)
        y(i) = max(1,abs(y0(i)+deltayaff));
        s(i) = max(1,abs(s0(i)+deltasaff));
    end

    X = x;

    k = 0;

    converged = false;

    while ~converged && k < maxIter

        k = k + 1;
        
        % Compute S^-1 Z
        D = diag(z./s);

        Hbar = H + C*D*C';
        
        sysmat = [Hbar -A;-A' zeros(size(A,2),size(A,2))];

        % [Lsys,Dsys] = ldl(sysmat);

        %% Affine direction
        
        rLbar = rL - C*D*(rC - rsz./z);
        
        sysrhs = -[rLbar;rA];
        
        % [deltaxyaff] = Lsys'\(Dsys\(Lsys\sysrhs));
        [deltaxyaff] = sysmat\sysrhs;

        deltaxaff = deltaxyaff(1:length(x));
        deltayaff = deltaxyaff((length(x)+1):end);

        deltazaff = -D*C'*deltaxaff + D*(rC-rsz./z);
        deltasaff = -rsz./z - deltazaff.*s./z;

        alphaaff = max([-z./deltazaff;-s./deltasaff]);

        %% Duality gap and centering parameter
        
        % Choose to use centering correction
        correctorstep = true;
        
        if correctorstep
            muaff = (z + alphaaff*deltazaff)'*(s + alphaaff*deltasaff)/mc;
    
            sigma = (muaff/mu)^3;
    
            %%% Affine-centering-correction direction %%%
    
            rszbar = rsz + deltasaff.*deltazaff - sigma*mu*ones(length(rsz),1);
    
            rLbar = rL - C*D*(rC - rszbar./z);
    
            sysrhs = -[rLbar;rA];
    
            % [deltaxy] = Lsys'\(Dsys\(Lsys\sysrhs));
            [deltaxy] = sysmat\sysrhs;
    
            deltax = deltaxy(1:length(x));
            deltay = deltaxy((length(x)+1):end);
    
            deltaz = -D*C'*deltax + D*(rC - rszbar./z);
            deltas = -rszbar./z - deltaz.*s./z;
    
            alpha = max([-z./deltaz;-s./deltas]);
            
        end
        %%% Update iteration %%%
        
        if ~correctorstep
            alpha = alphaaff;
            deltax = deltaxaff;
            deltay = deltayaff;
            deltaz = deltazaff;
            deltas = deltasaff;
        end

        eta = 0.5;
        alphabar = eta*alpha;
        
        x = x + alphabar*deltax;
        y = y + alphabar*deltay;
        z = z + alphabar*deltaz;
        s = s + alphabar*deltas;
        
        if isempty(A) % If no equality constraints
            rL = H*x+g-C*z;
            rA = [];
        else
            rL = H*x+g-A*y-C*z;
            rA = b-A'*x;
        end

        rC = s+d-C'*x;
        rsz = s.*z;

        mu = (z'*s)/mc;

        X = [X x];

        converged = norm(blkdiag(rL,rA,rC,mu),'inf') < tol;

    end
end