function [X] = InteriorPointMethodConvexQP(x0,y0,z0,s0,H,g,A,b,C,d)

    % Assume form
    % min f(x) = 1/2 x' H x + g' x
    % s.t.
    % A' x = b
    % C' x >= d
    %
    % We assume C is non-empty, whereas A may be empty, i.e. A = [];
    
    % Initialize
    x = x0;
    y = y0;
    z = z0;
    s = s0;

    X = x;

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

    mu = (z'*s)/mc;

    % while FirstOrderKKT()
    for i = 1:10
        
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
        deltasaff = -rsz./z - deltazaff.*s./z;

        alphaaff = max([-z./deltazaff;-s./deltasaff]);

        %%% Duality gap and centering parameter %%% 
        % This part is not working

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

        %%% Update iteration %%%

        % alpha = alphaaff;
        % deltax = deltaxaff;
        % deltay = deltayaff;
        % deltaz = deltazaff;
        % deltas = deltasaff;

        eta = 0.995;
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

    end
end