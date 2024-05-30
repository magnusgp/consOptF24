clear
clc
close all

% Sparsity
s = 0.8;

Beta = linspace(0.1,0.9,4);
N = round(linspace(10,350,6));

times = zeros(length(N),length(Beta),3);

options =  optimset('Display','off');

warning('off')

for i = 1:length(N)
    for j = 1:length(Beta)
    
        disp("Iteration")
        disp(N(i))
        
        n = N(i);
        beta = Beta(j);

        % Calculate m
        m = round(beta * n);
    
        % Generate system matrix
        A = randn(m,n);
    
        x0 = rand(n,1);
    
        x0(randperm(n-m)) = 0;
        
        b = A*x0;
    
        g = abs(randn(n,1));
    
        Bset = find(x0 > 0);
        Nset = setdiff(1:length(x0),Bset);
    
        % disp(max(x0(Nset)))
        % disp(rank(A(:,Bset)))
        % disp(size(A(:,Bset)))
    
        %% Linprog
        disp("Linprog")
    
        lb = zeros(n,1);
        options = optimoptions('linprog','Display','off');
        timefun = @() linprog(g, [], [], A, b, lb, [], options);
        times(i,j,1) = timeit(timefun);
    
        [xlinprog, ~, exitflag] = linprog(g, [], [], A, b, lb, [], options);
        disp(exitflag)
    
        %% Active set
        disp("AS")
    
        maxiter = 3000;
        timefun = @() lpsolverActiveSet(g,A,b,x0,Bset,Nset,maxiter);
        times(i,j,2) = timeit(timefun);
    
        [xAS,~,itAS] = lpsolverActiveSet(g,A,b,x0,Bset,Nset,maxiter);
    
        disp(itAS)
    
        %% Interior point
        disp("IP")
    
        maxIter = 1000;
        tol = 1.0e-9;
        timefun = @() lpsolverInteriorPoint(g,A,b,tol);
        times(i,j,3) = timeit(timefun);
        [xIP,lambda,s,itIP,info,~,~,~] = lpsolverInteriorPoint(g,A,b,tol);
    
        disp(itIP)
    
        disp("Differences")
        disp(norm(xlinprog-xAS,'inf'))
        disp(norm(xlinprog-xIP,'inf'))
    end
end

%%

figure;
subplot(1,2,1)
plot(N,squeeze(mean(times(:,:,1),2)),'-o',LineWidth=1.5)
hold on
plot(N,squeeze(mean(times(:,:,2),2)),'-o',LineWidth=1.5)
hold on
plot(N,squeeze(mean(times(:,:,3),2)),'-o',LineWidth=1.5)
legend("Linprog","Active set","Interior point",Location="northwest")
xlabel("Problem size (N)")
ylabel("CPU time (sec.)")
grid on

subplot(1,2,2)
plot(Beta,squeeze(mean(times(:,:,1),1)), 'o-',LineWidth=1.5);
hold on
plot(Beta,squeeze(mean(times(:,:,2),1)), 'o-',LineWidth=1.5);
hold on
plot(Beta,squeeze(mean(times(:,:,3),1)), 'o-',LineWidth=1.5);
legend("Linprog","Active set","Interior point",Location="northwest")
xlabel("Amount of constraints (beta)")
ylabel("CPU time (sec.)")
grid on

sgtitle("CPU time of quadprog, active set and interior point algorithms")

% if solver1
%     sysmat = [zeros(n) A' eye(n);
%               A zeros(m) zeros(m,n);
%               S zeros(n,m) X];
% 
%     rhs = [-rc;-rb;-rxs];
% 
%     sol = sysmat\rhs;
% 
%     dx = sol(1:length(x));
%     % dlambda = sol((length(x)+1):(length(x)+length(lambda)));
%     ds = sol((length(x)+length(lambda)+1):end);
% else 
%     D = diag(x./s);
% 
%     sysmat = [D A';
%               A zeros(m,m)];
% 
%     rhs = [-rc+Xinv*rxs;-rb];
% 
%     sol = sysmat\rhs;
% 
%     dx = sol(1:length(x));
%     % dlambda = sol((length(x)+1):(length(x)+length(lambda)));
% 
%     ds = -Xinv*rxs - Xinv*S*dx;
% end
% 
% if solver1
%     rhs = [-rc;-rb;-rxs - dx.*ds + sigma*mu];
% 
%     sol = sysmat\rhs;
% 
%     dx = sol(1:length(x));
%     dlambda = sol((length(x)+1):(length(x)+length(lambda)));
%     ds = sol((length(x)+length(lambda)+1):end);
% else 
%     temp = -x.*s - dx.*ds + sigma*mu;
% 
%     rhs = [-rc+Xinv*temp;-rb];
% 
%     sol = sysmat\rhs;
% 
%     dx = sol(1:length(x));
%     dlambda = sol((length(x)+1):(length(x)+length(lambda)));
% 
%     ds = -Xinv*temp - Xinv*S*dx;
% end