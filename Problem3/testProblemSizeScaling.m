clear
clc
close all

beta = 0.8;
alpha = 10;

% Sparsity
s = 0.9;

N = 10:50:400;
times = zeros(length(N),3);

options =  optimset('Display','off');

warning('off')

for i = 1:length(N)

    disp("Iteration")
    disp(N(i))
    
    n = N(i);
    
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
    times(i,1) = timeit(timefun);

    [xlinprog, ~, exitflag] = linprog(g, [], [], A, b, lb, [], options);
    disp(exitflag)

    %% Active set
    disp("AS")

    maxiter = 1000;
    timefun = @() lpsolverActiveSet(g,A,b,x0,Bset,Nset,maxiter);
    times(i,2) = timeit(timefun);

    [xAS,~,itAS] = lpsolverActiveSet(g,A,b,x0,Bset,Nset,maxiter);

    disp(itAS)

    %% Interior point
    disp("IP")

    maxIter = 1000;
    tol = 1.0e-9;
    timefun = @() lpsolverInteriorPoint(g,A,b,tol);
    times(i,3) = timeit(timefun);
    [xIP,lambda,s,itIP,info,~,~,~] = lpsolverInteriorPoint(g,A,b,tol);

    disp(itIP)

    disp("Differences")
    disp(norm(xlinprog-xAS,'inf'))
    disp(norm(xlinprog-xIP,'inf'))

end

%%

figure;
plot(N,times(:,1),'-o',LineWidth=1.5)
hold on
plot(N,times(:,2),'-o',LineWidth=1.5)
hold on
plot(N,times(:,3),'-o',LineWidth=1.5)
legend("Linprog","Active set","Interior point",Location="northwest")
xlabel("Problem size (N)")
ylabel("CPU time (sec.)")
grid on
sgtitle("CPU time of linprog, active set and interior point algorithms")
