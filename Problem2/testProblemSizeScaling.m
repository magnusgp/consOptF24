clear
clc
close all

beta = 0.8;
alpha = 100;

tol = 1e+05;

N = 10:10:100;
times = zeros(length(N),3);

for i = 1:length(N)
    
    n = N(i);
    
    % Calculate m
    m = round(beta * n);
    
    % Sparsity
    s = 0.5;
    
    % Generate sparse random matrices A and M
    C = sprandn(n, m, s, 0.5);
    
    M = sprandn(n, n, s, 0.5);
    
    % Generate H
    H = M * M' + alpha * eye(n);
    
    while cond(H) > 10
    
        M = sprandn(n, n, s, 0.5);
        
        % Generate H
        H = M * M' + alpha * eye(n);
    
    end
    
    % Generate x and lambda
    x0 = randn(n, 1);
    lambda_init = randn(m, 1);
    
    % Generate g and b
    g = H * x0 + C * lambda_init;
    d = C' * x0;
    
    y0 = [];
    z0 = ones(size(C,2),1);
    s0 = z0;
    
    timefun = @() quadprog(H,g,-C',-d,[],[]);
    times(i,1) = timeit(timefun);
    % [xquadprog,fval,exitflag,output] = quadprog(H,g,-C',-d,[],[]);
    
    timefun = @() qpsolverActiveSet(H,g,C,-d,x0);
    times(i,2) = timeit(timefun);
    % [xAS,lambdaAS,XAS,Wset,itAS] = qpsolverActiveSet(H,g,C,-d,x0);
    
    predictorCorrector = true;
    maxIter = 100;
    tol = 1.0e-5;
    timefun = @() qpsolverInteriorPoint(x0,y0,z0,s0,H,g,[],[],C,d,maxIter,tol,predictorCorrector);
    times(i,3) = timeit(timefun);
    % [xIPPC,lambdaIPPC,XIPPC,itIPPC] = qpsolverInteriorPoint(x0,y0,z0,s0,H,g,[],[],C,d,maxIter,tol,predictorCorrector);
    
    % disp(norm(xquadprog-xAS,'inf'))
    % disp(norm(xquadprog-xIPPC,'inf'))
    % disp(norm(xAS-xIPPC,'inf'))

end

%%

figure;
plot(N,times(:,1),'-o')
hold on
plot(N,times(:,2),'-o')
hold on
plot(N,times(:,3),'-o')
legend("Quadprog","Active set","Interior point")