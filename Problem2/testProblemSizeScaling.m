clear
clc
close all

% Diagonal scaling of H
alpha = 10;

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
        
        % Generate sparse random matrices A and M
        C = full(sprandn(n, m, s, 1));
        
        M = full(sprandn(n, n, s, 1));
        
        % Generate H
        H = M * M' + alpha * eye(n);
        
        while cond(H) > 5
    
            disp(cond(H))
        
            M = full(sprandn(n, n, s, 1));
            
            % Generate H
            H = M * M' + alpha * eye(n);
        
        end
        
        % Generate x and lambda
        x0 = randn(n, 1);
        lambda_init = randn(m, 1);
        
        % Generate g and b
        g = H * x0 + C * lambda_init;
        d = C' * x0;
        
        disp("Quadprog")
    
        timefun = @() quadprog(H,g,-C',-d,[],[],[],[],x0,options);
        times(i,j,1) = timeit(timefun);
        [xquadprog,~,exitflag,~] = quadprog(H,g,-C',-d,[],[]);
        disp(exitflag)
    
        disp("AS")
    
        tol = 1.0e-8;
        timefun = @() qpsolverActiveSet(H,g,C,-d,x0,tol);
        times(i,j,2) = timeit(timefun);
        % [xAS,lambdaAS,XAS,Wset,itAS] = qpsolverActiveSet(H,g,C,-d,x0,tol);
        
        disp("IPPC")
    
        predictorCorrector = true;
        maxIter = 200;
        tol = 1.0e-8;
        timefun = @() qpsolverInteriorPoint(x0,H,g,[],[],C,d,maxIter,tol,predictorCorrector);
        times(i,j,3) = timeit(timefun);
        % [xIPPC,lambdaIPPC,XIPPC,itIPPC] = qpsolverInteriorPoint(x0,H,g,[],[],C,d,maxIter,tol,predictorCorrector);
    
        % disp("IP")
        % 
        % predictorCorrector = false;
        % maxIter = 200;
        % tol = 1.0e-8;
        % timefun = @() qpsolverInteriorPoint(x0,y0,z0,s0,H,g,[],[],C,d,maxIter,tol,predictorCorrector);
        % times(i,4) = timeit(timefun);
        % [xIP,lambdaIP,XIP,itIP] = qpsolverInteriorPoint(x0,y0,z0,s0,H,g,[],[],C,d,maxIter,tol,predictorCorrector);
        % 
        % disp(itIP)
         
        % disp(norm(xquadprog-xAS,'inf'))
        % disp(norm(xquadprog-xIPPC,'inf'))

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
% hold on
% plot(N,times(:,4),'-o',LineWidth=1.5)
legend("Quadprog","AS","IP",Location="northwest")
xlabel("Problem size (N)")
ylabel("CPU time (sec.)")
grid on

subplot(1,2,2)
plot(Beta, squeeze(mean(times(:,:,1),1)), 'o-',LineWidth=1.5);
hold on
plot(Beta, squeeze(mean(times(:,:,2),1)), 'o-',LineWidth=1.5);
hold on
plot(Beta, squeeze(mean(times(:,:,3),1)), 'o-',LineWidth=1.5);
legend("Quadprog","AS","IP",Location="northwest")
xlabel("Amount of constraints (beta)")
ylabel("CPU time (sec.)")
grid on

sgtitle("CPU time of quadprog, active set and interior point algorithms")
