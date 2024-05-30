clear;
clc;
close all;
restoredefaultpath;

sysmat = load('QP_Test.mat');

H = sysmat.H;
g = sysmat.g;

lb = sysmat.l;
ub = sysmat.u;

A = [-sysmat.C'; sysmat.C'];
b = [-sysmat.dl; sysmat.du];

Aeq = [];
beq = [];

tic;
[xquadprog,fval,exitflag,output] = quadprog(H,g,A,b,Aeq,beq,lb,ub);
QuadprogTime = toc;

PlotSolutionQP(xquadprog)
sgtitle("Quadprog")

% sprintf("Quadprog \niter: %.0f \ntime: %.2f \n", output.iterations, QuadprogTime)

x0 = lb;

A = [eye(length(lb));
     -eye(length(ub));
     sysmat.C';
     -sysmat.C']';

b = -[lb;
     -ub;
     sysmat.dl;
     -sysmat.du];

tol = 1.0e-8;

tic;
[xAS,lambdaAS,XAS,Wset,itAS] = qpsolverActiveSet(H,g,A,b,x0,tol);
ASTime = toc;

PlotSolutionQP(xAS)
sgtitle("Active set")

% sprintf("Active set \niter: %.0f \ntime: %.2f \n", itAS, ASTime)

x0 = 50*ones(length(g),1);
y0 = [];
z0 = 10*ones(size(A,2),1);
s0 = z0;

maxIter = 100;
tol = 1.0e-8;

tic;
predictorCorrector = true;
[xIPPC,lambdaIPPC,XIPPC,itIPPC,rszres,rLres,rCres] = qpsolverInteriorPoint(x0,H,g,[],[],A,-b,maxIter,tol,predictorCorrector);
IPPCTime = toc;

PlotSolutionQP(xIPPC)
sgtitle("Interior point w. predictor corrector")

% sprintf("Interior point w. predictor corrector \niter: %.0f \ntime: %.2f \n", itIPPC, IPPCTime)

tic;
predictorCorrector = false;
[xIP,lambdaIP,XIP,itIP] = qpsolverInteriorPoint(x0,H,g,[],[],A,-b,maxIter,tol,predictorCorrector);
IPTime = toc;

PlotSolutionQP(xIP)
sgtitle("Interior point w.o. predictor corrector")

% sprintf("Interior point w.o. predictor corrector \niter: %.0f \ntime: %.2f \n", itIP, IPTime)

%%

figure
subplot(1,3,1)
plot(vecnorm(rszres,2),'-o',LineWidth=1.5)
xlabel("Iterations")
ylabel("$||r_{sz}||_2$")
grid on
title(sprintf("$r_{sz}$"))
subplot(1,3,2)
plot(vecnorm(rLres,2),'-o',LineWidth=1.5)
xlabel("Iterations")
ylabel("$||r_{L}||_2$")
grid on
title("$r_L$")
subplot(1,3,3)
plot(vecnorm(rCres,2),'-o',LineWidth=1.5)
xlabel("Iterations")
ylabel("$||r_{C}||_2$")
grid on
title("$r_C$")
sgtitle("Residuals w.r.t. iterations")

disp("Max error, iter and time of IP:")
disp(max(abs(xIPPC-xquadprog)))
disp(itIPPC)
disp(IPPCTime)

disp("Max error, iter and time of AS:")
disp(max(abs(xAS-xquadprog)))
disp(itAS)
disp(ASTime)

disp("Iter and time of quadprog")
disp(output.iterations)
disp(QuadprogTime)