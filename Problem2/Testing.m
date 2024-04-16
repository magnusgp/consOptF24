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

figure
tiledlayout(3,2)
PlotSolutionQP(xquadprog)

sprintf("Quadprog \niter: %f \ntime: %f \n", output.iterations, QuadprogTime)

x0 = lb;

A = [eye(length(lb));
     -eye(length(ub));
     sysmat.C';
     -sysmat.C']';

b = -[lb;
     -ub;
     sysmat.dl;
     -sysmat.du];

tic;
[xactiveset,lambdaopt,Wset,it] = qpsolverActiveSet(H,g,A,b,x0);
ActiveSetTime = toc;

PlotSolutionQP(xactiveset)

sprintf("Active set \niter: %f \ntime: %f \n", it, ActiveSetTime)

eps = 0.5;

x0 = x0 + eps;
y0 = [];
z0 = ones(size(A,2),1);
s0 = z0;

maxIter = 100;
tol = 1.0e-8;

[xinteriorpoint] = InteriorPointMethodConvexQP(x0,y0,z0,s0,H,g,[],[],A,b,maxIter,tol);

PlotSolutionQP(xinteriorpoint(:,end))
