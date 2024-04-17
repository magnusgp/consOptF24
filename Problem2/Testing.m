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

sprintf("Quadprog \niter: %.0f \ntime: %.2f \n", output.iterations, QuadprogTime)

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
[xAS,lambdaAS,XAS,Wset,itAS] = qpsolverActiveSet(H,g,A,b,x0);
ASTime = toc;

PlotSolutionQP(xAS)

sprintf("Active set \niter: %.0f \ntime: %.2f \n", itAS, ASTime)

x0 = lb;
y0 = [];
z0 = ones(size(A,2),1);
s0 = z0;

maxIter = 200;
tol = 1.0e-5;

tic;
[xIP,lambdaIP,XIP,itIP] = qpsolverInteriorPoint(x0,y0,z0,s0,H,g,[],[],A,-b,maxIter,tol);
IPTime = toc;

PlotSolutionQP(xIP)

sprintf("Interior point \niter: %.0f \ntime: %.2f \n", itIP, IPTime)
