clear;
clc;
close all;

sysmat = load('QP_Test.mat');

H = sysmat.H;
g = sysmat.g;

lb = sysmat.l;
ub = sysmat.u;

A = [-sysmat.C'; sysmat.C'];
b = [-sysmat.dl; sysmat.du];

Aeq = [];
beq = [];

[xquadprog,fval,exitflag,output] = quadprog(H,g,A,b,Aeq,beq,lb,ub);

PlotSolutionQP(xquadprog)

x0 = lb;

A = [eye(length(lb));
     -eye(length(ub));
     sysmat.C';
     -sysmat.C']';

b = -[lb;
     -ub;
     sysmat.dl;
     -sysmat.du];

[xactiveset,lambdaopt,Wset,it] = qpsolverActiveSet(H,g,A,b,x0);

PlotSolutionQP(xactiveset)
