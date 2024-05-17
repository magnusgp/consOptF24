%
% Test LP solver
%

clear
clc
close all

n = 200;
m = 50;

A = randn(m,n);

x = zeros(n,1);
x(1:m,1) = abs(rand(m,1));

s = zeros(n,1);
s((m+1):n,1) = abs(rand(n-m,1));

lambda = rand(m,1);

g = A'*lambda + s;
b = A*x;

x0 = ones(n,1);
lambda0 = zeros(m,1);
s0 = ones(n,1);

% [xlp,info,mulp,lambdalp,iter] = LPippd(g,A,b,ones(n,1));
tol = 1.0e-8;
[xout,lambdaout,sout,iter,info,rcres,rbres,mures] = lpsolverInteriorPoint(g,A,b,x0,lambda0,s0,tol);

disp(iter)
disp(max(abs(xout-x)))
disp(max(abs(sout-s)))
disp(max(abs(lambdaout-lambda)))

