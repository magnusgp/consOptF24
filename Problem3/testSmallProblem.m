clear
clc
close all
restoredefaultpath

% Want to solve
% min f(x) = g'x
% s.t.
% A'x = b
% x >= 0

A = [2 1];
b = [10];

g = [5;1];

x0AS = [5;0];
Bset = [1];
Nset = [2];

x0IP = [1;1];

% A = [1 1 1 0;
%      2 0.5 0 1];
% b = [5;8];
% 
% g = [-3;-2;0;0];
% 
% x0AS = [0;0;1;1];
% Bset = [3,4];
% Nset = [1,2];
% 
% x0IP = [1,1,1,1];

lb = zeros(length(g),1);

f = @(x) g'*x;

xlinprog = linprog(g,[],[],A,b,lb,[]);

disp("linprog:")
disp(xlinprog)

% Must start away from boundary
s0 = ones(length(x0IP),1);
lambda0 = zeros(size(A,1),1);

maxIter = 100;
tol = 10^-6;

[xIP,lambda,s,iterIP,info,rcres,rbres,mures,XIP] = lpsolverInteriorPoint(g,A,b,x0IP,lambda0,s0,tol);

disp("IP:")
disp(xIP)
disp(iterIP)

maxiter = 100;
[xAS,XAS,iterAS] = lpsolverActiveSet(g,A,b,x0AS,Bset,Nset,maxiter);

disp("AS:")
disp(xAS)
disp(iterAS)

%% Plot

N = 150;

x1start = 0;
x1end = 5;
x2start = 0;
x2end = 10;

x1 = x1start:(x1end-x1start)/N:x1end;
x2 = x2start:(x2end-x2start)/N:x2end;

[X1,X2] = meshgrid(x1,x2);

falt = @(x1,x2) g'*[x1;x2];

F = reshape(arrayfun(falt,X1(:),X2(:)),length(x1),length(x2));

v = -50:1:50;

% Inequalities
c1 = @(xx1,xx2) A*[xx1;xx2]==b;

C1 = reshape(arrayfun(c1,X1(:),X2(:)),length(x1),length(x2));

colormap jet
contour(X1,X2,F,v,'linewidth',2, 'HandleVisibility','off');
colorbar

hold on

plot(X1(C1),X2(C1),'-',LineWidth=3,Color='k')

plot(XIP(1,:),XIP(2,:),'-o','color','r','MarkerSize',5,'MarkerFaceColor','r','linewidth',1.5)
plot(XAS(1,:),XAS(2,:),'-o','color','m','MarkerSize',5,'MarkerFaceColor','m','linewidth',1.5)
legend("Ax=b","Interior point","Active set",Location="northeast",fontsize=10)