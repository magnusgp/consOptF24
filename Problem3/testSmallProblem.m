clear
clc
close all

% Want to solve
% min f(x) = g'x
% s.t.
% A'x = b
% x >= 0

A = [2 1];
b = [7];

g = [5;1];

f = @(x) g'*x;

x0 = [4;5];

xout = linprog(g,[],[],A,b,[0;0],[]);

disp("linprog:")
disp(xout)

%%

% Feasible starting point
x0 = [4;5];
y0 = [];
z0 = [1;1];
s0 = [1;1];

maxIter = 100;
tol = 10^-6;



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

plot(X1(C1),X2(C1),'-',LineWidth=1.5,Color='k')

% plot(XIP(1,:),XIP(2,:),'-o','color','r','MarkerSize',3,'MarkerFaceColor','r','linewidth',1.5)
% plot(XIPPC(1,:),XIPPC(2,:),'-o','color','b','MarkerSize',3,'MarkerFaceColor','b','linewidth',1.5)
% plot(XAS(1,:),XAS(2,:),'-o','color','m','MarkerSize',3,'MarkerFaceColor','m','linewidth',1.5)
% legend("Interior point w.o. corrector","Interior point w. corrector","Active set",Location="southwest")