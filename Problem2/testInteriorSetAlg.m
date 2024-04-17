clear
clc
close all

addpath('../');

Aineq = [2 1;6 1];
% Aeq = [-3 2];
Aeq = [];

bineq = [10;1];
% beq = -2;
beq = [];

H = [3 1;1 4];
g = [5;1];

f = @(x) 1/2*x'*H*x+g'*x;

x0 = [4;5];

% xout = fmincon(f,x0,A,b,Aeq,beq);
xout = fmincon(f,x0,-Aineq',-bineq,Aeq,beq);

disp("Fmincon:")
disp(xout)

% Need to find feasible x0 and active constraints at x0

% Feasible starting point
x0 = [4;5];
y0 = [];
z0 = [1;1];
s0 = [1;1];

% Assume form
% min f(x) = 1/2 x' H x + g' x
% s.t.
% A' x = b
% C' x >= d

A = Aeq;
b = beq;

C = Aineq;
d = bineq;

H = [3 1;1 4];
g = [5;1];

maxIter = 100;
tol = 10^-6;

[xopt,lambdaopt,X,it] = qpsolverInteriorPoint(x0,y0,z0,s0,H,g,A,b,C,d,maxIter,tol);

% disp("Interior set:")
% disp(X(:,(end-4):end))

%% Plot

N = 150;

x1start = -5;
x1end = 5;
x2start = -3;
x2end = 5;

x1 = x1start:(x1end-x1start)/N:x1end;
x2 = x2start:(x2end-x2start)/N:x2end;

[X1,X2] = meshgrid(x1,x2);

falt = @(x1,x2) 1/2*[x1;x2]'*H*[x1;x2]+g'*[x1;x2];

F = reshape(arrayfun(falt,X1(:),X2(:)),length(x1),length(x2));

v = 0:2:100;

% Inequalities
c1 = @(xx1,xx2) Aineq(:,1)'*[xx1;xx2]<=bineq(1);
c2 = @(xx1,xx2) Aineq(:,2)'*[xx1;xx2]<=bineq(2);

% Equalities
% c3 = @(xx1,xx2) Aeq*[xx1;xx2]==beq;

C1 = reshape(arrayfun(c1,X1(:),X2(:)),length(x1),length(x2));
C2 = reshape(arrayfun(c2,X1(:),X2(:)),length(x1),length(x2));
% C3 = reshape(arrayfun(c3,X1(:),X2(:)),length(x1),length(x2));

colormap jet
contour(X1,X2,F,v,'linewidth',2);
colorbar

hold on

opacity = 0.1;

mkrcol = 'k';
mkrsize = 10;

scatter([X1(C1)],[X2(C1)],mkrsize,"square","filled",'MarkerFaceColor',mkrcol,'MarkerEdgeColor',mkrcol,'MarkerFaceAlpha',opacity,'MarkerEdgeAlpha',opacity)
scatter([X1(C2)],[X2(C2)],mkrsize,"square","filled",'MarkerFaceColor',mkrcol,'MarkerEdgeColor',mkrcol,'MarkerFaceAlpha',opacity,'MarkerEdgeAlpha',opacity)

% hold on
% plot(X1(C3),X2(C3),mkrcol,'LineWidth',1.5)

hold on

plot(X(1,:),X(2,:),'-o','color','r','MarkerFaceColor','r','linewidth',2)
