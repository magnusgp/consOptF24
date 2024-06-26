clear
clc
close all

% Want to solve
% min f(x) = 1/2 x'Hx + g'x
% s.t.
% A'x = b
% C'x >= d

C = [2 1;6 1];
d = [10;1];

H = [3 1;1 4];
g = [5;1];

f = @(x) 1/2*x'*H*x+g'*x;

x0 = [4;5];

xout = fmincon(f,x0,-C',-d,[],[]);

disp("Fmincon:")
disp(xout)

% Feasible starting point
y0 = [];
z0 = [1;1];
s0 = [1;1];

maxIter = 100;
tol = 10^-6;

predictorCorrector = true;
[xIPPC,lambdaIPPC,XIPPC,itIPPC] = qpsolverInteriorPoint(x0,H,g,[],[],C,d,maxIter,tol,predictorCorrector);

predictorCorrector = false;
[xIP,lambdaIP,XIP,itIP] = qpsolverInteriorPoint(x0,H,g,[],[],C,d,maxIter,tol,predictorCorrector);

[xAS,lambdaAS,XAS,Wset,itAS] = qpsolverActiveSet(H,g,C,-d,x0,tol);

sprintf('IPPC iter: %.0f \nIP iter: %.0f \nAS iter: %.0f',itIPPC,itIP,itAS)

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
c1 = @(xx1,xx2) C(:,1)'*[xx1;xx2]<=d(1);
c2 = @(xx1,xx2) C(:,2)'*[xx1;xx2]<=d(2);

C1 = reshape(arrayfun(c1,X1(:),X2(:)),length(x1),length(x2));
C2 = reshape(arrayfun(c2,X1(:),X2(:)),length(x1),length(x2));

colormap jet
contour(X1,X2,F,v,'linewidth',2, 'HandleVisibility','off');
colorbar

hold on

opacity = 0.1;

mkrcol = 'k';
mkrsize = 10;

scatter([X1(C1)],[X2(C1)],mkrsize,"square","filled",'MarkerFaceColor',mkrcol,'MarkerEdgeColor',mkrcol,'MarkerFaceAlpha',opacity,'MarkerEdgeAlpha',opacity, 'HandleVisibility','off')
scatter([X1(C2)],[X2(C2)],mkrsize,"square","filled",'MarkerFaceColor',mkrcol,'MarkerEdgeColor',mkrcol,'MarkerFaceAlpha',opacity,'MarkerEdgeAlpha',opacity, 'HandleVisibility','off')

hold on

plot(XIP(1,:),XIP(2,:),'-o','color','r','MarkerSize',3,'MarkerFaceColor','r','linewidth',1.5)
plot(XIPPC(1,:),XIPPC(2,:),'-o','color','b','MarkerSize',3,'MarkerFaceColor','b','linewidth',1.5)
plot(XAS(1,:),XAS(2,:),'-o','color','m','MarkerSize',3,'MarkerFaceColor','m','linewidth',1.5)
legend("Interior point w.o. corrector","Interior point w. corrector","Active set",Location="southwest")

sgtitle("Comparing IP with and without corrector and AS")