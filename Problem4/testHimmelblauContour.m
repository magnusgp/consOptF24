clear
clc
close all

addpath('../Problem2')

f = @(x) (x(1).^2+x(2)-11).^2 + (x(1) + x(2).^2 - 7).^2;

% options = optimoptions("fmincon",...
%     "Algorithm","interior-point",...
%     "SubproblemAlgorithm","cg");

[xout1,fval1,exitflag1,output1] = fmincon(f,[-4;2],-[-4 10],0,[],[],[],[],@mycon);
[xout2,fval2,exitflag2,output2] = fmincon(f,[0;1],-[-4 10],0,[],[],[],[],@mycon);

disp(output1)
disp(output2)

sol1 = [-3.654605171;2.737718272]';
sol2 = [-0.298347614;2.895620843]';
sol3 = [-1.424243078;0.3314960330];

df = @(x) [4*x(1)^3 + (4*x(2) - 42)*x(1) + 2*x(2)^2 - 14;
           4*x(2)^3 + (4*x(1) - 26)*x(2) + 2*x(1)^2 - 22];

cE = @(x) (x(1) + 2).^2 - x(2);
cI = @(x) -4*x(1) + 10*x(2);

AE = @(x) [2*x(1)+4 -1];
AI = @(x) [-4 10];

ddL = @(x,y,z) [12*x(1)^2 + 4*x(2) - 2*y - 42, 4*x(1) + 4*x(2); 
                  4*x(1) + 4*x(2), 12*x(2)^2 + 4*x(1) - 26];

N = 150;

x1start = -5;
x1end = 5;
x2start = -5;
x2end = 5;

x1 = x1start:(x1end-x1start)/N:x1end;
x2 = x2start:(x2end-x2start)/N:x2end;

[X1,X2] = meshgrid(x1,x2);

falt = @(xx1,xx2) (xx1.^2+xx2-11).^2 + (xx1 + xx2.^2 - 7).^2;

F = reshape(arrayfun(falt,X1(:),X2(:)),length(x1),length(x2));

v = 0:20:500;

% Equalities
c1 = @(xx1,xx2) (xx1 + 2).^2 - xx2;
c2 = @(xx1,xx2) -4*xx1 + 10*xx2 <= 0;

% C1 = reshape(arrayfun(c1,X1(:),X2(:)),length(x(1)),length(x(2)));
C2 = reshape(arrayfun(c2,X1(:),X2(:)),length(x1),length(x2));

figure

colormap jet
contour(X1,X2,F,v,'linewidth',1, 'HandleVisibility','off');
colorbar

hold on

opacity = 0.1;
mkrcol = 'k';
mkrsize = 9;
scatter([X1(C2)],[X2(C2)],mkrsize,"square","filled",'MarkerFaceColor',mkrcol,'MarkerEdgeColor',mkrcol,'MarkerFaceAlpha',opacity,'MarkerEdgeAlpha',opacity, 'HandleVisibility','off')

hold on

fimplicit(c1,LineWidth = 2,Color='k')
xlabel("$x_1$")
ylabel("$x_2$")
xticks(-5:1:5)

hold on
plot(sol1(1),sol1(2),'*',MarkerSize=10,LineWidth=2,Color='r')
hold on
plot(sol2(1),sol2(2),'*',MarkerSize=10,LineWidth=2,Color='g')
hold on
plot(sol3(1),sol3(2),'*',MarkerSize=10,LineWidth=2,Color='b')

legend("","$x^*_1$: local minimum","$x^*_2$: local minimum","$x^*_3$: saddle",Location="southeast")

sgtitle("Contour plot with stationary points")