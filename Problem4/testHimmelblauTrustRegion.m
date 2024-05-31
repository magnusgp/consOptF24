clear
clc
close all

addpath('../Problem2')

f = @(x) (x(1).^2+x(2)-11).^2 + (x(1) + x(2).^2 - 7).^2;

options = optimoptions("fmincon",...
    "Algorithm","interior-point",...
    "SubproblemAlgorithm","cg");

xout1 = fmincon(f,[-4;2],-[-4 10],0,[],[],[],[],@mycon,options);
xout2 = fmincon(f,[0;1],-[-4 10],0,[],[],[],[],@mycon,options);
xout3 = fmincon(f,[4;2],-[-4 10],0,[],[],[],[],@mycon,options);

sol1 = [-3.654605171;2.737718272]';
sol2 = [-0.298347614;2.895620843]';

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

y0 = [1];
z0 = [1];
mu = 50;
Delta0 = 1;

x01 = [-4;2];
x02 = [4;2];
x03 = [0;1];
x04 = [0;-3];

[x,Xout1,iter1] = nlpsolverSQPTrustRegion(mu,Delta0,x01,y0,z0,cE,cI,AE,AI,f,df,ddL,@qpsolverInteriorPoint,true);
[x,Xout2,iter2] = nlpsolverSQPTrustRegion(mu,Delta0,x02,y0,z0,cE,cI,AE,AI,f,df,ddL,@qpsolverInteriorPoint,true);
[x,Xout3,iter3] = nlpsolverSQPTrustRegion(mu,Delta0,x03,y0,z0,cE,cI,AE,AI,f,df,ddL,@qpsolverInteriorPoint,true);
[x,Xout4,iter4] = nlpsolverSQPTrustRegion(mu,Delta0,x04,y0,z0,cE,cI,AE,AI,f,df,ddL,@qpsolverInteriorPoint,true);

figure
subplot(1,2,1)

colormap jet
contour(X1,X2,F,v,'linewidth',1, 'HandleVisibility','off');
xticks(x1start:1:x1end)

hold on

opacity = 0.1;
mkrcol = 'k';
mkrsize = 5;
scatter([X1(C2)],[X2(C2)],mkrsize,"square","filled",'MarkerFaceColor',mkrcol,'MarkerEdgeColor',mkrcol,'MarkerFaceAlpha',opacity,'MarkerEdgeAlpha',opacity, 'HandleVisibility','off')

hold on

fimplicit(c1,LineWidth = 2,Color='k')
xlabel("$x_1$")
ylabel("$x_2$")

hold on
plot(Xout1(1,:),Xout1(2,:),'-o',LineWidth=1.5,Color='m')
hold on
plot(Xout1(1,end),Xout1(2,end),'*',LineWidth=2,MarkerSize=14,Color='k')
hold on
plot(Xout2(1,:),Xout2(2,:),'-o',LineWidth=1.5,Color='m')
hold on
plot(Xout2(1,end),Xout2(2,end),'*',LineWidth=2,MarkerSize=14,Color='k')
hold on
plot(Xout3(1,:),Xout3(2,:),'-o',LineWidth=1.5,Color='m')
hold on
plot(Xout3(1,end),Xout3(2,end),'*',LineWidth=2,MarkerSize=14,Color='k')
hold on
plot(Xout4(1,:),Xout4(2,:),'-o',LineWidth=1.5,Color='m')
hold on
plot(Xout4(1,end),Xout4(2,end),'*',LineWidth=2,MarkerSize=14,Color='k')
title("With BFGS")

% disp(Xout1(:,end)-sol1')
% disp(Xout3(:,end)-sol2')

[x,Xout1,iter1] = nlpsolverSQPTrustRegion(mu,Delta0,x01,y0,z0,cE,cI,AE,AI,f,df,ddL,@qpsolverInteriorPoint,false);
[x,Xout2,iter2] = nlpsolverSQPTrustRegion(mu,Delta0,x02,y0,z0,cE,cI,AE,AI,f,df,ddL,@qpsolverInteriorPoint,false);
[x,Xout3,iter3] = nlpsolverSQPTrustRegion(mu,Delta0,x03,y0,z0,cE,cI,AE,AI,f,df,ddL,@qpsolverInteriorPoint,false);
[x,Xout4,iter4] = nlpsolverSQPTrustRegion(mu,Delta0,x04,y0,z0,cE,cI,AE,AI,f,df,ddL,@qpsolverInteriorPoint,false);

subplot(1,2,2)

colormap jet
contour(X1,X2,F,v,'linewidth',1, 'HandleVisibility','off');
xticks(x1start:1:x1end)
colorbar

hold on

opacity = 0.1;
mkrcol = 'k';
mkrsize = 5;
scatter([X1(C2)],[X2(C2)],mkrsize,"square","filled",'MarkerFaceColor',mkrcol,'MarkerEdgeColor',mkrcol,'MarkerFaceAlpha',opacity,'MarkerEdgeAlpha',opacity, 'HandleVisibility','off')

hold on

fimplicit(c1,LineWidth = 2,Color='k')
xlabel("$x_1$")
ylabel("$x_2$")

hold on
plot(Xout1(1,:),Xout1(2,:),'-o',LineWidth=1.5,Color='m')
hold on
plot(Xout1(1,end),Xout1(2,end),'*',LineWidth=2,MarkerSize=14,Color='k')
hold on
plot(Xout2(1,:),Xout2(2,:),'-o',LineWidth=1.5,Color='m')
hold on
plot(Xout2(1,end),Xout2(2,end),'*',LineWidth=2,MarkerSize=14,Color='k')
hold on
plot(Xout3(1,:),Xout3(2,:),'-o',LineWidth=1.5,Color='m')
hold on
plot(Xout3(1,end),Xout3(2,end),'*',LineWidth=2,MarkerSize=14,Color='k')
hold on
plot(Xout4(1,:),Xout4(2,:),'-o',LineWidth=1.5,Color='m')
hold on
plot(Xout4(1,end),Xout4(2,end),'*',LineWidth=2,MarkerSize=14,Color='k')
title("Without BFGS")

sgtitle("SQP with line search with and without BFGS")
