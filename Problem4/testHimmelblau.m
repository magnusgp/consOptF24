
f = @(x) (x(1).^2+x(2)-11).^2 + (x(1) + x(2).^2 - 7).^2;

df = @(x) [4*x(1).^3 + (4*x(2) - 42)*x(1) + 2*x(2).^2 - 14;
           4*x(2).^3 + (4*x(1) - 26)*x(2) + 2*x(1).^2 - 22];

cE = @(x) (x(1) + 2).^2 - x(2);
cI = @(x) -4*x(1) + 10*x(2);

AE = @(x) [2*x(1)+4;
           -1];
AI = @(x) [-4;
           10];

dL = @(x) [4*x(1).^3 + (-2*y + 4*x(2) - 42)*x(1) + 2*x(2).^2 - 4*y + 4*z - 14;
           4*x(2).^3 + (4*x(1) - 26)*x(2) + 2*x(1).^2 + y - 10*z - 22];

ddL = @(x) [12*x1^2 + 4*x2 - 2*y - 42 4*x1 + 4*x2; 4*x1 + 4*x2 12*x2^2 + 4*x1 - 26];

% [x] = nlpsolverLineSearch(x0,s0,cE,cI,AE,AI,f,df,dL,ddL);

%%

N = 150;

x1start = -5;
x1end = 5;
x2start = -5;
x2end = 5;

x(1) = x1start:(x1end-x1start)/N:x1end;
x(2) = x2start:(x2end-x2start)/N:x2end;

[X1,X2] = meshgrid(x(1),x(2));

falt = @(x(1),x(2)) (x(1).^2+x(2)-11).^2 + (x(1) + x(2).^2 - 7).^2;

F = reshape(arrayfun(falt,X1(:),X2(:)),length(x(1)),length(x(2)));

v = 0:20:500;

% Equalities
c1 = @(xx1,xx2) (xx1 + 2).^2 - xx2;
c2 = @(xx1,xx2) -4*xx1 + 10*xx2 <= 0;

% C1 = reshape(arrayfun(c1,X1(:),X2(:)),length(x(1)),length(x(2)));
C2 = reshape(arrayfun(c2,X1(:),X2(:)),length(x(1)),length(x(2)));

colormap jet
contour(X1,X2,F,v,'linewidth',1.5, 'HandleVisibility','off');
colorbar

hold on

opacity = 0.1;
mkrcol = 'k';
mkrsize = 10;
scatter([X1(C2)],[X2(C2)],mkrsize,"square","filled",'MarkerFaceColor',mkrcol,'MarkerEdgeColor',mkrcol,'MarkerFaceAlpha',opacity,'MarkerEdgeAlpha',opacity, 'HandleVisibility','off')

hold on

fimplicit(c1,LineWidth = 2,Color='k')