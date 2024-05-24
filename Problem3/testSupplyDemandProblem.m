clear
clc
close all
set(0,'defaultTextInterpreter','latex');
set(groot,'defaultAxesFontSize',16)

% Load system

% Load the problem data from the file LP_Test.mat
load('LP_Test.mat'); % This will load C, U, Pd_max, Pg_max

% The objective function coefficients need to be constructed based on your U and C vectors
% Since we are maximizing, we negate the objective function for linprog
g = [-U; C]; 

% The equality constraint (power balancing) seems to be the sum of p_d - sum of p_g = 0
% Assuming the number of p_d variables is equal to the length of U and the same for p_g and C
A = [ones(1, length(U)), -ones(1, length(C))];
b = 0;

% Lower bounds are zeros for both p_d and p_g, upper bounds are Pd_max and Pg_max respectively
lb = [zeros(length(U), 1); zeros(length(C), 1)]; 
ub = [Pd_max; Pg_max];

% Put system on normal form

[m, n] = size(A);

Abar = [A -A zeros(m,n) zeros(m,n);
        eye(n) -eye(n) -eye(n) zeros(n);
        eye(n) -eye(n) zeros(n) eye(n)];

bbar = [b;lb;ub];

gbar = [g;-g;zeros(n,1);zeros(n,1)];

[M, N] = size(Abar);

%% Linprog

% Solve the linear program using the dual-simplex algorithm
options = optimoptions('linprog','Algorithm','dual-simplex');
tic;
[xlinprog, fval, exitflag, output, lambda] = linprog(g, [], [], A, b, lb, ub, options);
t = toc;
% Extract the solution for p_d and p_g from the result xlinprog
p_d = xlinprog(1:length(U));
p_g = xlinprog(length(U)+1:end);

%% Interior point

% Solving

tol = 1.0e-9;
tic;
[xIPtemp,lambda,s,iterIP,info,rcres,rbres,mures] = lpsolverInteriorPoint(gbar,Abar,bbar,tol);
tIP = toc;

xIP = xIPtemp(1:n)-xIPtemp((n+1):(2*n));

% Extract the solution for p_d and p_g from the result xlinprog
p_dIS = xIP(1:length(U));
p_gIS = xIP(length(U)+1:end);

%% Active set

% Find basic feasible initial point
% Phase 1

tic;
A1 = [Abar ones(M,1) -eye(M) zeros(M);
      -Abar ones(M,1) zeros(M) -eye(M)];

b1 = [bbar;-bbar];

g1 = [zeros(N,1);1;zeros(2*M,1)];

t0 = max(abs(bbar));

x01 = [zeros(N,1);t0;t0-bbar;t0+bbar];

Bset1 = find(x01 > 0);
Nset1 = setdiff(1:length(x01),Bset1);

maxiter = 1000;
[sol1,~,iter1] = lpsolverActiveSet(g1,A1,b1,x01,Bset1,Nset1,maxiter);

% Phase 2

x0 = sol1(1:N);

Bset = ((N+1-M):N)';
Nset = setdiff(1:length(x0),Bset);

[xAStemp,~,iterAS] = lpsolverActiveSet(gbar,Abar,bbar,x0,Bset,Nset,maxiter);
tAS = toc;

xAS = xAStemp(1:n)-xAStemp((n+1):(2*n));

% Extract the solution for p_d and p_g from the result xlinprog
p_dAS = xAS(1:length(U));
p_gAS = xAS(length(U)+1:end);

%% Plotting

figure
subplot(1,3,1)
plotSupplyDemand(p_d,p_g,'Supply-demand curves with linprog \newline')

subplot(1,3,2)
plotSupplyDemand(p_dAS,p_gAS,'Supply-demand curves with AS \newline')

subplot(1,3,3)
plotSupplyDemand(p_dIS,p_gIS,'Supply-demand curves with IP \newline')

% Residuals from IP
figure
subplot(1,3,1)
plot(vecnorm(rcres,2),'-o',LineWidth=1.5)
xlabel("Iterations")
ylabel("$||r_{c}||_2$")
grid on
title(sprintf("$r_{c}$"))

subplot(1,3,2)
plot(vecnorm(rbres,2,1),'-o',LineWidth=1.5)
xlabel("Iterations")
ylabel("$||r_{b}||_2$")
grid on
title("$r_b$",Interpreter="latex")

subplot(1,3,3)
plot(mures,'-o',LineWidth=1.5)
xlabel("Iterations")
ylabel("$\mu$")
grid on
title("$\mu$")
sgtitle("Residuals w.r.t. iterations")

disp("Max error, iter and time of IP:")
disp(max(abs(xIP-xlinprog)))
disp(iterIP)
disp(tIP)

disp("Max error, iter and time of AS:")
disp(max(abs(xAS-xlinprog)))
disp(iterAS)
disp(tAS)

disp("Iter and time of linprog")
disp(output.iterations)
disp(t)
