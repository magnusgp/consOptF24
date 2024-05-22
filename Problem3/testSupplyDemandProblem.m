clear
clc
close all

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

%% Linprog

% Solve the linear program using the dual-simplex algorithm
options = optimoptions('linprog','Algorithm','dual-simplex');
tic;
[xlinprog, fval, exitflag, output, lambda] = linprog(g, [], [], A, b, lb, ub, options);
t = toc;
% Extract the solution for p_d and p_g from the result xlinprog
p_d = xlinprog(1:length(U));
p_g = xlinprog(length(U)+1:end);

% get number of iterations 
num_iter = output.iterations;

% Display the results
disp('The optimal demand prices, p_d, are:');
disp(p_d);
disp('The optimal generation prices, p_g, are:');
disp(p_g);

% display cpu time and number of iterations
disp('CPU time: ');
disp(t);
disp('Number of iterations: ');
disp(num_iter);

% Post-processing steps for supply-demand curve and market clearing price can be added here
% sum p_d and p_g to check if market is cleared
assert(sum(p_d) == sum(p_g), 'Market not cleared');
% display 
disp('sum of p_d: ');
disp(sum(p_d));
disp('sum of p_g: ');
disp(sum(p_g));

plotSupplyDemand(p_d,p_g)

%% Interior point

% Finding initial point

xtilde = Abar'*inv(Abar*Abar')*bbar;
lambdatilde = inv(Abar*Abar')*Abar*gbar;
stilde = gbar - Abar'*lambdatilde;

deltax = max([-(3/2)*min(xtilde),0]);
deltas = max([-(3/2)*min(stilde),0]);

e = ones(length(xtilde),1);

xIP = xtilde + deltax*e;
shat = stilde + deltas*e;

deltaxhat = 0.5*(xIP'*shat)/(e'*shat);
deltashat = 0.5*(xIP'*shat)/(e'*xIP);

x0 = xIP + deltaxhat*e;
lambda0 = lambdatilde;
s0 = shat + deltashat*e;

% Solving

tol = 1.0e-9;
[xIPtemp,lambda,s,iterIP,info,rcres,rbres,mures] = lpsolverInteriorPoint(gbar,Abar,bbar,x0,lambda0,s0,tol);

xIP = xIPtemp(1:n)-xIPtemp((n+1):(2*n));

disp("Max error:")
disp(max(abs(xIP-xlinprog)))

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
title("$r_b$")

subplot(1,3,3)
plot(mures,'-o',LineWidth=1.5)
xlabel("Iterations")
ylabel("$\mu$")
grid on
title("$\mu$")
sgtitle("Residuals w.r.t. iterations")

% Extract the solution for p_d and p_g from the result xlinprog
p_d = xIP(1:length(U));
p_g = xIP(length(U)+1:end);

plotSupplyDemand(p_d,p_g)

%% Active set

% Find basic feasible initial point

E = zeros(length(bbar));

for i = 1:length(bbar)
    if bbar(i) >= 0
        E(i,i) = 1;
    else
        E(i,i) = -1;
    end
end

g1 = [zeros(length(gbar),1);ones(length(bbar),1)];

x01 = [zeros(length(gbar),1);abs(bbar)];

A1 = [Abar E];

idx = 1:length(x01);
Nset = idx(1:length(gbar));
Bset = idx((length(gbar)+1):end);

[sol1,~,iter] = lpsolverActiveSet(g1,A1,bbar,x01,Bset,Nset);

options = optimoptions('linprog','Algorithm','dual-simplex');
sol1linprog = linprog(g1, [], [], A1, bbar, zeros(length(g1),1), [], options);

disp("here")
disp(g1'*sol1)
disp(g1'*sol1linprog)

g2 = [gbar;zeros(length(bbar),1)];

x02 = sol1;

A2 = [Abar eye(length(bbar))];

idx = 1:length(x02);
Nset = idx(1:length(gbar));
Bset = idx((length(gbar)+1):end);

[xAStemp,XAStemp,iterAS] = lpsolverActiveSet(g2,A2,bbar,x02,Bset,Nset);

xAStemp = xAStemp(1:length(gbar));

xAS = xAStemp(1:n)-xAStemp((n+1):(2*n));

disp("Max error:")
disp(max(abs(xAS-xlinprog)))

%%

% Extract the solution for p_d and p_g from the result xlinprog
p_d = xAS(1:length(U));
p_g = xAS(length(U)+1:end);

plotSupplyDemand(p_d,p_g)