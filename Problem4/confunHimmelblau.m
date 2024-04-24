function [c,ceq] = confunHimmelblau(x,p)
c = zeros(2,1);
ceq = zeros(0,1);

% Inequality constraints c(x) <= 0
tmp = x(1)+2;
c(1,1) = -(tmp*tmp - x(2));
c(2,1) = -(-4*x(1) + 10*x(2));