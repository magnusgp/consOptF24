function [c,ceq] = mycon(x)
    c = [];                      % Compute nonlinear inequalities at x.
    ceq = (x(1) + 2).^2 - x(2);  % Compute nonlinear equalities at x.
end