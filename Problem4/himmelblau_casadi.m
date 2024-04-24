import casadi.*

% Symbols/expressions
x1 = MX.sym('x1');
x2 = MX.sym('x2');

t1 = x1*x1 + x2 - 11;
t2 = x1 + x2*x2 - 7;
f = t1*t1 + t2*t2;

c1 = (x1+2)^2 - x2;
c2 = -4*x1 + 10*x2;
g = [c1; c2];

nlp = struct; % NLP declaration
nlp.x = [x1;x2]; % decision vars
nlp.f = f; % objective
nlp.g = g; % constraints

% Create solver instance
F = nlpsol('F','ipopt',nlp);

% Solve the problem using a guess
res = F('x0',[0.0 0.0],'ubg',1e8,'lbg',0,'lbx',[-5;-5],'ubx',[5;5]);

res.x
