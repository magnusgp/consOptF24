function [F, dF] = quadraticFdF(z,G,c,A,b)

% H = matrix for quadratic term
% 
% g = vector for linear term
% 
% A = matrix for constraints
% 
% b = RHS of constrains

sysmat = [G -A';A zeros(size(A,1),size(A,1))];

F = sysmat*z+[c;-b];

dF = sysmat;
