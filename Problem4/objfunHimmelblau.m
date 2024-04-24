function [f,dfdx] = objfungradHimmelblau(x,p)
tmp1 = x(1)*x(1)+x(2)-11;
tmp2 = x(1)+x(2)*x(2)-7;
f = tmp1*tmp1 + tmp2*tmp2;

% compute the gradient of f
if nargout > 1
    dfdx = zeros(2,1);
    dfdx(1,1) = 4*tmp1*x(1) + 2*tmp2;
    dfdx(2,1) = 2*tmp1 + 4*tmp2*x(2);
end
