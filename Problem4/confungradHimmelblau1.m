function [c,ceq,dcdx,dceqdx] = confungradHimmelblau1(x,p)
c = zeros(1,1);
ceq = zeros(0,1);

% Inequality constraints c(x) <= 0
tmp = x(1)+2;
c(1,1) = -(tmp*tmp - x(2));

% Compute constraint gradients
if nargout > 2
    dcdx = zeros(2,1);
    dceqdx = zeros(2,0);
    
    dcdx(1,1) = -2*tmp; % dc1dx1
    dcdx(2,1) = 1.0; % dc1dx2
end
