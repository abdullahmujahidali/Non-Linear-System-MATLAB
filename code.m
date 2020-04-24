
% system of equations
function f = F(p)
x = p(1); y = p(2);
%
f1 = y + x^2 -x -0.75;
f2 = y + 5*x*y - x^2;
%
f = [f1;f2];
end

%________________________


% Jacobi
function j = J(p)
x = p(1); y = p(2);
% partial derivatives
f1dx = 2*x - 1;
f1dy = 1;
%
f2dx = 5*y - 2*x;
f2dy = 5*x + 1;
%% Jacobian matrix
%
j = [f1dx f1dy;
f2dx f2dy];
end

%______________________

function [solution,iter] = NewtonRaphson(F,J,p0,tol)

p = p0;
e = Inf; % intitial error
%%
iter = 0;
while e>tol
% newton's
P = p - (J(p)\F(p));
e = norm(P-p); %error
p = P;
iter = iter + 1;
end

solution = p;
end

%________________________
