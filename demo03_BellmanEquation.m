
clear, clc

%% Model parameters
beta = 0.95;
r = 0.02;
y = 1;
R = 1+r;

%% Basis parameters
a = 0.0;
b = 4;
n = 12;
B = basisChebyshev(n,a,b,'gaussian','assets');

%% Asset nodes and initial guess V(A) = A
A0 = B.nodes;  % Initial nodes
Vnodes = A0;
Phi = B.Interpolation;
coef = Phi\Vnodes;

% Optimal assets for given nodes
A1 = A0;  %just to allocate memory


%% Value function iteration
tic
for it = 1:2000
    
    % Keep track of old value function
    coef_old = coef;
    
    % Solving maximization problem for each A0 node
    for k = 1:n
        f = @(A1) - log(A0(k) + y-A1/R) - beta*B.Evaluate(coef,A1);   % minimization problem
        [A1(k), fval] = fminbnd(f,0,R*A0(k));
        Vnodes(k) = -fval;
    end
    
    % Update value function
    coef = Phi\Vnodes;
    
    % Check for convergence, exit loop if solution found
    dcoef = norm(coef-coef_old);
    fprintf('%5d  %8.3e %8.2f\n',it,dcoef,toc)
    
    if dcoef < 1.0e-8
        fprintf('Solution found after %d iterations\n',it)
        break 
    end
end


%% Value function iteration, but solving FOC
coef2 = Phi\A0;

tic
for it = 1:2000
    
    % Keep track of old value function
    coef_old = coef2;
    
    % Solving maximization problem for each A0 node
    for k = 1:n
        f2 = @(A1) R*beta * (A0(k) + y -A1/R) * B.Evaluate(coef2,A1,1) - 1;   % minimization problem
        A1(k) = fzero(f2,A0(k));
        Vnodes(k) = log(A0(k) + y - A1(k)/R) + beta * B.Evaluate(coef2,A1(k));
    end
    
    % Update value function
    coef2 = Phi\Vnodes;
    
    % Check for convergence, exit loop if solution found
    dcoef = norm(coef2-coef_old);
    fprintf('%5d  %8.3e %8.2f\n',it,dcoef,toc)
    
    if dcoef < 1.0e-8
        fprintf('Solution found after %d iterations\n',it)
        break 
    end
end









%% Plot solution and approximation
AA = linspace(a,b,121)';

V = @(A0) log(A0)/(1-beta) + (beta*log(R) + beta*log(beta) + (1-beta)*log(1-beta)) / (1-beta)^2;
Vapprox = @(A0) B.Evaluate(coef,A0);

subplot(2,1,1)
plot(AA,[V(AA)  Vapprox(AA)])

subplot(2,1,2)
plot(AA,V(AA) - Vapprox(AA))




