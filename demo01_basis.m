%% DEMO: basis
%
% This script provides a few examples on using the <basis.m basis> class.
%
%
% Implemented by:
%
% * *Randall Romero-Aguilar*
% * |randall.romero@outlook.com|
% * Last updated: September 23, 2014

%% EXAMPLE 1: Using |basisChebyshev| to approximate a 1-D Chebyshev basis
%
% *PROBLEM*: Interpolating the function y = f(x) = 1 + sin(2*x) on the domain [0,pi], using 5 Gaussian nodes.
%
% * First, create the function |f|
f = @(x) 1 + sin(2*x);

%%
% * and the Chebyshev basis |B|. If the type of nodes is unspecified, Gaussian is computed by default
B = basisChebyshev(5,0,pi)

%% Interpolation matrix and nodes
% * Obtain the interpolation matrix |Phi|, evaluated at the basis nodes.
Phi = B.Interpolation

%%
% * The basis nodes are:
xnodes = B.nodes

%% Fitting a function
% * To set the interpolation coefficients |c|:
c = Phi\f(xnodes)

%%%
% or simply
%
%   c = B.Interpolation\f(B.nodes)


%% 
% Next plot the function |f| and its approximation. To evaluate the function defined by the basis |B| and coefficients
% |c| at values |xx| we use the Evaluate method:
%
%   y_approx = B.Evaluate(c,xx)
%
% We also plot the residuals, showing the residuals at the interpolating nodes (zero by construction)

xx = linspace(B.a,B.b,121)';
figure
subplot(2,1,1), plot(xx,[f(xx), B.Evaluate(c,xx)])
axis tight
legend('f = 1 + sin(2x)','approx.')

subplot(2,1,2), plot(xx, B.Evaluate(c,xx) - f(xx),'r')
title('Residuals using 5 nodes'), axis tight
hold on, plot(B.nodes,f(B.nodes) - B.Evaluate(c),'bo')

%% Adjusting the number of nodes
% To increase accuracy, we increase the number of nodes in |B| to 25.
B.n = 25;
c2 = B.Interpolation\f(B.nodes);
figure
subplot(2,1,1), plot(xx,[f(xx), B.Evaluate(c2,xx)])
axis tight
legend('f = 1 + sin(2x)','approx.')

subplot(2,1,2), plot(xx, B.Evaluate(c2,xx) - f(xx),'r')
title('Residuals using 25 nodes'), axis tight
