%% DEMO: basisChebyshev
%
% This script provides examples on using the <basisChebyshev.m basisChebyshev> class.
%
% Implemented by:
%
% * *Randall Romero-Aguilar*
% * |randall.romero@outlook.com|
% * Last updated: September 23, 2014


%% EXAMPLE 2: Using |basis| to approximate a 2-D Chebyshev basis
%
% *PROBLEM*: Approximate the Cobb-Douglas production function y = f(K,L) = K^0.4 L^0.6 for capital K in [0,2] and labor
% L in [0,1].
%
% * First, create the production function
y = @(K,L) K.^0.4 .* L.^0.6;

%% Using a tensor-product basis 
% Create a 2-D Chebyshev basis |T| by taking the tensor product of the individual dimensions
expand.method = 'tensor';
opts.nodetype = 'lobatto';

T = basis('cheb',[32 20],[0 0],[2 1],{'capital','labor'},opts,expand)

%%
% The coefficients for the function are then:
c_tensor = T.Interpolation\y(T.nodes(:,1),T.nodes(:,2));

%%
% Plot the functions
ngrid = 201;
kk = linspace(T.a(1),T.b(1),ngrid)';
ll = linspace(T.a(2),T.b(2),ngrid)';
[kk,ll] = meshgrid(kk,ll);

yy = y(kk,ll);
yyapprox = reshape(T.Evaluate(c_tensor,[kk(:),ll(:)]),ngrid,ngrid);

figure
subplot(2,1,1), surf(kk,ll,yy,'EdgeColor','none')
hold on
surf(kk,ll,yyapprox,'EdgeColor','none')
axis tight
xlabel(T.varname(1)), ylabel(T.varname(2)), zlabel('Production')


subplot(2,1,2), surf(kk,ll,yyapprox - yy,'EdgeColor','none')
title('Residuals using 640 nodes, tensor basis'), axis tight
xlabel(T.varname(1)), ylabel(T.varname(2)), zlabel('Residual')

%% Using a Smolyak basis
% Now repeat the exercise with Smolyak nodes, with degree and node parameters set to q = 3. Call the basis |S|
%
expand.method = 'smolyak';
expand.degreeParam = 6;
expand.nodeParam = 6;

S = basis('cheb',[80 50],[0 0],[2 1],{'capital','labor'},opts,expand)

%%
% The coefficients for the function are then:
c_smolyak = S.Interpolation\y(S.nodes(:,1),S.nodes(:,2));

%%
% Plot the functions
yyapprox = reshape(S.Evaluate(c_smolyak,[kk(:),ll(:)]),ngrid,ngrid);

figure
subplot(2,1,1), surf(kk,ll,yy,'EdgeColor','none')
hold on
surf(kk,ll,yyapprox,'EdgeColor','none')
axis tight
xlabel(S.varname(1)), ylabel(S.varname(2)), zlabel('Production')


subplot(2,1,2), surf(kk,ll,yyapprox - yy,'EdgeColor','none')
title('Residuals using 321 nodes, tensor basis'), axis tight
xlabel(S.varname(1)), ylabel(S.varname(2)), zlabel('Residual')

%%
% Notice that the residuals in Smolyak and tensor bases are of similar magnitude, despite the fact that Smolyak is using
% around half the number of nodes and bases functions (i.e., its interpolation matrix is around one-forth the size of
% the tensor interpolating matrix).
%
% Also notice the warning about the number of Smolyak nodes: The number of nodes for each dimension must be n = 2^k+1
% for some positive integer k. If user calls the constructor with different number of nodes, the constructor increases n
% to the next admissible value (from 80 to 129 for capital, and from 50 to 65 for labor).