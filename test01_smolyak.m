%% Smolyak nodes: replicating Judd et al 2014 plots
%
% This file replicates figures 1 and 4 of Judd et al 2014 paper on JEDC 44, pp 92-123.
% The purpose of this replication is simply to test the Matlab implementation of the
% Smolyak algorithm, which is included in the class <basis.m basis>
%
% Copyright (C) 2014 Randall Romero-Aguilar
%
% Licensed under the MIT license, see LICENSE.txt

%% Figure 1: Isotropic grids

opts.method = 'smolyak';
opts.nodetype = 'lobatto';
a = [-1 -1];
b = [1 1];
n = 9;

r =[-1 1];

figure('Position',[100 100 1050 200])
for mu=0:3
    opts.nodeParam = mu;
    B = basis(n,a,b,opts);
    subplot(1,5,mu+1)
    plot(B.nodes(:,1),B.nodes(:,2),'.','MarkerSize',9)
    title(['\mu = ',num2str(mu)])
    xlim(r),ylim(r)
end

B = basis(5,a,b);
subplot(1,5,5)
plot(B.nodes(:,1),B.nodes(:,2),'.','MarkerSize',9)
xlim(r),ylim(r)
title('5 x 5')

%% Figure 4 Anisotropic grids
n = [9 3];
muC = [1 0; 2 1; 3 1];

figure('Position',[100 100 750 200])
for k=1:3
    mu = muC(k,:);
    opts.nodeParam = mu;
    B = basis(n,a,b,opts);
    subplot(1,3,k)
    plot(B.nodes(:,1),B.nodes(:,2),'.','MarkerSize',9)
    title(['(\mu_1,\mu_2) = (',num2str(mu(1)),', ',num2str(mu(2)),')'])
    xlim(r),ylim(r)
end


%% Reference
% Judd, Maliar, Maliar and Valero 2014 Smolyak method for solving dynamic economic models:
% Lagrange interpolation, anisotropic grid and addaptive domain. Journal of Economic
% Dynamics & Control 44, pp. 92-13