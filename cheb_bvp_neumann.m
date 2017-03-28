% CHEB_BVP_NEUMANN solves BVP with nonhomogeneous Dirichlet 
% and Neumann boundary conditions.

clear all
close all

set(0,'DefaultLineLineWidth',1.5, ...
      'DefaultAxesLineWidth',1.5, ...
      'DefaultAxesFontSize',14, ...
      'DefaultTextFontSize',14, ...
      'DefaultTextInterpreter', 'latex', ...
      'DefaultAxesTickLabelInterpreter','latex');
  
n = 20;

[D, x] = cheb(n);
D2 = D*D;
A = [D(1,:); D2(2:n,:); zeros(1,n) 1];
b = [1; 0.25*exp(2*(x(2:n) + 1)); 1];
u = A\b;

x = 0.5*(1+x);
xx = linspace(0, 1, 100);
u_exact = (exp(4*xx) - 4*xx*(exp(4) - 8) + 15)/16;

plot(xx, u_exact, '-k', x, u, 'ok')
xlabel('$x$')
ylabel('$u$')
axis square
grid on

print -deps ../Figures/cheb_bvp_neumann.eps
