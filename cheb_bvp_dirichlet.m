% CHEB_BVP_DIRICHLET olves BVP with homogeneous Dirichlet boundary conditions

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
A = D2(2:n,2:n);
b = exp(4*x(2:n));
u = [0; A\b; 0];

xx = linspace(-1, 1, 100);
u_exact = (exp(4*xx) - xx*sinh(4) - cosh(4))/16;

plot(xx, u_exact, '-k', x, u, 'ok')
xlabel('$x$')
ylabel('$u$')
axis square
grid on

print -deps ../Figures/cheb_bvp_dirichlet.eps
