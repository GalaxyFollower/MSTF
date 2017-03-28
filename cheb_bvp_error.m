% Solves BVP with nonhomogeneous Dirichlet and Neumann boundary
% conditions and plots error.

clear all
close all

set(0,'DefaultLineLineWidth',1.5, ...
      'DefaultAxesLineWidth',1.5, ...
      'DefaultAxesFontSize',14, ...
      'DefaultTextFontSize',14, ...
      'DefaultTextInterpreter', 'latex', ...
      'DefaultAxesTickLabelInterpreter','latex');

n = 2;
for k = 1:8

    [D, x] = cheb(n);
    D2 = D*D;
    A = D2(2:n,2:n);
    b = exp(4*x(2:n));
    u = [0; A\b; 0];

    u_exact = (exp(4*x) - x*sinh(4) - cosh(4))/16;
    err = norm(u - u_exact)
    
    loglog(n, err, 'ok')
    hold on
    
    n = n*2;

end

xlabel('$n$')
ylabel('Error')
axis square
grid on
hold off
print -deps ../Figures/cheb_bvp_error.eps

