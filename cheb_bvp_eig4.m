% CHEB_BVP_EIG4 solves fourth-order eigenvalue problem.

clear all
close all

set(0,'DefaultLineLineWidth',1.5, ...
      'DefaultAxesLineWidth',1.5, ...
      'DefaultAxesFontSize',14, ...
      'DefaultTextFontSize',14, ...
      'DefaultTextInterpreter', 'latex', ...
      'DefaultAxesTickLabelInterpreter','latex');
  
n = 21;

% Setup and solve eigenvalue problem

[D, x] = cheb(n);
D2 = D*D; D3 = D*D2; D4 = D*D3;
A = (diag(1-x.^2)*D4 - 8*diag(x)*D3 - 12*D2)*diag(1./(1-x.^2));
A = A(2:n,2:n);
[v, lam] = eig(A);

% Add boundary values of u in eigenvectors

u = [ zeros(1,n-1); v; zeros(1,n-1) ];

% Sort eigenvalues and eigenvectors in order of inincreasing
% magnitude

[beta, k] = sort(2*diag(lam).^0.25);
u = u(:,k);

% Estimate error of eigenvalues

err = abs(cos(beta).*cosh(beta)-1);
semilogy(err,'o')
axis([0 floor(n/2) 1e-15 1e20])
axis square
xlabel('Mode')
ylabel('Error')

print -deps ../Figures/cheb_bvp_eig4_val.eps

% Plot first four eigenvectors and compare with exact solution.
% Normalise so the amplitude is the same.

x = 0.5*(1+x)
beta_exact = fsolve(@(z) cosh(z).*cos(z) - 1, beta(1:4))
xx = linspace(0, 1, 100);
[B, X] = meshgrid(beta_exact(1:4), xx);
u_exact = - (sinh(B.*X) - sin(B.*X))./(sinh(B) - sin(B)) ...
        + (cosh(B.*X) - cos(B.*X))./(cosh(B) - cos(B));

for k = 1:4
    subplot(2,2,k)
    plot(xx,u_exact(:,k)/max(abs(u_exact(:,k))),'-k', ...
          x,u(:,k)/max(abs(u(:,k))),'ok')
    xlabel('$x$')
    ylabel('$u$')
    title(sprintf('%s%4.2f','$\beta=$ ',beta(k)))
    axis square
end

print -deps ../Figures/cheb_bvp_eig4_vec.eps










