clear all
close all

set(0,'DefaultLineLineWidth',1.5, ...
      'DefaultAxesLineWidth',1.5, ...
      'DefaultAxesFontSize',14, ...
      'DefaultTextFontSize',14, ...
      'DefaultTextInterpreter', 'latex', ...
      'DefaultAxesTickLabelInterpreter','latex');
  
n = 20;
xj = linspace(-1,1,n+1)';
uj = 1./(1+16*xj.^2);
x = linspace(-1,1,10*n)';
wj = baryfit(xj,uj);
u = baryval(x,xj,uj,wj);

plot(x,u,'k-',xj,uj,'ko')
xlabel('$x$')
ylabel('$u$')
axis([-1 1 -0.5 1.5])
axis square
grid on
% print -deps ../Figures/runge_uniform.eps

xj = cos(pi/n*(0:n)');
uj = 1./(1+16*xj.^2);
x = linspace(-1,1,10*n)';
wj = baryfit(xj,uj);
u = baryval(x,xj,uj,wj);

figure
plot(x,u,'k-',xj,uj,'ko')
xlabel('$x$')
ylabel('$u$')
axis([-1 1 -0.5 1.5])
axis square
grid on
% print -deps ../Figures/runge_cheb.eps