function wj = baryfit(xj, fj)

n = length(xj);
[Xj, X] = meshgrid(xj,xj);
wj = 1./prod(X - Xj + eye(n),2); % column

