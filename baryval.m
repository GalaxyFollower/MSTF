function y = baryval(x, xj, fj, wj)

[Xj, X] = meshgrid(xj, x);
dX = X - Xj;
[i, j] = find(dX == 0);
dX(i,j) = 1;
dX = 1./dX;
y = wj.*fj(:);
y = (dX*y)./(dX*wj);
y(i) = fj(j);
