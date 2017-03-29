clear all;
close all;

fsode = @(x,y,beta) [ stuff ];
fsbc = @(ya,yb) [ stuff ];
fsinit = @(x) [ stuff ];

a = 0;
b = 8;
beta = 0;

solinit = bvpinit(linspace(a,b,10), fsinit);
sol = bvp4c(@(x,y) fsode(x,y,beta), fsbc, solinit);

x = linspace(a, b, 100);
y = deval(sol, x);
plot(x, y(2,:))

