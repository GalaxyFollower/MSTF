% FFT_BURGERS solves Burgers equation using Fourier pseudospectral method.

N = 32;
t = linspace(0, 2, 25);
nu = 0.1;
f = @sin;

dx = 2*pi/N;
x = dx*(0:N-1)';

ik = 1i*[0:N/2-1 0 -N/2+1:-1]';
k2 = [0:N/2-1 -N/2:-1]'.^2;

burgers = @(t, uh) -fft( real(ifft(uh)) .* real(ifft(ik.*uh)) ) - nu*k2.*uh;

u0 = f(x);
uh0 = fft(u0);

[~, uh] = ode45(burgers, t, uh0);
u = real(ifft(uh.')); % Careful! ' by itself does a complex conjugate transpose!

plot(x, u(:,end))




