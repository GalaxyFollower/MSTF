% FFT_ADVDIFF solves advection-diffusion equation using 
% Fourier pseudospectral method.

N = 16;
c = 10;
nu = 1;
t = 1;

x = 2*pi*(0:N-1)'/N;
f = exp(sin(x));

ik = 1i*[0:N/2-1 0 -N/2+1:-1]';
k2 = [0:N/2-1 -N/2:-1]'.^2;

fh = fft(f);
uh = fh.*exp(-(c*ik + nu*k2)*t);
u = real(ifft(uh));

plot(x,u) 