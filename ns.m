% NS solves Navier-Stokes equations in a triply 2*pi periodic domain.

close all
clear all

% Coefficients for third-order low-storage Runge-Kutta scheme.

a = [ 0.0 -5.0/9.0  -153.0/128.0 ];
b = [ 1.0/3.0 15.0/16.0 8.0/15.0 ];

% Output parameters

pathname = '/home/tmattner/Teaching/MSTF/Notes/Run';
jobname = 'test'

% Simulation parameters.

n = 32;       % grid size
kinvis = 0.2; % kinematic viscosity
cfl = 0.5;    % time-step restriction
tdump = 0.1;  % time between data dumps
ndump = 100;  % number of dumps
idump = 0;    % initial dump number

% Grid

dx = 2*pi/n;
x = dx*[0:n-1];
[x, y, z] = ndgrid(x, x, x);

% Wavenumbers

dk = 1;
k2 = [0:n/2-1 -n/2:-1].^2;
[kx2, ky2, kz2] = ndgrid(k2, k2, k2);
k2 = kx2 + ky2 + kz2;

ik = 1i*[0:n/2-1 0 -n/2+1:-1];
[ikx, iky, ikz] = ndgrid(ik, ik, ik);

% Initial conditions. Read from file, if it exists.

filename = sprintf('%s/%s%3.3i%s', pathname, jobname, idump, '.h5')

if (exist(filename) == 2)
    
    u = h5read(filename, '/u');
    v = h5read(filename, '/v');
    w = h5read(filename, '/w');
    t = h5read(filename, '/t');
    
else

    u = sin(z) + cos(y);
    v = sin(x) + cos(z);
    w = sin(y) + cos(x);       
    t = 0;

end

% Initialise variables.

ush = fftn(u);
vsh = fftn(v);
wsh = fftn(w);
uh = ush;
vh = vsh;
wh = wsh;
ts = t;

while (t < tdump*ndump)
    
    % Time-step is limited by dt*u/dx < cfl for stability and
    % dt*kinvis/dx^2 < cfl for accuracy. 
    
    dt = cfl*min([dx^2/kinvis dx/max(abs(u(:)) + abs(v(:)) + abs(w(:)))]);
    dt = min([dt, tdump*(idump + 1) - t, tdump ]); 
    tic
    
    % Runge-Kutta step.
    
    for i = 1:length(a)
        
        % Nonlinear terms.
        
        u = ifftn(uh);
        v = ifftn(vh);
        w = ifftn(wh);
        
        dudx = ifftn(ikx.*uh);
        dudy = ifftn(iky.*uh);
        dudz = ifftn(ikz.*uh);
        dvdx = ifftn(ikx.*vh);
        dvdy = ifftn(iky.*vh);
        dvdz = ifftn(ikz.*vh);
        dwdx = ifftn(ikx.*wh);
        dwdy = ifftn(iky.*wh);
        dwdz = ifftn(ikz.*wh);
        
        hxh = 0.5*fftn(u.*dudx + v.*dudy + w.*dudz);
        hyh = 0.5*fftn(u.*dvdx + v.*dvdy + w.*dvdz);
        hzh = 0.5*fftn(u.*dwdx + v.*dwdy + w.*dwdz);
        
        hxh = hxh + 0.5*(ikx.*fftn(u.*u) + iky.*fftn(u.*v) + ikz.*fftn(u.*w));
        hyh = hyh + 0.5*(ikx.*fftn(v.*u) + iky.*fftn(v.*v) + ikz.*fftn(v.*w));
        hzh = hzh + 0.5*(ikx.*fftn(w.*u) + iky.*fftn(w.*v) + ikz.*fftn(w.*w));
        
        % Runge-Kutta substep.
        
        ts = a(i)*ts + dt;
        visc = exp(-kinvis*b(i)*ts.*k2);
        ush = visc.*(a(i)*ush + dt*(-hxh));
        vsh = visc.*(a(i)*vsh + dt*(-hyh));
        wsh = visc.*(a(i)*wsh + dt*(-hzh));
                
        t = b(i)*ts + t;        
        uh = b(i)*ush + visc.*uh;
        vh = b(i)*vsh + visc.*vh;
        wh = b(i)*wsh + visc.*wh;
        
        % Pressure projection.
        
        ph = -(ikx.*uh + iky.*vh + ikz.*wh)./k2;
        ph(1,1,1) = 0;
        uh = uh - ikx.*ph;
        vh = vh - iky.*ph;
        wh = wh - ikz.*ph;
        
    end
    
    u = ifftn(uh);
    v = ifftn(vh);
    w = ifftn(wh);
    
    toc
    
    % Dump data. Will fail if files already exist.
    
    if (abs(t - tdump*(idump + 1)) < 1e-8)
        
        idump = idump + 1;
        
        filename = sprintf('%s/%s%3.3i%s', pathname, jobname, idump, '.h5')
        h5create(filename, '/u', size(u))
        h5create(filename, '/v', size(v))
        h5create(filename, '/w', size(w))
        h5create(filename, '/t', size(t))
        h5create(filename, '/kinvis', size(kinvis))

        h5write(filename, '/u', gather(u))
        h5write(filename, '/v', gather(v))
        h5write(filename, '/w', gather(w))
        h5write(filename, '/t', gather(t))
        h5write(filename, '/kinvis', gather(kinvis))
        
    end
    
    % Monitor solution, compute stats, plot (but not for large jobs), etc
    
    err = (abs(u - exp(-kinvis*t)*(sin(z) + cos(y))) ...
        +  abs(v - exp(-kinvis*t)*(sin(x) + cos(z))) ...
        +  abs(w - exp(-kinvis*t)*(sin(y) + cos(x))));
    errmax = max(err(:))

    tke = 0.5*sum(u(:).^2 + v(:).^2 + w(:).^2)/n^3 
    
end


