% Creates a sequence of images 

close all
clear all

set(groot, 'DefaultLineLineWidth', 1, ...
    'DefaultAxesLineWidth', 1, ...
    'DefaultAxesFontSize', 10, ...
    'DefaultTextFontSize', 10, ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultAxesTickLabelInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultLegendFontSize', 10, ...
    'DefaultColorbarTickLabelInterpreter', 'latex', ...
    'DefaultColorbarFontSize', 10);
  
% Data

pathname = '/home/tmattner/Teaching/MSTF/2017/Assignments/Assign3/Data/';
jobname = 'turb_64_re200_'

% Simulation parameters.

n = 64;       % grid size
ndump = 100;  % number of dumps
idump = 1;    % initial dump number

% Grid

dx = 2*pi/n;
x = dx*[0:n-1];
[x, y, z] = ndgrid(x, x, x);

% Permute for plotting

x = permute(x, [2 1 3]);
y = permute(y, [2 1 3]);
z = permute(z, [2 1 3]);

% Wavenumbers

dk = 1;
k2 = [0:n/2-1 -n/2:-1].^2;
[kx2, ky2, kz2] = ndgrid(k2, k2, k2);
k2 = kx2 + ky2 + kz2;

ik = 1i*[0:n/2-1 0 -n/2+1:-1];
[ikx, iky, ikz] = ndgrid(ik, ik, ik);

% Read data.

while (idump <= ndump)
    
    % Read data

    filename = sprintf('%s/%s%3.3i%s', pathname, jobname, idump, '.h5')
    u = h5read(filename, '/u');
    v = h5read(filename, '/v');
    w = h5read(filename, '/w');
    kinvis = h5read(filename, '/kinvis');
    
    uh = fftn(u);
    vh = fftn(v);
    wh = fftn(w);
    
    dudx = ifftn(ikx.*uh, 'symmetric');
    dudy = ifftn(iky.*uh, 'symmetric');
    dudz = ifftn(ikz.*uh, 'symmetric');
    dvdx = ifftn(ikx.*vh, 'symmetric');
    dvdy = ifftn(iky.*vh, 'symmetric');
    dvdz = ifftn(ikz.*vh, 'symmetric');
    dwdx = ifftn(ikx.*wh, 'symmetric');
    dwdy = ifftn(iky.*wh, 'symmetric');
    dwdz = ifftn(ikz.*wh, 'symmetric');
    
    % Enstrophy
        
    wx = dwdy - dvdz;
    wy = dudz - dwdx;
    wz = dvdx - dudy;
    w2 = wx.^2 + wy.^2 + wz.^2;

    % Plot enstrophy
    
    w2 = permute(w2, [2 1 3]);
    w2mean = mean(w2(:));
    clf
    h = slice(x, y, z, w2, [2*pi-dx], [2*pi-dx], [dx]);
    hold on
    for k = 1:length(h)
        h(k).FaceColor = 'interp';
        h(k).EdgeColor = 'none';
    end
    colorbar
    title(colorbar, '$\langle \omega_i \omega_i \rangle$', ...
        'interpreter', 'latex', 'fontsize', 10)
    caxis([0 4*w2mean])
    isosurface(x, y, z, w2, 6*w2mean);
    isocaps(x, y, z, w2, 6*w2mean);
    camlight headlight
    axis([0 2*pi 0 2*pi 0 2*pi])
    axis vis3d
    box on
    xticks([0 2 4 6])
    yticks([0 2 4 6])
    zticks([0 2 4 6])
    xlabel('$y$')
    ylabel('$x$')
    zlabel('$z$')
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    hold off
    
    drawnow      
    filename = sprintf('%s/%s%3.3i%s', pathname, jobname, idump, '_w2.jpg');
    print(filename, '-djpeg', '-r150')
    
    idump = idump + 1
    
end
