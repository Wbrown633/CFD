% Solves the 2D heat equation with an explicit finite difference scheme
% Code modified from template provided by: 
% http://geodynamics.usc.edu/~becker/teaching/557/problem_sets/problem_set_
% fd_2dheat.pdf
% Many thanks! 

clear
%Physical parameters
L = .5; % Width of plate [m]
H = .5; % Height of plate [m]
Tbot = 20; % Temperature of bottom wall [C]
Twall = 10; % Temperature of left wall [C]
Tinit = 10; % Initial Temperature guess [C]

alpha = 9.7e-5; % Thermal diffusivity of Al [m2/s]

% Numerical parameters
nx = 100; % # gridpoints in x-direction
ny = 100; % # gridpoints in z-direction
nt = 20000; % Number of timesteps to compute
dx = L/(nx-1); % Spacing of grid in x-direction
dy = H/(ny-1); % Spacing of grid in z-direction
[x2d,z2d] = meshgrid(-L/2:dx:L/2, -H:dy:0); % create grid
% Compute stable timestep
dt = min([dx,dy])^2/alpha/4;
% Setup initial linear temperature profile
T = zeros(ny,nx);

T(:,:) = Tinit;
time = 0;
for n=1:nt
    % Compute new temperature
    Tnew = zeros(ny,nx);
    Tnew(:,:) = Tinit;
    sx = alpha*dt/dx^2;
    sy = alpha*dt/dy^2;
    for j=2:nx-1
        for i=2:ny-1
            Tnew(i,j) = T(i,j) + sx*(T(i,j+1) -2*T(i,j)+ T(i,j-1)) +sy*(T(i+1,j)-2*T(i,j) + T(i-1,j));
        end
    end
% Set boundary conditions
Tnew(1,:) = Tbot;
Tnew(:,1) = Twall;
for i=2:ny-1
    Tnew(i,nx) = Tnew(i,nx-1);
end

for j=2:nx-1
    Tnew(ny,j) = Tnew(ny-1,j);
end 

T = Tnew;
time = time+dt;
%Plot solution every 50 timesteps
if (mod(n,50)==0)
    figure(1), clf
    pcolor(x2d/1e3,z2d/1e3,Tnew); shading interp; colorbar
    hold on
    contour(x2d/1e3,z2d/1e3,Tnew,[100:100:1500],'k');
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('Temperature [C]')
    title(['Temperature Diffusion in Aluminum Plate '])
    drawnow
end
end

csvwrite('CFDHomeworkExplicit.csv',T);
