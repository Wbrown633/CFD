% Solves the 2D heat equation with an explicit finite difference scheme
% Code modified from template provided by: 
% http://geodynamics.usc.edu/~becker/teaching/557/problem_sets/problem_set_
% fd_2dheat.pdf
% Many thanks! 

clear
% NOTE: This iteration was designed to be used with equal horizontal and 
% vertical spacing and size. dx = dy and L = H 
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
converged = 1e-6; % Convergence criteria
itteration_difference = 1; % Initialize check variable to run while loop
itterations = 0; % Counter to track how many times we've run
nt = 20000; % Maximum iterations to be run [shouldn't be reached]  
dx = L/(nx-1); % Spacing of grid in x-direction
dy = H/(ny-1); % Spacing of grid in z-direction
[x2d,z2d] = meshgrid(-L/2:dx:L/2, -H:dy:0); % create grid
% Compute stable timestep
dt = min([dx,dy])^2/alpha/4;
% Setup initial linear temperature profile
T = zeros(ny,nx);

T(:,:) = Tinit;
time = 0;
while itteration_difference > converged && itterations <= nt;

% Create New matrix 
Tnew = zeros(ny,nx);
Tnew(:,:) = Tinit;

% Set Wall boundary conditions
Tnew(1,:) = Tbot;
Tnew(:,1) = Twall;

% Compute new temperature

    for i=2:ny-1
        for j=2:nx-1
            Tnew(i,j) = (T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1))/4;
        end
    end
    
% Set Insulation Boundary Conditions
for i=2:ny-1
    Tnew(i,nx) = Tnew(i,nx-1);
end

for j=2:nx-1
    Tnew(ny,j) = Tnew(ny-1,j);
end
% Subtract new matrix from the previous itteration, the absolute value of 
% the max value of the resulting 

differenceMatrix = Tnew - T;
abs(differenceMatrix);
itteration_difference = max(max(differenceMatrix));

% We need two 'max' calls because the first returns the largest value in 
% each column, the max of this array gives us the largest overall value



T = Tnew;
time = time+dt;
itterations = itterations + 1; 

%Plot solution every 50 timesteps
if (mod(itterations,50)==0)
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

csvwrite('CFDHomeworkIdeal.csv',T);