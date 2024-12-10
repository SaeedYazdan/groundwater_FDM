% Parameters
Lx = 1000;   % Length of the domain in the x-direction (m)
Ly = 1000;   % Length of the domain in the y-direction (m)
Nx = 50;     % Number of grid points in x-direction
Ny = 50;     % Number of grid points in y-direction
dx = Lx / (Nx - 1);  % Grid spacing in x-direction (m)
dy = Ly / (Ny - 1);  % Grid spacing in y-direction (m)
dt = 1;      % Time step (s)
T = 1000;    % Total simulation time (s)
Kx = 1e-4;   % Hydraulic conductivity in x-direction (m/s)
Ky = 1e-4;   % Hydraulic conductivity in y-direction (m/s)
Ss = 1e-5;   % Specific storage (1/m)
Q_well = -0.1; % Pumping rate of the well (m^3/s) (negative means extraction)

% Well location (arbitrary grid point)
well_x = round(Nx / 2);  % Well location in x-direction
well_y = round(Ny / 3);  % Well location in y-direction

% Initialize topography (elevation)
[X, Y] = meshgrid(linspace(0, Lx, Nx), linspace(0, Ly, Ny));

% Define a simple topography (hill in the center, valley at the edges)
topography = 10 * exp(-((X - Lx/2).^2 + (Y - Ly/2).^2) / (2 * (Lx/4)^2)) ...
             - 5 * (X / Lx);  % Gaussian hill in the center + slope from left to right

% Initialize hydraulic head (h) array based on topography
h = topography + 100;  % Assume initial head is 100 + topography

% Boundary conditions: fixed head at boundaries
h(:, 1) = 95;  % Left boundary (Dirichlet condition)
h(:, end) = 90; % Right boundary (Dirichlet condition)
h(1, :) = 95;   % Top boundary (Dirichlet condition)
h(end, :) = 95; % Bottom boundary (Dirichlet condition)

% Create a copy of the hydraulic head to store the updated values
h_new = h;

% Calculate the pumping rate at the well location in terms of head
source_term = Q_well * dt / (Ss * dx * dy);  % Source/sink term due to well

% Finite difference method: iterate over time steps
for t = 0:dt:T
    % Loop over interior points (excluding boundaries)
    for i = 2:Nx-1
        for j = 2:Ny-1
            % Calculate derivatives using central differences
            dhdx2 = (h(j, i+1) - 2*h(j, i) + h(j, i-1)) / dx^2;
            dhdy2 = (h(j+1, i) - 2*h(j, i) + h(j-1, i)) / dy^2;
            
            % Update hydraulic head using the finite difference formula
            h_new(j, i) = h(j, i) + dt * (Kx * dhdx2 + Ky * dhdy2) / Ss;
        end
    end
    
    % Apply the well as a source or sink at the specified location
    h_new(well_y, well_x) = h_new(well_y, well_x) + source_term;
    
    % Update the hydraulic head for the next time step
    h = h_new;
end

% Plot the final hydraulic head distribution
figure;
contourf(X, Y, h, 20);  % 2D contour plot of hydraulic head
colorbar;
title('Hydraulic Head Distribution with Topography and Well');
xlabel('X (m)');
ylabel('Y (m)');

% Plot the topography
figure;
surf(X, Y, topography); % 3D plot of the topography
title('Topography (Surface Elevation)');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Elevation (m)');
