% Parameters
nx = 50;            % Number of grid points in x-direction
ny = 50;            % Number of grid points in y-direction
Lx = 1000;          % Length of the domain in x (m)
Ly = 1000;          % Length of the domain in y (m)
K = 1e-2;           % Hydraulic conductivity (m/s)

% Discretization
dx = Lx / (nx - 1); % Grid spacing in x
dy = Ly / (ny - 1); % Grid spacing in y
x = linspace(0, Lx, nx);
y = linspace(0, Ly, ny);
h = ones(ny, nx);   % Initialize head matrix

% Boundary Conditions
h(:, 1) = 10;       % Left boundary fixed head
h(:, end) = 10;      % Right boundary fixed head

% Well Locations
h(1, 17:28) = -10;

% Finite Difference Solution
for iter = 1:1000
    h_old = h;
    for i = 2:nx-1
        for j = 2:ny-1
                % Standard finite difference calculation
                h(j, i) = (h_old(j, i+1) + h_old(j, i-1) + h_old(j+1, i) + h_old(j-1, i)) / 4;
        end
    end
    if max(max(abs(h - h_old))) < 1e-5
        break;
    end
end

% Visualization
figure;
contourf(x, y, h, 50);
colorbar;
title('Groundwater Head Distribution');
xlabel('Distance (m)');
ylabel('Distance (m)');
hold on;
legend('Head');
