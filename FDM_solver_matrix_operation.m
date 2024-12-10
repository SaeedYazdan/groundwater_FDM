% Parameters
Nx = 20;           % Number of grid points in x-direction
Ny = 20;           % Number of grid points in y-direction
h_top = 10;        % Top boundary condition

% Build the coefficient matrix
A = buildCoeffMatrix(Nx, Ny);

% Build the right-hand side vector
b = buildRHSVector(Nx, Ny, h_top);

% Solve the linear system and get the full solution with boundary conditions
h2D = solveLinearSystem(A, b, Nx, Ny, h_top);

% Plot the solution
plotSolution(h2D, Nx, Ny);

% Function Definitions

function A = buildCoeffMatrix(Nx, Ny)
    % Constructs the coefficient matrix A for the linear system Ax = b
    
    % Create the diagonal, upper, and lower components of D
    Ddiag = -4 * eye(Nx - 1);
    Dupper = diag(ones(1, Nx - 2), 1);
    Dlower = diag(ones(1, Nx - 2), -1);
    
    % Create D matrix
    D = Ddiag + Dupper + Dlower;
    
    % Create a block matrix where the diagonals are D
    Ds = repmat({D}, 1, Nx - 1);
    A = blkdiag(Ds{:});
    
    % Create identity diagonals
    I = ones((Nx - 1) * (Nx - 2), 1);
    Iupper = diag(I, Nx - 1);
    Ilower = diag(I, -(Nx - 1));
    
    % Add identity diagonals to A
    A = A + Iupper + Ilower;
end

function b = buildRHSVector(Nx, Ny, h_top)
    % Constructs the right-hand side vector b for Ax = b
    b = zeros((Nx - 1)^2, 1);
    b(end - Nx + 2:end) = -h_top;
end

function h2D = solveLinearSystem(A, b, Nx, Ny, h_top)
    % Solves the system Ax = b and adds boundary conditions
    
    % Solve for h vector and reshape to 2D
    h = A \ b;
    h = reshape(h, Nx - 1, Ny - 1);
    
    % Create a 2D array including boundary conditions
    h2D = zeros(Nx + 1, Ny + 1);
    h2D(1, :) = h_top;               % Top boundary condition
    h2D(2:end-1, 2:end-1) = flipud(h);  % Solution in interior
end

function plotSolution(h2D, Nx, Ny)
    % Plots the solution on a 2D meshgrid
    
    % Create mesh grid
    x = linspace(0, 1, Nx + 1);
    y = linspace(1, 0, Ny + 1);
    [X, Y] = meshgrid(x, y);
    
    % Plot solution
    figure;
    contourf(X, Y, h2D, 10);
    colorbar;
    xlabel('x [m]');
    ylabel('y [m]');
    title('2D Steady State Flow Solution');
end
