% MATLAB Code for Fully Coupled Hydrogeophysical Inversion of Synthetic Salt Tracer Experiments

% Clear workspace and command window
clear; clc;

% Define simulation parameters
Lx = 2.8; % Length of the domain (m)
Ly = 0.05; % Width of the domain (m)
Lz = 0.6; % Height of the domain (m)

nx = 280; % Number of cells in x direction
nz = 60; % Number of cells in z direction
Dx = Lx / nx; % Grid spacing in x direction (m)
Dz = Lz / nz; % Grid spacing in z direction (m)

% Hydraulic and transport parameters
hydraulic_gradient = 0.01; % Hydraulic gradient
c_in = 0.4; % Concentration in inflow (g/L)
porosity = 0.4; % Porosity
alpha_l = 0.01; % Longitudinal dispersivity (m)
alpha_t = 0.001; % Transverse dispersivity (m)
D_m = 1e-9; % Pore diffusion coefficient (m^2/s)

% Geoelectrical parameters
I = 2; % Current (mA)
s0 = 0.03; % Base electrical conductivity (S/m)
lambda = 0.06; % Linear dependence of electrical conductivity on concentration (S*m^2/kg)

% Time parameters
t_inj = 2 * 3600; % Duration of injection (s)
dt = 300; % Time discretization in transient calculation (s)

% Geostatistical values of log hydraulic conductivity field
s_Y2 = 2; % Variance
l_x = 0.4; % Correlation length in x direction (m)
l_y = 0.05; % Correlation length in y direction (m)
b_Y = log(1e-7); % Prior mean of log hydraulic conductivity (m/s)
R_bb = 10; % Uncertainty of prior mean

% Measurement errors
s_h = 1e-3; % Hydraulic head measurement error (m)
s_R = 0.01; % Electrical resistance measurement error (1%)

% Initialize hydraulic conductivity field
K_true = exp(b_Y + randn(nx, nz) * sqrt(s_Y2)); % True hydraulic conductivity field

% Simulate tracer injection and measure electrical potential
% Initialize arrays for measurements
D_phi = zeros(nx, nz); % Electrical potential difference
concentration = zeros(nx, nz); % Concentration field

% Simulate tracer injection
for t = 1:t_inj/dt
    % Update concentration field based on hydraulic conductivity
    concentration = update_concentration(concentration, K_true, hydraulic_gradient, c_in, dt, porosity, alpha_l, alpha_t);
    
    % Calculate electrical potential difference based on concentration
    D_phi = calculate_electrical_potential(concentration, s0, lambda);
end

% Add measurement noise
D_phi_noisy = D_phi + randn(size(D_phi)) * s_R;

% Compute temporal moments of potential perturbation
[mean_arrival_time, zeroth_moment] = compute_temporal_moments(D_phi_noisy, dt);

% Inversion process to estimate hydraulic conductivity
% Initialize estimated hydraulic conductivity field
K_estimated = ones(nx, nz) * exp(b_Y); % Initial guess

% Perform inversion using geostatistical method
for iter = 1:20 % Number of iterations
    % Compute sensitivity matrix
    H = compute_sensitivity_matrix(K_estimated, mean_arrival_time, zeroth_moment);
    
    % Update estimated hydraulic conductivity
    K_estimated = update_hydraulic_conductivity(K_estimated, H, mean_arrival_time, zeroth_moment);
    
    % Check for convergence (optional)
    if check_convergence(K_estimated, K_true)
        break;
    end
end

% Plot results
figure;
subplot(1, 2, 1);
imagesc(K_true); colorbar; title('True Hydraulic Conductivity');
xlabel('X Direction'); ylabel('Z Direction');

subplot(1, 2, 2);
imagesc(K_estimated); colorbar; title('Estimated Hydraulic Conductivity');
xlabel('X Direction'); ylabel('Z Direction');

% Functions used in the simulation
function concentration = update_concentration(concentration, K, hydraulic_gradient, c_in, dt, porosity, alpha_l, alpha_t)
    % Update the concentration field based on the hydraulic conductivity
    % and the advection-dispersion equation
    % This is a simplified representation of the transport process
    % using a finite difference approach.
    
    % Calculate the flux based on Darcy's law
    q = -K .* hydraulic_gradient; % Specific discharge
    % Update concentration using advection-dispersion equation
    concentration = concentration + (q * c_in * dt) / porosity; % Simplified update
end

function D_phi = calculate_electrical_potential(concentration, s0, lambda)
    % Calculate the electrical potential difference based on the concentration
    % and the relationship between concentration and electrical conductivity
    s_prime = s0 + lambda * concentration; % Perturbed conductivity
    D_phi = -log(s_prime); % Simplified potential calculation
end

function [mean_arrival_time, zeroth_moment] = compute_temporal_moments(D_phi, dt)
    % Compute the zeroth and first temporal moments of the electrical potential
    zeroth_moment = sum(D_phi(:)); % Total potential
    mean_arrival_time = sum((1:length(D_phi(:)))' * D_phi(:)) / zeroth_moment; % Mean arrival time
end

function H = compute_sensitivity_matrix(K, mean_arrival_time, zeroth_moment)
    % Compute the sensitivity matrix for the inversion process
    % This function calculates how sensitive the measurements are to changes in K
    H = zeros(length(K(:)), 1); % Initialize sensitivity matrix
    for i = 1:length(K(:))
        % Perturb K slightly and recalculate mean arrival time
        K_perturbed = K;
        K_perturbed(i) = K_perturbed(i) * 1.01; % 1% perturbation
        new_mean_arrival_time = compute_temporal_moments(calculate_electrical_potential(update_concentration(zeros(size(K)), K_perturbed, 0.01, 0.4, 300, 0.4, 0.01, 0.001), 0.03, 0.06), 300);
        H(i) = (new_mean_arrival_time - mean_arrival_time) / (0.01 * K(i)); % Sensitivity calculation
    end
end

function K_updated = update_hydraulic_conductivity(K, H, mean_arrival_time, zeroth_moment)
    % Update the hydraulic conductivity based on the inversion process
    % This function uses a simple gradient descent approach
    alpha = 0.01; % Learning rate
    K_updated = K - alpha * H; % Update rule
end

function converged = check_convergence(K_estimated, K_true)
    % Check for convergence of the estimated hydraulic conductivity
    tolerance = 1e-3; % Convergence tolerance
    error = norm(K_estimated - K_true); % Calculate error
    converged = error < tolerance; % Check if within tolerance
end