% Clear all variables and close all figures
clc
clearvars
close all

% Parameters
N = 3;                  % Number of astronauts
V = 30e3;               % Volume in liters (30 m^3 = 30,000 liters)
t_final = 4 * 24 * 60;  % Total time in minutes (4 days)
dt = 0.01;              % Time step in minutes
c_initial = 0;          % Initial CO₂ concentration (assumed zero)
c_target = 0.04;        % Maximum allowable CO₂ concentration (4%)

% Define the analytical solution function for c(t)
analytical_c = @(k, t) ((0.3 * N) / (V * k)) * (1 - exp(-k * t));

% Define the function f(k) = c(t_final) - c_target
f = @(k) analytical_c(k, t_final) - c_target;

% Provide an initial guess for k
k_initial_guess = 1e-4;  % Initial guess for k (1/min)

% Use fzero to find the value of k that keeps c(t_final) = c_target
[k_solution, ~, exitflag] = fzero(f, k_initial_guess);

% Check if the root-finding was successful
if exitflag <= 0
    error('Root finding did not converge. Try a different initial guess.');
else
    fprintf('The rate constant k that keeps the astronauts alive for 4 days is %.6e L/min in the linear case.\n', k_solution);
end

% Time vector for plotting
t = 0:dt:t_final;

% Compute the analytical solution using the found k
c_analytical = analytical_c(k_solution, t);

% Initialize arrays for the numerical solution
num_steps = length(t);
c_numerical = zeros(1, num_steps);
c_numerical(1) = c_initial;

% Euler's Method loop to compute the numerical solution
for i = 2:num_steps
    dcdt = (0.3 * N) / V - k_solution * c_numerical(i - 1);
    c_numerical(i) = c_numerical(i - 1) + dt * dcdt;
end

% Plot both the numerical and analytical solutions
figure;
plot(t, c_numerical, 'r-', 'DisplayName', 'Numerical Solution (Euler)');
hold on;
plot(t, c_analytical, 'b--', 'DisplayName', 'Analytical Solution');
xlabel('Time (minutes)');
ylabel('CO₂ Concentration (fraction)');
legend;
title('CO₂ Concentration Over Time (Linear Case)');
grid on;

% Clear all variables and close all figures
clearvars


% Parameters
N = 3;                  % Number of astronauts
V = 30e3;               % Volume in liters (30 m^3 = 30,000 liters)
t_final = 4 * 24 * 60;  % Total time in minutes (4 days)
dt = 0.01;              % Time step in minutes
c_initial = 0;          % Initial CO₂ concentration (assumed zero)
c_target = 0.04;        % Maximum allowable CO₂ concentration (4%)

% Since an analytical solution may not be readily available for the nonlinear case,
% we'll use numerical methods to find k.

% Define the function to compute c(t_final) numerically for a given k
function c_end = compute_c_end_nonlinear(k, N, V, t_final, dt, c_initial)
    % Initialize variables
    c = c_initial;  % Initial concentration
    t = 0;
    while t < t_final
        dcdt = (0.3 * N) / V - k * c^2;
        c = c + dt * dcdt;
        t = t + dt;
        % Ensure concentration remains non-negative
        if c < 0
            c = 0;
        end
    end
    c_end = c;
end

% Define the function f(k) = c(t_final) - c_target
f = @(k) compute_c_end_nonlinear(k, N, V, t_final, dt, c_initial) - c_target;

% Provide an initial guess for k
k_initial_guess = 1e-2;  % Initial guess for k (1/min)

% Use fzero to find the value of k that keeps c(t_final) = c_target
[k_solution, fval, exitflag] = fzero(f, k_initial_guess);

% Check if the root-finding was successful
if exitflag <= 0
    error('Root finding did not converge. Try a different initial guess.');
else
    fprintf('The rate constant k that keeps the astronauts alive for 4 days is %.6e L/min in the non-linear case.\n', k_solution);
end

% Time vector for plotting
t = 0:dt:t_final;
num_steps = length(t);

% Initialize array for the numerical solution
c_numerical = zeros(1, num_steps);
c_numerical(1) = c_initial;

% Euler's Method loop to compute the numerical solution
for i = 2:num_steps
    dcdt = (0.3 * N) / V - k_solution * c_numerical(i - 1)^2;
    c_numerical(i) = c_numerical(i - 1) + dt * dcdt;
    % Ensure concentration remains non-negative
    if c_numerical(i) < 0
        c_numerical(i) = 0;
    end
end

% Plot the numerical solution
figure;
plot(t, c_numerical, 'r-', 'DisplayName', 'Numerical Solution (Euler)');
xlabel('Time (minutes)');
ylabel('CO₂ Concentration (fraction)');
legend;
title('CO₂ Concentration Over Time (Nonlinear Case)');
grid on;
