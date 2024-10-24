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

% Linear Case
% Define the analytical solution function for c(t)
analytical_c = @(k, t) ((0.3 * N) / (V * k)) * (1 - exp(-k * t));

% Define the function f(k) = c(t_final) - c_target for the linear case
f_linear = @(k) analytical_c(k, t_final) - c_target;

% Provide an initial guess for k
k_initial_guess_linear = 1e-4;  % Initial guess for k (1/min)

% Use fzero to find the value of k for the linear case
[k_solution_linear, ~, exitflag_linear] = fzero(f_linear, k_initial_guess_linear);

% Check if the root-finding was successful for the linear case
if exitflag_linear <= 0
    error('Root finding for the linear case did not converge. Try a different initial guess.');
else
    fprintf('The rate constant k that keeps the astronauts alive for 4 days is %.6e L/min in the linear case.\n', k_solution_linear);
end

% Time vector for plotting
t = 0:dt:t_final;

% Compute the analytical solution using the found k (linear case)
c_analytical_linear = analytical_c(k_solution_linear, t);

% Initialize arrays for the numerical solution (linear case)
num_steps = length(t);
c_numerical_linear = zeros(1, num_steps);
c_numerical_linear(1) = c_initial;

% Euler's Method loop to compute the numerical solution (linear case)
for i = 2:num_steps
    dcdt = (0.3 * N) / V - k_solution_linear * c_numerical_linear(i - 1);
    c_numerical_linear(i) = c_numerical_linear(i - 1) + dt * dcdt;
end

% Plot both the numerical and analytical solutions (linear case)
figure;
plot(t, c_numerical_linear, 'r-', 'DisplayName', 'Numerical Solution (Euler)');
hold on;
plot(t, c_analytical_linear, 'b--', 'DisplayName', 'Analytical Solution');
xlabel('Time (minutes)');
ylabel('CO₂ Concentration (fraction)');
legend;
title('CO₂ Concentration Over Time (Linear Case)');
grid on;

% Second-Order Nonlinear Case
% Fixed second-order analytical solution function
analytical_c_second_order = @(k, t) sqrt((0.3 * N) / (k * V)) * ...
    ((exp(2 * t .* k * sqrt((0.3 * N) / (k * V))) - 1) ./ ...
    (exp(2 * t .* k * sqrt((0.3 * N) / (k * V))) + 1));

% Function to compute c(t_final) numerically for the second-order nonlinear case
function c_end = compute_c_end_nonlinear(k, N, V, t_final, dt, c_initial)
    % Initialize variables for the nonlinear case
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
    c_end = c;  % Return the final concentration
end

% Define the function f(k) = c(t_final) - c_target for the second-order case
f_second_order = @(k) compute_c_end_nonlinear(k, N, V, t_final, dt, c_initial) - c_target;

% Provide an initial guess for k (second-order case)
k_initial_guess_second_order = 1e-2;  % Initial guess for k (1/min)

% Use fzero to find the value of k for the second-order case
[k_solution_second_order, fval, exitflag_second_order] = fzero(f_second_order, k_initial_guess_second_order);

% Check if the root-finding was successful for the second-order case
if exitflag_second_order <= 0
    error('Root finding for the second-order case did not converge. Try a different initial guess.');
else
    fprintf('The rate constant k that keeps the astronauts alive for 4 days is %.6e L/min in the second-order nonlinear case.\n', k_solution_second_order);
end

% Compute the analytical solution using the found k (second-order case)
c_analytical_second_order = analytical_c_second_order(k_solution_second_order, t);

% Initialize array for the numerical solution (second-order case)
c_numerical_second_order = zeros(1, num_steps);
c_numerical_second_order(1) = c_initial;

% Euler's Method loop to compute the numerical solution (second-order case)
for i = 2:num_steps
    dcdt = (0.3 * N) / V - k_solution_second_order * c_numerical_second_order(i - 1)^2;
    c_numerical_second_order(i) = c_numerical_second_order(i - 1) + dt * dcdt;
    % Ensure concentration remains non-negative
    if c_numerical_second_order(i) < 0
        c_numerical_second_order(i) = 0;
    end
end

% Plot the numerical and analytical solutions (second-order case)
figure;
plot(t, c_numerical_second_order, 'r-', 'DisplayName', 'Numerical Solution (Euler)');
hold on;
plot(t, c_analytical_second_order, 'b--', 'DisplayName', 'Analytical Solution');
xlabel('Time (minutes)');
ylabel('CO₂ Concentration (fraction)');
legend;
title('CO₂ Concentration Over Time (Second-Order Nonlinear Case)');
grid on;