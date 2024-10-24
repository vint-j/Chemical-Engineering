% Clear variables and command window
clear;
clc;

% Constants and initial conditions
sphere_diameter_m = 0.03;                     % Diameter of the sphere (m)
sphere_volume_m3 = pi * (sphere_diameter_m^3) / 6;  % Volume of the sphere (m^3)
gravity_mps2 = 9.8;                           % Gravitational acceleration (m/s^2)
fluid_density_kgpm3 = 1261;                   % Density of the fluid (kg/m^3)
sphere_density_kgpm3 = 2000;                  % Density of the sphere (kg/m^3)
sphere_mass_kg = sphere_volume_m3 * sphere_density_kgpm3;  % Mass of the sphere (kg)
spring_constant = 1;                          % Spring constant (dimensionless)
dynamic_viscosity_pa_s = 1.4;                 % Dynamic viscosity of the fluid (Pa·s)

% Time array for the analytical solution (linspace gives 100 points between 0 and 15 seconds)
time_array_s = linspace(0, 15, 100);

% Calculate constants for the velocity equation
A_constant = sphere_volume_m3 * gravity_mps2 * (sphere_density_kgpm3 - fluid_density_kgpm3) / sphere_mass_kg;
B_constant = spring_constant * sphere_diameter_m * dynamic_viscosity_pa_s / sphere_mass_kg;

% Analytical velocity solution over time
analytical_velocity_mps = (A_constant / B_constant) * (1 - exp(-B_constant * time_array_s));

% Plot the analytical velocity solution
figure;
plot(time_array_s, analytical_velocity_mps, 'b', 'DisplayName', 'Analytical Solution')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Analytical and Euler Method Velocity of the Sphere')
hold on;  % Hold the plot to add the numerical solution later

% Euler Method Variables
velocity_array_euler(1) = 0.0;  % Initial velocity (m/s)
time_array_euler(1) = 0.0;      % Initial time (s)
time_step_dt = 0.01;            % Time step for Euler method (s)

% Numerical solution using Euler's method
for i = 2:100 / time_step_dt
    time_array_euler(i) = (i - 1) * time_step_dt;  % Update time at each step
    dvdt = my_odeweek2(time_array_euler(i - 1), velocity_array_euler(i - 1));  % Call ODE function for velocity change
    velocity_array_euler(i) = velocity_array_euler(i - 1) + time_step_dt * dvdt;  % Euler update for velocity
end

% Plot the Euler method velocity solution
plot(time_array_euler, velocity_array_euler, 'r', 'DisplayName', 'Euler Method')
xlim([0 5])  % Set x-axis limit for better visualization
legend('show')

% Characteristic time calculation
characteristic_time_s = pi / 6 * sphere_diameter_m^2 * (sphere_density_kgpm3 - fluid_density_kgpm3) / dynamic_viscosity_pa_s;
fprintf('Characteristic time (s): %.4f\n', characteristic_time_s);

% Inverse of the characteristic time
inverse_characteristic_time = 1 / characteristic_time_s;
fprintf('Inverse of the characteristic time: %.4f\n', inverse_characteristic_time);

% ===========================
% Function: my_odeweek2
% ===========================

% This function calculates the rate of change of velocity (dv/dt) of the sphere
% as it falls through a fluid under gravity, with drag effects.
function dvdt = my_odeweek2(~, velocity_current)

    % Define the constants
    sphere_diameter_m = 0.03;  % Diameter of the sphere (m)
    sphere_volume_m3 = pi * (sphere_diameter_m^3) / 6;  % Volume of the sphere (m^3)
    gravity_mps2 = 9.8;  % Gravitational acceleration (m/s^2)
    fluid_density_kgpm3 = 1261;  % Fluid density (kg/m^3)
    sphere_density_kgpm3 = 2000;  % Sphere density (kg/m^3)
    sphere_mass_kg = sphere_volume_m3 * sphere_density_kgpm3;  % Sphere mass (kg)
    spring_constant = 1;  % Spring constant or proportionality constant
    dynamic_viscosity_pa_s = 1.4;  % Dynamic viscosity of the fluid (Pa·s)

    % Compute constants for velocity rate equation
    A_constant = sphere_volume_m3 * gravity_mps2 * (sphere_density_kgpm3 - fluid_density_kgpm3) / sphere_mass_kg;  % Driving force
    B_constant = spring_constant * sphere_diameter_m * dynamic_viscosity_pa_s / sphere_mass_kg;  % Drag effect

    % Calculate the rate of change of velocity (dv/dt)
    dvdt = A_constant - (B_constant * velocity_current);  % Rate of velocity change due to gravity and drag
end