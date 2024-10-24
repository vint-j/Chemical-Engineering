% Clear all variables and command window
clear;
clc;

% Initial conditions for height and time
height_array(1) = 5;  % Initial height of the fluid (m)
time_array(1) = 0.0;  % Initial time (s)
time_step_dt = 0.1;   % Time step for Euler method (s)

% Euler method loop to solve ODE numerically
for i = 2:10000
    % Update time at each step
    time_array(i) = (i-1) * time_step_dt;
    
    % Call ODE function to calculate rate of change of height (dh/dt)
    dhdt = my_odeweek3d(time_array(i-1), height_array(i-1));
    
    % Update the height using Euler's method
    height_array(i) = height_array(i-1) + time_step_dt * dhdt;
    
    % Stop the loop if the height becomes negative
    if height_array(i) < 0
        break;
    end
end

% Plot the Euler method solution for height over time
figure;
plot(time_array, height_array)
ylim([0 5])  % Set y-axis limits to keep height between 0 and 5
xlim([0 time_array(end)])  % Set x-axis limit up to the last time step
xlabel("Time (s)")
ylabel("Height (m)")
title('Height of Fluid Over Time (Euler Method)')

% Constants for the analytical solution
rho_fluid = 1000;  % Density of the fluid (kg/m^3)
g = 9.8;  % Gravitational acceleration (m/s^2)
initial_height_ho = 5;  % Initial height of the fluid (m)
Dp_pipe_diameter = 0.1;  % Diameter of the pipe (m)
Dt_tank_diameter = 1;  % Diameter of the tank (m)
dynamic_viscosity_mu = 0.001;  % Dynamic viscosity of the fluid (Pa·s)

% Constants related to flow and proportionality
k0_constant = 1;  % Proportionality constant
k1_constant = (4/pi) * k0_constant;
k2_constant = (4/pi) * (1/k1_constant)^0.5;
k_proportionality = k2_constant / 2;

% Constant for the analytical solution equation
a_constant = (Dp_pipe_diameter^2 / Dt_tank_diameter^2) * k_proportionality * g^0.5;

% Time array for the analytical solution (linspace gives 100 points between 0 and 120 seconds)
time_array_analytical = linspace(0, 120, 100);

% Analytical solution for height over time
height_analytical = (initial_height_ho^0.5 - a_constant .* time_array_analytical).^2;

% Plot the analytical solution for height over time
figure;
plot(time_array_analytical, height_analytical)
ylim([0 5])  % Set y-axis limits to keep height between 0 and 5
xlim([0 time_array_analytical(end)])  % Set x-axis limit up to the last time step
xlabel("Time (s)")
ylabel("Height (m)")
title('Height of Fluid Over Time (Analytical Solution)')

% ===========================
% Function: my_odeweek3d
% ===========================

% This function calculates the rate of change of height (dh/dt) of fluid 
% draining from a tank through a pipe, accounting for gravitational effects.
function dhdt = my_odeweek3d(~, height_current)

    % Constants inside the ODE function
    rho_fluid = 1000;  % Density of the fluid (kg/m^3)
    g = 9.8;  % Gravitational acceleration (m/s^2)
    Dp_pipe_diameter = 0.1;  % Diameter of the pipe (m)
    Dt_tank_diameter = 1;  % Diameter of the tank (m)
    dynamic_viscosity_mu = 0.001;  % Dynamic viscosity of the fluid (Pa·s)

    % Constants related to flow and proportionality
    k0_constant = 1;  % Proportionality constant
    k1_constant = (4/pi) * k0_constant;
    k_proportionality = (4/pi) * (1/k1_constant)^0.5;

    % Calculate the rate of change of height (dh/dt)
    dhdt = -1 * (Dp_pipe_diameter^2 / Dt_tank_diameter^2) * (g * height_current)^0.5 * k_proportionality;
end