clc
clearvars
close all
format long

fprintf('Question 1-a)\n')
time = 0:1:4;
aperture_diameter_m = 1 / 100;
circumference_cm = [72.5, 65.5, 57.0, 44.5, 30.5];
circumference_m = circumference_cm ./ 100;
diameter_m = circumference_m ./ pi;
volume_of_balloon_m3 = (pi .* diameter_m.^3) ./ 6;
volumetric_flowrate_m3_per_s = -diff(volume_of_balloon_m3) ./ diff(time) % negative sign because it is deflating
aperture_area_m = pi .* aperture_diameter_m.^2 ./4;
average_volumetric_flowrate_m3_per_s = mean(volumetric_flowrate_m3_per_s);
velocity = volumetric_flowrate_m3_per_s ./ aperture_area_m;
average_velocity_mps = average_volumetric_flowrate_m3_per_s / aperture_area_m
fprintf('From the Beaufort scale, a gale is 17.2 - 20.7 m/s, as average velocity = %.2f m/s, this would be a gale, 8 on the scale.\n\n', average_velocity_mps)

fprintf('Question 1-b)\n\n')
density_air = 1.295; % kg / m3
viscosity_air = 1.81e-5; % kg / m*s
reynolds_number = (density_air * aperture_diameter_m * average_velocity_mps) / viscosity_air;
fprintf('Reynolds number = %.0f \n\n', reynolds_number) % rounding to nearest integer

fprintf('Question 1-f\n')

k = 1 / ((9*pi/16)^(1/3));

initial_volume = volume_of_balloon_m3(1)
vb = k .* (aperture_area_m).^0.5 .* average_velocity_mps ./ (initial_volume - volumetric_flowrate_m3_per_s.* time((1:end-1))).^(1/3)
plot(time(1:end-1),vb)

%for distance dsb/dt = vb
distance = trapz(time(1:end-1), vb)

fprintf('Question 1-g)\n')

% Constants
initial_radius_r0 = 0.5 * (0.725 / pi); % (m) Radius calculated from circumference
initial_volume_V0 = (4/3) * pi * (initial_radius_r0^3); % (m^3) Initial volume of the balloon
air_volume_flow_rate_Vfloair = 1.45 * 10^(-3); % (m^3/s) Volume flow rate of air
air_velocity = 18.49; % (m/s) Velocity of air
aperture_area_Aaperture = pi * (0.005^2); % (m^2) Area of aperture
air_density_rho = 1.225; % (kg/m^3) Air density
balloon_mass_Mballoon = 0.004; % (kg) Mass of the balloon
drag_coefficient_k = 0.5; % Drag coefficient
time_step_dt = 0.001; % (s) Time step for Euler method
velocity_threshold = 0.01; % Threshold for stopping velocity (m/s)

% Preallocate time, velocity, and acceleration arrays
time_array_t = 0:time_step_dt:6; % Time array from 0 to 6 seconds
balloon_velocity_v = zeros(size(time_array_t)); % Preallocate velocity of balloon array
balloon_acceleration_a = zeros(size(time_array_t)); % Preallocate acceleration of balloon array

% Initial conditions
balloon_velocity_v(1) = 0; % Initial velocity of balloon
balloon_acceleration_a(1) = 0; % Initial acceleration is zero (at t=0)

% Variable to store the time when the balloon stops
min_velocity_time = NaN; 
min_velocity = Inf;  % Initialize as very large value

% Euler method loop to calculate balloon velocity and acceleration over time
for i = 2:length(time_array_t)
    time_current = time_array_t(i-1);
    
    % Calculate current balloon cross-sectional area
    balloon_cross_section_area_Aballoon = pi * ((3 * (initial_volume_V0 - air_volume_flow_rate_Vfloair * time_current) / (4 * pi))^(1/3))^2;
    
    % Compute the rate of change of velocity (acceleration = dv/dt)
    dvdt = (air_density_rho * (air_velocity^2 * aperture_area_Aaperture * drag_coefficient_k - ...
        balloon_velocity_v(i-1)^2 * balloon_cross_section_area_Aballoon + ...
        balloon_velocity_v(i-1) * air_velocity * aperture_area_Aaperture)) ...
        / (balloon_mass_Mballoon + air_density_rho * (initial_volume_V0 - air_volume_flow_rate_Vfloair * time_current));
    
    % Store the acceleration at this time step
    balloon_acceleration_a(i) = dvdt;
    
    % Update the velocity using Euler's method
    balloon_velocity_v(i) = balloon_velocity_v(i-1) + time_step_dt * dvdt;
    
    % Track the minimum velocity and corresponding time between 4 and 5 seconds
    if time_current >= 4 && time_current <= 5
        if balloon_velocity_v(i) < min_velocity
            min_velocity = balloon_velocity_v(i);
            min_velocity_time = time_current;
        end
    end
end

% Now clamp velocities and accelerations after the minimum velocity is found
for i = 1:length(time_array_t)
    if time_array_t(i) > min_velocity_time
        balloon_velocity_v(i) = 0; % Clamp velocity to zero after the minimum point
        balloon_acceleration_a(i) = 0; % Clamp acceleration to zero
    end
end

% Plot velocity of the balloon over time
figure;
plot(time_array_t, balloon_velocity_v, "DisplayName", "Velocity of Balloon (Euler)")
grid on
legend show
xlabel("Time (s)")
ylabel("Velocity (m/s)")
title("Velocity Profile of Balloon (Euler Method)")

% Plot acceleration of the balloon over time
figure;
plot(time_array_t, balloon_acceleration_a, "DisplayName", "Acceleration of Balloon (Euler)")
grid on
legend show
xlabel("Time (s)")
ylabel("Acceleration (m/s^2)")
title("Acceleration Profile of Balloon (Euler Method)")

% Display the time when the minimum velocity is reached
if ~isnan(min_velocity_time)
    fprintf('The balloon reaches minimum velocity at approximately t = %.2f seconds\n', min_velocity_time);
else
    fprintf('The balloon did not reach minimum velocity within the simulated time range.\n');
end