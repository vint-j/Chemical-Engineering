import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.integrate import quad


os.system('cls' if os.name == 'nt' else 'clear') # clear the terminal

print('TASK 1')
print()
# data arrays
particle_diameter = np.array([10, 20, 40, 80, 160, 320])
volume_percentage = np.array([0, 0.01, 0.09, 0.9, 9, 90])

# calculate cumulative and differential frequencies
cumulative_frequency_Qv = np.cumsum(volume_percentage / 100)
differential_frequency_qv = np.diff(cumulative_frequency_Qv, prepend=0) / np.diff(particle_diameter, prepend=particle_diameter[0])

# calculate volume of a single particle
diameter_of_sphere = particle_diameter
volume_of_sphere = (4/3) * np.pi * (diameter_of_sphere / 2) ** 3  # volume of a sphere in cubic microns

# calculate number fraction based on volume percentage (without total volume)
N = (volume_percentage / 100) / volume_of_sphere  # number of particles per volume percentage
total_N = np.sum(N)
print(f'total number of particles = {total_N}')
number_fraction = N / total_N

# calculate number cumulative frequency and differential frequency
number_cumulative_frequency_Qn = np.cumsum(number_fraction)
number_differential_frequency_qn = np.diff(number_cumulative_frequency_Qn, prepend=0) / np.diff(particle_diameter, prepend=particle_diameter[0])

# calculate weighted mean diameters
number_weighted_mean_diameter = np.sum(number_fraction * particle_diameter) / np.sum(number_fraction)
volume_weighted_mean_diameter = np.sum(number_fraction * particle_diameter**4) / np.sum(number_fraction * particle_diameter**3)
print(np.sum(number_fraction * particle_diameter**4), np.sum(number_fraction * particle_diameter**3))

# estimate volume-weighted median diameter with linear interpolation
median_index = np.where(cumulative_frequency_Qv >= 0.5)[0][0] # basically here we are asking where is the cumulative frequency Qv > 0.5 (50%) and returning the index it happens at
if median_index > 0: # if this index isnt the first value then we proceed
    x1, x2 = particle_diameter[median_index - 1], particle_diameter[median_index] # we set x1 and x2 to the particle diameter at index of the median -1 and the median respectively
    y1, y2 = cumulative_frequency_Qv[median_index - 1], cumulative_frequency_Qv[median_index] # similarly we assign y1 and y2 to the cumulative frequency at these indexes
    volume_weighted_median_diameter = x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1) # then we use the linear interpolation formula x1 + (0.5 - y1) * gradient
else:
    volume_weighted_median_diameter = particle_diameter[median_index] # this is just a backup incase the median index is 0, it returns the first value of the particle diameter

# create pandas dataframe
data = {
    'Grid Points': np.arange(len(particle_diameter)),
    'Particle Diameter': particle_diameter,
    'Volume Percentage (%)': volume_percentage,
    'Cumulative Frequency Qv': cumulative_frequency_Qv,
    'Differential Frequency qv (1 / micron)': differential_frequency_qv,
    'Volume of one particle (micron^3)': volume_of_sphere,
    'Number of particles': N,
    'Number fraction': number_fraction,
    'Cumulative number fraction': number_cumulative_frequency_Qn,
    'Differential number frequency qn (1 / micron)': number_differential_frequency_qn
}
df = pd.DataFrame(data)

print(df)
print("Number Weighted Mean Diameter:", number_weighted_mean_diameter)
print("Volume Weighted Mean Diameter:", volume_weighted_mean_diameter)
print("Volume Weighted Median Diameter (Interpolated):", volume_weighted_median_diameter)

print()
print('TASK 2')
print()

# calculate volume equivalent diameters for needle-shaped particles
# given: square cross-section of 5 micron x 5 micron, length varying between 20 and 200 microns
cross_section_side = 5  # microns
lengths = np.array([20, 200])  # microns

# calculate volume of needle-shaped particles
volume_needle = cross_section_side**2 * lengths  # volume of needle = cross-sectional area * length

# calculate equivalent sphere diameters for the given lengths
volume_equivalent_diameter = (6 * volume_needle / np.pi) ** (1/3)  # rearrange volume of sphere formula to solve for diameter

print("Volume Equivalent Diameters for Lengths of 20 and 200 microns:", volume_equivalent_diameter)

# calculate cumulative frequency distribution Qn(x) for the particle length and tabulate values for lengths: 10, 100, 1000 microns
particle_lengths = np.array([10, 100, 1000])  # microns

# calculate cumulative frequency Qn for given lengths
Qn_values = []
for length in particle_lengths:
    if length <= 20:
        Qn_values.append(0)
    elif 20 < length <= 200:
        Qn_values.append(0.00556 * (length - 20)) # we subtract the 20 to basically get the region we are interested in since it scales linearly starting at 20, so say it was 100, we would do 100-20 which would be 80 which is the effective range we want to look at 
    else:
        Qn_values.append(0.00556 * (200 - 20))  # beyond 200 microns, Qn remains constant

# normalize Qn to ensure it is between 0 and 1
Qn_values = np.array(Qn_values) / Qn_values[-1]

# tabulate Qn values
Qn_data = {
    'Particle Length (micron)': particle_lengths,
    'Cumulative Frequency Qn': Qn_values
}
Qn_df = pd.DataFrame(Qn_data)
print(Qn_df)


# # UNCOMMENT THIS TO GET THE GRAPH
# # sketch the graph of Qn(x)

# plt.figure(figsize=(4, 4))
# plt.plot(particle_lengths, Qn_values, marker='o', linestyle='-', color='b')
# for i, txt in enumerate(Qn_values):
#     plt.text(particle_lengths[i], Qn_values[i], f'{txt:.2f}', fontsize=9, ha='right', va='bottom')
# plt.xlabel('Particle Length (micron)')
# plt.ylabel('Cumulative Frequency Qn')
# plt.title('Cumulative Frequency Distribution Qn(x)')
# plt.minorticks_on()
# plt.grid(True, which='both', linestyle='--', linewidth=0.5)
# plt.show()

# Define the limits of integration
lower_limit = 20  # microns
upper_limit = 200  # microns

# Define the number fraction distribution as a function of particle length
def number_fraction_length_function(x):
    return 0.00556  # Assuming a uniform distribution similar to the original example

# Calculate the number weighted average length using integration
numerator_number_weighted, _ = quad(lambda x: number_fraction_length_function(x) * x, lower_limit, upper_limit)
denominator_number_weighted, _ = quad(number_fraction_length_function, lower_limit, upper_limit)
number_weighted_average_length = numerator_number_weighted / denominator_number_weighted

# Calculate the length weighted average length using integration
numerator_length_weighted, _ = quad(lambda x: number_fraction_length_function(x) * x**2, lower_limit, upper_limit)
denominator_length_weighted, _ = quad(lambda x: number_fraction_length_function(x) * x, lower_limit, upper_limit)
length_weighted_average_length = numerator_length_weighted / denominator_length_weighted

print("Number Weighted Average Length:", number_weighted_average_length)
print("Length Weighted Average Length:", length_weighted_average_length)

print()
print('TASK 3')
print()


from mpmath import mp

# Set precision for more accurate numerical calculations
mp.dps = 50  # Increase the decimal precision

# (a) Calculate the differential frequency distribution qn(x)
# Define the cumulative frequency distribution function Qn(x)
def Qn(x):
    if x <= 200:
        return 0
    elif 200 < x <= 600:
        return 0.0025 * (x - 200)
    else:
        return 1

# Define the differential frequency distribution qn(x) as the derivative of Qn(x)
def qn(x):
    if 200 < x <= 600:
        return 0.0025
    else:
        return 0

# Generate data for plotting qn(x)
x_values = np.linspace(100, 700, 500)
qn_values = [qn(x) for x in x_values]

# Plot qn(x)
# plt.plot(x_values, qn_values, label="qn(x)")
# plt.xlabel("Diameter (x) [μm]")
# plt.ylabel("Differential Frequency Distribution qn(x)")
# plt.title("Differential Frequency Distribution for Particle Diameter")
# plt.grid()
# plt.legend()
# plt.show()

# (b) Calculate the number weighted average diameter and the volume weighted average diameter
# Number weighted average diameter
numerator_number_weighted, _ = quad(lambda x: qn(x) * x, 200, 600, epsabs=1e-10, epsrel=1e-10)
denominator_number_weighted, _ = quad(qn, 200, 600, epsabs=1e-10, epsrel=1e-10)
number_weighted_average_diameter = numerator_number_weighted / denominator_number_weighted

# Volume weighted average diameter
numerator_volume_weighted, _ = quad(lambda x: qn(x) * (x ** 3), 200, 600, epsabs=1e-10, epsrel=1e-10)
denominator_volume_weighted, _ = quad(lambda x: qn(x) * (x ** 2), 200, 600, epsabs=1e-10, epsrel=1e-10)
volume_weighted_average_diameter = numerator_volume_weighted / denominator_volume_weighted

print("Number Weighted Average Diameter:", number_weighted_average_diameter, "μm")
print("Volume Weighted Average Diameter:", volume_weighted_average_diameter, "μm")

# (c) calculate the total volume and total surface area per 1 kg of particles
# given data
density = 2000  # kg/m^3
thickness = 5  # microns
thickness_m = thickness * 1e-6  # Convert to meters

# total volume per 1 kg of particles
volume_per_particle = lambda x: (mp.pi * (x * 1e-6 / 2) ** 2) * thickness_m  # volume of a cylinder (m^3)
total_volume, _ = quad(lambda x: qn(x) * volume_per_particle(x), 200, 600, epsabs=1e-10, epsrel=1e-10)
number_of_particles_per_kg = 1 / (total_volume * density)

print("Total Volume per 1 kg of Particles:", total_volume * number_of_particles_per_kg, "m^3")

# total surface area per 1 kg of particles
surface_area_per_particle = lambda x: (2 * mp.pi * (x * 1e-6 / 2) ** 2) + (2 * mp.pi * (x * 1e-6 / 2) * thickness_m)  # surface area of a cylinder (m^2)
total_surface_area, _ = quad(lambda x: qn(x) * surface_area_per_particle(x), 200, 600, epsabs=1e-10, epsrel=1e-10)

print("Total Surface Area per 1 kg of Particles:", total_surface_area * number_of_particles_per_kg, "m^2")


print()
print('TASK 4')
print()

