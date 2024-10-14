import os
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# constants
R = 8.314
C_0 = 1

# the constants used (17.34 and 42.04 as well as the activation energies) come from page 218 of chemical reaction engineering 3rd edition by octave levenspiel

# rate constants as functions of temperature (in celsius)
def k1(T):
    return np.exp(17.34 - 48900 / (R * (273.15 + T))) / 60  # Convert from per minute to per second

def k2(T):
    return np.exp(42.04 - 124200 / (R * (273.15 + T))) / 60  # Convert from per minute to per second

# rate equation
def rate_equation(T, X, rate):
    return k1(T) * C_0 * (1 - X) - k2(T) * C_0 * X - rate

# function to calculate equilibrium conversion X for a given temperature
def equilibrium_conversion(T):
    return fsolve(lambda X: rate_equation(T, X, 0), 0.5)[0]

# function to calculate conversion X for a given temperature and rate
def conversion_at_rate(T, rate):
    return fsolve(lambda X: rate_equation(T, X, rate), 0.5)[0]

# plotting function
def plot_conversion(rates, save_plot=False):
    temperatures = np.linspace(-10, 300, 500)
    
    # calculate equilibrium conversion
    X_eq = [equilibrium_conversion(T) for T in temperatures]
    
    # plotting
    plt.figure(figsize=(10, 6))
    
    # plot equilibrium conversion with thicker line
    plt.plot(temperatures, X_eq, linewidth=3, color='orange', label='Equilibrium Conversion (Rate = 0)', zorder=5)
    
    # plot conversion for each rate
    colors = plt.cm.viridis(np.linspace(0, 1, len(rates)))
    for rate, color in zip(rates, colors):
        X_rate = [conversion_at_rate(T, rate) for T in temperatures]
        plt.plot(temperatures, X_rate, linewidth=2, color=color, label=f'Conversion (Rate = {rate})', zorder=4)
    
    plt.xlim(-10, 300)
    plt.ylim(0, 1)
    plt.xlabel('Temperature (Celsius)')
    plt.ylabel('Conversion, X')
    plt.title('Conversion vs Temperature')
    plt.legend(title=f'$C_{{A0}} = {C_0}$ mol / L')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.minorticks_on()
    
    if save_plot:
        plt.savefig('conversion_vs_temperature.png')
        print("Plot saved as 'conversion_vs_temperature.png'")
    
    plt.show()

# main function to run the script
def main():
    os.system('cls' if os.name == 'nt' else 'clear')
    time.sleep(1)  # delay to ensure terminal is cleared before prompting for input
    try:
        rate_input = input("Please enter the rate values separated by commas: ")
        rates = [float(rate.strip()) for rate in rate_input.split(',')]
    except ValueError:
        print("All rates must be numbers.")
        return
    
    plot_conversion(rates)
    
    time.sleep(1)  # delay before prompting to save the plot
    save_input = input("Do you want to save the plot? (Y/N): ").strip().lower()
    if save_input == 'y':
        plot_conversion(rates, save_plot=True)
    else:
        print("Plot not saved.")

if __name__ == "__main__":
    main()