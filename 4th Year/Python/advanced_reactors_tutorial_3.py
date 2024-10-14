import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import matplotlib.patheffects as path_effects

# constants
R = 8.314 # ideal gas constant
C_0 = 1 # initial concentration mol / L

# the constants used (17.34 and 42.04 as well as the activation energies) come from page 218 of chemical reaction engineering 3rd edition by octave levenspiel
# the tutorial question is based off of this example (9.2) and these values are valid for a reversible exothermic reaction that is first order and of the form A <=> R

# based on the graph that is produced when using the rates given in the question, the best thing to do for getting a decent conversion within a reasonable time would be to start off at a higher temperature
# which would allow the conversion to reach a local maximum of the iso-rate line with a higher value, then cool towards the target local maximum of an iso-rate line with a lower rate meaning higher conversion
# this is basically jumping across each peak of the iso-rate lines to ensure that we are utilising the higher reaction rate until it gets to its maximum conversion

# if we dont do the cooling, then we will start to approach the equilibrium line and conversion and rate will become unfeasible
# the approach of starting off hot, then cooling is better than starting off cold and then heating on the iso-rate line with the desired conversion, because it means we can get the benefit of higher rates at the start
# when comparing the different ca0 value graphs, we see that by increasing ca0 we are able to achieve a higher conversion from each iso-rate line, it is more pronounced the increase for the higher rate lines


# rate constants as functions of temperature (in celsius)
def k1(T):
    return np.exp(17.34 - 48900 / (R * (273.15 + T)))

def k2(T):
    return np.exp(42.04 - 124200 / (R * (273.15 + T)))

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
def plot_conversion(rates, save_plot=False, rate_input=''):
    temperatures = np.linspace(-10, 110, 500) if rate_input == 'levenplot' else np.linspace(-10, 180, 500)
    
    # calculate equilibrium conversion
    X_eq = [equilibrium_conversion(T) for T in temperatures]
    
    # plotting
    plt.figure(figsize=(6, 5))
    
    # plot equilibrium conversion with the same thickness as other lines, in red
    plt.plot(temperatures, X_eq, linewidth=2, color='black', label='Equilibrium Conversion (Rate = 0)', zorder=5)
    
    # plot conversion for each rate
    colors = plt.cm.viridis(np.linspace(0, 1, len(rates))) # this is just to make the other lines look nice, it linearly picks a colour from a dark blue to yellow gradient based on how many iso-rate lines you want to plot
    for rate, color in zip(rates, colors):
        X_rate = [conversion_at_rate(T, rate) for T in temperatures]
        plt.plot(temperatures, X_rate, linewidth=2, color=color, zorder=4)
        # Annotate the rate value at the intersection with y = -x + 1 (diagonal line)
        diagonal_line = -temperatures / 110 + 1
        intersection_index = np.argmin(abs(np.array(X_rate) - diagonal_line))
        text = plt.text(
            temperatures[intersection_index], 
            X_rate[intersection_index], 
            f'{rate:.3f}', 
            fontsize=8, 
            color=color, 
            bbox=dict(facecolor='lightgrey', alpha=0.5, pad=1, boxstyle="round,pad=0.3"), 
            zorder=6
        )
        # Add black or white outline to text for better readability based on color brightness
        outline_color = 'white' if np.mean(color[:3]) < 0.35 else 'black'
        text.set_path_effects([path_effects.Stroke(linewidth=0.5, foreground=outline_color), path_effects.Normal()])
    
    plt.xlim(-10, 110) if rate_input == 'levenplot' else plt.xlim(-10, 180)
    plt.ylim(0, 1)
    plt.xlabel('Temperature (Celsius)')
    plt.ylabel('Conversion, X')
    plt.title('Conversion vs Temperature')
    plt.legend(title=f'$C_{{A0}} = {C_0}$ mol / L', loc='upper right', fontsize=6)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.minorticks_on()
    
    # Adjust layout
    plt.tight_layout()
    
    if save_plot:
        filename = 'conversion_vs_temperature_levenplot.png' if rate_input == 'levenplot' else 'conversion_vs_temperature.png'
        plt.savefig(filename, dpi=300)
        print(f"Plot saved as '{filename}'")
    
    plt.show()

# main function to run the script
def main():
    os.system('cls' if os.name == 'nt' else 'clear')
    try:
        rate_input = input("Please enter the rate values separated by commas, or enter 'levenplot' to see the example from levenspiels book: ")
        if rate_input.lower() == 'levenplot':
            rates = [0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5]
        else:
            rates = [float(rate.strip()) for rate in rate_input.split(',')]
    except ValueError:
        print("All rates must be numbers.")
        return
    
    plot_conversion(rates, rate_input='levenplot' if rate_input.lower() == 'levenplot' else '')
    
    save_input = input("Do you want to save the plot? (Y/N): ").strip().lower() # strip removes any leading or trailing whitespace from user input, lower converts to lowercase, meaning the selection is not case sensitive, any value passed other than "y" or "Y" will result in the plot not saving
    if save_input == 'y':
        plot_conversion(rates, save_plot=True, rate_input='levenplot' if rate_input.lower() == 'levenplot' else '')
    else:
        print("Plot not saved.")

if __name__ == "__main__": # this just checks to see if this script is being run directly as it should be, if the module name matches then it is being ran directly and will then run the main loop
    main()