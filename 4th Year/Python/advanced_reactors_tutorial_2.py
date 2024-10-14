import math
import matplotlib.pyplot as plt
import numpy as np 

conversion_target = 0.75
step = 0.001
t2_celcius = np.arange(0, 100 + step, step)  
t2 = [temp + 273.15 for temp in t2_celcius]  
t1 = 298.15
r = 8.3145
dh0 = -75300

def get_ke2(t2_value):
    ke2 = 300 * (math.exp((-dh0 / r) * ((t1 - t2_value) / (t1 * t2_value))))
    return ke2

ke2 = [get_ke2(temp) for temp in t2]
xae = [x / (1 + x) for x in ke2]

result = [t for t, x in zip(t2_celcius, xae) if x > conversion_target] # return t (temperature) for t and x in the zip function of the temps and the conversion, if x > conversion target

if result:
    print(f'The maximum temperature that gives a conversion of {conversion_target * 100:.2f}% is {result[-1]:.2f} °C')
else:
    print(f'No temperature could be found which gave the target conversion of {conversion_target * 100:.2f}%')

# UNCOMMENT THE CODE BELOW TO SEE THE PLOT

plt.figure(figsize=(6,6),facecolor='none')
plt.plot(t2_celcius, xae, color='red')
plt.xlabel(r'Temperature ($\degree$C)')
plt.ylabel('Conversion')
plt.minorticks_on()
plt.grid(which='both')
plt.gca().set_facecolor('none')
plt.show()
