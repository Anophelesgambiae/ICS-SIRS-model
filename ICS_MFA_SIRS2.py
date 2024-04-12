import random
import numpy as np
import matplotlib.pyplot as plt

# Parameters:
beta = 2.0  # mean infectivity per infected
p1 = 1/5  # probability that a infected from pathogen 1 becomes recovered
p2 = 1/5.2  # probability that a infected from pathogen 2 becomes recovered
omega = 1/6 # lost of immunity rate
S0 = 0.95  # initial ratio of susceptible population
I0 = 0.05  # initial ratio of infected population
h = 10 # The scaling of infectivity (in the Hill function it is the half_max_value)  
timesteps = 500 # number of timesteps to run
N_population = 10000 # Total population

S_data = np.zeros(timesteps+1)
I1_data = np.zeros(timesteps+1)
I2_data = np.zeros(timesteps+1)
R_data = np.zeros(timesteps+1)
time_data = np.arange(timesteps+1) 

S_data[0] = S0
I1_data[0] = I0/2
I2_data[0] = I0/2

'''
Return a value on the Hill function. The Hill function brings a saturaded relation between x and the result.
As example the result can be the enzym activity and x be the concentration of the substrate.
The property's of the Hill function:
- if x = h then result = 0.5
- if x converge to endless then result converge to max_value 
- the Hillcoëficient determines how steep the function becomes near x = h.
'''
def Hill_function(x, max_value = 1, half_max_value = 0.5, Hillcoëficient = 1):

    result = max_value * pow(x, Hillcoëficient) / (pow(x,Hillcoëficient) + pow(half_max_value, Hillcoëficient))
    return result


for t in range(0,timesteps):

    S_data[t+1] = S_data[t] - Hill_function(8*beta*(I1_data[t]+I2_data[t]), half_max_value = h)*S_data[t] + omega*R_data[t]

    I1_data[t+1] = I1_data[t] + (8*beta*I1_data[t])/(8*beta*(I1_data[t]+I2_data[t]) + h)*S_data[t]
    I2_data[t+1] = I2_data[t] + (8*beta*I2_data[t])/(8*beta*(I1_data[t]+I2_data[t]) + h)*S_data[t]

    # For low I_1 in the calculation it is possible to have negative numbers. 
    # That's why we set for negative numbers the population to zero.
    if I1_data[t+1] < 0:
        I1_data[t+1] = 0

    R_data[t+1] = R_data[t] - omega*R_data[t]
    
    # Because we work with a probability that a infected becomes recovered, we need now absolute values 
    # of the population and not ratios. That's why we define a total population N, and for each infected person we calculate of
    # he/she becomes recovered. If so, we add a ratio on R_data and remove a ratio on I_data.
    for i in range(0, round(N_population*I1_data[t])):
        if random.random() < p1:
            I1_data[t+1] -= 1/N_population
            R_data[t+1] += 1/N_population
    
    for j in range(0, round(N_population*I2_data[t])):
        if random.random() < p2:
            I2_data[t+1] -= 1/N_population
            R_data[t+1] += 1/N_population 
               

# Plot all four populations
plt.plot(time_data, S_data, label = 'S', color = 'gray', linewidth = 2.5, alpha = 0.8)
plt.plot(time_data, I1_data, label = 'I\u2081', color = 'lightcoral', linewidth = 2.5, alpha = 0.8)
plt.plot(time_data, I2_data, label = 'I\u2082', color = 'darkred', linewidth = 2.5, alpha = 0.8)
plt.plot(time_data, R_data, label = 'R', color = 'green', linewidth = 2.5, alpha = 0.8)
plt.xlabel('time')
plt.ylabel('population ratio')
plt.legend(markerscale = 1.5)
plt.show()