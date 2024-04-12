import random
import numpy as np
import matplotlib.pyplot as plt

# Parameters:
beta = 1.5  # mean infectivity per infected
p = 1/5  # probability that a infected becomes recovered
omega = 1/6 # lost of immunity rate
S0 = 0.95  # initial ratio of susceptible population
I0 = 0.05  # intial ratio of infected population
h = 10 # The scaling of infectivity (in the Hill function it is the half_max_value)  
timesteps = 50 # number of timesteps to run
N_population = 10000 # Total population

S_data = np.zeros(timesteps+1)
I_data = np.zeros(timesteps+1)
R_data = np.zeros(timesteps+1)
time_data = np.arange(timesteps+1) 

S_data[0] = S0
I_data[0] = I0

'''
Return a value on the Hill function. The Hill function brings a saturaded relation between x and the result.
As example the result can be the enzym activity and x be the concentration of the substrate.
The property's of the Hill function:
- if x = h then result = 0.5
- if x converge to endless then result converge to max_value 
- the Hillcoëficient determines how steep the function becomes near x = h.
'''
def Hill_function(x, max_value = 1, half_max_value = 0.5, Hillcoëficient = 1):

    result = max_value * pow(x, Hillcoëficient) / (pow(x, Hillcoëficient) + pow(half_max_value, Hillcoëficient))
    return result

for t in range(0,timesteps):

    S_data[t+1] = S_data[t] - Hill_function(8*beta*I_data[t], half_max_value = h)*S_data[t] + omega*R_data[t]
    I_data[t+1] = I_data[t] + Hill_function(8*beta*I_data[t], half_max_value = h)*S_data[t]
    R_data[t+1] = R_data[t] - omega*R_data[t]
    
    # Because we work with a probability that a infected becomes recovered, we need now absolute values 
    # of the population and not ratios. That's why we define a total population N, and for each infected person we calculate of
    # he/she becomes recovered. If so, we add a ratio on R_data and remove a ratio on I_data.
    for i in range(0, round(N_population*I_data[t])):
        if random.random() < p:
            I_data[t+1] -= 1/N_population
            R_data[t+1] += 1/N_population         

# Plot all three populations
plt.plot(time_data, S_data, label = 'S', color = 'gray', linewidth = 2, alpha = 0.8)
plt.plot(time_data, I_data, label = 'I', color = 'red', linewidth = 2, alpha = 0.8)
plt.plot(time_data, R_data, label = 'R', color = 'green', linewidth = 2, alpha = 0.8)
plt.xlabel('time')
plt.ylabel('population ratio')
plt.legend(markerscale = 1.5)
plt.show()
   