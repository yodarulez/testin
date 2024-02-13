# 1D HEAT CONDUCTION code

import math
import matplotlib.pyplot as plt
import numpy as np
import csv

# Data

# physical data
coeff_a = 10**-5 # thermal diffusivity a = k/ρCv, k: thermal conductivity, Cv specific heat capacity @ constant volume
total_length = 1 # length of the rod
total_time = 3000 # total calculation time
T_ic = 0 # Initial condition
T_bc = 100 # Boundary condition

# discretization data
grid_points = 11 # number of discrete points in the rod
coeff_s = 0.5 # s <= 0.5 for stability in FTCS scheme
no_of_terms = 10 # number of terms to be considered in the analytical solution

# calculated data
delta_x = total_length / (grid_points - 1) # cell size
delta_t = round((coeff_s * delta_x ** 2) / coeff_a, 2) # time-step
time_steps = int(total_time / delta_t) # total number of time steps

# Matrices initialization/declaration
T_calc = [[0 for i in range(grid_points)] for n in range(time_steps + 1)]
T_anal = [[0 for i in range(grid_points)] for n in range(time_steps + 1)]  
x = [0 for i in range(grid_points)]
t = [0 for n in range(time_steps + 1)]


# Initial & boundary conditions
for n in range (time_steps + 1):
    for i in range (grid_points):
        if n == 0: # initial conditions
            T_calc[n][i] = T_ic
        if i == 0 or i == grid_points - 1: # boundary conditions
            if n == 0 and T_ic != T_bc: # average between boundary and initial values for the corner points
                T_calc[n][i] = 0.5 * (T_ic + T_bc)
            else:
                T_calc[n][i] = T_bc 
        
# Numerical scheme
def FTCS (): 
    for n in range (time_steps):
        for i in range (1, grid_points - 1):
            T_calc[n + 1][i] = coeff_s * T_calc[n][i - 1] + (1 - 2 * coeff_s) * T_calc[n][i] + coeff_s * T_calc[n][i + 1]
                                    
FTCS()

# Analytical solution
for n in range(time_steps + 1):
    t[n] = delta_t * n # t time elapsed since initial state
    for i in range(grid_points):
        x[i] = delta_x * i # x distance from the left boundary
        sum1 = 0
        for m in range (1, no_of_terms):
            term = (2 * m - 1) * math.pi
            sum1 = sum1 + (400 / term) * math.sin(term * x[i]) * math.exp(-coeff_a * term * term * t[n]) 
            anal_solution = T_bc - sum1
            T_anal[n][i] = round(anal_solution, 3)
            
# RMS error / convergence criterion

for n in range(time_steps + 1):
    sum2 = 0
    for i in range(grid_points):
        sum2 = sum2 + (T_calc[n][i] - T_anal[n][i])**2
        rms = math.sqrt(sum2/grid_points)
    print (rms)

# Plots
# Calculated temperature
for n in range(time_steps + 1):
    plt.plot(x, T_calc[n], label='T_calculated' if n == 0 else None, color='blue', linestyle='dashed', marker='x',markerfacecolor='blue', markersize=12)
plt.xlabel('Points')
plt.ylabel('Temperature')
plt.title('1D Heat Conduction')
plt.legend()
plt.show()

#fuckah
#cocksuckah


# # Analytical temperature
# for n in range(time_steps):
#     plt.plot(x, T_anal[n], label='T_analytical', color='red', linestyle='dashed', marker='o',markerfacecolor='red', markersize=12)
# plt.xlabel('Points')
# plt.ylabel('Temperature')
# plt.title('1D Heat Conduction')
# plt.legend()
# plt.show()

# # CSV file
# file = open(r'C:\Users\Nick\cfdcode\results', 'w')
# writer = csv.writer(file)
# writer.writerow(T_calc[time_steps])
# file.close()

    
    
