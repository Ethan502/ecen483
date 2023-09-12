# mass-spring-damper Parameter File
import numpy as np

# Physical parameters of the arm known to the controller
m = 0.5  # mass kg
k = 11  # spring constant Kg/s^2
b = 0.45 # damping coefficient Kg/s

# parameters for animation
length = 5.0
width = 1.0

# Initial Conditions
z0 = 0.25  # initial position of mass, m
zdot0 = 2  # initial velocity of mass m/s

# Simulation Parameters
t_start = 0.0 # Start time of simulation
t_end =  30.0 # End time of simulation
Ts = 0.01  # sample time for simulation
t_plot = 0.1 # the plotting and animation is updated at this rate

# dirty derivative parameters
# sigma =  # cutoff freq for dirty derivative
# beta =   # dirty derivative gain

# saturation limits
# F_max =   # Max force, N

