import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import HummingbirdDynamics
from ctrlLonPID import ctrlLonPID

# instantiate pendulum, controller, and reference classes
hummingbird = HummingbirdDynamics(alpha=0.1)
controller = ctrlLonPID()
theta_ref = SignalGenerator(amplitude=0.5, frequency=0.1)

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()

t = P.t_start  # time starts at t_start
y = hummingbird.h()
while t < P.t_end:  # main simulation loop

    # Propagate dynamics at rate Ts
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        r = np.array([[theta_ref.square(t)], [0.]])
        u, y_ref = controller.update(r, y)
        y = hummingbird.update(u)  # Propagate the dynamics
        t = t + P.Ts  # advance time by Ts

    # update animation and data plots at rate t_plot
    state = hummingbird.state
    force = (u[0][0] + u[1][0])*P.km
    torque = (u[0][0] - u[1][0])*P.d
    animation.update(t,state)
    dataPlot.update(t,state,y_ref,force,torque)

    # the pause causes figure to be displayed during simulation
    plt.pause(0.01)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
