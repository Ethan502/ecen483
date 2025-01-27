import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import HummingbirdDynamics

# instantiate reference input classes
phi_ref = SignalGenerator(amplitude=1.5, frequency=0.05)
theta_ref = SignalGenerator(amplitude=0.5, frequency=0.05)
psi_ref = SignalGenerator(amplitude=0.5, frequency=.05)
force_1 = SignalGenerator(amplitude=10,frequency=1)
force_2 = SignalGenerator(amplitude=10,frequency=1)

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()
hummingbird = HummingbirdDynamics()

t = P.t_start  # time starts at t_start
while t < P.t_end:
    t_next = t + P.t_plot

    while t < t_next:
        phi = phi_ref.sin(t)
        theta = 0#theta_ref.sin(t)
        psi = 0#psi_ref.sin(t)
        u_1 = force_1.sin(t)
        u_2 = force_2.sin(t)
        u = np.array([[u_1],[u_2]])
        state = hummingbird.state
        y = hummingbird.update(u)
        t = t + P.Ts
        ref = np.array([[0], [0], [0]])
        force = 0
        torque = 0
    
    animation.update(t,state)
    dataPlot.update(t,state,ref,force,torque)
    plt.pause(0.1)

    animation.update(t, state)
    dataPlot.update(t, state, ref, force, torque)

    t = t + P.t_plot  # advance time by t_plot
    plt.pause(0.5)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
