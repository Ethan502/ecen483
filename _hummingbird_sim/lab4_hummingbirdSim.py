import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import HummingbirdDynamics
from ctrlEquilibrium import CtrlEquilibrium

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
control = CtrlEquilibrium()

t = P.t_start  # time starts at t_start
while t < P.t_end:
    t_next = t + P.t_plot

    while t < t_next:
        u = control.update()
        y = hummingbird.update(u)
        t = t + P.Ts

    ref = np.array([[0], [0], [0]])
    force = (u[0][0] + u[1][0])*P.km
    torque = control.tau
    state = hummingbird.state
    animation.update(t,state)
    dataPlot.update(t,state,ref,force,torque)
    plt.pause(0.1)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
