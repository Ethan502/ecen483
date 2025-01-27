import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics
from ctrlStateFeedbackIntegrator import ctrlStateFeedbackIntegrator

mass = massDynamics(alpha=0.2)
controller = ctrlStateFeedbackIntegrator()
reference = signalGenerator(amplitude=1,frequency=0.05)
disturbance = signalGenerator(amplitude=0.25)

dataPlot = dataPlotter()
animation = massAnimation()

t = P.t_start
y = mass.h()

while t < P.t_end:
    t_next_plot = t + P.t_plot

    while t < t_next_plot:
        r = reference.square(t)
        d = disturbance.step(t)
        n = 0.0
        u = controller.update(r,mass.state)
        y = mass.update(u[0] + d)
        t = t + P.Ts

    # update animation and data plots
    animation.update(mass.state)
    dataPlot.update(t,r,mass.state,u)

    # pause the figure
    plt.pause(0.01)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()


