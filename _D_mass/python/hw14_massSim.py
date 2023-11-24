import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics
from ctrlObserver import ctrlObserver
from dataPlotterObserver import dataPlotterObserver

mass = massDynamics(alpha=0.0)
controller = ctrlObserver()
reference = signalGenerator(amplitude=1,frequency=0.05)
disturbance = signalGenerator(amplitude=0.25)
noise = signalGenerator(amplitude=0.01)

dataPlot = dataPlotter()
animation = massAnimation()
dataPlotObserver = dataPlotterObserver()

t = P.t_start
y = mass.h()

while t < P.t_end:
    t_next_plot = t + P.t_plot

    while t < t_next_plot:
        r = reference.square(t)
        d = disturbance.step(t)
        n = noise.random(t)
        u, x_hat = controller.update(r,y + n)
        y = mass.update(u[0] + d)
        t = t + P.Ts

    # update animation and data plots
    animation.update(mass.state)
    dataPlot.update(t,r,mass.state,u)
    dataPlotObserver.update(t,mass.state,x_hat)

    # pause the figure
    plt.pause(0.001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()


