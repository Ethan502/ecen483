import matplotlib.pyplot as plt
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics

mass = massDynamics()
reference = signalGenerator(amplitude=0.01,frequency=0.1)
force = signalGenerator(amplitude=10, frequency=1)

dataPlot = dataPlotter()
animation = massAnimation()

t = P.t_start
while t < P.t_end:
    t_next = t + P.t_plot

    while t < t_next:
        r = reference.square(t)
        F = force.sin(t)

        y = mass.update(F)
        t = t + P.Ts
    
    animation.update(mass.state)
    dataPlot.update(t,r,mass.state,F)

    plt.pause(0.0001)

print('Press key to close')
plt.waitforbuttonpress()
plt.close()