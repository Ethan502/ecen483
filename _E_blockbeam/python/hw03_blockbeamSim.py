import matplotlib.pyplot as plt
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockbeamDynamics import blockbeamDynamics

blockbeam = blockbeamDynamics()
reference = signalGenerator(amplitude=0.5, frequency=0.02)
force = signalGenerator(amplitude=0.5, frequency=1)

dataPlot = dataPlotter()
animation = blockbeamAnimation()

t = P.t_start
while t < P.t_end:
    t_next = t + P.t_plot
    while t < t_next:
        r = reference.square(t)
        F = force.sin(t)
        y = blockbeam.update(F)
        t = t + P.Ts

    animation.update(blockbeam.state)
    dataPlot.update(t,r,blockbeam.state,F)
    plt.pause(0.0001)

print('Press key to close')
plt.waitforbuttonpress()
plt.close()