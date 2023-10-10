import matplotlib.pyplot as plt
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockbeamDynamics import blockbeamDynamics
from ctrlPD import ctrlPD

blockbeam = blockbeamDynamics()
reference = signalGenerator(amplitude=0.15, frequency=0.1, y_offset=0.25)
force = signalGenerator(amplitude=0.5,frequency=1,y_offset=11.5)

dataPlot = dataPlotter()
animation = blockbeamAnimation()
controller = ctrlPD()

t = P.t_start
while t < P.t_end:
    t_next = t + P.t_plot
    while t < t_next:
        r = reference.square(t)
        state = blockbeam.state
        u = controller.update(r,state)
        y = blockbeam.update(u)
        t = t + P.Ts

    animation.update(blockbeam.state)
    dataPlot.update(t,r,blockbeam.state,u)
    plt.pause(0.0001)

print('Press key to close')
plt.waitforbuttonpress()
plt.close()