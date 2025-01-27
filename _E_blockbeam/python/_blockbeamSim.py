import matplotlib.pyplot as plt
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockbeamDynamics import blockbeamDynamics
from testPID import ctrlPID

blockbeam = blockbeamDynamics()
reference = signalGenerator(amplitude=0.15, frequency=0.1, y_offset=0.25)
force = signalGenerator(amplitude=0.5,frequency=1,y_offset=11.5)

dataPlot = dataPlotter()
animation = blockbeamAnimation()
controller = ctrlPID()

t = P.t_start
while t < P.t_end:
    t_next = t + P.t_plot
    while t < t_next:
        r = reference.square(t)
        z = blockbeam.state[0][0]
        theta = blockbeam.state[1][0]
        u = controller.update(r,z,theta)
        y = blockbeam.update(u)
        t = t + P.Ts

    animation.update(blockbeam.state)
    dataPlot.update(t,r,blockbeam.state,u)
    plt.pause(0.01)

print('Press key to close')
plt.waitforbuttonpress()
plt.close()