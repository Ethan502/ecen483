import matplotlib.pyplot as plt
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockbeamDynamics import blockbeamDynamics
from ctrlStateFeedback import ctrlStateFeedback

blockbeam = blockbeamDynamics(alpha=0.2)
reference = signalGenerator(amplitude=0.15, frequency=0.05, y_offset=0.25)
disturbance = signalGenerator(amplitude=0.1)

dataPlot = dataPlotter()
animation = blockbeamAnimation()
controller = ctrlStateFeedback()

t = P.t_start
y = blockbeam.h()
while t < P.t_end:
    t_next = t + P.t_plot
    while t < t_next:
        r = reference.square(t)
        d = disturbance.square(t)
        u = controller.update(r,blockbeam.state)
        y = blockbeam.update(u)
        t = t + P.Ts
    # Update the animation and dataplot
    animation.update(blockbeam.state)
    dataPlot.update(t,r,blockbeam.state,u)
    plt.pause(0.01)

print('Press key to close')
plt.waitforbuttonpress()
plt.close()