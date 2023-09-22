import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from dataPlotter import dataPlotter
from blockbeamAnimation import blockbeamAnimation
from signalGenerator import signalGenerator

reference = signalGenerator(amplitude=0.5,frequency=0.1)
theta_ref = signalGenerator(amplitude=np.pi/8,frequency=0.1)
z_ref = signalGenerator(amplitude=0.05,frequency=0.5, y_offset=0.2)
f_ref = signalGenerator(amplitude=5,frequency=0.5)

dataPlot = dataPlotter()
animation = blockbeamAnimation()

t = P.t_start
while t < P.t_end:
    r = reference.square(t)
    theta = theta_ref.square(t)
    z = z_ref.sin(t)
    f = f_ref.sawtooth(t)

    state = np.array([[z],[theta],[0.0],[0.0]])
    animation.update(state)
    dataPlot.update(t,r,state,f)

    t = t + P.t_plot
    plt.pause(0.05)

print('Press key to close')
plt.waitforbuttonpress()
plt.close()
    