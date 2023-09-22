import matplotlib.pyplot as plt
import VTOLParam as P
import numpy as np
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import VTOLDynamics

vtol = VTOLDynamics()
reference = signalGenerator(amplitude=0.5, frequency=0.1)
force_1 = signalGenerator(amplitude=10,frequency=1)
force_2 = signalGenerator(amplitude=10,frequency=1)

dataPlot = dataPlotter()
animation = VTOLAnimation()

t = P.t_start
while t < P.t_end:
    t_next = t + P.t_plot

    while t < t_next:
        r = reference.square(t)
        u_1 = force_1.sin(t)
        u_2 = force_2.sin(t)
        u = np.array([[u_1],[u_2]])
        y = vtol.update(u)
        t = t + P.Ts
    
    animation.update(vtol.state)
    dataPlot.update(t,r,vtol.state,u)
    plt.pause(0.0001)

print('Press key to close')
plt.waitforbuttonpress()
plt.close()
