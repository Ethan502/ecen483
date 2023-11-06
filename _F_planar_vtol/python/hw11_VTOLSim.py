import matplotlib.pyplot as plt
import VTOLParam as P
import numpy as np
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import VTOLDynamics
from ctrlStateFeedback import ctrlStateFeedback


vtol = VTOLDynamics(alpha=0.0)
z_reference = signalGenerator(amplitude=2.5, frequency=0.08, y_offset=3)
h_reference = signalGenerator(amplitude=3, frequency=0.03,y_offset=5)
dataPlot = dataPlotter()
animation = VTOLAnimation()
controller = ctrlStateFeedback()

t = P.t_start
y = vtol.h()
while t < P.t_end:
    t_next = t + P.t_plot
    while t < t_next:
        zref = z_reference.square(t)
        href = h_reference.square(t)
        n = np.array([[0.0],[0.0],[0.0]])
        d = np.array([[0.0],[0.0]])
        u = controller.update(href,zref,vtol.state)
        y = vtol.update(u + d)
        t = t + P.Ts
    
    animation.update(vtol.state)
    dataPlot.update(t,vtol.state,zref,href,u[0][0],u[1][0])
    plt.pause(0.01)

print('Press key to close')
plt.waitforbuttonpress()
plt.close()
