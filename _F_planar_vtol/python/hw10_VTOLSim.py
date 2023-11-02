import matplotlib.pyplot as plt
import VTOLParam as P
import numpy as np
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import VTOLDynamics
from ctrlPID import ctrlPID


vtol = VTOLDynamics()
z_reference = signalGenerator(amplitude=2.5, frequency=0.08, y_offset=3)
h_reference = signalGenerator(amplitude=3, frequency=0.03,y_offset=5)
dataPlot = dataPlotter()
animation = VTOLAnimation()
controller = ctrlPID()

t = P.t_start
while t < P.t_end:
    t_next = t + P.t_plot
    while t < t_next:
        zref = z_reference.square(t)
        href = h_reference.square(t)
        n = np.array([[0.0],[0.0],[0.0]])
        d = np.array([[0.0],[0.0]])
        z = vtol.state[0][0]
        h = vtol.state[1][0]
        theta = vtol.state[2][0]
        u = controller.update(href,zref,z,h,theta)
        y = vtol.update(u + d)
        t = t + P.Ts
    
    animation.update(vtol.state)
    dataPlot.update(t,vtol.state,zref,href,u[0][0],u[1][0])
    plt.pause(0.01)

print('Press key to close')
plt.waitforbuttonpress()
plt.close()
