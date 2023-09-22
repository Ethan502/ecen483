import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from signalGenerator import signalGenerator

reference = signalGenerator(amplitude=0.5, frequency=0.1)
thetaRef = signalGenerator(amplitude=np.pi/8, frequency=0.5)
zRef = signalGenerator(amplitude=4,frequency=0.1,y_offset=5)
alt_ref = signalGenerator(amplitude=2,frequency=0.1,y_offset=2)
ctl_ref = signalGenerator(amplitude=5, frequency=.5)


dataPlot = dataPlotter()
animation = VTOLAnimation()

t = P.t_start
while t < P.t_end:
    r = reference.square(t)
    theta = thetaRef.sin(t)
    z = zRef.sin(t)
    alt = alt_ref.square(t)
    ctl = ctl_ref.sawtooth(t)

    state = np.array([[z],[alt],[theta],[0.0]])
    animation.update(x=state)
    dataPlot.update(t=t,states=state,z_ref=P.z0,h_ref=P.h0,torque=3,force=P.fmax)

    t = t + P.t_plot
    plt.pause(0.001)

print('Press key to close')
plt.waitforbuttonpress()
plt.close()