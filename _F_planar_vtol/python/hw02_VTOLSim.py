import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from signalGenerator import signalGenerator

z_reference = signalGenerator(amplitude=0.5, frequency=0.1)
thetaRef = signalGenerator(amplitude=np.pi/8, frequency=0.5)
zRef = signalGenerator(amplitude=4,frequency=0.1,y_offset=5)
alt_ref = signalGenerator(amplitude=2,frequency=0.1,y_offset=2)
ctl_ref = signalGenerator(amplitude=5, frequency=.5)

# Added this line
h_reference = signalGenerator(amplitude=0.5,frequency=0.1,y_offset=1.0)



dataPlot = dataPlotter()
animation = VTOLAnimation()

t = P.t_start
while t < P.t_end:

    # added this line
    h_r = h_reference.square(t)

    z_ref = z_reference.sin(t) #changed this one to a sin wave
    theta = thetaRef.sin(t)
    z = zRef.sin(t)
    alt = alt_ref.square(t)
    f = ctl_ref.sawtooth(t) #just changed this variable name
    tau = ctl_ref.sawtooth(t) # added this line
    motor_thrusts = P.mixing * np.array([[f],[tau]]) # added this line

    state = np.array([[z],[alt],[theta],[0.0],[0.0],[0.0]]) # added some zeros
    animation.update(x=state)


    dataPlot.update(t=t,states=state,z_ref=z_ref,h_ref=h_r,torque=3,force=P.fmax)
    
    # this line is what is in the solution, but the solution said that the values for
    # my above commented out line can be 0, so it is still correct
    #dataPlot.update(t,state,z_ref,h_r,motor_thrusts,torque=3)

    t = t + P.t_plot
    plt.pause(0.001)

print('Press key to close')
plt.waitforbuttonpress()
plt.close()