import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter

#instantiate the reference input classes
reference = signalGenerator(amplitude=0.5,frequency=0.1)
distance = signalGenerator(amplitude=1,frequency=0.5,y_offset=0.2)
tauRef = signalGenerator(amplitude=5, frequency=.5)

dataPlot = dataPlotter()
animation = massAnimation()

t = P.t_start
while t < P.t_start:
    #set the variables
    r = reference.square(t)
    d = distance.sin(t)
    tau = tauRef.sawtooth(t)
    # update the animation
    state = np.array([[d],[0.0]])
    animation.update(state)
    dataPlot.update(t,r,state,)
