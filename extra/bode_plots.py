import control as cnt
import numpy as np
import matplotlib.pyplot as plt

p = [1,0.1,0.6]

g = cnt.tf([1],[1,1,0.0])
print(g)
mag,phase,omega = cnt.bode(g)
plt.show()