import numpy as np
import massParam as P

class ctrlPD:
    def __init__(self):
        tr = 0.7
        zeta = 0.7
        wn = 2.2/tr

        alpha0 = wn**2
        alpha1 = 2 * zeta * wn

        a1 = P.b / P.m
        a0 = P.k / P.m
        b0 = 1 / P.m

        self.kd = (alpha1 - 0.1)/0.2
        self.kp = (alpha0 - 0.6)/0.2

        print('kp: ', self.kp)
        print('kd: ', self.kd)

    def update(self, z_r, state):
        z = state[0][0]
        zdot = state[1][0]

        force = self.kp * (z_r - z) - self.kd * zdot
        force = saturate(force,P.F_max)
        return force

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u     
