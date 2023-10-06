import numpy as np
import massParam as P

class ctrlPD:
    def __init__(self):
        tr = 2
        zeta = 0.7
        wn = 2.2/tr

        alpha0 = wn**2
        alpha1 = 2 * zeta * wn

        a1 = P.b / P.m
        a0 = P.k / P.m
        b0 = 1 / P.m

        self.kd = (alpha1 - 0.1)/0.2
        self.kp = (alpha0 - 0.6)/0.2

