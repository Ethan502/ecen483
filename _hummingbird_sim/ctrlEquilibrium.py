import numpy as np
import hummingbirdParam as P

class CtrlEquilibrium:
    def __init__(self):
        self.tau = 0
        self.F = ((P.m1*P.ell1 + P.m2*P.ell2)*P.g)/P.ellT

    def update(self):
        u_l = (1/(2*P.km))*(self.F + self.tau/P.d)
        u_r = (1/(2*P.km))*(self.F - self.tau/P.d)
        return np.array([[u_l], [u_r]])