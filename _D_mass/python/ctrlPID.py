import numpy as np
import massParam as P

class ctrlPID:
    def __init__(self):
        tr = 0.2
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

        # Init variables for the PID control
        self.ki = 0.1
        self.z_delayed = P.z0
        self.zdot_delayed = P.zdot0
        self.error_sum = 0.0
        self.error_sum_prev = 0.0
        self.error_z_prev = 0.0

    def update(self, z_r, z):
        zdot = P.beta * self.zdot_delayed + (1-P.beta)/P.Ts * (z - self.z_delayed)
        error_z = (z_r - z)
        self.error_sum = self.error_sum + P.Ts/2 * (error_z + self.error_z_prev) 


        force = self.kp * (z_r - z) - self.kd * zdot + self.ki * self.error_sum
        force = saturate(force,P.F_max)
        return force

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u     
