import numpy as np
import blockbeamParam as P

class ctrlPD:
    def __init__(self):
        tr_th = 0.1
        zeta = 0.707
        M = 10
        wn_th = 2.2/tr_th
        
        b0_th = P.length/((P.m2*(P.length**2))/3 + P.m1 * (P.ze**2))

        alpha1_th = 2 * zeta * wn_th 
        alpha0_th = wn_th**2

        self.kd_th = alpha1_th / b0_th
        self.kp_th = alpha0_th / b0_th

        print("Kp_th: ", self.kp_th)
        print("Kd_th: ", self.kd_th)

        tr_z = M * tr_th
        b0_z = P.g
        wn_z = 2.2/tr_z

        alpha1_z = 2 * zeta * wn_z
        alpha0_z = wn_z**2

        self.kp_z = -alpha0_z / b0_z
        self.kd_z = -alpha1_z / b0_z

        print("Kp_z: ", self.kp_z)
        print("Kd_z: ", self.kd_z)

        k_DC = (b0_th*self.kp_th)/(b0_th*self.kp_th)


    def update(self, z_r, state):
        z = state[0][0]
        theta = state[1][0]
        zdot = state[2][0]
        thetadot = state[3][0]

        theta_ref = self.kp_z * (z_r - z) - self.kd_z * zdot

        F = self.kp_th * (theta_ref - theta) - self.kd_th * thetadot
        Fe = P.m1*P.g*z / P.length + 0.5*P.m2*P.g
        return saturate(F + Fe,P.Fmax)

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u    
