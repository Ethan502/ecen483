import numpy as np
import VTOLParam as P

class ctrlPD:
    def __init__(self):
        tr_h = 3
        zeta = 0.707
        wn_h = 2.2/tr_h

        b0_h = 1/(P.mc + 2*P.mr)
        alpha1_h = 2 * zeta * wn_h
        alpha0_h = wn_h**2

        self.kd_h = alpha1_h / b0_h
        self.kp_h = alpha0_h / b0_h

        print("Kp_h: ", self.kp_h)
        print("Kd_h: ", self.kd_h)



        # Inner Loop ----------------------------

        tr_th = 0.1
        wn_th = 2.2 / tr_th
        b0_th = 1/(P.Jc + 2 * P.mr * (P.d**2))
        alpha1_th = 2 * zeta * wn_th
        alpha0_th = wn_th**2

        self.kp_th = alpha0_th / b0_th
        self.kd_th = alpha1_th / b0_th

        print("Kp_th: ", self.kp_th)
        print("Kd_th: ", self.kd_th)

        DC_th = 1

        # ----------------------------------------

        # Outer Loop -----------------------------
        M = 10
        tr_z = M * tr_th
        wn_z = 2.2 / tr_z
        b0_z = -P.g
        alpha1_z = 2 * zeta * wn_z
        alpha0_z = wn_z**2

        a1 = P.mu/(P.mc + 2 * P.mr)

        self.kp_z = alpha0_z / b0_z
        self.kd_z = (alpha1_z - a1) / b0_z

        print("Kp_z: ", self.kp_z)
        print("Kd_z: ", self.kd_z)
        # -----------------------------------------

    def update(self,h_r,z_r,state):
        z = state[0][0]
        h = state[1][0]
        theta = state[2][0]
        zdot = state[3][0]
        hdot = state[4][0]
        thetadot = state[5][0]

        F_tilde = self.kp_h * (h_r - h) - self.kd_h * hdot

        theta_ref = self.kp_z * (z_r - z) - self.kd_z * zdot
        theta_ref = saturate(theta_ref, 10 * np.pi/180)
        tau = self.kp_th * (theta_ref - theta) - self.kd_th * thetadot
        tau = saturate(tau,P.fmax)
        return np.array([[saturate(F_tilde + P.Fe,2*P.fmax)], [tau]])

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u
        



