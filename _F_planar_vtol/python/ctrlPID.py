import numpy as np
import VTOLParam as P

class ctrlPID:
    def __init__(self):
        tr_h = 1.5
        zeta = 0.707
        wn_h = 2.2/tr_h

        b0_h = 1/(P.mc + 2*P.mr)
        alpha1_h = 2 * zeta * wn_h
        alpha0_h = wn_h**2

        self.kd_h = alpha1_h / b0_h
        self.kp_h = alpha0_h / b0_h

        print("Kp_h: ", self.kp_h)
        print("Kd_h: ", self.kd_h)

        self.ki_h = 0.5
        self.h_delayed = P.h0
        self.hdot_delayed = P.hdot0
        self.error_h_sum = 0.0
        self.error_h_sum_prev = 0.0
        self.error_h_prev = 0.0


        # Inner Loop ----------------------------

        tr_th = 0.5
        wn_th = 2.2 / tr_th
        b0_th = 1/(P.Jc + 2 * P.mr * (P.d**2))
        alpha1_th = 2 * zeta * wn_th
        alpha0_th = wn_th**2

        self.kp_th = alpha0_th / b0_th
        self.kd_th = alpha1_th / b0_th

        print("Kp_th: ", self.kp_th)
        print("Kd_th: ", self.kd_th)

        self.theta_delayed = P.theta0
        self.thetadot_delayed = P.thetadot0

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

        self.ki_z = 0.5
        self.z_delayed = P.z0
        self.zdot_delayed = P.zdot0
        self.error_z_sum = 0.0
        self.error_z_sum_prev = 0.0
        self.error_z_prev = 0.0

        print("Kp_z: ", self.kp_z)
        print("Kd_z: ", self.kd_z)
        # -----------------------------------------

    def update(self,h_r,z_r,h,z,theta):
        #z = state[0][0]
        #h = state[1][0]
        #theta = state[2][0]
        #zdot = state[3][0]
        #hdot = state[4][0]
        #thetadot = state[5][0]

        thetadot = P.beta * self.thetadot_delayed + (1-P.beta)/P.Ts * (theta-self.theta_delayed)

        hdot = P.beta * self.hdot_delayed + (1-P.beta)/P.Ts * (h-self.h_delayed)
        error_h = (h_r - h)
        self.error_h_sum += P.Ts/2 * (error_h + self.error_h_prev)

        zdot = P.beta * self.zdot_delayed + (1-P.beta)/P.Ts * (z - self.z_delayed)
        error_z = (z_r - z)
        self.error_z_sum += P.Ts/2 * (error_z + self.error_z_prev)

        self.h_delayed = h
        self.z_delayed = z
        self.theta_delayed = theta
        self.hdot_delayed = hdot
        self.zdot_delayed = zdot
        self.theta_delayed = thetadot
        self.error_h_prev = error_h
        self.error_z_prev = error_z

        # PID for the h
        F_tilde = self.kp_h * (h_r - h) - self.kd_h * hdot + self.ki_h * self.error_h_sum

        # PID for the z
        theta_ref = self.kp_z * (z_r - z) - self.kd_z * zdot + self.ki_z * self.error_z_sum
        theta_ref = saturate(theta_ref, 10 * np.pi/180)
        tau = self.kp_th * (theta_ref - theta) - self.kd_th * thetadot
        tau = saturate(tau,P.fmax)
        return np.array([[saturate(F_tilde + P.Fe,2*P.fmax)], [tau]])

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u
        



