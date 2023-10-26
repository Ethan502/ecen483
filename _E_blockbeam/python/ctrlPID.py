import numpy as np
import blockbeamParam as P

class ctrlPID:
    def __init__(self):

        # Inner Loop ___________________________________________________________________
        tr_th = 1
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

        # _________________________Outer Loop __________________________________

        self.z_delayed = P.z0
        self.theta_delayed = P.theta0
        self.zdot_delayed = P.zdot0
        self.thetadot_delayed = P.thetadot0
        self.error_sum = 0.0
        self.error_sum_prev = 0.0
        self.error_z_prev = 0.0

        tr_z = M * tr_th
        b0_z = P.g
        wn_z = 2.2/tr_z

        alpha1_z = 2 * zeta * wn_z
        alpha0_z = wn_z**2

        self.kp_z = -alpha0_z / b0_z
        self.kd_z = -alpha1_z / b0_z
        self.ki_z = 0.0

        print("Kp_z: ", self.kp_z)
        print("Kd_z: ", self.kd_z)

        k_DC = (b0_th*self.kp_th)/(b0_th*self.kp_th)


    def update(self, z_r, state):
        z = state[0][0]
        theta = state[1][0]
        #zdot = state[2][0]
        #thetadot = state[3][0]

        # Calculate Derivatives
        zdot = P.beta * self.zdot_delayed + (1-P.beta)/P.Ts * (z - self.z_delayed)
        thetadot = P.beta * self.theta_delayed + (1-P.beta)/P.Ts * (theta - self.theta_delayed)

        # Integrate Error
        error_z = (z_r - z)
        if np.abs(zdot) < 0.05:
            self.error_sum = self.error_sum + P.Ts/2 * (error_z + self.error_z_prev)

        # Do the PID on the theta_ref
        theta_ref = self.kp_z * (z_r - z) - self.kd_z * zdot + self.ki_z * self.error_sum

        F = self.kp_th * (theta_ref - theta) - self.kd_th * thetadot
        Fe = P.m1*P.g*z / P.length + 0.5*P.m2*P.g

        F_sat = saturate(F + Fe,P.Fmax)

        self.z_delayed = z
        self.theta_delayed = theta
        self.zdot_delayed = zdot
        self.thetadot_delayed = thetadot
        self.error_z_prev = error_z

        return F_sat

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u    
