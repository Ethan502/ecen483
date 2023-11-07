import numpy as np
import VTOLParam as P

class ctrlPID:
    def __init__(self):
        tr_h = 2.0
        zeta = 0.707
        wn_h = 2.2/tr_h

        b0_h = 1/(P.mc + 2*P.mr)
        alpha1_h = 2 * zeta * wn_h
        alpha0_h = wn_h**2

        self.kd_h = alpha1_h / b0_h
        self.kp_h = alpha0_h / b0_h

        print("Kp_h: ", self.kp_h)
        print("Kd_h: ", self.kd_h)

        self.ki_h = 0.75
        self.h_d1 = P.h0
        self.hdot = P.hdot0
        self.integrator_h = 0.0
        self.error_h_d1 = 0.0

        # Inner Loop ----------------------------

        tr_th = 0.3
        wn_th = 2.2 / tr_th
        b0_th = 1/(P.Jc + 2 * P.mr * (P.d**2))
        alpha1_th = 2 * zeta * wn_th
        alpha0_th = wn_th**2

        self.kp_th = alpha0_th / b0_th
        self.kd_th = alpha1_th / b0_th

        print("Kp_th: ", self.kp_th)
        print("Kd_th: ", self.kd_th)

        self.theta_d1 = P.theta0
        self.thetadot = P.thetadot0

        DC_th = 1

        # Outer Loop -----------------------------
        M = 10
        tr_z = M/2 * tr_th
        wn_z = 2.2 / tr_z
        b0_z = -P.g
        alpha1_z = 2 * zeta * wn_z
        alpha0_z = wn_z**2

        a1 = P.mu/(P.mc + 2 * P.mr)

        self.kp_z = alpha0_z / b0_z
        self.kd_z = (alpha1_z - a1) / b0_z

        self.ki_z = 0.5
        self.z_d1 = P.z0
        self.zdot = P.zdot0
        self.integrator_z = 0.0
        self.error_z_d1 = 0.0

        print("Kp_z: ", self.kp_z)
        print("Kd_z: ", self.kd_z)
        # -----------------------------------------
        # Saturation values
        self.Fe = (P.mc + 2.0 * P.mr) * P.g
        self.theta_max = 10.0 * np.pi / 180.0

    def update(self,h_r,z_r,z,h,theta):
        # z = y[0][0]
        # h = y[1][0]
        # theta = y[2][0]
        # _____________________Altitude Control______________________
        # Find the error in h
        error_h = h_r - h
        # Derive hdot
        self.hdot = P.beta * self.hdot + (1-P.beta) * ((h-self.h_d1) / P.Ts)
        # Integrate the error
        if np.abs(self.hdot) < 0.03:
            self.integrator_h = self.integrator_h + (P.Ts/2) * (error_h + self.error_h_d1)
        # PID control on height for force
        F_tilde = self.kp_h * error_h + self.ki_h * self.integrator_h - self.kd_h * self.hdot
        F_unsat = self.Fe + F_tilde
        F = saturate(F_unsat,2*P.fmax)

        # _____________________Position Control______________________
        # Error in z
        error_z = z_r - z
        # Derive z
        self.zdot = P.beta * self.zdot + (1-P.beta) * ((z - self.z_d1) / P.Ts)
        # Integrate error in z
        self.integrator_z = self.integrator_z + (P.Ts/2) * (error_z - self.error_z_d1)
        # PID control for z
        theta_r_unsat = self.kp_z * error_z + self.ki_z * self.integrator_z - self.kd_z * self.zdot
        # saturate theta
        theta_r = saturate(theta_r_unsat,1.5*self.theta_max)
        # Anti-windup on the integrator
        # if self.ki_z != 0.0:
        #     self.integrator_z = self.integrator_z + P.Ts / self.ki_z * (theta_r - theta_r_unsat)
        
        # ___________________Pitch Control______________________________
        # Error in theta
        error_th = theta_r - theta
        # Derive theta
        self.thetadot = P.beta * self.thetadot + (1-P.beta) * ((theta-self.theta_d1) / P.Ts)
        # PD control on theta to get tau
        tau_unsat = self.kp_th * error_th - self.kd_th * self.thetadot
        # Saturate theta
        tau = saturate(tau_unsat, 2*P.fmax*P.d)

        # Update time-related vars
        self.error_h_d1 = error_h
        self.h_d1 = h
        self.error_z_d1 = error_z
        self.z_d1 = z
        self.theta_d1 = theta

        # Return the thrusts
        motor_thrusts = P.mixing @ np.array([[F],[tau]])
        return motor_thrusts


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u
        



