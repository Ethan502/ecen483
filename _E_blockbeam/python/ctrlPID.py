import numpy as np
import blockbeamParam as P

class ctrlPID:
    def __init__(self):

        # Inner Loop ___________________________________________________________________
        tr_th = 0.15
        zeta = 0.707
        M = 10
        wn_th = 2.2/tr_th

        self.theta_max = 10.0 * np.pi / 180.0
        
        b0_th = P.length/((P.m2*(P.length**2))/3.0 + P.m1 * (P.ze**2))

        alpha1_th = 2 * zeta * wn_th 
        alpha0_th = wn_th**2

        self.kd_th = alpha1_th / b0_th
        self.kp_th = alpha0_th / b0_th

        print("Kp_th: ", self.kp_th)
        print("Kd_th: ", self.kd_th)

        # _________________________Outer Loop __________________________________
        tr_z = M * tr_th
        b0_z = P.g
        wn_z = 2.2/tr_z

        alpha1_z = 2 * zeta * wn_z
        alpha0_z = wn_z**2

        self.kp_z = -alpha0_z / b0_z
        self.kd_z = -alpha1_z / b0_z
        self.ki_z = -0.12

        print("Kp_z: ", self.kp_z)
        print("Kd_z: ", self.kd_z)

        k_DC = (b0_th*self.kp_th)/(b0_th*self.kp_th)

        #_______________________Variables for math___________________________
        self.thetadot = 0
        self.theta_d1 = 0.0
        self.zdot = 0.0
        self.z_d1 = 0.0
        self.integrator_z = 0
        self.error_z_d1 = 0


    def update(self,z_r,y):
        z = y[0][0]
        theta = y[1][0]

        # ___________________ Outer Loop (z-control)_____________________________
        # Get the error in z
        error_z = (z_r - z)
        # Calculate Derivatives
        self.zdot = P.beta * self.zdot + (1-P.beta)/P.Ts * (z - self.z_d1)
        # Integrate Error
        if np.abs(self.zdot) < 0.04:
            self.integrator_z = self.integrator_z + (P.Ts/2) * (error_z + self.error_z_d1)
        # Do the PID on the theta_ref
        theta_ref = self.kp_z * error_z - self.kd_z * self.zdot + self.ki_z * self.integrator_z
        # Saturate the theta
        theta_ref = saturate(theta_ref,self.theta_max)

        # ___________________Inner Loop (theta-control)___________________________
        error_th = theta_ref - theta
        # Differentiate theta
        self.thetadot = P.beta * self.thetadot + (1-P.beta)/P.Ts * (theta - self.theta_d1)
        # PD control on theta
        F_tilde = self.kp_th * error_th - self.kd_th * self.thetadot
        # Feedback linearizing force
        Fe = P.m1*P.g*(z / P.length) + 0.5*P.m2*P.g
        # Saturate the force
        F = saturate(F_tilde + Fe,P.Fmax)
        # Update variables
        self.error_z_d1 = error_z
        self.z_d1 = z
        self.theta_d1 = theta
        # Return the computed force
        return F

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u    
