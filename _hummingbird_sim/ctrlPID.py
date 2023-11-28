import numpy as np
import hummingbirdParam as P


class ctrlPID:
    def __init__(self):

        #######################################################################################
        #                               Longitude Control                                     #
        #######################################################################################
        # tuning parameters
        tr_pitch = 0.3
        zeta_pitch = 0.707
        self.ki_pitch = -0.013
        # gain calculation
        b_theta = P.ellT/(P.m1 * P.ell1**2 + P.m2 * P.ell2**2 + P.J1y + P.J2y)
        #print('b_theta: ', b_theta)
        wn_pitch = 2.2 / tr_pitch
        self.kp_pitch = (wn_pitch**2)/b_theta
        self.kd_pitch = (2.0 * zeta_pitch * wn_pitch)/b_theta
        # print gains to terminal
        print('kp_pitch: ', self.kp_pitch)
        print('ki_pitch: ', self.ki_pitch)
        print('kd_pitch: ', self.kd_pitch) 
        # sample rate of the controller
        self.Ts = P.Ts
        # dirty derivative parameters
        sigma = 0.05  # cutoff freq for dirty derivative
        self.beta = (2 * sigma - self.Ts) / (2 * sigma + self.Ts)
        # delayed variables
        self.theta_d1 = P.theta0
        self.theta_dot = P.thetadot0
        self.integrator_theta = 0.0
        self.error_theta_d1 = 0.0  # pitch error delayed by 1

        #######################################################################################
        #                               Latitude Control                                      #
        #######################################################################################
        
        # ---------------------- INNER LOOP (ROLL)-----------------------------------
        tr_phi = 0.1
        zeta_phi = 0.707
        # gain calculations
        wn_phi = 2.2 / tr_phi
        self.kp_phi = (wn_phi**2)*P.J1x
        self.kd_phi = (2 * zeta_phi * wn_phi)*P.J1x
        print('kp_phi: ', self.kp_phi)
        print('kd_phi: ', self.kd_phi)
        # delayed variables
        self.phi_d1 = P.phi0
        self.phi_dot = P.phidot0

        # ---------------------- OUTER LOOP (YAW)--------------------------------------
        M = 10
        tr_psi = tr_phi * M
        zeta_psi = 0.707
        self.ki_psi = 0.001
        # gain calculations
        b_psi = (P.ellT * P.Fe)/(P.JT + P.J1z)
        wn_psi = 2.2 / tr_psi
        self.kp_psi = (wn_psi**2)/b_psi
        self.kd_psi = (2 * zeta_psi * wn_psi)/b_psi
        print('kp_psi: ', self.kp_psi)
        print('ki_psi: ', self.ki_psi)
        print('kd_psi: ', self.kd_psi)

        # delayed variables
        self.psi_d1 = P.psi0
        self.psi_dot = P.psidot0
        self.integrator_psi = 0.0
        self.error_psi_d1 = 0.0




    def update(self, r, y):
        theta_ref = r[0][0]
        psi_ref = r[1][0]
        phi = y[0][0]
        theta = y[1][0]
        psi = y[2][0]
        # -------------------------------Pitch control------------------------------------------
        error_theta = theta_ref - theta
        self.theta_dot = self.beta * self.theta_dot + (1-self.beta)/P.Ts * (theta - self.theta_d1)
        if np.abs(self.theta_dot) < 0.02:
            self.integrator_theta += (P.Ts/2) * (error_theta + self.error_theta_d1)
        force_fl = (P.m1*P.ell1 + P.m2*P.ell2)*(P.g/P.ellT)*np.cos(theta)
        force_unsat = force_fl + (error_theta * self.kp_pitch - self.kd_pitch * self.theta_dot + self.ki_pitch * self.integrator_theta)
        force = saturate(force_unsat, -P.force_max, P.force_max)
        # -------------------------------Roll/yaw control (outer)-------------------------------------------
        error_psi = psi_ref - psi
        self.psi_dot = self.beta * self.psi_dot + (1-self.beta)/P.Ts * (psi - self.psi_d1)
        if np.abs(self.psi_dot) < 0.03:
            self.integrator_psi += (P.Ts/2.0) * (error_psi + self.error_psi_d1)
        phi_ref = self.kp_psi * error_psi - self.kd_psi * self.psi_dot + self.ki_psi * self.integrator_psi
        #--------------------------------Roll/yaw control (inner)-------------------------------------------
        error_phi = phi_ref - phi
        self.phi_dot = self.beta * self.phi_dot + (1-self.beta)/P.Ts * (phi - self.phi_d1)
        torque_unsat = self.kp_phi * error_phi - self.kd_phi * self.phi_dot
        torque = saturate(torque_unsat,-P.torque_max, P.torque_max)
        
        
        # convert force and torque to pwm signals
        pwm = np.array([[force + torque / P.d],               # u_left
                      [force - torque / P.d]]) / (2 * P.km)   # r_right          
        pwm = saturate(pwm, 0, 1)

        # update all delayed variables
        self.theta_d1 = theta
        self.error_theta_d1 = error_theta
        self.phi_d1 = phi
        self.psi_d1 = psi
        self.error_psi_d1 = error_psi
        # return pwm plus reference signals
        return pwm, np.array([[phi_ref], [theta_ref], [psi_ref]])


def saturate(u, low_limit, up_limit):
    if isinstance(u, float) is True:
        if u > up_limit:
            u = up_limit
        if u < low_limit:
            u = low_limit
    else:
        for i in range(0, u.shape[0]):
            if u[i][0] > up_limit:
                u[i][0] = up_limit
            if u[i][0] < low_limit:
                u[i][0] = low_limit
    return u




