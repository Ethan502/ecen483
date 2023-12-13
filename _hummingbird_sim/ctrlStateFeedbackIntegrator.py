import numpy as np
import control as cnt
import hummingbirdParam as P


class ctrlStateFeedbackIntegrator:
    def __init__(self):
        #--------------------------------------------------
        # State Feedback Control Design
        #--------------------------------------------------
        # tuning parameters
        tr_th = 1.0
        wn_th = 2.2 / tr_th  
        M = 2  
        zeta_th = 0.707 
        pi_th = 1 
        tr_phi = 0.1
        tr_psi = tr_phi * M
        wn_psi = 2.2 / tr_psi     
        zeta_psi = 0.707
        wn_phi =  2.2 / tr_phi       
        zeta_phi = 0.707
        pi_psi = 0.5
        # hard code Ackerman's formula
        alpha1_lon = 2 * zeta_th * wn_th + pi_th
        alpha2_lon = 2 * zeta_th * wn_th * pi_th + wn_th**2
        alpha3_lon = wn_th**2*pi_th
        self.k_th = alpha2_lon / P.b_theta
        self.k_thdot = alpha1_lon / P.b_theta
        self.ki_lon = alpha3_lon / P.b_theta
        alpha1_lat = 2*wn_phi*zeta_phi+2*wn_psi*zeta_psi+pi_psi
        alpha2_lat = wn_phi**2 + 4*wn_phi*wn_psi*zeta_phi*zeta_psi + 2*wn_phi*pi_psi*zeta_phi+wn_psi**2+2*wn_psi*pi_psi*zeta_psi
        alpha3_lat = 2*wn_phi**2*wn_psi*zeta_psi+wn_phi**2*pi_psi+2*wn_phi*wn_psi**2*zeta_phi+4*wn_phi*wn_psi*pi_psi*zeta_phi*zeta_psi+wn_psi**2*pi_psi
        alpha4_lat = wn_phi**2*wn_psi**2+2*wn_phi**2*wn_psi*pi_psi*zeta_psi+2*wn_phi*wn_psi**2*pi_psi*zeta_phi
        alpha5_lat = wn_phi**2*wn_psi**2*pi_psi
        b1 = 1/P.J1x
        a1 = P.ellT*P.Fe/(P.JT+P.J1z)
        self.k_phi = alpha2_lat/b1
        self.k_psi = alpha4_lat/(a1*b1)
        self.k_phidot = alpha1_lat / b1
        self.k_psidot = alpha3_lat/(a1 * b1)
        self.ki_lat = alpha5_lat / (b1*a1)
        # print gains to terminal
        print('K_lon: [', self.k_th, ',', self.k_thdot, ']')
        print('ki_lon: ', self.ki_lon)         
        print('K_lat: [', self.k_phi, ',', self.k_psi, ',', self.k_phidot, ',', self.k_psidot, ']')
        print('ki_lat: ', self.ki_lat)        
        #--------------------------------------------------
        # saturation limits
        theta_max = 30.0 * np.pi / 180.0  # Max theta, rads
        #--------------------------------------------------
        self.Ts = P.Ts
        sigma = 0.05  # cutoff freq for dirty derivative
        self.beta = (2 * sigma - self.Ts) / (2 * sigma + self.Ts)
        self.phi_d1 = 0.
        self.phi_dot = 0.
        self.theta_d1 = 0.
        self.theta_dot = 0.
        self.psi_d1 = 0.
        self.psi_dot = 0.     
        self.phidot_d1 = 0.0
        self.psidot_d1 = 0.0
        self.thetadot_d1 = 0.0   
        # variables to implement integrator
        self.integrator_th = 0.0  
        self.error_th_d1 = 0.0  
        self.integrator_psi = 0.0  
        self.error_psi_d1 = 0.0 

    def update(self, r, y):
        theta_ref = r[0][0]
        psi_ref = r[1][0]
        phi = y[0][0]
        theta = y[1][0]
        psi = y[2][0]
        force_equilibrium = P.Fe * np.cos(theta) 
        # update differentiators
        self.phi_dot = self.beta * self.phidot_d1 + (1-self.beta) * ((phi-self.phi_d1) / P.Ts)
        self.phi_d1 = phi
        self.theta_dot = self.beta * self.thetadot_d1 + (1-self.beta) * ((theta-self.theta_d1) / P.Ts)
        self.theta_d1 = theta
        self.psi_dot = self.beta * self.psidot_d1 + (1-self.beta) * ((psi-self.psi_d1) / P.Ts) 
        self.psi_d1 = psi
        # integrate error
        error_th = theta_ref - theta
        error_psi = psi_ref - psi
        self.integrator_th = self.integrator_th + (self.Ts / 2.0) * (error_th + self.error_th_d1)
        self.integrator_psi = self.integrator_psi + (self.Ts / 2.0) * (error_psi + self.error_psi_d1)
        self.error_th_d1 = error_th
        self.error_psi_d1 = error_psi

        # longitudinal control
        force_unsat = -np.array([self.k_th,self.k_thdot]) @ np.array([theta,self.theta_dot]) + self.ki_lon * self.integrator_th
        force = saturate(force_unsat, -P.force_max, P.force_max)
        # lateral control
        torque_unsat = -np.array([self.k_phi,self.k_psi,self.k_phidot,self.k_psidot]) @ np.array([phi,psi,self.phi_dot,self.psi_dot]) + self.ki_lat * self.integrator_psi
        torque = saturate(torque_unsat, -P.torque_max, P.torque_max)
        # convert force and torque to pwm signals
        pwm = np.array([[force + torque / P.d],               # u_left
                      [force - torque / P.d]]) / (2 * P.km)   # r_right          
        pwm = saturate(pwm, 0, 1)

        self.phidot_d1 = self.phi_dot
        self.thetadot_d1 = self.theta_dot
        self.psidot_d1 = self.psi_dot


        return pwm, np.array([[0], [theta_ref], [psi_ref]])


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
