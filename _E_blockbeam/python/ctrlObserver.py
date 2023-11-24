import numpy as np
import control as cnt
import blockbeamParam as P
from scipy import signal

class ctrlObserver:
    def __init__(self):
        tr_th = 0.2
        zeta_th = 0.707
        wn_th = 2.2 / tr_th

        M = 10
        tr_z = M * tr_th
        zeta_z = 0.707
        wn_z = 2.2 / tr_z

        z_e = P.length/2.0

        pI = np.array([-2.0])

        tr_z_obs = tr_z / 10.0
        tr_th_obs = tr_th / 10.0

        self.A = np.array([
            [0.0,0.0,1.0,0.0],
            [0.0,0.0,0.0,1.0],
            [0.0,-1*P.g,0.0,0.0],
            [-1*(P.m1*P.g)/(((P.m2*P.length**2)/3)+P.m1*z_e**2),0.0,0.0,0.0]
            ])
        self.B = np.array([[0.0],
                      [0.0],
                      [0.0],
                      [P.length/(((P.m2*P.length**2)/3)+P.m1*z_e**2)]
                    ])

        self.C = np.array([[1.0,0.0,0.0,0.0],
                      [0.0,1.0,0.0,0.0]])
        
        self.Cr = np.array([[1.0, 0.0, 0.0, 0.0]])

        # Augment the matrices to add the integrator
        A1 = np.vstack((
                np.hstack((self.A, np.zeros((4,1)))),
                np.hstack((-self.Cr, np.zeros((1,1))))
                ))
        B1 = np.vstack((self.B, np.zeros((1,1))))

        # Get the desired poles
        des_char_poly = np.convolve(np.convolve(
            [1,2*zeta_z*wn_z,wn_z**2],
            [1,2*zeta_th*wn_th,wn_th**2]
            ),
            np.poly(pI))
        des_poles = np.roots(des_char_poly)

        # Compute the gains of the system
        if np.linalg.matrix_rank(cnt.ctrb(self.A, self.B)) != 4:
            print("The system is not controllable")
        else:
            K1 = cnt.acker(A1, B1, des_poles)
            self.K = K1[0][0:4]
            self.ki = K1[0][4]
        print('K: ', self.K)
        print('ki: ', self.ki)

        wn_z_obs = 2.2 / tr_z_obs
        wn_th_obs = 2.2 / tr_th_obs
        des_obs_char_poly = np.convolve([1, 2*zeta_z*wn_z_obs, wn_z_obs**2],
                                        [1, 2*zeta_th*wn_th_obs, wn_th_obs**2])
        des_obs_poles = np.roots(des_obs_char_poly)
        # Compute the observer gains if the system in controllable
        if np.linalg.matrix_rank(cnt.ctrb(self.A.T, self.C.T)) != 4:
            print("The system is not observable")
        else:
            self.L = signal.place_poles(self.A.T, self.C.T, 
                                        des_obs_poles).gain_matrix.T
        print('K: ', self.K)
        print('ki ', self.ki)
        print('L^T: ', self.L.T)
        #-------------------------------------------------------------------------------
        # Variables for the integral calculation
        self.integrator_z = 0.0
        self.error_z_d1 = 0.0
        self.x_hat = np.array([
            [0.0],  # initial estimate for z_hat
            [0.0],  # initial estimate for theta_hat
            [0.0],  # initial estimate for z_hat_dot
            [0.0]])  # initial estimate for theta_hat_dot
        self.F_d1 = 0.0  
    
    def update(self, z_r, y):
        x_hat = self.update_observer(y)
        z_hat = self.Cr @ x_hat
        # Integrate the error
        error_z = z_r - z_hat
        self.integrator_z = self.integrator_z \
            + (P.Ts / 2.0) * (error_z + self.error_z_d1)
        self.error_z_d1 = error_z
        # Linearized state
        x_e = np.array([[P.ze],[0.0],[0.0],[0.0]])
        x_tilde = x_hat - x_e
        # Compute state feedback controller
        F_e = (P.m1*P.g / P.length) * z_hat + P.m2*P.g/2.0
        F_tilde = -self.K @ x_tilde - self.ki * self.integrator_z
        F_unsat = F_e + F_tilde
        F = saturate(F_unsat[0][0], P.Fmax)
        self.anti_windup(F,F_unsat)
        self.F_d1 = F
        return F, x_hat

    def update_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f(self.x_hat, y_m)
        F2 = self.observer_f(self.x_hat + P.Ts / 2 * F1, y_m)
        F3 = self.observer_f(self.x_hat + P.Ts / 2 * F2, y_m)
        F4 = self.observer_f(self.x_hat + P.Ts * F3, y_m)
        self.x_hat += P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        return self.x_hat

    def observer_f(self, x_hat, y):
        # xhatdot = A*xhat + B*u + L(y-C*xhat)
        xhat_dot = self.A @ x_hat \
                   + self.B * self.F_d1 \
                   + self.L @ (y-self.C @ x_hat)
        return xhat_dot
    
    def anti_windup(self, F, F_unsat):
        if self.ki != 0.0:
            self.integrator_z = self.integrator_z + P.Ts/self.ki*(F-F_unsat)

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u

