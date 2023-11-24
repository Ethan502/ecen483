import numpy as np
import control as cnt
import massParam as P

class ctrlDisturbanceObserver:
    def __init__(self):
        tr = 1.0
        zeta = 0.707
        w_n = 2.2/tr
        wn_obs = 10
        dist_obs_pole = 5.5 # something to tune
        pI = np.array([-0.62])

        A = np.array([[0.0,1.0],
                      [-P.k/P.m,-1.0*P.b/P.m]])
        B = np.array([[0.0],
                      [1.0/P.m]])
        C = np.array([[1.0,0.0]])

        # Find the desired poles
        des_char_poly = np.convolve([1,2*zeta*w_n,w_n**2],np.poly(pI))
        des_poles = np.roots(des_char_poly)

        A1 = np.vstack((np.hstack((A, np.zeros((np.size(A,1),1)))), 
                        np.hstack((-C, np.array([[0.0]]))) ))
        B1 = np.vstack( (self.B, 0.0) )

        # Find the gains, but only if the system in controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1,B1)) != np.size(A1,1):
            print("The system is not controllable")
        else:
            K1 = cnt.place(A1,B1,des_poles)
            self.K = K1[0][0:2]
            self.ki = K1[0][2]

        # Observer Design
        # Augmented Matrices
        A2 = np.concatenate((
                            np.concatenate((A,B),axis=1),
                            np.zeros((1,3))),
                            axis=0)
        B2 = np.concatenate((B, np.zeros((1,1))),axis=0)
        C2 = np.concatenate((C,np.zeros((1,1))),axis=1)

        # Get some poles
        des_char_est = np.array([1.0, 2.0*zeta*wn_obs, wn_obs**2])
        des_obs_char_poly = np.convolve(des_char_est,[1,dist_obs_pole])

        des_obs_poles = np.poles(des_obs_char_poly)
        # Compute the gains if the system is controllable
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A2.T, C2.T)) != 3:
            print("The system is not observerable")
        else:
            L2 = cnt.acker(A2.T, C2.T, des_obs_poles).T
        print('K: ', self.K)
        print('ki ', self.ki)
        print('L^T: ', L2.T)
        # --------------------------------------------------------------------------
        # Variables for the integrator
        self.integrator = 0.0
        self.error_d1 = 0.0
        self.obs_state = np.array([
            [0.0],  # F_hat_0
            [0.0],  # Fdot_hat_0
            [0.0],  # estimate of the disturbance
        ])
        self.F_d1 = 0.0 # control force, delayed for 1 sample
        self.L = L2
        self.A = A2
        self.B = B1
        self.C = C2


    def update(self,z_r,y):
        # update the observer
        x_hat, d_hat = self.update_observer(y)
        z_hat = x_hat[0][0]
        # integrate the error
        error = z_r - z_hat
        self.integrator = self.integrator + (P.Ts / 2.0) * (error + self.error_d1)
        self.error_d1 = error
        # compute the feedback linearizing force
        F_fl = P.k * z_hat
        # Compute the state feedback controller
        F_tilde = -self.K @ x_hat - self.ki * self.integrator - d_hat

        # Saturate the total force and return
        F = saturate(F_fl + F_tilde,P.F_max)
        self.F_d1 = F
        return F, x_hat
    
    def update_observer(self,y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f(self.x_hat, y_m)
        F2 = self.observer_f(self.x_hat + P.Ts / 2 * F1, y_m)
        F3 = self.observer_f(self.x_hat + P.Ts / 2 * F2, y_m)
        F4 = self.observer_f(self.x_hat + P.Ts * F3, y_m)
        self.obs_state += P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        x_hat = self.obs_state[0:2]
        d_hat = self.obs_state[2][0]
        return x_hat, d_hat


    def observer_f(self, x_hat, y_m):
        F_fl = P.k * y_m[0][0]
        # xhatdot = A*(xhat-xe) + B*(u-ue) + L(y-C*xhat)
        xhat_dot = self.A @ x_hat \
                 + self.B * (self.F_d1 - F_fl)\
                 + self.L * (y_m - self.C @ x_hat)
        return xhat_dot


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u