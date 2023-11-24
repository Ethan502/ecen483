import numpy as np
import control as cnt
import massParam as P

class ctrlStateFeedbackIntegrator:
    def __init__(self):
        tr = 1.0
        zeta = 0.707
        w_n = 2.2/tr
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
        B1 = np.vstack( (B, 0.0) )

        print(A1)
        print(B1)

        # Find the gains, but only if the system in controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1,B1)) != np.size(A1,1):
            print("The system is not controllable")
        else:
            K1 = cnt.place(A1,B1,des_poles)
            self.K = K1[0][0:2]
            self.ki = K1[0][2]
        print('K: ', self.K)
        print('ki: ', self.ki)
        print('Desired poles: ', des_poles)

        # variables for the integrator
        self.integrator = 0.0
        self.error_d1 = 0.0

    def update(self,z_r,x):
        z = x[0][0]
        # find the error in z
        error = z_r - z
        # integrate the error
        self.integrator = self.integrator + (P.Ts / 2.0) * (error + self.error_d1)
        self.error_d1 = error
        # compute the feedback linearizing force
        F_fl = P.k * z
        # Use PID to find the F tilde
        F_tilde = -self.K @ x - self.ki * self.integrator
        # Saturate the total force and return
        F = saturate(F_fl + F_tilde,P.F_max)
        return F


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u