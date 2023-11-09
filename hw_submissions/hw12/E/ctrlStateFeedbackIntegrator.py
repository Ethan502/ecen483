import numpy as np
import control as cnt
import blockbeamParam as P

class ctrlStateFeedbackIntegrator:
    def __init__(self):
        tr_th = 0.15
        zeta_th = 0.707
        wn_th = 2.2 / tr_th

        M = 10
        tr_z = M * tr_th
        zeta_z = 0.707
        wn_z = 2.2 / tr_z

        z_e = P.length/2.0

        pI = np.array([-2.0])

        A = np.array([
            [0.0,0.0,1.0,0.0],
            [0.0,0.0,0.0,1.0],
            [0.0,-1*P.g,0.0,0.0],
            [-1*(P.m1*P.g)/(((P.m2*P.length**2)/3)+P.m1*z_e**2),0.0,0.0,0.0]
            ])
        B = np.array([[0.0],
                      [0.0],
                      [0.0],
                      [P.length/(((P.m2*P.length**2)/3)+P.m1*z_e**2)]
                    ])

        C = np.array([[1.0,0.0,0.0,0.0],
                      [0.0,1.0,0.0,0.0]])
        
        Cr = np.array([[1.0, 0.0, 0.0, 0.0]])

        # Augment the matrices to add the integrator
        A1 = np.vstack((
                np.hstack((A, np.zeros((4,1)))),
                np.hstack((-Cr, np.zeros((1,1))))
                ))
        B1 = np.vstack((B, np.zeros((1,1))))

        # Get the desired poles
        des_char_poly = np.convolve(np.convolve(
            [1,2*zeta_th*wn_th,wn_th**2],
            [1,2*zeta_z*wn_z,wn_z**2]
            ),
            np.poly(pI))
        des_poles = np.roots(des_char_poly)

        # Compute the gains of the system
        if np.linalg.matrix_rank(cnt.ctrb(A, B)) != 4:
            print("The system is not controllable")
        else:
            K1 = cnt.acker(A1, B1, des_poles)
            self.K = K1[0][0:4]
            self.ki = K1[0][4]
        print('K: ', self.K)
        print('ki: ', self.ki)

        # Variables for the integral calculation
        self.integrator_z = 0.0
        self.error_z_d1 = 0.0

    def update(self,z_r,x):
        z = x[0][0]
        # Get the error in z
        error_z = z_r - z
        # Integrate the error
        self.integrator_z = self.integrator_z + (P.Ts / 2.0) * (error_z + self.error_z_d1)
        self.error_z_d1 = error_z
        # Find the equilibrium force
        F_e = (P.m1*P.g / P.length) * z + P.m2*P.g/2.0
        # Compute the state feedback controller
        F_tilde = -self.K @ x - self.ki * self.integrator_z
        F_unsat = F_e + F_tilde
        # Saturate and return
        F = saturate(F_unsat[0],P.Fmax)
        return F

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u