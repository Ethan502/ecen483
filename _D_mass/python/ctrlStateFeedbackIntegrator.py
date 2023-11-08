import numpy as np
import control as cnt
import massParam as P

class ctrlStateFeedback:
    def __init__(self):
        tr = 2
        zeta = 0.707
        w_n = 2.2/tr
        pI = np.array([-0.5])

        A = np.array([[0.0,1.0],
                      [-P.k/P.m,-1.0*P.b/P.m]])
        B = np.array([[0.0],
                      [1.0/P.m]])
        C = np.array([[1.0,0.0]])

        # Find the desired poles
        des_char_poly = np.convolve([1,2*zeta*w_n,w_n**2],np.poly(pI))
        des_poles = np.roots(des_char_poly)

        # Augment the matrices
        A1 = np.vstack((
            np.hstack((A, np.zeros((np.size(A,1),1)))),
            np.hstack((-C,0))
            ))
        B1 = np.vstack((B,0))

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

    def update(self,z_r,x):
        z = x[0][0]
        # Find the F_fl with feedback linearization
        F_fl = P.k * z
        # Compute the state feedback controller
        F_tilde = -self.K @ x + self.kr * z_r
        # Find the total force
        F_unsat = F_fl + F_tilde[0][0]
        # Saturate and return
        F = saturate(F_unsat,P.F_max)
        return F

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u